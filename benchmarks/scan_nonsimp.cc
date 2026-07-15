// Scan a sample of non-simplicial-iDT isomers from
// claude-projects/delaunay/non_simplical.spirals and report:
//   - input iDT topology (n_self, n_multi)
//   - solver outcome: succeeded vs FAILED-by-which-invariant
//   - final T(0) topology (n_self, n_multi)
//
// Usage:
//   scan_nonsimp <spiral_file> --all [--out <failures.tsv>]
//   scan_nonsimp <spiral_file> --sample <n> [--seed <s>] [--out <failures.tsv>]
//
// The --out file accumulates one TSV row per failing isomer:
//   N\tinput_self\tinput_multi\tfinal_self\tfinal_multi\tkappa\tflips\tsteps\tfail_type\tname
// flushed after every failure so a kill -9 keeps prior data.

#include "fullerenes/triangulation.hh"
#include "fullerenes/delaunay.hh"
#include "fullerenes/delaunay_alexandrov.hh"
#include "fullerenes/spiral.hh"

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <string>
#include <fstream>
#include <vector>
#include <random>
#include <map>

using namespace std;

struct EdgeStats {
  int n_self_loops = 0;
  int n_multi_pairs = 0;
};

static EdgeStats inspect(const DelaunayTriangulation& T) {
  EdgeStats s;
  map<pair<int,int>, int> by_pair;
  for (int h = 0; h < T.nh; h += 2) {
    if (!T.alive(h)) continue;
    int u = T.he_origin[h], v = T.dest(h);
    if (u == v) s.n_self_loops++;
    else        by_pair[{min(u,v), max(u,v)}]++;
  }
  for (auto& [k, c] : by_pair) if (c > 1) s.n_multi_pairs++;
  return s;
}

struct Outcome {
  enum Kind { OK, FAIL_KAPPA, FAIL_T0_NONSIMPLICIAL, FAIL_TBAR_NONSIMPLE, FAIL_BUILD };
  Kind kind;
  EdgeStats input;
  EdgeStats final;
  double final_kappa = 0;
  int n_flips = 0;
  int n_steps = 0;
  int N = 0;            // number of carbon atoms
  string name;
};

static bool valid_spiral_string(const string& s) {
  // Cheap sanity check before handing to the parser, which may call
  // terminate()/abort() on malformed input.  Expected shape:
  //   "C<N>-[<scheme>:<comma-list>]-fullerene" or with -[jumps] segment.
  if (s.size() < 12 || s[0] != 'C') return false;
  // Must contain the closing "]-fullerene".
  if (s.find("]-fullerene") == string::npos) return false;
  // Must have an opening '[' and matching close.
  size_t lb = s.find('['), rb = s.rfind(']');
  if (lb == string::npos || rb == string::npos || lb >= rb) return false;
  return true;
}

static Outcome run_one(const string& spiral_name) {
  Outcome o;
  o.name = spiral_name;
  if (!valid_spiral_string(spiral_name)) {
    o.kind = Outcome::FAIL_BUILD;
    return o;
  }
  try {
    spiral_nomenclature sn(spiral_name);
    Triangulation T(sn);
    o.N = (int)T.N;       // dual nv = N/2 + 2
    auto D = DelaunayTriangulation::compute(T);
    o.input = inspect(D);

    AlexandrovSolver solver;
    solver.D = D;
    auto coords = solver.solve();
    o.final = inspect(solver.D);
    o.final_kappa = solver.stats_final_kappa;
    o.n_flips = solver.stats_flips;
    o.n_steps = solver.stats_steps;

    // solve() now always returns positions when reconstruct succeeded;
    // status enum is the source of truth.
    using S = AlexandrovSolver::ValidationStatus;
    switch (solver.stats_status) {
      case S::OK:                       o.kind = Outcome::OK; break;
      case S::FAIL_KAPPA_NOT_CONVERGED: o.kind = Outcome::FAIL_KAPPA; break;
      case S::FAIL_NOT_SIMPLE:
        o.kind = solver.stats_t0_simplicial ? Outcome::FAIL_TBAR_NONSIMPLE
                                              : Outcome::FAIL_T0_NONSIMPLICIAL;
        break;
      case S::FAIL_RECONSTRUCT:
      case S::FAIL_VOLUME_DEGENERATE:
      case S::FAIL_SELF_INTERSECTING:
      case S::FAIL_NOT_CONVEX:
        o.kind = Outcome::FAIL_T0_NONSIMPLICIAL;  // legacy bucket; new finer
                                                   // labels live in scan_diag_max
        break;
    }
  } catch (...) {
    o.kind = Outcome::FAIL_BUILD;
  }
  return o;
}

static const char* kind_str(Outcome::Kind k) {
  switch (k) {
    case Outcome::OK:                    return "OK";
    case Outcome::FAIL_KAPPA:            return "FAIL_KAPPA";
    case Outcome::FAIL_T0_NONSIMPLICIAL: return "FAIL_T0_NONSIMPLICIAL";
    case Outcome::FAIL_TBAR_NONSIMPLE:   return "FAIL_TBAR_NONSIMPLE";
    case Outcome::FAIL_BUILD:            return "FAIL_BUILD";
  }
  return "UNKNOWN";
}

int main(int argc, char** argv) {
  string spiral_file, out_file;
  int n_sample = -1;       // -1 = all
  int seed = 42;

  for (int i = 1; i < argc; i++) {
    string a = argv[i];
    if (a == "--all")            n_sample = -1;
    else if (a == "--sample" && i + 1 < argc) n_sample = atoi(argv[++i]);
    else if (a == "--seed"   && i + 1 < argc) seed     = atoi(argv[++i]);
    else if (a == "--out"    && i + 1 < argc) out_file = argv[++i];
    else if (spiral_file.empty()) spiral_file = a;
    else { fprintf(stderr, "unknown arg: %s\n", argv[i]); return 1; }
  }
  if (spiral_file.empty()) {
    fprintf(stderr,
      "Usage:\n"
      "  %s <spiral_file> --all [--out <failures.tsv>]\n"
      "  %s <spiral_file> --sample <n> [--seed <s>] [--out <failures.tsv>]\n",
      argv[0], argv[0]);
    return 1;
  }

  ifstream in(spiral_file);
  if (!in) { fprintf(stderr, "cannot open %s\n", spiral_file.c_str()); return 1; }

  vector<string> all;
  string line;
  while (getline(in, line)) if (!line.empty()) all.push_back(line);
  fprintf(stderr, "loaded %zu isomers from %s\n", all.size(), spiral_file.c_str());

  vector<int> idx;
  if (n_sample < 0 || n_sample >= (int)all.size()) {
    idx.resize(all.size());
    for (size_t i = 0; i < all.size(); i++) idx[i] = (int)i;
  } else {
    mt19937 rng(seed);
    vector<int> all_idx(all.size());
    for (size_t i = 0; i < all.size(); i++) all_idx[i] = (int)i;
    shuffle(all_idx.begin(), all_idx.end(), rng);
    idx.assign(all_idx.begin(), all_idx.begin() + n_sample);
  }

  FILE* outfp = nullptr;
  if (!out_file.empty()) {
    outfp = fopen(out_file.c_str(), "w");
    if (!outfp) { fprintf(stderr, "cannot open %s for writing\n", out_file.c_str()); return 1; }
    // setvbuf MUST be called before any I/O operation on the stream (POSIX).
    setvbuf(outfp, nullptr, _IOLBF, 0);   // line-buffered: survive kill -9
    fprintf(outfp, "#N\tin_self\tin_multi\tout_self\tout_multi\tkappa\tflips\tsteps\tfail_type\tname\n");
  }

  int n_ok = 0, n_kappa = 0, n_t0 = 0, n_tbar = 0, n_build = 0;

  for (size_t i = 0; i < idx.size(); i++) {
    if (i % 1000 == 0)
      fprintf(stderr, "  #%zu / %zu  (ok=%d, t0=%d, tbar=%d, kappa=%d, build=%d)\n",
              i, idx.size(), n_ok, n_t0, n_tbar, n_kappa, n_build);
    auto o = run_one(all[idx[i]]);
    switch (o.kind) {
      case Outcome::OK:                     n_ok++;     break;
      case Outcome::FAIL_KAPPA:             n_kappa++;  break;
      case Outcome::FAIL_T0_NONSIMPLICIAL:  n_t0++;     break;
      case Outcome::FAIL_TBAR_NONSIMPLE:    n_tbar++;   break;
      case Outcome::FAIL_BUILD:             n_build++;  break;
    }
    if (outfp && o.kind != Outcome::OK) {
      fprintf(outfp, "%d\t%d\t%d\t%d\t%d\t%.6e\t%d\t%d\t%s\t%s\n",
              o.N, o.input.n_self_loops, o.input.n_multi_pairs,
              o.final.n_self_loops, o.final.n_multi_pairs,
              o.final_kappa, o.n_flips, o.n_steps,
              kind_str(o.kind), o.name.c_str());
    }
  }

  printf("\n=== Summary ===\n");
  printf("Total scanned:               %d\n", (int)idx.size());
  printf("OK (solver succeeded):       %d (%.2f%%)\n", n_ok, 100.0*n_ok/idx.size());
  printf("FAIL: T(0) non-simplicial:   %d (%.2f%%)  [I-1 violation, PALC misconverged]\n",
         n_t0, 100.0*n_t0/idx.size());
  printf("FAIL: T̄(0) non-simple poly:  %d (%.2f%%)\n", n_tbar, 100.0*n_tbar/idx.size());
  printf("FAIL: |κ| > 0.01:            %d (%.2f%%)\n", n_kappa, 100.0*n_kappa/idx.size());
  printf("FAIL: build error:           %d\n", n_build);
  if (outfp) {
    fclose(outfp);
    printf("\nFailures (one TSV row each) written to: %s\n", out_file.c_str());
  }
  return 0;
}
