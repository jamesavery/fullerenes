// Per-isomer diagnostic-summary scanner: runs the Alexandrov solver with
// record_diag=true on each input spiral, post-processes the per-step
// trajectory, and emits one TSV row per isomer with summary scalars.
//
// Designed to characterize the distribution of "near-π edge counts"
// across the full 1.03M non-simplicial-iDT database.  OpenMP-parallel,
// line-buffered, resumable.
//
// Usage:
//   scan_diag_max <spiral_file> --out summary.tsv [--start <idx>] [--limit <n>]
//
// TSV columns:
//   N  max_n_near_pi_001  max_n_near_pi_0001  final_n_near_pi_001
//   final_n_near_pi_0001  final_min_h_sq  kind  name
//
// `kind` mirrors scan_nonsimp: OK / FAIL_KAPPA / FAIL_T0_NONSIMPLICIAL /
// FAIL_TBAR_NONSIMPLE / FAIL_BUILD.

#include "fullerenes/triangulation.hh"
#include "fullerenes/delaunay.hh"
#include "fullerenes/delaunay_alexandrov.hh"
#include "fullerenes/spiral.hh"

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <string>
#include <vector>
#include <atomic>
#include <mutex>

#ifdef _OPENMP
#include <omp.h>
#endif

using namespace std;

struct Summary {
  int N = 0;
  int max_n_near_pi_001 = 0;
  int max_n_near_pi_0001 = 0;
  int final_n_near_pi_001 = 0;
  int final_n_near_pi_0001 = 0;
  double final_min_h_sq = 0;
  double volume = 0;
  // Volume normalized by typical-edge-cubed: V / (mean_edge_length^3).
  // Removes overall-scale dependence so volumes across different N
  // are directly comparable.
  double volume_normalized = 0;
  string kind;
  string name;
};

static bool valid_spiral_string(const string& s) {
  if (s.size() < 12 || s[0] != 'C') return false;
  if (s.find("]-fullerene") == string::npos) return false;
  size_t lb = s.find('['), rb = s.rfind(']');
  if (lb == string::npos || rb == string::npos || lb >= rb) return false;
  return true;
}

static Summary run_one(const string& spiral_name) {
  Summary s;
  s.name = spiral_name;
  if (!valid_spiral_string(spiral_name)) { s.kind = "FAIL_BUILD"; return s; }
  try {
    spiral_nomenclature sn(spiral_name);
    Triangulation T(sn);
    s.N = (int)T.N;
    auto D = DelaunayTriangulation::compute(T);

    AlexandrovSolver solver;
    solver.D = D;
    solver.record_diag = true;
    auto coords = solver.solve();

    // Post-process diag_trace.
    int max1 = 0, max2 = 0;
    for (auto& e : solver.diag_trace) {
      if (e.n_near_pi_001 > max1)  max1 = e.n_near_pi_001;
      if (e.n_near_pi_0001 > max2) max2 = e.n_near_pi_0001;
    }
    s.max_n_near_pi_001 = max1;
    s.max_n_near_pi_0001 = max2;
    if (!solver.diag_trace.empty()) {
      auto& last = solver.diag_trace.back();
      s.final_n_near_pi_001 = last.n_near_pi_001;
      s.final_n_near_pi_0001 = last.n_near_pi_0001;
      s.final_min_h_sq = last.min_h_sq;
    }

    // Classify outcome from the solver's ValidationStatus enum.  Note
    // coords are now ALWAYS populated when reconstruction succeeded —
    // solver.valid() (== ValidationStatus::OK) is the OK predicate, NOT
    // !coords.empty().
    s.kind = AlexandrovSolver::status_str(solver.stats_status);

    // Compute reconstructed polytope volume.  solve() now ALWAYS returns
    // positions when reconstruction succeeded (even if validation
    // failed) — so we can compute volume for failed cases too.
    // Normalize by mean-edge-length³ to remove overall scale dependence.
    auto& pos = coords;
    if (!pos.empty()) {
      double vol6 = 0;
      for (int f = 0; f < solver.D.nf; f++) {
        if (solver.D.f_he[f] < 0) continue;
        int ha = solver.D.f_he[f];
        int hb = solver.D.he_next[ha];
        int hc = solver.D.he_next[hb];
        auto& a = pos[solver.D.he_origin[ha]];
        auto& b = pos[solver.D.he_origin[hb]];
        auto& c = pos[solver.D.he_origin[hc]];
        vol6 += a.dot(b.cross(c));
      }
      s.volume = std::abs(vol6) / 6.0;
      // Mean edge length over alive non-bigon edges.
      double sum_l = 0; int n_e = 0;
      for (int h = 0; h < solver.D.nh; h += 2) {
        if (!solver.D.alive(h)) continue;
        sum_l += solver.D.he_length[h];
        n_e++;
      }
      double mean_l = (n_e > 0) ? sum_l / n_e : 1.0;
      s.volume_normalized = s.volume / (mean_l * mean_l * mean_l);
    }
  } catch (...) {
    s.kind = "FAIL_BUILD";
  }
  return s;
}

int main(int argc, char** argv) {
  string spiral_file, out_file;
  int start = 0, limit = -1;

  for (int i = 1; i < argc; i++) {
    string a = argv[i];
    if      (a == "--out"   && i+1 < argc) out_file = argv[++i];
    else if (a == "--start" && i+1 < argc) start    = atoi(argv[++i]);
    else if (a == "--limit" && i+1 < argc) limit    = atoi(argv[++i]);
    else if (spiral_file.empty()) spiral_file = a;
    else { fprintf(stderr, "unknown arg: %s\n", argv[i]); return 1; }
  }
  if (spiral_file.empty() || out_file.empty()) {
    fprintf(stderr,
      "Usage: %s <spiral_file> --out summary.tsv [--start <idx>] [--limit <n>]\n",
      argv[0]);
    return 1;
  }

  ifstream in(spiral_file);
  if (!in) { fprintf(stderr, "cannot open %s\n", spiral_file.c_str()); return 1; }
  vector<string> all;
  string line;
  while (getline(in, line)) if (!line.empty()) all.push_back(line);
  fprintf(stderr, "loaded %zu isomers\n", all.size());

  int total = (int)all.size();
  if (start < 0) start = 0;
  if (start >= total) { fprintf(stderr, "start out of range\n"); return 1; }
  int end = (limit > 0 && start + limit < total) ? start + limit : total;

  FILE* fp = fopen(out_file.c_str(), "w");
  if (!fp) { fprintf(stderr, "cannot open %s for writing\n", out_file.c_str()); return 1; }
  setvbuf(fp, nullptr, _IOLBF, 0);
  fprintf(fp,
    "#N\tmax_n_near_pi_001\tmax_n_near_pi_0001\tfinal_n_near_pi_001"
    "\tfinal_n_near_pi_0001\tfinal_min_h_sq\tvolume\tvolume_norm\tkind\tname\n");

  std::mutex out_mtx;
  std::atomic<int> done{0};

  #pragma omp parallel for schedule(dynamic, 64)
  for (int i = start; i < end; i++) {
    auto s = run_one(all[i]);
    {
      std::lock_guard<std::mutex> lk(out_mtx);
      fprintf(fp, "%d\t%d\t%d\t%d\t%d\t%.6e\t%.6e\t%.6e\t%s\t%s\n",
              s.N, s.max_n_near_pi_001, s.max_n_near_pi_0001,
              s.final_n_near_pi_001, s.final_n_near_pi_0001,
              s.final_min_h_sq, s.volume, s.volume_normalized,
              s.kind.c_str(), s.name.c_str());
    }
    int d = ++done;
    if (d % 5000 == 0)
      fprintf(stderr, "  %d / %d\n", d, end - start);
  }

  fclose(fp);
  return 0;
}
