// Scan a buckygen size for iDTs that contain self-loops or multi-edges.
// Reports the first few of each type, and whether the solver succeeds on them.
//
// Usage:  scan_topo N [max_isomers_to_scan]

#include "fullerenes/triangulation.hh"
#include "fullerenes/delaunay.hh"
#include "fullerenes/delaunay_alexandrov.hh"
#include "fullerenes/buckygen-wrapper.hh"

#include <cstdio>
#include <cmath>
#include <vector>
#include <map>

using namespace std;

struct EdgeStats {
  int n_self_loops = 0;
  int n_multi_pairs = 0;
  bool has_unequal_multi = false;  // any multi-edge with ℓ₁ ≠ ℓ₂
};

static EdgeStats inspect(const DelaunayTriangulation& T) {
  EdgeStats s;
  map<pair<int,int>, vector<double>> by_pair;
  for (int h = 0; h < T.nh; h += 2) {
    if (!T.alive(h)) continue;
    int u = T.he_origin[h], v = T.dest(h);
    if (u == v) s.n_self_loops++;
    else        by_pair[{min(u,v), max(u,v)}].push_back(T.he_length[h]);
  }
  for (auto& [k, ls] : by_pair) {
    if (ls.size() > 1) {
      s.n_multi_pairs++;
      double mn = *min_element(ls.begin(), ls.end());
      double mx = *max_element(ls.begin(), ls.end());
      if ((mx - mn) / mx > 1e-6) s.has_unequal_multi = true;
    }
  }
  return s;
}

int main(int argc, char** argv) {
  if (argc < 2) { fprintf(stderr, "Usage: %s N [max]\n", argv[0]); return 1; }
  int N = atoi(argv[1]);
  int max_scan = (argc >= 3) ? atoi(argv[2]) : 100000;

  auto Q = BuckyGen::start(N, false, false);
  Triangulation T;
  int idx = 0;
  int n_self = 0, n_multi_eq = 0, n_multi_ne = 0;
  int report_self = 5, report_multi_eq = 5, report_multi_ne = 5;

  while (idx < max_scan && BuckyGen::next_fullerene(Q, T)) {
    if (idx % 100 == 0) fprintf(stderr, "scanning #%d (so far self=%d, multi-eq=%d, multi-ne=%d)\n", idx, n_self, n_multi_eq, n_multi_ne);
    DelaunayTriangulation D;
    try {
      D = DelaunayTriangulation::compute(T);
    } catch (...) { idx++; continue; }
    auto s = inspect(D);
    bool note = false;
    if (s.n_self_loops > 0 && report_self > 0) {
      printf("C%d #%d: SELF-LOOP iDT (n_self=%d, n_multi=%d)\n", N, idx, s.n_self_loops, s.n_multi_pairs);
      report_self--; note = true;
    }
    if (s.n_multi_pairs > 0 && s.has_unequal_multi && report_multi_ne > 0) {
      printf("C%d #%d: MULTI-EDGE-UNEQUAL iDT (n_multi=%d, n_self=%d)\n", N, idx, s.n_multi_pairs, s.n_self_loops);
      report_multi_ne--; note = true;
    }
    if (s.n_multi_pairs > 0 && !s.has_unequal_multi && s.n_self_loops == 0 && report_multi_eq > 0) {
      printf("C%d #%d: MULTI-EDGE-EQUAL iDT (n_multi=%d)\n", N, idx, s.n_multi_pairs);
      report_multi_eq--; note = true;
    }
    if (s.n_self_loops > 0) n_self++;
    if (s.n_multi_pairs > 0 && s.has_unequal_multi) n_multi_ne++;
    else if (s.n_multi_pairs > 0) n_multi_eq++;

    // Also: log the FIRST self-loop encountered (always) so background scans don't miss any.
    if (s.n_self_loops > 0 && n_self == 1) {
      printf(">>> FIRST SELF-LOOP iDT at C%d #%d\n", N, idx);
    }

    // Run solver on these noteworthy cases
    if (note) {
      AlexandrovSolver solver;
      solver.D = D;
      auto coords = solver.solve();
      auto sf = inspect(solver.D);
      printf("    solver: %s  final|κ|=%.2e  flips=%d  T(0) self=%d multi=%d\n",
             coords.empty() ? "FAIL(reconstruct)" : "OK",
             solver.stats_final_kappa, solver.stats_flips,
             sf.n_self_loops, sf.n_multi_pairs);
    }
    idx++;
  }
  BuckyGen::stop(Q);
  printf("\nC%d scanned %d isomers: self-loop=%d, multi-eq=%d, multi-ne=%d\n",
         N, idx, n_self, n_multi_eq, n_multi_ne);
  return 0;
}
