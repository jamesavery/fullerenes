// Probe a single failing isomer: dump multi-edge θ_e values and Cayley-Menger
// diagnostics at the converged (T(0), r(0)).
//
// Usage: probe_multi <spiral_name>
// e.g.:  probe_multi 'C134-[GS:1,2,12,17,36,43,50,56,57,62,67,69]-fullerene'

#include "fullerenes/triangulation.hh"
#include "fullerenes/delaunay.hh"
#include "fullerenes/delaunay_alexandrov12.hh"
#include "fullerenes/spiral.hh"

#include <cstdio>
#include <cmath>
#include <vector>
#include <map>
#include <set>

using namespace std;

int main(int argc, char** argv) {
  if (argc < 2) { fprintf(stderr, "Usage: %s <spiral_name>\n", argv[0]); return 1; }
  string name = argv[1];
  spiral_nomenclature sn(name);
  Triangulation T(sn);
  auto D = DelaunayTriangulation::compute(T);

  AlexandrovSolver solver;
  solver.D = D;
  solver.verbose = true;
  auto coords = solver.solve();

  printf("\n=== Final state ===\n");
  printf("dual N = %d, |κ| = %.3e\n", solver.D.nv, solver.stats_final_kappa);
  printf("flips during solve: %d\n", solver.stats_flips);
  printf("T(0) simplicial: %s\n", solver.stats_t0_simplicial ? "yes" : "no");
  printf("T̄(0) simple polygonal: %s (n_cells=%d)\n",
         solver.stats_tbar_simple_polygonal ? "yes" : "no",
         solver.stats_tbar_n_cells);
  printf("returned %zu coords\n", coords.size());

  // Identify multi-edges and self-loops.
  printf("\n=== Multi-edges and self-loops in T(0) at κ=0 ===\n");
  map<pair<int,int>, vector<int>> by_pair;
  vector<int> self_loops;
  for (int h = 0; h < solver.D.nh; h += 2) {
    if (!solver.D.alive(h)) continue;
    int u = solver.D.he_origin[h], v = solver.D.dest(h);
    if (u == v) self_loops.push_back(h);
    else by_pair[{min(u,v), max(u,v)}].push_back(h);
  }

  printf("self-loops: %zu\n", self_loops.size());
  for (int h : self_loops) {
    printf("  h=%d v=%d ℓ=%.6f\n", h, solver.D.he_origin[h], solver.D.he_length[h]);
  }

  // Inessential mask
  auto bi_mask = AlexandrovSolver::inessential_edges(solver.D, solver.r);

  printf("\nmulti-edge pairs (θ at r_before_extrap vs final r):\n");
  bool have_pre = !solver.r_before_extrap.empty();
  for (auto& [key, hs] : by_pair) {
    if (hs.size() <= 1) continue;
    printf("  (%d,%d): %zu copies\n", key.first, key.second, hs.size());
    for (int h : hs) {
      bool is_bigon = (solver.D.he_face[h] == solver.D.he_face[h ^ 1]);
      double theta_pre   = have_pre && !is_bigon
        ? AlexandrovSolver::theta(solver.D, solver.r_before_extrap, h)
        : NAN;
      double theta_final = is_bigon ? NAN : AlexandrovSolver::theta(solver.D, solver.r, h);
      printf("    h=%d  ℓ=%.6f  θ_pre=%.6f  θ_final=%.6f%s\n",
             h, solver.D.he_length[h],
             theta_pre, theta_final,
             is_bigon ? " [bigon]" : "");
    }
  }
  printf("\n  stats_extrap_kappa = %.3e  (max|κ| just after extrapolation)\n",
         solver.stats_extrap_kappa);

  // Tesselation comparison
  vector<int> labels(solver.D.nv);
  for (int v = 0; v < solver.D.nv; v++) labels[v] = v;
  auto coc_tess = solver.D.canonical_tesselation(labels);
  auto bi_tess  = AlexandrovSolver::polytope_tesselation(solver.D, solver.r, labels);
  int n_coc = 0, n_bi = 0;
  auto coc_mask = solver.D.cocircular_edges();
  for (int h = 0; h < solver.D.nh; h += 2) {
    if (!solver.D.alive(h)) continue;
    if (coc_mask[h]) n_coc++;
    if (bi_mask[h])  n_bi++;
  }
  printf("\n=== Tesselations ===\n");
  printf("canonical (cocircular): %d cells, %d cocircular edges\n",
         coc_tess.n_cells(), n_coc);
  printf("B-I (inessential):      %d cells, %d inessential edges\n",
         bi_tess.n_cells(), n_bi);

  // Print T̄ cells (B-I) with vertex labels and edge length²
  printf("\n=== T̄(0) cells (B-I tesselation) ===\n");
  int idx = 0;
  for (auto& cell : bi_tess.cells) {
    printf("  cell %d (n=%zu):", idx, cell.size());
    for (auto& [v, L] : cell) printf(" (%d, ℓ²=%lld)", v, L);
    printf("\n");
    idx++;
  }
  return 0;
}
