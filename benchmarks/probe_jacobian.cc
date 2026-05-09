// Probe the Jacobian eigenvalue spectrum at the converged-or-stalled
// state of a single isomer.  Reports the spectrum of J = ∂κ/∂r, plus
// per-vertex r values, edge-length-derived stats, and key
// AlexandrovSolver outputs.
//
// Usage: probe_jacobian <spiral_name> [--gauge-snap]

#include "fullerenes/triangulation.hh"
#include "fullerenes/delaunay.hh"
#include "fullerenes/delaunay_alexandrov12.hh"
#include "fullerenes/spiral.hh"

#include <cstdio>
#include <string>

using namespace std;

int main(int argc, char** argv) {
  if (argc < 2) { fprintf(stderr, "Usage: %s <spiral> [--gauge-snap] [--gauge-project] [--tsvd] [--rcond <r>]\n", argv[0]); return 1; }
  string name = argv[1];
  bool gauge_snap = false, gauge_project = false, tsvd = false;
  double rcond = 5e-3;
  for (int i = 2; i < argc; i++) {
    string a = argv[i];
    if      (a == "--gauge-snap")    gauge_snap    = true;
    else if (a == "--gauge-project") gauge_project = true;
    else if (a == "--tsvd")          tsvd          = true;
    else if (a == "--rcond" && i+1 < argc) rcond = std::stod(argv[++i]);
  }
  spiral_nomenclature sn(name);
  Triangulation T(sn);
  auto D = DelaunayTriangulation::compute(T);

  AlexandrovSolver solver;
  solver.D = D;
  solver.palc_gauge_snap    = gauge_snap;
  solver.palc_gauge_project = gauge_project;
  solver.palc_tsvd          = tsvd;
  solver.palc_tsvd_rcond    = rcond;
  solver.solve();
  printf("# gauge_snap=%s  gauge_project=%s  tsvd=%s  rcond=%.1e\n",
         gauge_snap ? "true" : "false",
         gauge_project ? "true" : "false",
         tsvd ? "true" : "false", rcond);

  printf("# %s\n", name.c_str());
  printf("# status=%s  kappa=%.6e  vol_norm=%.6e  flips=%d  steps=%d\n",
         AlexandrovSolver::status_str(solver.stats_status),
         solver.stats_final_kappa, solver.stats_volume_norm,
         solver.stats_flips, solver.stats_steps);
  printf("# tbar_n_cells=%d  t0_simplicial=%d  tbar_simple=%d\n",
         solver.stats_tbar_n_cells, solver.stats_t0_simplicial,
         solver.stats_tbar_simple_polygonal);

  // Final state
  printf("\n# r values (per vertex):\n");
  for (int v = 0; v < (int)solver.r.size(); v++)
    printf("  r[%2d] = %.6f\n", v, solver.r[v]);

  // kappa values
  auto kappa = AlexandrovSolver::kappa(solver.D, solver.r);
  printf("\n# kappa values (per vertex):\n");
  for (int v = 0; v < (int)kappa.size(); v++)
    printf("  kappa[%2d] = %+.4e\n", v, kappa[v]);

  // Jacobian spectrum
  auto eigs = AlexandrovSolver::jacobian_eigvals(solver.D, solver.r);
  printf("\n# Jacobian eigenvalues (J = ∂κ/∂r), ascending:\n");
  for (int i = 0; i < (int)eigs.size(); i++)
    printf("  λ[%2d] = %+.6e\n", i, eigs[i]);

  // Spread
  if (!eigs.empty()) {
    double min_abs = 1e30, max_abs = 0;
    int n_neg = 0;
    for (auto e : eigs) {
      double a = std::abs(e);
      if (a < min_abs) min_abs = a;
      if (a > max_abs) max_abs = a;
      if (e < 0) n_neg++;
    }
    printf("\n# spectrum summary: smallest |λ| = %.4e, largest |λ| = %.4e, "
           "cond = %.4e, neg eigvals = %d/%d (Lorentzian signature should be (1, n-1))\n",
           min_abs, max_abs, max_abs / min_abs, n_neg, (int)eigs.size());
  }

  // Sign(det J)
  int sgn = AlexandrovSolver::jacobian_det_sign(solver.D, solver.r);
  printf("# sign(det J) = %d\n", sgn);

  return 0;
}
