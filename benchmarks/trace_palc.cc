// Per-step PALC + Newton trajectory dump for one isomer (Task #28
// characterization run).  Solves with record_diag=true and emits the
// resulting DiagEntry sequence as TSV — one row per step.
//
// Usage:
//   trace_palc <spiral_name> [--out trace.tsv]
//
// TSV columns:
//   phase  step  t  ds  nit  kappa_max
//   theta_min_dist_to_pi  n_near_pi_01  n_near_pi_001  n_near_pi_0001
//   min_h_sq  r_cv  det_J_sign  n_flips_cum  n_non_bigon_alive

#include "fullerenes/triangulation.hh"
#include "fullerenes/delaunay.hh"
#include "fullerenes/delaunay_alexandrov.hh"
#include "fullerenes/spiral.hh"

#include <cstdio>
#include <cstring>
#include <string>

using namespace std;

int main(int argc, char** argv) {
  if (argc < 2) {
    fprintf(stderr, "Usage: %s <spiral_name> [--out file.tsv]\n", argv[0]);
    return 1;
  }
  string name = argv[1];
  string out_path;
  for (int i = 2; i < argc; i++) {
    if (string(argv[i]) == "--out" && i + 1 < argc) out_path = argv[++i];
  }

  spiral_nomenclature sn(name);
  Triangulation T(sn);
  auto D = DelaunayTriangulation::compute(T);

  AlexandrovSolver solver;
  solver.D = D;
  solver.record_diag = true;
  solver.solve();

  FILE* fp = out_path.empty() ? stdout : fopen(out_path.c_str(), "w");
  if (!fp) { fprintf(stderr, "cannot open %s\n", out_path.c_str()); return 1; }

  fprintf(fp, "phase\tstep\tt\tds\tnit\tkappa_max"
              "\ttheta_min_dist_to_pi\tn_near_pi_01\tn_near_pi_001\tn_near_pi_0001"
              "\tmin_h_sq\tr_cv\tdet_J_sign\tn_flips_cum\tn_non_bigon_alive\n");
  for (auto& e : solver.diag_trace) {
    fprintf(fp,
            "%c\t%d\t%.10g\t%.10g\t%d\t%.10g"
            "\t%.10g\t%d\t%d\t%d"
            "\t%.10g\t%.10g\t%d\t%d\t%d\n",
            e.phase, e.step, e.t, e.ds, e.nit, e.kappa_max,
            e.theta_min_dist_to_pi, e.n_near_pi_01, e.n_near_pi_001, e.n_near_pi_0001,
            e.min_h_sq, e.r_cv, e.det_J_sign, e.n_flips_cum, e.n_non_bigon_alive);
  }

  // Trailing summary line as a comment (final solver outcome).
  fprintf(fp, "# final_kappa=%.6e flips=%d steps=%d t0_simpl=%d "
              "tbar_n_cells=%d tbar_simple=%d\n",
          solver.stats_final_kappa, solver.stats_flips, solver.stats_steps,
          solver.stats_t0_simplicial, solver.stats_tbar_n_cells,
          solver.stats_tbar_simple_polygonal);

  if (!out_path.empty()) fclose(fp);
  return 0;
}
