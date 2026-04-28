// Run the Alexandrov solver on a single spiral name and dump the resulting
// 12-cone-point polytope as .mol2 (iDT edges as bonds).  solve() now
// always returns reconstructed positions, even when validation fails —
// so degenerate cases are still inspectable here.  Header records the
// solver's ValidationStatus.
//
// Usage:
//   spiral_mol2 <spiral_name> [--out file.mol2]
//
// Atom labels: vertex index 0..11 written as "C0" .. "C11" with element C.

#include "fullerenes/triangulation.hh"
#include "fullerenes/delaunay.hh"
#include "fullerenes/delaunay_alexandrov12.hh"
#include "fullerenes/spiral.hh"

#include <cstdio>
#include <string>
#include <vector>
#include <set>
#include <utility>

using namespace std;

int main(int argc, char** argv) {
  if (argc < 2) {
    fprintf(stderr, "Usage: %s <spiral_name> [--out file.mol2]\n", argv[0]);
    return 1;
  }
  string name = argv[1];
  string out_path;
  for (int i = 2; i < argc; i++) {
    if (string(argv[i]) == "--out" && i+1 < argc) out_path = argv[++i];
  }

  spiral_nomenclature sn(name);
  Triangulation T(sn);
  auto D = DelaunayTriangulation::compute(T);

  AlexandrovSolver solver;
  solver.D = D;
  auto pos = solver.solve();   // always returns positions if reconstruct succeeded
  if (pos.empty()) {
    fprintf(stderr, "reconstruction failed for %s (status=%s)\n",
            name.c_str(),
            AlexandrovSolver::status_str(solver.stats_status));
    return 1;
  }

  // Compute mean edge length for vol_norm reporting.
  double sum_l = 0; int n_e = 0;
  for (int h = 0; h < solver.D.nh; h += 2) {
    if (!solver.D.alive(h)) continue;
    sum_l += solver.D.he_length[h]; n_e++;
  }
  double mean_l = (n_e > 0) ? sum_l / n_e : 1.0;

  // Volume.
  double vol6 = 0;
  for (int f = 0; f < solver.D.nf; f++) {
    if (solver.D.f_he[f] < 0) continue;
    int ha = solver.D.f_he[f];
    int hb = solver.D.he_next[ha];
    int hc = solver.D.he_next[hb];
    vol6 += pos[solver.D.he_origin[ha]].dot(
              pos[solver.D.he_origin[hb]].cross(pos[solver.D.he_origin[hc]]));
  }
  double vol = std::abs(vol6) / 6.0;
  double vol_norm = vol / (mean_l * mean_l * mean_l);

  // Collect unique undirected edges from alive half-edges.
  set<pair<int,int>> edges;
  for (int h = 0; h < solver.D.nh; h += 2) {
    if (!solver.D.alive(h)) continue;
    int u = solver.D.he_origin[h], v = solver.D.dest(h);
    if (u == v) continue;             // skip self-loops in mol2 (no bond)
    edges.insert({min(u,v), max(u,v)});
  }

  FILE* fp = out_path.empty() ? stdout : fopen(out_path.c_str(), "w");
  if (!fp) { fprintf(stderr, "cannot open %s\n", out_path.c_str()); return 1; }

  fprintf(fp, "# spiral: %s\n", name.c_str());
  fprintf(fp, "# nv=%d nh=%d  vol=%.6e  mean_edge_length=%.6e  vol_norm=%.6e\n",
          solver.D.nv, solver.D.nh, vol, mean_l, vol_norm);
  fprintf(fp, "# kappa=%.3e  flips=%d  steps=%d  t0_simpl=%d  tbar_n_cells=%d  tbar_simple=%d  vol_norm_solver=%.6e\n",
          solver.stats_final_kappa, solver.stats_flips, solver.stats_steps,
          solver.stats_t0_simplicial, solver.stats_tbar_n_cells,
          solver.stats_tbar_simple_polygonal, solver.stats_volume_norm);
  fprintf(fp, "# status=%s  convex=%d  no_self_intersect=%d\n",
          AlexandrovSolver::status_str(solver.stats_status),
          solver.stats_polytope_convex, solver.stats_polytope_no_self_intersect);
  fprintf(fp, "@<TRIPOS>MOLECULE\n%s\n%d %zu 0 0 0\nSMALL\nNO_CHARGES\n",
          name.c_str(), solver.D.nv, edges.size());
  fprintf(fp, "@<TRIPOS>ATOM\n");
  for (int v = 0; v < solver.D.nv; v++) {
    fprintf(fp, "%d C%d %.6f %.6f %.6f C\n",
            v + 1, v, pos[v][0], pos[v][1], pos[v][2]);
  }
  fprintf(fp, "@<TRIPOS>BOND\n");
  int bi = 1;
  for (auto& e : edges) {
    fprintf(fp, "%d %d %d 1\n", bi++, e.first + 1, e.second + 1);
  }

  if (!out_path.empty()) fclose(fp);
  return 0;
}
