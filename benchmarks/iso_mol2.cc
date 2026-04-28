// Run the Alexandrov solver on a fullerene specified by (N, isomer_idx)
// produced by buckygen, and dump the resulting 12-cone-point polytope
// as .mol2.  Convenience wrapper around spiral_mol2's logic — handles
// the buckygen → canonical-spiral-name lookup.
//
// Usage:
//   iso_mol2 <N> <idx> [--ipr] [--out file.mol2]
//
// Notes:
//   - <idx> is 0-based index into buckygen's enumeration of size-N
//     fullerenes.
//   - --ipr restricts to the IPR (isolated pentagon rule) sub-list,
//     in which case <idx> is into that smaller list.

#include "fullerenes/triangulation.hh"
#include "fullerenes/delaunay.hh"
#include "fullerenes/delaunay_alexandrov12.hh"
#include "fullerenes/spiral.hh"
#include "fullerenes/buckygen-wrapper.hh"

#include <cstdio>
#include <cstring>
#include <string>
#include <vector>
#include <set>
#include <utility>

using namespace std;

int main(int argc, char** argv) {
  if (argc < 3) {
    fprintf(stderr,
      "Usage: %s <N> <idx> [--ipr] [--out file.mol2]\n", argv[0]);
    return 1;
  }
  int N = atoi(argv[1]);
  int idx = atoi(argv[2]);
  bool ipr = false;
  string out_path;
  for (int i = 3; i < argc; i++) {
    string a = argv[i];
    if      (a == "--ipr") ipr = true;
    else if (a == "--out" && i+1 < argc) out_path = argv[++i];
  }

  // Enumerate via buckygen up to and including idx, take the dual T.
  auto Q = BuckyGen::start(N, ipr);
  Graph G;
  Triangulation T;
  for (int i = 0; i <= idx; i++) {
    if (!BuckyGen::next_fullerene(Q, G)) {
      fprintf(stderr, "buckygen exhausted before idx=%d (got %d isomers)\n",
              idx, i);
      BuckyGen::stop(Q);
      return 1;
    }
    if (i == idx) T = Triangulation(G);  // buckygen returns the dual triangulation
  }
  BuckyGen::stop(Q);

  spiral_nomenclature sn(T, spiral_nomenclature::FULLERENE,
                            spiral_nomenclature::TRIANGULATION,
                            /*rarest_special_start=*/true);
  string name = "C" + to_string(N) + "-" + sn.to_string();

  auto D = DelaunayTriangulation::compute(T);
  AlexandrovSolver solver;
  solver.D = D;
  auto pos = solver.solve();
  if (pos.empty()) {
    fprintf(stderr, "reconstruction failed for %s (status=%s)\n",
            name.c_str(),
            AlexandrovSolver::status_str(solver.stats_status));
    return 1;
  }

  // Compute mean edge length and volume (same convention as spiral_mol2).
  double sum_l = 0; int n_e = 0;
  for (int h = 0; h < solver.D.nh; h += 2) {
    if (!solver.D.alive(h)) continue;
    sum_l += solver.D.he_length[h]; n_e++;
  }
  double mean_l = (n_e > 0) ? sum_l / n_e : 1.0;
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
    if (u == v) continue;
    edges.insert({min(u,v), max(u,v)});
  }

  FILE* fp = out_path.empty() ? stdout : fopen(out_path.c_str(), "w");
  if (!fp) { fprintf(stderr, "cannot open %s\n", out_path.c_str()); return 1; }

  fprintf(fp, "# C%d isomer idx=%d %s\n", N, idx, ipr ? "(IPR)" : "");
  fprintf(fp, "# spiral: %s\n", name.c_str());
  fprintf(fp, "# nv=%d  vol=%.6e  mean_edge_length=%.6e  vol_norm=%.6e\n",
          solver.D.nv, vol, mean_l, vol_norm);
  fprintf(fp, "# kappa=%.3e  flips=%d  steps=%d  status=%s\n",
          solver.stats_final_kappa, solver.stats_flips, solver.stats_steps,
          AlexandrovSolver::status_str(solver.stats_status));
  fprintf(fp, "# tbar_n_cells=%d  tbar_simple=%d  convex=%d  no_self_intersect=%d\n",
          solver.stats_tbar_n_cells, solver.stats_tbar_simple_polygonal,
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
