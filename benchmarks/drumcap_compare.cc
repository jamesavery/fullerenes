// Compare the Alexandrov-solver fixed point ("drum-cap") with the
// Deltahedron::fromExtensionPathOptimized geometry ("true polytope")
// for a single fullerene isomer.  Useful for investigating Class-A
// drum-cap failures (see claude-projects/delaunay-geometry/failure-suite.md).
//
// Outputs three mol2 files in <out_prefix>.{bi,delt,cones}.mol2:
//   .bi.mol2     12 cone vertices and iDT edges from AlexandrovSolver
//   .delt.mol2   N/2+2 vertices and dual-graph edges from Deltahedron
//   .cones.mol2  12 degree-5 vertices extracted from .delt (no edges)
//
// Usage: drumcap_compare <spiral> <out_prefix>

#include "fullerenes/triangulation.hh"
#include "fullerenes/delaunay.hh"
#include "fullerenes/delaunay_alexandrov12.hh"
#include "fullerenes/spiral.hh"
#include "fullerenes/buckinverse.hh"
#include "fullerenes/deltahedron.hh"

#include <cstdio>
#include <fstream>
#include <string>
#include <vector>

using namespace std;
using namespace buckinverse;

static void write_mol2(const string& path, const string& name,
                        const vector<coord3d>& pos,
                        const vector<pair<int,int>>& edges) {
  ofstream out(path);
  if (!out) { fprintf(stderr, "cannot write %s\n", path.c_str()); exit(1); }
  out << "@<TRIPOS>MOLECULE\n" << name << "\n";
  out << pos.size() << " " << edges.size() << " 0 0 0\n";
  out << "SMALL\nNO_CHARGES\n\n";
  out << "@<TRIPOS>ATOM\n";
  for (int i = 0; i < (int)pos.size(); i++) {
    out << (i + 1) << " V" << i << " "
        << pos[i][0] << " " << pos[i][1] << " " << pos[i][2]
        << " C 1 RES 0.0\n";
  }
  out << "@<TRIPOS>BOND\n";
  for (int i = 0; i < (int)edges.size(); i++)
    out << (i + 1) << " " << (edges[i].first + 1) << " "
        << (edges[i].second + 1) << " 1\n";
}

static vector<pair<int,int>> idt_edges(const DelaunayTriangulation& T) {
  vector<pair<int,int>> es;
  for (int h = 0; h < T.nh; h += 2)
    if (T.alive(h))
      es.push_back({T.he_origin[h], T.dest(h)});
  return es;
}

// Edges of a triangulation / deltahedron via the indexed adjacency
// API (D.degree(v), D[v][i]).
template <class G>
static vector<pair<int,int>> dense_graph_edges(const G& T) {
  vector<pair<int,int>> es;
  for (int u = 0; u < T.N; u++) {
    int d = (int)T.degree(u);
    for (int i = 0; i < d; i++) {
      int v = T[u][i];
      if (u < v) es.push_back({u, v});
    }
  }
  return es;
}

int main(int argc, char** argv) {
  if (argc < 3) {
    fprintf(stderr, "Usage: %s <spiral> <out_prefix>\n", argv[0]);
    return 1;
  }
  string spiral = argv[1];
  string prefix = argv[2];

  // 1. Parse the spiral and build the dual triangulation.
  spiral_nomenclature sn(spiral);
  Triangulation T(sn);
  fprintf(stderr, "Built dual triangulation: %d vertices\n", T.N);

  // 2. Run the Alexandrov solver — produces the (drum-cap) fixed point.
  auto D = DelaunayTriangulation::compute(T);
  AlexandrovSolver solver;
  solver.D         = D;
  solver.palc_tsvd = true;
  auto bi_coords = solver.solve();
  fprintf(stderr, "AlexandrovSolver: status=%s  kappa=%.3e  vol_norm=%.3e\n",
          AlexandrovSolver::status_str(solver.stats_status),
          solver.stats_final_kappa, solver.stats_volume_norm);
  write_mol2(prefix + ".bi.mol2", spiral + ":bi",
             bi_coords, idt_edges(solver.D));
  fprintf(stderr, "wrote %s.bi.mol2 (12 cone vertices, %zu iDT edges)\n",
          prefix.c_str(), idt_edges(solver.D).size());

  // 3. Run the Deltahedron extension-path optimiser — produces the
  // full N/2+2 vertex 3D embedding.
  ReducibleDual rd(T);
  ExtensionPath ep = rd.reduceToExtensionPath(20);
  fprintf(stderr, "ExtensionPath: seed=%d steps=%zu\n",
          (int)ep.seed, ep.steps.size());
  Deltahedron Delt = Deltahedron::fromExtensionPathOptimized(ep);
  vector<coord3d> delt_coords(Delt.points.begin(), Delt.points.end());
  auto delt_es = dense_graph_edges(Delt);
  write_mol2(prefix + ".delt.mol2", spiral + ":delt",
             delt_coords, delt_es);
  fprintf(stderr, "wrote %s.delt.mol2 (%zu vertices, %zu edges)\n",
          prefix.c_str(), delt_coords.size(), delt_es.size());

  // 4. Extract degree-5 vertices (= pentagons = cone points) from
  // the Deltahedron — these are the vertices of the candidate true
  // Alexandrov polytope.
  vector<coord3d> cone_coords;
  for (int v = 0; v < Delt.N; v++)
    if (Delt.degree(v) == 5)
      cone_coords.push_back(delt_coords[v]);
  fprintf(stderr, "Extracted %zu degree-5 vertices (expected 12)\n",
          cone_coords.size());
  write_mol2(prefix + ".cones.mol2", spiral + ":cones",
             cone_coords, {});
  fprintf(stderr, "wrote %s.cones.mol2 (%zu cone vertices)\n",
          prefix.c_str(), cone_coords.size());

  // 5. Quantitative comparison: stats on the two 12-vertex clouds.
  auto stats = [&](const char* label, const vector<coord3d>& pos) {
    if (pos.empty()) return;
    coord3d cen(0, 0, 0);
    for (const auto& p : pos) cen = cen + p;
    cen = cen * (1.0 / pos.size());
    double rsum2 = 0;
    coord3d ext_min(1e30, 1e30, 1e30), ext_max(-1e30, -1e30, -1e30);
    for (const auto& p : pos) {
      coord3d d = p - cen;
      rsum2 += d.dot(d);
      for (int k = 0; k < 3; k++) {
        if (p[k] < ext_min[k]) ext_min[k] = p[k];
        if (p[k] > ext_max[k]) ext_max[k] = p[k];
      }
    }
    double r_rms = std::sqrt(rsum2 / pos.size());
    coord3d span = ext_max - ext_min;
    fprintf(stderr,
            "  %-6s  centroid=(%+.3f, %+.3f, %+.3f)  RMS_radius=%.3f  "
            "bbox=(%.2f, %.2f, %.2f)  flat?=%s\n",
            label, cen[0], cen[1], cen[2], r_rms,
            span[0], span[1], span[2],
            (std::min({span[0], span[1], span[2]}) < 0.05 * r_rms)
              ? "YES" : "no");
  };
  fprintf(stderr, "\n=== Geometry comparison ===\n");
  stats("BI",    bi_coords);
  stats("cones", cone_coords);

  return 0;
}
