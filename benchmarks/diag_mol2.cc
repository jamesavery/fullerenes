// Produce a .mol2 file from the (possibly failing) Alexandrov solve.
// Tries each face of T* as a seed for BFS reconstruction, emits whichever
// works.  If none work, emits a placeholder mol2 with the cone-point positions
// the BFS could recover before failing, plus all r* values printed to stderr.
//
// Usage: diag_mol2 <graph_file> <marker> <out.mol2>

#include "fullerenes/triangulation.hh"
#include "fullerenes/delaunay.hh"
#include "fullerenes/delaunay_alexandrov.hh"
#include "fullerenes/geometry.hh"

#include <cstdio>
#include <cmath>
#include <vector>
#include <string>
#include <fstream>
#include <sstream>

using namespace std;

static Graph parse_graph_dump(const string& path, const string& marker) {
  ifstream in(path);
  if (!in) { fprintf(stderr, "cannot open %s\n", path.c_str()); exit(1); }
  string line;
  vector<vector<node_t>> rows;
  int N = 0;
  bool in_graph = false, found = false;
  while (getline(in, line)) {
    if (!in_graph) {
      if (line.rfind("GRAPH ", 0) == 0) {
        if (marker.empty() || line.find(marker) != string::npos) {
          auto pos = line.find("N=");
          if (pos == string::npos) continue;
          N = atoi(line.c_str() + pos + 2);
          rows.assign(N, {});
          in_graph = true; found = true;
        }
      }
    } else {
      if (line.rfind("END_GRAPH", 0) == 0) break;
      auto colon = line.find(':');
      if (colon == string::npos) continue;
      int v = atoi(line.c_str());
      stringstream ss(line.substr(colon + 1));
      int u;
      while (ss >> u) rows[v].push_back(u);
    }
  }
  if (!found) { fprintf(stderr, "no graph found matching '%s'\n", marker.c_str()); exit(1); }
  Graph G(N, GRAPH_DMAX);
  for (int v = 0; v < N; v++) G.assign_row(v, rows[v]);
  return G;
}

// Robust BFS reconstruction: try each face as seed; succeed if all 12 vertices placed.
// If seed face fails, try next face.  Returns coords (possibly with NaNs for unplaced).
static vector<coord3d> robust_reconstruct(const DelaunayTriangulation& T, const vector<double>& r) {
  int n = T.nv;
  for (int seed_attempt = 0; seed_attempt < T.nf; seed_attempt++) {
    if (T.f_he[seed_attempt] < 0) continue;
    vector<coord3d> pos(n, coord3d(0,0,0));
    vector<bool> placed(n, false), face_done(T.nf, false);

    int f0 = seed_attempt;
    int h0 = T.f_he[f0], h1 = T.he_next[h0], h2 = T.he_next[h1];
    int i = T.he_origin[h0], j = T.he_origin[h1], k = T.he_origin[h2];
    if (i == j || j == k || i == k) continue;  // degenerate seed face

    auto gram = [](double r_u, double r_v, double L_uv) {
      return (r_u*r_u + r_v*r_v - L_uv*L_uv) / 2;
    };

    pos[i] = coord3d(r[i], 0, 0);
    double g_ij = gram(r[i], r[j], T.he_length[h0]);
    double jx   = g_ij / r[i];
    double jy2  = r[j]*r[j] - jx*jx;
    if (jy2 < -1e-8 * r[j]*r[j]) continue;
    pos[j] = coord3d(jx, sqrt(max(0.0, jy2)), 0);

    double g_ik = gram(r[i], r[k], T.he_length[h2]);
    double g_jk = gram(r[j], r[k], T.he_length[h1]);
    double kx   = g_ik / r[i];
    double ky   = (pos[j][1] > 1e-15) ? (g_jk - kx*pos[j][0]) / pos[j][1] : 0;
    double kz2  = r[k]*r[k] - kx*kx - ky*ky;
    if (kz2 < -1e-8 * r[k]*r[k]) continue;
    pos[k] = coord3d(kx, ky, sqrt(max(0.0, kz2)));
    placed[i] = placed[j] = placed[k] = true;
    face_done[f0] = true;

    vector<int> queue = {f0};
    int head = 0;
    while (head < (int)queue.size()) {
      int f = queue[head++];
      int hf = T.f_he[f];
      for (int s = 0; s < 3; s++, hf = T.he_next[hf]) {
        int ht = hf ^ 1, fa = T.he_face[ht];
        if (fa < 0 || face_done[fa]) continue;
        int u = T.he_origin[hf], v = T.dest(hf), w = T.he_origin[T.prev(ht)];
        if (u == v || u == w || v == w) {
          // degenerate face; skip placement but propagate
          face_done[fa] = true;
          queue.push_back(fa);
          continue;
        }
        if (!placed[w]) {
          double g_wu = gram(r[w], r[u], T.he_length[T.he_next[ht]]);
          double g_wv = gram(r[w], r[v], T.he_length[T.prev(ht)]);
          coord3d pu = pos[u], pv = pos[v];
          double uu = pu.dot(pu), uv = pu.dot(pv), vv = pv.dot(pv);
          double det = uu*vv - uv*uv;
          if (fabs(det) < 1e-20) { face_done[fa] = true; queue.push_back(fa); continue; }
          double a = (g_wu*vv - g_wv*uv) / det;
          double b = (g_wv*uu - g_wu*uv) / det;
          coord3d proj = pu*a + pv*b;
          coord3d nrm  = pu.cross(pv);
          double nl = nrm.norm();
          double gamma_sq = r[w]*r[w] - proj.dot(proj);
          if (gamma_sq < -1e-8 * r[w]*r[w] || nl < 1e-15) {
            face_done[fa] = true; queue.push_back(fa); continue;
          }
          double gamma = sqrt(max(0.0, gamma_sq));
          double scale = gamma / nl;
          // Side: against the face's previously-placed third vertex
          int hf2 = T.f_he[f];
          int old_w = -1;
          for (int s2 = 0; s2 < 3; s2++, hf2 = T.he_next[hf2])
            if (T.he_origin[hf2] != u && T.he_origin[hf2] != v) { old_w = T.he_origin[hf2]; break; }
          coord3d p_old = (old_w >= 0) ? pos[old_w] : coord3d(0,0,0);
          double side = (p_old - proj).dot(nrm);
          pos[w] = (side > 0) ? proj - nrm*scale : proj + nrm*scale;
          placed[w] = true;
        }
        face_done[fa] = true;
        queue.push_back(fa);
      }
    }

    int n_placed = 0; for (bool b : placed) if (b) n_placed++;
    if (n_placed == n) {
      fprintf(stderr, "Reconstruction OK from seed face %d (%d vertices placed)\n", f0, n_placed);
      return pos;
    }
    // partial — keep trying other seed faces
  }
  fprintf(stderr, "All seed faces failed; returning empty.\n");
  return {};
}

static void write_mol2(const string& path, const vector<coord3d>& pos,
                        const DelaunayTriangulation& T, const string& molname) {
  ofstream out(path);
  if (!out) { fprintf(stderr, "cannot write %s\n", path.c_str()); exit(1); }
  // Build edge list (one per undirected edge).
  vector<pair<int,int>> edges;
  for (int h = 0; h < T.nh; h += 2) {
    if (!T.alive(h)) continue;
    int u = T.he_origin[h], v = T.dest(h);
    edges.push_back({u, v});
  }
  out << "@<TRIPOS>MOLECULE\n" << molname << "\n";
  out << pos.size() << " " << edges.size() << " 0 0 0\n";
  out << "SMALL\nNO_CHARGES\n\n";
  out << "@<TRIPOS>ATOM\n";
  for (int i = 0; i < (int)pos.size(); i++) {
    out << (i+1) << " V" << i << " " << pos[i][0] << " " << pos[i][1] << " " << pos[i][2]
        << " C 1 RES 0.0\n";
  }
  out << "@<TRIPOS>BOND\n";
  for (int i = 0; i < (int)edges.size(); i++)
    out << (i+1) << " " << (edges[i].first+1) << " " << (edges[i].second+1) << " 1\n";
}

int main(int argc, char** argv) {
  if (argc < 4) { fprintf(stderr, "Usage: %s graph_file marker out.mol2\n", argv[0]); return 1; }
  string path = argv[1], marker = argv[2], outpath = argv[3];

  Graph G = parse_graph_dump(path, marker);
  Triangulation T_dual(G);
  auto D = DelaunayTriangulation::compute(T_dual);

  AlexandrovSolver solver;
  solver.D = D;
  solver.verbose = true;
  auto coords = solver.solve();

  fprintf(stderr, "After solve: |κ|=%.3e, %zu coords from solve\n",
          solver.stats_final_kappa, coords.size());
  fprintf(stderr, "Radii r*:\n");
  for (size_t i = 0; i < solver.r.size(); i++) fprintf(stderr, "  r[%zu] = %.6f\n", i, solver.r[i]);

  // If solver returned empty, try robust reconstruct.
  if (coords.empty()) {
    fprintf(stderr, "Solver returned empty; attempting robust reconstruction.\n");
    coords = robust_reconstruct(solver.D, solver.r);
  }

  if (coords.empty()) {
    fprintf(stderr, "No reconstruction possible.  Emitting placeholder mol2 with vertices on a circle.\n");
    coords.assign(solver.D.nv, coord3d(0,0,0));
    for (int i = 0; i < solver.D.nv; i++) {
      double t = 2*M_PI*i/solver.D.nv;
      coords[i] = coord3d(cos(t), sin(t), 0) * solver.r[i];
    }
  }

  write_mol2(outpath, coords, solver.D, marker);
  fprintf(stderr, "Wrote %s\n", outpath.c_str());
  return 0;
}
