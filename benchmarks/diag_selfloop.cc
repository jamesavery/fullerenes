// Self-loop / multi-edge diagnostic for the Alexandrov solver.
//
// Loads a graph in the GRAPH FAIL/GFAIL adjacency-list format produced by
// bench_alexandrov, computes the iDT, reports self-loop and multi-edge counts
// at the input iDT, runs the solver, and reports the final state.
//
// Usage:
//   diag_selfloop <graph_file> [start_marker]
//   diag_selfloop --buckygen N idx
//
// Where <graph_file> is a text file with a section like:
//   GRAPH FAIL Cn #idx N=80
//     0: 75 74 38 20 64 79
//     1: 14 22 15 33 53 43
//     ...
//   END_GRAPH
//
// start_marker is optional; if given, finds the first GRAPH section whose
// header contains start_marker (e.g. "C156 #4229114").

#include "fullerenes/triangulation.hh"
#include "fullerenes/delaunay.hh"
#include "fullerenes/delaunay_alexandrov12.hh"
#include "fullerenes/buckygen-wrapper.hh"

#include <cstdio>
#include <cmath>
#include <vector>
#include <set>
#include <map>
#include <string>
#include <fstream>
#include <sstream>

using namespace std;

// ---------- Δ-complex inspection helpers ----------

struct EdgeStats {
  int n_self_loops = 0;
  int n_multi_pairs = 0;          // distinct (u,v) pairs with >1 edge between them
  int n_multi_total_edges = 0;    // total edges minus distinct-vertex-pair count
  vector<pair<int,int>> self_loops;             // (vertex, half-edge h)
  map<pair<int,int>, vector<int>> by_pair;      // (min(u,v),max(u,v)) -> list of half-edges (one per edge)
};

static EdgeStats inspect(const DelaunayTriangulation& T) {
  EdgeStats s;
  for (int h = 0; h < T.nh; h += 2) {
    if (!T.alive(h)) continue;
    int u = T.he_origin[h], v = T.dest(h);
    if (u == v) {
      s.n_self_loops++;
      s.self_loops.push_back({u, h});
    } else {
      auto key = make_pair(min(u,v), max(u,v));
      s.by_pair[key].push_back(h);
    }
  }
  for (auto& [k, hs] : s.by_pair) {
    if (hs.size() > 1) {
      s.n_multi_pairs++;
      s.n_multi_total_edges += (int)hs.size() - 1;  // extra copies beyond the first
    }
  }
  return s;
}

static void report(const char* tag, const DelaunayTriangulation& T) {
  auto s = inspect(T);
  printf("[%s] nv=%d  edges(unique pairs)=%zu  self-loops=%d  multi-edge pairs=%d (extra copies=%d)\n",
         tag, T.nv, (size_t)s.by_pair.size() + s.n_self_loops, s.n_self_loops, s.n_multi_pairs, s.n_multi_total_edges);
  if (s.n_self_loops > 0) {
    printf("  self-loops at vertices:");
    for (auto [v, h] : s.self_loops) printf(" v=%d(h=%d,L=%.6f)", v, h, T.he_length[h]);
    printf("\n");
  }
  if (s.n_multi_pairs > 0) {
    printf("  multi-edge pairs:\n");
    for (auto& [k, hs] : s.by_pair) {
      if (hs.size() <= 1) continue;
      printf("    (%d,%d): %zu edges, lengths=", k.first, k.second, hs.size());
      for (int h : hs) printf(" %.6f", T.he_length[h]);
      printf("\n");
    }
  }
}

// ---------- Graph parsing ----------

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
          in_graph = true;
          found = true;
          printf("Loading: %s\n", line.c_str());
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
  if (!found) { fprintf(stderr, "no graph found matching '%s' in %s\n", marker.c_str(), path.c_str()); exit(1); }
  Graph G(N, GRAPH_DMAX);
  for (int v = 0; v < N; v++) G.assign_row(v, rows[v]);
  return G;
}

// ---------- Buckygen path ----------

static Triangulation from_buckygen(int N, int target_idx) {
  auto Q = BuckyGen::start(N, false, false);
  Triangulation T;
  int idx = 0;
  while (BuckyGen::next_fullerene(Q, T)) {
    if (idx == target_idx) {
      BuckyGen::stop(Q);
      return T;
    }
    idx++;
  }
  BuckyGen::stop(Q);
  fprintf(stderr, "buckygen exhausted before index %d (N=%d)\n", target_idx, N);
  exit(1);
}

// ---------- main ----------

int main(int argc, char** argv) {
  Triangulation T_dual;
  if (argc >= 4 && string(argv[1]) == "--buckygen") {
    int N = atoi(argv[2]), idx = atoi(argv[3]);
    T_dual = from_buckygen(N, idx);
    printf("Loaded buckygen #%d N=%d (dual nv=%d)\n", idx, N, T_dual.N);
  } else if (argc >= 2) {
    string path = argv[1], marker = (argc >= 3 ? argv[2] : "");
    Graph G = parse_graph_dump(path, marker);
    T_dual = Triangulation(G);
    printf("Loaded dual graph: %d vertices\n", T_dual.N);
  } else {
    fprintf(stderr, "Usage:\n  %s <graph_file> [marker]\n  %s --buckygen N idx\n", argv[0], argv[0]);
    return 1;
  }

  // Compute iDT.
  auto D_input = DelaunayTriangulation::compute(T_dual);
  printf("\n=== Input iDT ===\n");
  report("input iDT", D_input);

  // Run solver
  AlexandrovSolver solver;
  solver.D = D_input;
  solver.verbose = true;
  auto coords = solver.solve();
  printf("\n=== After solver ===\n");
  printf("steps=%d flips=%d newton=%d  final|κ|=%.3e (%s)\n",
         solver.stats_steps, solver.stats_flips, solver.stats_newton_total,
         solver.stats_final_kappa,
         coords.empty() ? "FAILED" : "succeeded");
  printf("  T(0) simplicial:        %s\n",
         solver.stats_t0_simplicial ? "yes" : "NO (invariant I-1 violated)");
  printf("  T̄(0) simple polygonal:  %s  (n_cells=%d)\n",
         solver.stats_tbar_simple_polygonal ? "yes" : "NO (invariant I-1 violated)",
         solver.stats_tbar_n_cells);
  report("final T(0)", solver.D);

  // Compare: were self-loops/multi-edges in input still present at the end?
  auto si = inspect(D_input);
  auto sf = inspect(solver.D);
  printf("\n=== Δ during solve ===\n");
  printf("  self-loops:   input=%d  final=%d  (delta %+d)\n",
         si.n_self_loops, sf.n_self_loops, sf.n_self_loops - si.n_self_loops);
  printf("  multi pairs:  input=%d  final=%d  (delta %+d)\n",
         si.n_multi_pairs, sf.n_multi_pairs, sf.n_multi_pairs - si.n_multi_pairs);
  printf("  total flips during solve: %d\n", solver.stats_flips);

  // Compare canonical (cocircular) vs B-I (inessential) tesselations of T(0).
  // For a non-degenerate convex polytope these should coincide on flat-2-face
  // diagonals; divergence flags an interesting case.
  vector<int> labels(solver.D.nv);
  for (int v = 0; v < solver.D.nv; v++) labels[v] = v;
  auto coc_tess = solver.D.canonical_tesselation(labels);
  auto bi_tess  = AlexandrovSolver::polytope_tesselation(solver.D, solver.r, labels);
  auto bi_mask  = AlexandrovSolver::inessential_edges(solver.D, solver.r);
  auto coc_mask = solver.D.cocircular_edges();
  int n_coc = 0, n_bi = 0, n_both = 0, n_only_coc = 0, n_only_bi = 0;
  for (int h = 0; h < solver.D.nh; h += 2) {
    if (!solver.D.alive(h)) continue;
    bool c = coc_mask[h], b = bi_mask[h];
    if (c) n_coc++;
    if (b) n_bi++;
    if (c && b) n_both++;
    if (c && !b) n_only_coc++;
    if (!c && b) n_only_bi++;
  }
  printf("\n=== Tesselations of T(0) ===\n");
  printf("  canonical (cocircular):  %d cells, %d tight edges\n",
         coc_tess.n_cells(), n_coc);
  printf("  B-I (inessential, θ=π):  %d cells, %d tight edges\n",
         bi_tess.n_cells(), n_bi);
  printf("  agree: %d  only-cocircular: %d  only-inessential: %d\n",
         n_both, n_only_coc, n_only_bi);
  printf("  fingerprint match: %s\n",
         coc_tess.fingerprint() == bi_tess.fingerprint() ? "YES" : "NO");

  // For each self-loop or multi-edge in the FINAL T(0): check whether ℓ ≈ chord (3D).
  if (!coords.empty()) {
    printf("\n=== Phantom-edge check (ℓ_intrinsic vs 3D chord) at final T(0) ===\n");
    for (int h = 0; h < solver.D.nh; h += 2) {
      if (!solver.D.alive(h)) continue;
      int u = solver.D.he_origin[h], v = solver.D.dest(h);
      double L = solver.D.he_length[h];
      double chord = (u == v) ? 0.0 : (coords[u] - coords[v]).norm();
      double rel_err = (L > 1e-15) ? fabs(L - chord) / L : 0.0;
      bool is_self_loop = (u == v);
      // List self-loops and multi-edge copies in particular
      auto key = make_pair(min(u,v), max(u,v));
      bool is_multi = !is_self_loop && sf.by_pair.count(key) && sf.by_pair[key].size() > 1;
      if (is_self_loop || is_multi || rel_err > 1e-6)
        printf("  h=%d edge (%d,%d) ℓ=%.6f chord=%.6f rel_err=%.3e %s%s\n",
               h, u, v, L, chord, rel_err,
               is_self_loop ? "[SELF-LOOP]" : "",
               is_multi ? "[MULTI]" : "");
    }
  }
  return 0;
}
