// Correctness tests for GraphView::hamiltonian_cycle_count (src/c++/hamiltonian.cc).
//
// Three independent oracles:
//
//   * small graphs whose Hamiltonian-cycle count is known in closed form
//     (K_n, K_a,b, the cube, the Petersen graph) -- the "any degree"
//     generalisation, which no fullerene exercises;
//   * an unpruned enumeration on random graphs of every density, which the
//     pruned search must match exactly;
//   * the legacy Fortran program's own data: the enumerated min/max counts
//     per size (hamilton.f PathStatistic tables) and the per-isomer ncycham
//     field of the isomer database.
//
// All three count UNDIRECTED cycles once -- the convention the function
// implements.  Only the third needs the database, and skips without it.

#include "fullerene-test-main.hh"

#include "fullerenes/fullerenegraph.hh"
#include "fullerenes/graph.hh"
#include "fullerenes/isomerdb.hh"

#include <algorithm>
#include <cstdint>
#include <cstdio>
#include <functional>
#include <optional>
#include <random>
#include <span>
#include <utility>
#include <vector>

namespace {

// hamilton.f PathStatistic data tables (@ref hamilton.f:508-522): the exact
// enumerated extremes of the Hamiltonian-cycle count over all isomers of a
// given size, general and IPR.
const int ham_min[] = {30,0,34,24,18,20,40,28,42,24,68,44,120,76,152,80,
  262,66,440,173,618,288,1062,197,1750,320,2688,1182,4230,1596,
  7110,2400,10814,1980,17905,1280,29944,7930,46231,13307,72168,
  20754,119540,40912,184445,5120};
const int ham_max[] = {30,0,34,24,43,32,76,66,128,96,280,150,327,260,512,
  410,806,642,1746,1068,3040,1802,3340,3096,6018,4818,10428,7832,
  15926,12226,35200,20856,39067,33427,76063,51586,117106,90221,
  209692,156288,417280,249148,686286,421194,1104223,743346};
const int ham_ipr_min[] = {1090,0,0,0,0,2790,3852,4794,6078,6988,9004,11226,
  14748,17853,22661,29277,36949,44730,60070,71950,93986,35907,
  149920,180243,237580,244254,383218,457235,630059,723505,1038971,
  1368498};
const int ham_ipr_max[] = {1090,0,0,0,0,2790,3852,4794,6643,8244,10970,13614,
  18260,21756,28652,36852,47054,59118,78044,95694,131690,161148,
  207165,257746,351976,426750,571622,699908,1013844,1151918,
  1590875,1888558};

// One value per even size, starting at N_first: the shape all four tables
// have, so the index arithmetic is written once.
struct SizeTable {
  int N_first;
  std::span<const int> value;
  std::optional<int> at(int N) const {
    if (N < N_first || (N - N_first) % 2 != 0) return std::nullopt;
    const size_t k = size_t(N - N_first) / 2;
    return k < value.size() ? std::optional<int>(value[k]) : std::nullopt;
  }
};

struct CountRange { int64_t min, max; };

// The enumerated extremes for this size and corpus, where the Fortran tabulated
// them (general C20-C110, IPR C60-C122).
std::optional<CountRange> enumerated_range(int N, IsomerDB::Corpus corpus) {
  const bool IPR = corpus == IsomerDB::Corpus::IPR;
  const SizeTable lo{IPR ? 60 : 20, IPR ? std::span<const int>(ham_ipr_min) : std::span<const int>(ham_min)};
  const SizeTable hi{IPR ? 60 : 20, IPR ? std::span<const int>(ham_ipr_max) : std::span<const int>(ham_max)};
  const std::optional<int> mn = lo.at(N), mx = hi.at(N);
  if (!mn || !mx) return std::nullopt;
  return CountRange{*mn, *mx};
}

// A simple undirected graph from an edge list, in vector order.  The row width
// is the edge list's own maximum degree: these cases are not fullerenes and
// GRAPH_DMAX does not bound them.
Graph graph_of(int N, const std::vector<std::pair<int,int>>& edges) {
  std::vector<std::vector<node_t>> adj(N);
  for (auto [u, v] : edges) { adj[u].push_back(v); adj[v].push_back(u); }
  int width = 1;
  for (const auto& row : adj) width = std::max(width, int(row.size()));
  Graph g(size_t(N), width);
  for (int u = 0; u < N; u++) for (node_t v : adj[u]) if (u < v) g.insert_edge({u, v});
  return g;
}

Graph complete(int n) {
  std::vector<std::pair<int,int>> e;
  for (int u = 0; u < n; u++) for (int v = u + 1; v < n; v++) e.push_back({u, v});
  return graph_of(n, e);
}
Graph complete_bipartite(int a, int b) {
  std::vector<std::pair<int,int>> e;
  for (int u = 0; u < a; u++) for (int v = 0; v < b; v++) e.push_back({u, a + v});
  return graph_of(a + b, e);
}

// Unpruned reference: every directed Hamiltonian cycle from vertex 0, halved.
// N < 3 is excluded because the walk below indexes vertex 0, not because the
// function under test special-cases it -- the oracle must not inherit that.
int64_t brute_force_count(const Graph& g) {
  const int N = g.N;
  if (N < 1) return 0;
  std::vector<char> used(N, 0);
  int64_t directed = 0;
  std::function<void(node_t,int)> walk = [&](node_t u, int depth) {
    if (depth == N) { if (g.edge_exists({u, 0})) directed++; return; }
    for (node_t v : g.nbrs(u)) if (!used[v]) { used[v] = 1; walk(v, depth + 1); used[v] = 0; }
  };
  used[0] = 1; walk(0, 1);
  return N < 3 ? 0 : directed / 2;   // a 1- or 2-cycle is not a cycle
}

}  // namespace

// Closed-form counts, and the degenerate shapes: N < 3, a path, a star, a
// pendant vertex, a disconnected graph.
TEST(HamiltonianCount, KnownCounts) {
  const struct { const char* name; Graph g; int64_t expected; } cases[] = {
    {"K4",             complete(4),            3},
    {"K5",             complete(5),           12},
    {"K6",             complete(6),           60},
    {"K7",             complete(7),          360},
    {"K3,3",           complete_bipartite(3, 3), 6},
    {"K4,4",           complete_bipartite(4, 4), 72},
    {"K3,4",           complete_bipartite(3, 4),  0},   // odd bipartition
    {"cube Q3",        graph_of(8, {{0,1},{1,2},{2,3},{3,0},{4,5},{5,6},{6,7},{7,4},{0,4},{1,5},{2,6},{3,7}}), 6},
    {"Petersen",       graph_of(10, {{0,1},{1,2},{2,3},{3,4},{4,0},{0,5},{1,6},{2,7},{3,8},{4,9},{5,7},{7,9},{9,6},{6,8},{8,5}}), 0},
    {"empty",          graph_of(0, {}),        0},
    {"one vertex",     graph_of(1, {}),        0},
    {"one edge",       graph_of(2, {{0,1}}),   0},
    {"triangle",       graph_of(3, {{0,1},{1,2},{2,0}}), 1},
    {"path",           graph_of(4, {{0,1},{1,2},{2,3}}), 0},
    {"star",           graph_of(4, {{0,1},{0,2},{0,3}}), 0},
    {"pendant vertex", graph_of(5, {{0,1},{1,2},{2,3},{3,0},{0,2},{1,3},{4,0}}), 0},
    {"two triangles",  graph_of(6, {{0,1},{1,2},{2,0},{3,4},{4,5},{5,3}}), 0},
  };
  for (const auto& c : cases) {
    SCOPED_TRACE(c.name);
    EXPECT_EQ(c.g.hamiltonian_cycle_count(), c.expected);
  }
}

// Brute-force cross-check on small random graphs of every density.  The edge
// order is shuffled so the neighbour rotations vary, which is what the
// canonical-direction rule depends on; the vertex labelling is not permuted,
// so this does not exercise independence from the search's fixed start.
TEST(HamiltonianCount, MatchesBruteForceOnRandomGraphs) {
  std::mt19937 rng(20260819);
  int n_graphs = 0;
  for (int N = 3; N <= 9; N++)
    for (int trial = 0; trial < 300; trial++) {
      const double p = 0.15 + 0.85 * (trial % 10) / 9.0;
      std::vector<std::pair<int,int>> e;
      for (int u = 0; u < N; u++) for (int v = u + 1; v < N; v++)
        if (std::uniform_real_distribution<double>(0, 1)(rng) < p) e.push_back({u, v});
      std::shuffle(e.begin(), e.end(), rng);
      Graph g = graph_of(N, e);
      ASSERT_EQ(g.hamiltonian_cycle_count(), brute_force_count(g)) << "N=" << N << " trial " << trial;
      n_graphs++;
    }
  printf("[          ] %d random graphs agree with brute force\n", n_graphs);
}

// The dodecahedron, built from its spiral, so this needs no database.
TEST(HamiltonianCount, C20DodecahedronHas30) {
  FullereneGraph dodecahedron(20, {0,1,2,3,4,5,6,7,8,9,10,11});
  ASSERT_EQ(dodecahedron.N, 20);
  EXPECT_EQ(dodecahedron.hamiltonian_cycle_count(), 30);
}

// Every isomer of every requested size, in both corpora that carry ground
// truth: per-isomer equality with the database's ncycham, and min/max equality
// with the Fortran's enumerated tables.
TEST(HamiltonianCount, MatchesEnumeratedTablesAndDatabase) {
  ASSERT_FALSE(fullerene_test::sizes().empty()) << "--sizes= would make this sweep assert nothing";
  size_t n_percycle_checked = 0, n_sizes_swept = 0;
  for (int N : fullerene_test::sizes())
    for (IsomerDB::Corpus corpus : {IsomerDB::Corpus::All, IsomerDB::Corpus::IPR}) {
      const bool IPR = corpus == IsomerDB::Corpus::IPR;
      const char* tag = IPR ? " IPR" : "";
      if (!IsomerDB::is_installed(N, corpus)) continue;

      IsomerDB db = IsomerDB::readPDB(N, IPR);
      const std::optional<CountRange> table = enumerated_range(N, corpus);
      // A size with neither per-isomer counts nor a table would assert
      // nothing; refuse it before spending the enumeration, not after.
      ASSERT_TRUE(db.with_ncycham || table) << "C" << N << tag << ": no ground truth for this size";
      ASSERT_GT(db.entries.size(), 0u) << "empty database for C" << N << tag;

      int64_t mn = INT64_MAX, mx = 0;
      for (size_t i = 0; i < db.entries.size(); i++) {
        const int64_t c = IsomerDB::makeIsomer(N, db.entries[i]).hamiltonian_cycle_count();
        mn = std::min(mn, c);
        mx = std::max(mx, c);
        if (db.with_ncycham) {
          EXPECT_EQ(c, db.entries[i].ncycham) << "C" << N << tag << " isomer " << i + 1;
          n_percycle_checked++;
        }
      }
      if (table) {
        EXPECT_EQ(mn, table->min) << "min count mismatch at C" << N << tag;
        EXPECT_EQ(mx, table->max) << "max count mismatch at C" << N << tag;
      }
      n_sizes_swept++;
      printf("[          ] C%d%s: %zu isomers, counts in [%lld, %lld]\n",
             N, tag, db.entries.size(), (long long)mn, (long long)mx);
    }
  if (n_sizes_swept == 0) GTEST_SKIP() << "none of the requested sizes is installed";
  printf("[          ] per-isomer ncycham comparisons: %zu\n", n_percycle_checked);
}

int main(int argc, char** argv) {
  return fullerene_test::run(argc, argv, {20, 24, 26, 28, 30, 32, 34, 36, 38, 40}, "hamiltonian");
}
