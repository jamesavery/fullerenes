#include <gtest/gtest.h>
#include <string>
#include <algorithm>
#include "fullerenes/triangulation.hh"
#include "fullerenes/spiral.hh"

using namespace std;

// C28 dual spiral
static const vector<int> C28_spiral = {5,5,5,6,5,6,5,6,5,5,5,5,5,5,5,6};

// Representative GC(k,l) pairs covering small to medium sizes.
// N_carbon = 2*(V0-2)*T = 28*T where T = k^2 + kl + l^2, V0 = 16.
static const pair<int,int> test_cases[] = {
  {1, 0},   // T=1,    N=28
  {2, 0},   // T=4,    N=112
  {2, 1},   // T=7,    N=196
  {3, 0},   // T=9,    N=252
  {3, 1},   // T=13,   N=364
  {3, 2},   // T=19,   N=532
  {4, 0},   // T=16,   N=448
  {4, 1},   // T=21,   N=588
  {4, 3},   // T=37,   N=1036
  {5, 0},   // T=25,   N=700
  {5, 2},   // T=39,   N=1092
  {6, 0},   // T=36,   N=1008
  {6, 1},   // T=43,   N=1204
  {7, 0},   // T=49,   N=1372
  {7, 3},   // T=67,   N=1876
  {8, 0},   // T=64,   N=1792
  {10, 0},  // T=100,  N=2800
  {10, 3},  // T=139,  N=3892
  {12, 5},  // T=199,  N=5572
  {15, 0},  // T=225,  N=6300
  {15, 7},  // T=274,  N=7672
  {20, 0},  // T=400,  N=11200
  {25, 0},  // T=625,  N=17500
  {26, 0},  // T=676,  N=18928
};
static const int n_cases = sizeof(test_cases) / sizeof(test_cases[0]);

// Validate that a Triangulation is a valid fullerene dual
static void check_fullerene_dual(const Triangulation& T) {
  EXPECT_TRUE(T.is_consistently_oriented()) << "Orientation is inconsistent";

  int Nf = T.N;
  int expected_triangles = 2 * (Nf - 2);
  EXPECT_EQ((int)T.triangles.size(), expected_triangles)
    << "Triangle count violates Euler formula";

  int deg5 = 0, deg6 = 0, other = 0;
  for(int u = 0; u < Nf; u++) {
    int d = T.degree(u);
    if(d == 5) deg5++;
    else if(d == 6) deg6++;
    else other++;
  }
  EXPECT_EQ(deg5, 12) << "Must have exactly 12 degree-5 nodes";
  EXPECT_EQ(other, 0) << "All nodes should be degree 5 or 6";
}

// Check that nb_list is a cyclic rotation of expected_list
static bool is_cyclic_rotation(const vector<node_t>& a, const vector<node_t>& b) {
  if(a.size() != b.size()) return false;
  int n = a.size();
  for(int offset = 0; offset < n; offset++){
    bool match = true;
    for(int j = 0; j < n; j++){
      if(a[j] != b[(j + offset) % n]){ match = false; break; }
    }
    if(match) return true;
  }
  return false;
}

// Test: graph -> spiral -> graph round-trip using the permutation.
// get_spiral returns perm where perm[i] = G-node at spiral position i.
// G' (wound-up) has node i at spiral position i.
// So G'.neighbours[i] should be the relabeled version of G.neighbours[perm[i]],
// i.e. G' = inv_perm(G).
TEST(SpiralRoundtrip, PermutationIsomorphism) {
  Triangulation C28dual(C28_spiral);

  for(int i = 0; i < n_cases; i++) {
    auto [k, l] = test_cases[i];
    int T = k*k + k*l + l*l;
    int N_carbon = 28 * T;
    SCOPED_TRACE("GC(" + to_string(k) + "," + to_string(l) + ") N=" + to_string(N_carbon));

    // Generate graph via GC transform
    Triangulation G = C28dual.GCtransform(k, l);
    int N = G.N;

    // Compute canonical spiral with permutation
    vector<int> spiral_code;
    jumplist_t jumps;
    vector<vector<node_t>> permutations;
    G.get_spiral(spiral_code, jumps, permutations, true, true);

    ASSERT_FALSE(permutations.empty()) << "No spiral permutation found";
    const vector<node_t>& perm = permutations[0];

    ASSERT_EQ((int)perm.size(), N) << "Permutation size mismatch";

    // Build inverse permutation: inv_perm[perm[i]] = i
    vector<node_t> inv_perm(N);
    for(int u = 0; u < N; u++)
      inv_perm[perm[u]] = u;

    // Wind up the spiral to get G'
    Triangulation Gprime(spiral_code, jumps, true);

    ASSERT_EQ(G.N, Gprime.N)
      << "Node count mismatch: G has " << G.N << " vs G' has " << Gprime.N;

    // Verify G' = inv_perm(G):
    // For each node i in G', its neighbour list should be
    // {inv_perm[v] : v in G.neighbours[perm[i]]} in cyclic order.
    for(int u = 0; u < N; u++) {
      auto gprime_nb = Gprime.nbrs(u);
      auto g_nb = G.nbrs(perm[u]);

      ASSERT_EQ(gprime_nb.size(), g_nb.size())
        << "Degree mismatch at G'-node " << u << " (G-node " << perm[u] << ")";

      // Relabel G's neighbours through inv_perm
      vector<node_t> relabeled(g_nb.size());
      for(size_t j = 0; j < g_nb.size(); j++)
        relabeled[j] = inv_perm[g_nb[j]];

      // The canonical spiral may use a CW or CCW starting triple;
      // the windup always builds CW. So the round-trip may invert orientation.
      vector<node_t> relabeled_rev(relabeled.rbegin(), relabeled.rend());

      EXPECT_TRUE(is_cyclic_rotation(vector<node_t>(gprime_nb.begin(), gprime_nb.end()), relabeled) ||
                  is_cyclic_rotation(vector<node_t>(gprime_nb.begin(), gprime_nb.end()), relabeled_rev))
        << "Neighbour list mismatch at G'-node " << u << " (G-node " << perm[u] << ")";
    }
  }
}

// Test: graph -> spiral -> graph round-trip with CW_only=true.
// When restricted to CW starting triples, orientation is preserved through
// the round-trip, so we can check exact cyclic match (no reversal needed).
TEST(SpiralRoundtrip, PermutationIsomorphismCW) {
  Triangulation C28dual(C28_spiral);

  for(int i = 0; i < n_cases; i++) {
    auto [k, l] = test_cases[i];
    int T = k*k + k*l + l*l;
    int N_carbon = 28 * T;
    SCOPED_TRACE("GC(" + to_string(k) + "," + to_string(l) + ") N=" + to_string(N_carbon));

    // Generate graph via GC transform
    Triangulation G = C28dual.GCtransform(k, l);
    int N = G.N;

    // Compute canonical CW-only spiral with permutation
    vector<int> spiral_code;
    jumplist_t jumps;
    vector<vector<node_t>> permutations;
    G.get_spiral(spiral_code, jumps, permutations, true, true, true);  // CW_only=true

    ASSERT_FALSE(permutations.empty()) << "No spiral permutation found";
    const vector<node_t>& perm = permutations[0];

    ASSERT_EQ((int)perm.size(), N) << "Permutation size mismatch";

    // Build inverse permutation: inv_perm[perm[i]] = i
    vector<node_t> inv_perm(N);
    for(int u = 0; u < N; u++)
      inv_perm[perm[u]] = u;

    // Wind up the spiral to get G'
    Triangulation Gprime(spiral_code, jumps, true);

    ASSERT_EQ(G.N, Gprime.N)
      << "Node count mismatch: G has " << G.N << " vs G' has " << Gprime.N;

    // Verify G' = inv_perm(G) with EXACT cyclic match (orientation preserved):
    for(int u = 0; u < N; u++) {
      auto gprime_nb = Gprime.nbrs(u);
      auto g_nb = G.nbrs(perm[u]);

      ASSERT_EQ(gprime_nb.size(), g_nb.size())
        << "Degree mismatch at G'-node " << u << " (G-node " << perm[u] << ")";

      // Relabel G's neighbours through inv_perm
      vector<node_t> relabeled(g_nb.size());
      for(size_t j = 0; j < g_nb.size(); j++)
        relabeled[j] = inv_perm[g_nb[j]];

      // CW_only ensures orientation is preserved — exact cyclic match required
      EXPECT_TRUE(is_cyclic_rotation(vector<node_t>(gprime_nb.begin(), gprime_nb.end()), relabeled))
        << "Neighbour list mismatch at G'-node " << u << " (G-node " << perm[u] << ")";
    }
  }
}

// Test: spiral -> graph -> spiral round-trip.
// Wind up a spiral, compute canonical spiral of the result, verify it matches.
TEST(SpiralRoundtrip, SpiralToGraphToSpiral) {
  Triangulation C28dual(C28_spiral);

  for(int i = 0; i < n_cases; i++) {
    auto [k, l] = test_cases[i];
    int T = k*k + k*l + l*l;
    int N_carbon = 28 * T;
    SCOPED_TRACE("GC(" + to_string(k) + "," + to_string(l) + ") N=" + to_string(N_carbon));

    // Step 1: Generate graph and compute its canonical spiral
    Triangulation graph = C28dual.GCtransform(k, l);
    spiral_nomenclature sn1(graph, spiral_nomenclature::FULLERENE,
                            spiral_nomenclature::TRIANGULATION, true);

    // Step 2: Wind up from the spiral
    Triangulation from_spiral(sn1);

    // Step 3: Compute canonical spiral again from wound-up graph
    spiral_nomenclature sn2(from_spiral, spiral_nomenclature::FULLERENE,
                            spiral_nomenclature::TRIANGULATION, true);

    // Verify: canonical spiral data (jumps + spiral_code) matches.
    // (The to_string() may differ in construction_scheme prefix due to
    //  auto-detection, but the underlying spiral is the invariant.)
    EXPECT_EQ(sn1.spiral, sn2.spiral)
      << "Spiral data differs after spiral->graph->spiral round-trip";
  }
}
