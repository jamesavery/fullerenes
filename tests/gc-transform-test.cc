#include <gtest/gtest.h>
#include "fullerenes/triangulation.hh"
#include "fullerenes/unfold.hh"

using namespace std;

// Expected node count for GC(k,l) on a triangulation with V vertices, E edges, F faces:
//   T = k^2 + kl + l^2  (number of sub-triangles per original triangle)
//   F_new = F * T,  E_new = F_new * 3 / 2,  V_new = 2 + E_new - F_new  (Euler)
static int expected_gc_nodes(int V, int E, int F, int k, int l) {
  int T = k*k + k*l + l*l;
  int F_new = F * T;
  int E_new = F_new * 3 / 2;
  return 2 + E_new - F_new;
}

// Check that a triangulation has exactly 12 degree-5 nodes and the rest degree-6
static void check_degree_preservation(const Triangulation& result) {
  int deg5_count = 0, deg6_count = 0, other = 0;
  for(int u = 0; u < result.N; u++) {
    int d = result.neighbours[u].size();
    if(d == 5) deg5_count++;
    else if(d == 6) deg6_count++;
    else other++;
  }
  EXPECT_EQ(deg5_count, 12);
  EXPECT_EQ(other, 0) << "All nodes should be degree 5 or 6";
  EXPECT_EQ(deg6_count, result.N - 12);
}

// ===== C20 dual fixture =====

class GCTransformTest : public ::testing::Test {
protected:
  // C20 dual: icosahedron with 12 nodes, all degree 5
  Triangulation C20dual = Triangulation(vector<int>(12, 5));
  // Icosahedron: V=12, E=30, F=20
  static constexpr int V0 = 12, E0 = 30, F0 = 20;
};

// --- Halma (l==0) path tests ---

TEST_F(GCTransformTest, IdentityPreservesN) {
  Triangulation result = C20dual.GCtransform(1, 0);
  EXPECT_EQ(result.N, C20dual.N);
}

TEST_F(GCTransformTest, HalmaNodeCount_2_0) {
  Triangulation result = C20dual.GCtransform(2, 0);
  EXPECT_EQ(result.N, expected_gc_nodes(V0, E0, F0, 2, 0));
}

TEST_F(GCTransformTest, HalmaNodeCount_3_0) {
  Triangulation result = C20dual.GCtransform(3, 0);
  EXPECT_EQ(result.N, expected_gc_nodes(V0, E0, F0, 3, 0));
}

TEST_F(GCTransformTest, HalmaNodeCount_4_0) {
  Triangulation result = C20dual.GCtransform(4, 0);
  EXPECT_EQ(result.N, expected_gc_nodes(V0, E0, F0, 4, 0));
}

TEST_F(GCTransformTest, HalmaDegreePreservation_2_0) {
  Triangulation result = C20dual.GCtransform(2, 0);
  check_degree_preservation(result);
}

TEST_F(GCTransformTest, HalmaDegreePreservation_3_0) {
  Triangulation result = C20dual.GCtransform(3, 0);
  check_degree_preservation(result);
}

TEST_F(GCTransformTest, HalmaConnectivity) {
  Triangulation result = C20dual.GCtransform(2, 0);
  vector<vector<node_t>> components = result.connected_components();
  EXPECT_EQ(components.size(), 1u);
}

TEST_F(GCTransformTest, HalmaTriangleCount) {
  Triangulation result = C20dual.GCtransform(2, 0);
  // F_new = F * T = 20 * 4 = 80
  EXPECT_EQ((int)result.triangles.size(), F0 * 4);
}

// --- Chiral (l!=0, unfold/fold path) tests ---

TEST_F(GCTransformTest, ChiralNodeCount_1_1) {
  Triangulation result = C20dual.GCtransform(1, 1);
  EXPECT_EQ(result.N, expected_gc_nodes(V0, E0, F0, 1, 1));
}

TEST_F(GCTransformTest, ChiralNodeCount_2_1) {
  Triangulation result = C20dual.GCtransform(2, 1);
  EXPECT_EQ(result.N, expected_gc_nodes(V0, E0, F0, 2, 1));
}

TEST_F(GCTransformTest, ChiralNodeCount_3_1) {
  Triangulation result = C20dual.GCtransform(3, 1);
  EXPECT_EQ(result.N, expected_gc_nodes(V0, E0, F0, 3, 1));
}

TEST_F(GCTransformTest, ChiralDegreePreservation_1_1) {
  Triangulation result = C20dual.GCtransform(1, 1);
  check_degree_preservation(result);
}

TEST_F(GCTransformTest, ChiralDegreePreservation_2_1) {
  Triangulation result = C20dual.GCtransform(2, 1);
  check_degree_preservation(result);
}

TEST_F(GCTransformTest, ChiralConnectivity_2_1) {
  Triangulation result = C20dual.GCtransform(2, 1);
  vector<vector<node_t>> components = result.connected_components();
  EXPECT_EQ(components.size(), 1u);
}

TEST_F(GCTransformTest, ChiralDegreePreservation_3_1) {
  Triangulation result = C20dual.GCtransform(3, 1);
  check_degree_preservation(result);
}

TEST_F(GCTransformTest, ChiralNodeCount_3_2) {
  Triangulation result = C20dual.GCtransform(3, 2);
  EXPECT_EQ(result.N, expected_gc_nodes(V0, E0, F0, 3, 2));
}

TEST_F(GCTransformTest, ChiralDegreePreservation_3_2) {
  Triangulation result = C20dual.GCtransform(3, 2);
  check_degree_preservation(result);
}

TEST_F(GCTransformTest, ChiralNodeCount_4_1) {
  Triangulation result = C20dual.GCtransform(4, 1);
  EXPECT_EQ(result.N, expected_gc_nodes(V0, E0, F0, 4, 1));
}

TEST_F(GCTransformTest, ChiralDegreePreservation_4_1) {
  Triangulation result = C20dual.GCtransform(4, 1);
  check_degree_preservation(result);
}

TEST_F(GCTransformTest, ChiralNodeCount_4_2) {
  Triangulation result = C20dual.GCtransform(4, 2);
  EXPECT_EQ(result.N, expected_gc_nodes(V0, E0, F0, 4, 2));
}

TEST_F(GCTransformTest, ChiralDegreePreservation_4_2) {
  Triangulation result = C20dual.GCtransform(4, 2);
  check_degree_preservation(result);
}

TEST_F(GCTransformTest, ChiralNodeCount_4_3) {
  Triangulation result = C20dual.GCtransform(4, 3);
  EXPECT_EQ(result.N, expected_gc_nodes(V0, E0, F0, 4, 3));
}

TEST_F(GCTransformTest, ChiralDegreePreservation_4_3) {
  Triangulation result = C20dual.GCtransform(4, 3);
  check_degree_preservation(result);
}

// --- Unfold/fold infrastructure tests ---

TEST_F(GCTransformTest, UnfoldProducesOutline) {
  Unfolding u(C20dual);
  EXPECT_GT(u.outline.size(), 0u);
  for(const auto& p : u.outline)
    EXPECT_LT(p.second, C20dual.N);
}

TEST_F(GCTransformTest, GCDreduce) {
  Unfolding u(C20dual);
  Unfolding gcu(u * Eisenstein(2, 0));
  auto reduced = Unfolding::GCDreduce(gcu.outline);
  EXPECT_EQ(reduced.size(), gcu.outline.size());
  for(size_t i = 0; i < reduced.size(); i++)
    EXPECT_EQ(reduced[i].second, gcu.outline[i].second);
}

TEST_F(GCTransformTest, GCDreduce_Segments) {
  Unfolding u(C20dual);
  Unfolding gcu(u * Eisenstein(3, 0));
  auto reduced = Unfolding::GCDreduce(gcu.outline);
  // Each segment in the reduced outline should be 1/3 the length of the original
  for(size_t i = 0; i < reduced.size(); i++) {
    Eisenstein orig_seg = gcu.outline[(i+1)%gcu.outline.size()].first - gcu.outline[i].first;
    Eisenstein red_seg  = reduced[(i+1)%reduced.size()].first - reduced[i].first;
    EXPECT_EQ(orig_seg, red_seg * Eisenstein(3, 0));
  }
}

TEST_F(GCTransformTest, ScaledUnfoldPreservesOutlineSize) {
  Unfolding u(C20dual);
  Unfolding gcu(u * Eisenstein(2, 1));
  EXPECT_EQ(gcu.outline.size(), u.outline.size());
}

TEST_F(GCTransformTest, FoldUnfoldRoundtrip) {
  Unfolding u(C20dual);
  Unfolding scaled(u * Eisenstein(1, 0));
  Folding f(scaled);
  Triangulation result = f.fold();
  EXPECT_EQ(result.N, C20dual.N);
}

// ===== C28 dual fixture =====

class C28GCTest : public ::testing::Test {
protected:
  // C28 dual from spiral representation: 16 nodes
  Triangulation C28dual = Triangulation(vector<int>({5,5,5,6,5,6,5,6,5,5,5,5,5,5,5,6}));
  // C28 dual: V=16, E=42, F=28
  static constexpr int V0 = 16, E0 = 42, F0 = 28;
};

TEST_F(C28GCTest, Identity) {
  Triangulation result = C28dual.GCtransform(1, 0);
  EXPECT_EQ(result.N, C28dual.N);
}

TEST_F(C28GCTest, Halma_2_0_NodeCount) {
  Triangulation result = C28dual.GCtransform(2, 0);
  EXPECT_EQ(result.N, expected_gc_nodes(V0, E0, F0, 2, 0));
}

TEST_F(C28GCTest, Halma_2_0_Degrees) {
  Triangulation result = C28dual.GCtransform(2, 0);
  check_degree_preservation(result);
}

TEST_F(C28GCTest, Halma_3_0_NodeCount) {
  Triangulation result = C28dual.GCtransform(3, 0);
  EXPECT_EQ(result.N, expected_gc_nodes(V0, E0, F0, 3, 0));
}

TEST_F(C28GCTest, Halma_3_0_Degrees) {
  Triangulation result = C28dual.GCtransform(3, 0);
  check_degree_preservation(result);
}

// --- C28 chiral tests ---
TEST_F(C28GCTest, Chiral_1_1_NodeCount) {
  Triangulation result = C28dual.GCtransform(1, 1);
  EXPECT_EQ(result.N, expected_gc_nodes(V0, E0, F0, 1, 1));
}

TEST_F(C28GCTest, Chiral_2_1_NodeCount) {
  Triangulation result = C28dual.GCtransform(2, 1);
  EXPECT_EQ(result.N, expected_gc_nodes(V0, E0, F0, 2, 1));
}
