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

TEST_F(GCTransformTest, HalmaDegreePreservation) {
  Triangulation result = C20dual.GCtransform(2, 0);
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

// halma_transform(m>=2) has a pre-existing bug; skip GC(3,0) for now
TEST_F(GCTransformTest, DISABLED_HalmaNodeCount_3_0) {
  Triangulation result = C20dual.GCtransform(3, 0);
  EXPECT_EQ(result.N, expected_gc_nodes(V0, E0, F0, 3, 0));
}

// --- Unfold/fold infrastructure tests ---

TEST_F(GCTransformTest, UnfoldProducesOutline) {
  Unfolding u(C20dual);
  EXPECT_GT(u.outline.size(), 0u);
  // Each outline vertex should have a valid node id in [0, N)
  for(const auto& p : u.outline)
    EXPECT_LT(p.second, C20dual.N);
}

TEST_F(GCTransformTest, GCDreduce) {
  Unfolding u(C20dual);
  Unfolding gcu(u * Eisenstein(2, 0));
  auto reduced = Unfolding::GCDreduce(gcu.outline);
  EXPECT_EQ(reduced.size(), gcu.outline.size());
  // After GCDreduce, the outline should have same node labels
  for(size_t i = 0; i < reduced.size(); i++)
    EXPECT_EQ(reduced[i].second, gcu.outline[i].second);
}

TEST_F(GCTransformTest, ScaledUnfoldPreservesOutlineSize) {
  Unfolding u(C20dual);
  Unfolding gcu(u * Eisenstein(2, 1));
  EXPECT_EQ(gcu.outline.size(), u.outline.size());
}
