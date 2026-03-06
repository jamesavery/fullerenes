#include <gtest/gtest.h>
#include "fullerenes/delaunay.hh"
#include "fullerenes/isomerdb.hh"

// Build a fullerene dual triangulation from N and isomer index.
static Triangulation make_dual(int N, int idx, bool IPR = false) {
  IsomerDB db = IsomerDB::readPDB(N, IPR);
  FullereneGraph G = IsomerDB::makeIsomer(N, db.entries[idx]);
  // Get the dual triangulation from the cubic graph
  PlanarGraph PG(G);
  Triangulation T(PG.dual_graph());
  return T;
}

// Verify basic properties of the reduced triangulation.
static void verify_reduced(const FulleroidDelaunay& D, int expected_verts) {
  // Correct number of vertices
  EXPECT_EQ(D.N, expected_verts);

  // Edge lengths are symmetric
  EXPECT_TRUE(D.edge_lengths_are_symmetric());

  // All edge lengths are positive
  for (node_t u = 0; u < D.N; u++)
    for (node_t v : D.neighbours[u])
      EXPECT_GT(D.get_length(u, v), 0) << "Edge (" << u << "," << v << ") has non-positive length";

  // Triangulation is Delaunay
  EXPECT_TRUE(D.is_delaunay()) << "Result is not Delaunay";

  // Euler formula for genus 0: V - E + F = 2
  // With V vertices and all faces triangular: E = 3V - 6, F = 2V - 4
  int V = D.N;
  int E = 0;
  for (node_t u = 0; u < D.N; u++)
    E += D.neighbours[u].size();
  E /= 2;

  if (expected_verts == 12) {
    // For 12 vertices on genus-0: E=30, F=20
    EXPECT_EQ(E, 30) << "Expected 30 edges for 12 vertices on genus 0";
  }

  // V - E + F = 2, F = E - V + 2
  int F = E - V + 2;
  EXPECT_EQ(F, 2 * V - 4) << "Euler formula check failed";

  // Check orientation is consistent
  EXPECT_TRUE(D.is_consistently_oriented()) << "Orientation is broken";
}

// ============================================================================
// Tests
// ============================================================================

TEST(IntrinsicDelaunay, C20) {
  // C20 has 12 vertices all degree-5, no flat vertices to remove.
  Triangulation T = make_dual(20, 0, false);
  FulleroidDelaunay D(T);
  D.remove_flat_vertices();
  verify_reduced(D, 12);
}

TEST(IntrinsicDelaunay, C24) {
  Triangulation T = make_dual(24, 0, false);
  FulleroidDelaunay D(T);
  D.remove_flat_vertices();
  verify_reduced(D, 12);
}

TEST(IntrinsicDelaunay, C28) {
  // C28 has 2 isomers
  for (int idx = 0; idx < 2; idx++) {
    Triangulation T = make_dual(28, idx, false);
    FulleroidDelaunay D(T);
    D.remove_flat_vertices();
    verify_reduced(D, 12);
  }
}

TEST(IntrinsicDelaunay, C60_Ih) {
  // C60 IPR (buckminsterfullerene)
  Triangulation T = make_dual(60, 0, true);
  FulleroidDelaunay D(T);
  D.remove_flat_vertices();
  verify_reduced(D, 12);
}

TEST(IntrinsicDelaunay, C40_AllIsomers) {
  int n_isomers = IsomerDB::number_isomers(40, "Any", false);
  for (int idx = 0; idx < n_isomers; idx++) {
    Triangulation T = make_dual(40, idx, false);
    FulleroidDelaunay D(T);
    D.remove_flat_vertices();
    verify_reduced(D, 12);
  }
}

TEST(IntrinsicDelaunay, C60_Sample) {
  // Test a sample of C60 isomers
  int n_isomers = IsomerDB::number_isomers(60, "Any", false);
  int step = std::max(1, n_isomers / 50);  // ~50 samples
  int tested = 0;
  for (int idx = 0; idx < n_isomers; idx += step) {
    Triangulation T = make_dual(60, idx, false);
    FulleroidDelaunay D(T);
    D.remove_flat_vertices();
    verify_reduced(D, 12);
    tested++;
  }
  std::cout << "Tested " << tested << " / " << n_isomers << " C60 isomers" << std::endl;
}

TEST(IntrinsicDelaunay, AllEquilateral) {
  // For C20 (icosahedron), all edges should remain length 1
  // since all vertices are degree-5 (no flat vertices to remove).
  Triangulation T = make_dual(20, 0, false);
  FulleroidDelaunay D(T);
  D.remove_flat_vertices();

  for (node_t u = 0; u < D.N; u++)
    for (node_t v : D.neighbours[u])
      EXPECT_DOUBLE_EQ(D.get_length(u, v), 1.0);
}

TEST(IntrinsicDelaunay, EdgeLengthsPositive) {
  // After removing flat vertices from C60 Ih, all edge lengths should be > 0
  Triangulation T = make_dual(60, 0, true);
  FulleroidDelaunay D(T);
  D.remove_flat_vertices();

  double min_len = 1e30, max_len = 0;
  for (node_t u = 0; u < D.N; u++)
    for (node_t v : D.neighbours[u]) {
      double len = D.get_length(u, v);
      min_len = std::min(min_len, len);
      max_len = std::max(max_len, len);
    }

  std::cout << "C60 Ih reduced: edge lengths in [" << min_len << ", " << max_len << "]" << std::endl;
  EXPECT_GT(min_len, 0);
  EXPECT_LT(max_len, 100);  // sanity check
}
