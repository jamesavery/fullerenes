#include <gtest/gtest.h>
#include "fullerenes/triangulation.hh"
#include "fullerenes/unfold.hh"

#include <set>

using namespace std;

// Exact-label fidelity: same vertex count and identical undirected adjacency per id.
static bool same_neighbour_sets(const Triangulation& a, const Triangulation& b){
  if(a.N != b.N) return false;
  for(int u=0; u<a.N; u++){
    auto an = a.nbrs(u), bn = b.nbrs(u);
    if(set<node_t>(an.begin(),an.end()) != set<node_t>(bn.begin(),bn.end())) return false;
  }
  return true;
}

// Expected node count for GC(k,l) on a triangulation with V vertices, E edges, F faces:
//   T = k^2 + kl + l^2  (number of sub-triangles per original triangle)
//   F_new = F * T,  E_new = F_new * 3 / 2,  V_new = 2 + E_new - F_new  (Euler)
static int expected_gc_nodes(int V, int E, int F, int k, int l) {
  int T = k*k + k*l + l*l;
  int F_new = F * T;
  int E_new = F_new * 3 / 2;
  return 2 + E_new - F_new;
}

// Validate that a Triangulation is a valid fullerene dual:
//  1) Planar: orientation is consistent (every directed arc belongs to exactly one face)
//  2) Triangulation: every face is a triangle
//  3) Degree constraint: exactly 12 vertices of degree 5, rest degree 6
//  4) Euler relation: Nf = N_carbon/2 + 2, i.e. the number of triangles F = 2*(Nf-2)
static void check_fullerene_dual(const Triangulation& T) {
  // 1) Consistent orientation: every directed arc is part of exactly one face.
  //    This verifies planarity of the embedding.
  EXPECT_TRUE(T.is_consistently_oriented()) << "Orientation is inconsistent";

  // 2) Every face is a triangle: compute faces via next_on_face and check sizes.
  //    The triangles member is populated by the constructor, so verify it
  //    matches the expected count from Euler's formula.
  int Nf = T.N;
  int expected_triangles = 2 * (Nf - 2);
  EXPECT_EQ((int)T.triangles().size(), expected_triangles)
    << "Triangle count violates Euler formula: expected 2*(Nf-2) = " << expected_triangles;

  // Also verify each stored triangle is actually a valid face by tracing with next_on_face
  for(const auto& tri : T.triangles()) {
    node_t u = tri[0], v = tri[1], w = tri[2];
    EXPECT_EQ(T.next_on_face(u, v), w)
      << "Triangle {" << u << "," << v << "," << w << "}: next_on_face(" << u << "," << v << ") != " << w;
    EXPECT_EQ(T.next_on_face(v, w), u)
      << "Triangle {" << u << "," << v << "," << w << "}: next_on_face(" << v << "," << w << ") != " << u;
    EXPECT_EQ(T.next_on_face(w, u), v)
      << "Triangle {" << u << "," << v << "," << w << "}: next_on_face(" << w << "," << u << ") != " << v;
  }

  // 3) Degree constraint: exactly 12 degree-5 vertices, rest degree-6
  int deg5 = 0, deg6 = 0, other = 0;
  for(int u = 0; u < Nf; u++) {
    int d = T.degree(u);
    if(d == 5) deg5++;
    else if(d == 6) deg6++;
    else other++;
  }
  EXPECT_EQ(deg5, 12) << "Must have exactly 12 degree-5 nodes (pentagons)";
  EXPECT_EQ(other, 0) << "All nodes should be degree 5 or 6";
  EXPECT_EQ(deg6, Nf - 12);

  // 4) Nf = N_carbon/2 + 2, equivalently N_carbon = 2*(Nf - 2)
  //    The number of carbon atoms is the number of triangles (dual faces = cubic vertices).
  int N_carbon = expected_triangles;
  EXPECT_EQ(Nf, N_carbon / 2 + 2)
    << "Nf=" << Nf << " should equal N/2+2=" << N_carbon/2+2 << " where N=" << N_carbon;
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
  check_fullerene_dual(result);
}

TEST_F(GCTransformTest, HalmaDegreePreservation_3_0) {
  Triangulation result = C20dual.GCtransform(3, 0);
  check_fullerene_dual(result);
}

TEST_F(GCTransformTest, HalmaConnectivity) {
  Triangulation result = C20dual.GCtransform(2, 0);
  vector<vector<node_t>> components = result.connected_components();
  EXPECT_EQ(components.size(), 1u);
}

TEST_F(GCTransformTest, HalmaTriangleCount) {
  Triangulation result = C20dual.GCtransform(2, 0);
  // F_new = F * T = 20 * 4 = 80
  EXPECT_EQ((int)result.triangles().size(), F0 * 4);
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
  check_fullerene_dual(result);
}

TEST_F(GCTransformTest, ChiralDegreePreservation_2_1) {
  Triangulation result = C20dual.GCtransform(2, 1);
  check_fullerene_dual(result);
}

TEST_F(GCTransformTest, ChiralConnectivity_2_1) {
  Triangulation result = C20dual.GCtransform(2, 1);
  vector<vector<node_t>> components = result.connected_components();
  EXPECT_EQ(components.size(), 1u);
}

TEST_F(GCTransformTest, ChiralDegreePreservation_3_1) {
  Triangulation result = C20dual.GCtransform(3, 1);
  check_fullerene_dual(result);
}

TEST_F(GCTransformTest, ChiralNodeCount_3_2) {
  Triangulation result = C20dual.GCtransform(3, 2);
  EXPECT_EQ(result.N, expected_gc_nodes(V0, E0, F0, 3, 2));
}

TEST_F(GCTransformTest, ChiralDegreePreservation_3_2) {
  Triangulation result = C20dual.GCtransform(3, 2);
  check_fullerene_dual(result);
}

TEST_F(GCTransformTest, ChiralNodeCount_4_1) {
  Triangulation result = C20dual.GCtransform(4, 1);
  EXPECT_EQ(result.N, expected_gc_nodes(V0, E0, F0, 4, 1));
}

TEST_F(GCTransformTest, ChiralDegreePreservation_4_1) {
  Triangulation result = C20dual.GCtransform(4, 1);
  check_fullerene_dual(result);
}

TEST_F(GCTransformTest, ChiralNodeCount_4_2) {
  Triangulation result = C20dual.GCtransform(4, 2);
  EXPECT_EQ(result.N, expected_gc_nodes(V0, E0, F0, 4, 2));
}

TEST_F(GCTransformTest, ChiralDegreePreservation_4_2) {
  Triangulation result = C20dual.GCtransform(4, 2);
  check_fullerene_dual(result);
}

TEST_F(GCTransformTest, ChiralNodeCount_4_3) {
  Triangulation result = C20dual.GCtransform(4, 3);
  EXPECT_EQ(result.N, expected_gc_nodes(V0, E0, F0, 4, 3));
}

TEST_F(GCTransformTest, ChiralDegreePreservation_4_3) {
  Triangulation result = C20dual.GCtransform(4, 3);
  check_fullerene_dual(result);
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

// Exact label fidelity: fold(unfold(T)) reproduces T's neighbour sets per ORIGINAL id.
// DISABLED until assemble_cone handles all developments (fold-engine WIP); the
// restore logic is in place, but the fold itself is not yet valid on every path.
TEST_F(GCTransformTest, DISABLED_FoldUnfoldRoundtripRestoresLabels) {
  Unfolding u(C20dual);
  Folding f(u * Eisenstein(1, 0));
  Triangulation result = f.fold();
  EXPECT_EQ(result.N, C20dual.N);
  check_fullerene_dual(result);
  EXPECT_TRUE(same_neighbour_sets(result, C20dual))
    << "round-trip did not restore the original labelling";
}

// The Triangulation ctor relabels cones-first and stores the permutation.
TEST_F(GCTransformTest, ConePermThreading) {
  Unfolding u(C20dual);
  EXPECT_EQ(u.cone_perm.size(), (size_t)C20dual.N);
  // cones (degree != 6) occupy a prefix of the id range; no cone follows a hexagon.
  bool seen_hex = false;
  for(int v = 0; v < u.graph.N; v++){
    if(u.graph.degree(v) == 6) seen_hex = true;
    else EXPECT_FALSE(seen_hex) << "cone (deg!=6) after a hexagon at id " << v;
  }
}

// GCtransform goes through fold(): original cones keep ids 0..11 (their degrees),
// and GC-introduced vertices get fresh ids >= N_orig and are degree 6.
// DISABLED until assemble_cone handles chiral folds (fold-engine WIP).
TEST_F(GCTransformTest, DISABLED_ChiralPreservesOriginalLabels_2_1) {
  Triangulation result = C20dual.GCtransform(2, 1);
  check_fullerene_dual(result);
  for(int o = 0; o < C20dual.N; o++)
    EXPECT_EQ(result.degree(o), C20dual.degree(o)) << "original id " << o << " not preserved";
  for(int g = C20dual.N; g < result.N; g++)
    EXPECT_EQ(result.degree(g), 6) << "GC-introduced id " << g << " should be degree 6";
}

// Outline-only Unfolding has no permutation: fold() degrades gracefully (identity
// restore over the outline's own labels), still yielding a valid oriented dual.
// DISABLED until the fold engine produces valid duals on the outline-only path.
TEST_F(GCTransformTest, DISABLED_OutlineOnlyFoldGracefulDegradation) {
  Unfolding u(C20dual);
  Unfolding outline_only(u.outline);          // ctor 3: empty cone_perm
  EXPECT_EQ(outline_only.cone_perm.size(), 0u);
  Folding f(outline_only);
  Triangulation result = f.fold();
  EXPECT_TRUE(result.is_consistently_oriented());
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
  check_fullerene_dual(result);
}

TEST_F(C28GCTest, Halma_3_0_NodeCount) {
  Triangulation result = C28dual.GCtransform(3, 0);
  EXPECT_EQ(result.N, expected_gc_nodes(V0, E0, F0, 3, 0));
}

TEST_F(C28GCTest, Halma_3_0_Degrees) {
  Triangulation result = C28dual.GCtransform(3, 0);
  check_fullerene_dual(result);
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

// C28 has degree-6 vertices, so its cones-first permutation is non-trivial.
TEST_F(C28GCTest, ConePermThreading) {
  Unfolding u(C28dual);
  EXPECT_EQ(u.cone_perm.size(), (size_t)C28dual.N);
  // exactly the 12 cones occupy ids 0..11
  for(int v = 0; v < 12; v++) EXPECT_NE(u.graph.degree(v), 6) << "id " << v << " should be a cone";
  for(int v = 12; v < u.graph.N; v++) EXPECT_EQ(u.graph.degree(v), 6) << "id " << v << " should be a hexagon";
}

// Non-identity restore through fold(): fold(unfold(C28)) must reproduce C28's
// neighbour sets per ORIGINAL id (exercises cone_perm.inverse() on a real permutation).
// DISABLED until assemble_cone handles deg-6-bearing chiral folds (fold-engine WIP).
TEST_F(C28GCTest, DISABLED_FoldUnfoldRoundtripRestoresLabels) {
  Unfolding u(C28dual);
  Folding f(u * Eisenstein(1, 0));
  Triangulation result = f.fold();
  EXPECT_EQ(result.N, C28dual.N);
  check_fullerene_dual(result);
  EXPECT_TRUE(same_neighbour_sets(result, C28dual))
    << "round-trip did not restore the original C28 labelling";
}

// ===== Multi-isomer GC tests =====

struct IsomerSpec {
  const char* name;
  vector<int> spiral;
  // Derived: V = spiral.size(), E = 3V-6, F = 2V-4
};

static const IsomerSpec test_isomers[] = {
  {"C42-1", {5,5,5,5,6,6,5,6,6,6,6,5,6,6,5,6,5,5,6,5,6,5,5}},
  {"C42-2", {5,5,5,6,5,6,5,6,6,5,6,6,6,6,5,6,5,6,5,5,6,5,5}},
  {"C44-1", {5,5,5,6,5,6,6,6,6,5,6,6,6,5,5,6,5,6,5,5,6,5,6,5}},
  {"C44-2", {5,5,5,6,5,6,6,6,6,6,5,6,6,5,5,6,5,5,6,5,5,6,5,6}},
  {"C52",   {5,5,5,6,6,6,6,6,6,5,6,6,5,6,6,6,5,5,6,6,5,5,6,5,5,6,6,5}},
};

TEST(MultiIsomerGCTest, HalmaNodeCount) {
  for(const auto& iso : test_isomers) {
    SCOPED_TRACE(iso.name);
    int V = iso.spiral.size(), E = 3*V - 6, F = 2*V - 4;
    Triangulation dual(vector<int>(iso.spiral));

    EXPECT_EQ(dual.N, V);

    for(int k : {2, 3}) {
      SCOPED_TRACE("GC(" + to_string(k) + ",0)");
      Triangulation result = dual.GCtransform(k, 0);
      EXPECT_EQ(result.N, expected_gc_nodes(V, E, F, k, 0));
      check_fullerene_dual(result);
    }
  }
}

TEST(MultiIsomerGCTest, ChiralNodeCount) {
  for(const auto& iso : test_isomers) {
    SCOPED_TRACE(iso.name);
    int V = iso.spiral.size(), E = 3*V - 6, F = 2*V - 4;
    Triangulation dual(vector<int>(iso.spiral));

    for(auto [k, l] : vector<pair<int,int>>{{1,1}, {2,1}, {3,1}}) {
      SCOPED_TRACE("GC(" + to_string(k) + "," + to_string(l) + ")");
      Triangulation result = dual.GCtransform(k, l);
      EXPECT_EQ(result.N, expected_gc_nodes(V, E, F, k, l));
      check_fullerene_dual(result);
    }
  }
}
