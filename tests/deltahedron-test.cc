#include <gtest/gtest.h>
#include <cmath>
#include "fullerenes/deltahedron.hh"
#include "fullerenes/polyhedron.hh"
#include "fullerenes/isomerdb.hh"
#include "C44geometries.hh"

using namespace std;

// Build an icosahedron Deltahedron from C20's dual
static Deltahedron make_icosahedron() {
  Polyhedron P20 = Polyhedron::C20();
  Polyhedron ico = P20.dual();  // icosahedron: 12 vertices, 20 triangular faces
  return Deltahedron(ico);
}

// Expected node count for GC(k,l) on a triangulation:
//   T = k^2+kl+l^2, F_new = F*T, E_new = 3*F_new/2, V_new = 2 + E_new - F_new
static int expected_gc_nodes(int V, int E, int F, int k, int l = 0) {
  int T = k*k + k*l + l*l;
  int F_new = F * T;
  int E_new = F_new * 3 / 2;
  return 2 + E_new - F_new;
}

class DeltahedronTest : public ::testing::Test {
protected:
  Deltahedron ico = make_icosahedron();
  // Icosahedron: V=12, E=30, F=20
  static constexpr int V0 = 12, E0 = 30, F0 = 20;
};

// ===== Construction tests =====

TEST_F(DeltahedronTest, ConstructionFromPolyhedron) {
  EXPECT_EQ(ico.N, V0);
  EXPECT_EQ((int)ico.points.size(), V0);
  EXPECT_EQ((int)ico.triangles.size(), F0);
}

TEST_F(DeltahedronTest, ConstructionFromTriangulationAndPoints) {
  Triangulation T(vector<int>(12, 5));  // C20 dual from spiral
  // Use icosahedron points (may have different node ordering, just check sizes)
  Deltahedron D(T, ico.points);
  EXPECT_EQ(D.N, V0);
  EXPECT_EQ((int)D.points.size(), V0);
}

TEST_F(DeltahedronTest, ComputeDualFaces) {
  auto faces = ico.compute_dual_faces();
  EXPECT_EQ((int)faces.size(), F0);
  for(const auto& f : faces)
    EXPECT_EQ((int)f.size(), 3);
}

// ===== GC Identity test =====

TEST_F(DeltahedronTest, GC_1_0_Identity) {
  Deltahedron result = ico.GCtransform(1, 0);
  EXPECT_EQ(result.N, ico.N);
  EXPECT_EQ((int)result.points.size(), ico.N);

  // Points should be identical (k=1: barycentric weights are trivial)
  for(int u = 0; u < ico.N; u++){
    for(int d = 0; d < 3; d++){
      EXPECT_NEAR(result.points[u][d], ico.points[u][d], 1e-12)
        << "vertex " << u << " coord " << d;
    }
  }
}

// ===== Node count tests =====

TEST_F(DeltahedronTest, GC_2_0_NodeCount) {
  Deltahedron result = ico.GCtransform(2, 0);
  EXPECT_EQ(result.N, expected_gc_nodes(V0, E0, F0, 2));
  EXPECT_EQ((int)result.points.size(), result.N);
}

TEST_F(DeltahedronTest, GC_3_0_NodeCount) {
  Deltahedron result = ico.GCtransform(3, 0);
  EXPECT_EQ(result.N, expected_gc_nodes(V0, E0, F0, 3));
  EXPECT_EQ((int)result.points.size(), result.N);
}

TEST_F(DeltahedronTest, GC_4_0_NodeCount) {
  Deltahedron result = ico.GCtransform(4, 0);
  EXPECT_EQ(result.N, expected_gc_nodes(V0, E0, F0, 4));
  EXPECT_EQ((int)result.points.size(), result.N);
}

// ===== Old vertex scaling test =====
// After GC(k,0), the original vertices (which appear at grid corners)
// should have coordinates equal to k * original, since the barycentric
// formula at a corner gives weight k to that corner and 0 to the others.

TEST_F(DeltahedronTest, GC_2_0_OldVertexScaling) {
  int k = 2;
  Deltahedron result = ico.GCtransform(k, 0);

  // Original vertices 0..N-1 are preserved with the same node IDs
  for(int u = 0; u < ico.N; u++){
    for(int d = 0; d < 3; d++){
      EXPECT_NEAR(result.points[u][d], k * ico.points[u][d], 1e-12)
        << "k=" << k << " vertex " << u << " coord " << d;
    }
  }
}

TEST_F(DeltahedronTest, GC_3_0_OldVertexScaling) {
  int k = 3;
  Deltahedron result = ico.GCtransform(k, 0);

  for(int u = 0; u < ico.N; u++){
    for(int d = 0; d < 3; d++){
      EXPECT_NEAR(result.points[u][d], k * ico.points[u][d], 1e-12)
        << "k=" << k << " vertex " << u << " coord " << d;
    }
  }
}

// ===== Surface preservation: new vertices lie on parent triangle plane =====

TEST_F(DeltahedronTest, GC_2_0_SurfacePreservation) {
  Deltahedron result = ico.GCtransform(2, 0);

  // Every new vertex should have zero distance to some parent triangle
  for(int u = ico.N; u < result.N; u++){
    double min_dist = INFINITY;
    for(const auto& tri : ico.triangles){
      Tri3D T(ico.points[tri[0]], ico.points[tri[1]], ico.points[tri[2]]);
      double d = T.distance(result.points[u]);
      if(d < min_dist) min_dist = d;
    }
    // The point lies on the surface scaled by k=2, so check against
    // the scaled triangles instead
    min_dist = INFINITY;
    for(const auto& tri : ico.triangles){
      Tri3D T(ico.points[tri[0]]*2, ico.points[tri[1]]*2, ico.points[tri[2]]*2);
      double d = T.distance(result.points[u]);
      if(d < min_dist) min_dist = d;
    }
    EXPECT_NEAR(min_dist, 0.0, 1e-10)
      << "new vertex " << u << " not on any scaled parent triangle plane";
  }
}

// ===== Topology consistency: GCtransform topology matches halma_transform =====

TEST_F(DeltahedronTest, GC_2_0_TopologyMatchesHalma) {
  Deltahedron result = ico.GCtransform(2, 0);
  Triangulation topo = ico.halma_transform(1);  // halma_transform(k-1)

  EXPECT_EQ(result.N, topo.N);
  EXPECT_EQ(result.triangles.size(), topo.triangles.size());
  for(int u = 0; u < result.N; u++){
    EXPECT_EQ(result.neighbours[u].size(), topo.neighbours[u].size())
      << "degree mismatch at node " << u;
  }
}

TEST_F(DeltahedronTest, GC_3_0_TopologyMatchesHalma) {
  Deltahedron result = ico.GCtransform(3, 0);
  Triangulation topo = ico.halma_transform(2);

  EXPECT_EQ(result.N, topo.N);
  EXPECT_EQ(result.triangles.size(), topo.triangles.size());
}

// ===== Face-interior harmonic property =====
// Within a single face, the coordinate map P(a,b) = (n-b)*P0 + (b-a)*P1 + a*P2
// is affine. The 6 neighbours of an interior grid point are symmetric, so their
// average equals the center point exactly. This does NOT hold for edge vertices
// (neighbours span two different affine maps across a crease) or corner vertices.

TEST_F(DeltahedronTest, GC_3_0_FaceInteriorHarmonicProperty) {
  int k = 3;
  int m = k - 1;
  vector<map<edge_t,node_t>> face_grids;
  static_cast<const Triangulation&>(ico).halma_transform(m, &face_grids);
  Deltahedron D = ico.GCtransform(k, 0);

  // Identify face-interior vertices: grid points (a,b) with a>0, b<m+1, a<b
  vector<bool> is_face_interior(D.N, false);
  for(int i = 0; i < (int)ico.triangles.size(); i++){
    for(const auto& [ab, node_id] : face_grids[i]){
      int a = ab.first, b = ab.second;
      if(a > 0 && b < m+1 && a < b)
        is_face_interior[node_id] = true;
    }
  }

  int n_tested = 0;
  for(int u = 0; u < D.N; u++){
    if(!is_face_interior[u]) continue;
    n_tested++;
    ASSERT_EQ((int)D.neighbours[u].size(), 6) << "face-interior vertex " << u << " should have degree 6";

    coord3d avg;
    for(int v : D.neighbours[u]) avg += D.points[v];
    avg /= 6.0;

    for(int d = 0; d < 3; d++){
      EXPECT_NEAR(avg[d], D.points[u][d], 1e-12)
        << "vertex " << u << " coord " << d;
    }
  }
  EXPECT_EQ(n_tested, 20) << "expected 20 face-interior vertices for GC(3,0) on icosahedron";
}

TEST_F(DeltahedronTest, GC_5_0_FaceInteriorHarmonicProperty) {
  int k = 5;
  int m = k - 1;
  vector<map<edge_t,node_t>> face_grids;
  static_cast<const Triangulation&>(ico).halma_transform(m, &face_grids);
  Deltahedron D = ico.GCtransform(k, 0);

  vector<bool> is_face_interior(D.N, false);
  for(int i = 0; i < (int)ico.triangles.size(); i++){
    for(const auto& [ab, node_id] : face_grids[i]){
      int a = ab.first, b = ab.second;
      if(a > 0 && b < m+1 && a < b)
        is_face_interior[node_id] = true;
    }
  }

  int n_tested = 0;
  for(int u = 0; u < D.N; u++){
    if(!is_face_interior[u]) continue;
    n_tested++;

    coord3d avg;
    for(int v : D.neighbours[u]) avg += D.points[v];
    avg /= D.neighbours[u].size();

    for(int d = 0; d < 3; d++){
      EXPECT_NEAR(avg[d], D.points[u][d], 1e-12)
        << "vertex " << u << " coord " << d;
    }
  }
  // (k-2)(k-1)/2 interior vertices per face * 20 faces = 3*4/2 * 20 = 120
  EXPECT_EQ(n_tested, 120);
}

// ===== No NaN or zero-norm coordinates =====

TEST_F(DeltahedronTest, GC_3_0_NoNaNPoints) {
  Deltahedron result = ico.GCtransform(3, 0);
  for(int u = 0; u < result.N; u++){
    for(int d = 0; d < 3; d++){
      EXPECT_FALSE(std::isnan(result.points[u][d]))
        << "NaN at vertex " << u << " coord " << d;
    }
    EXPECT_GT(result.points[u].norm(), 1e-10)
      << "zero-norm point at vertex " << u;
  }
}

// ===== Smoothing test =====

TEST_F(DeltahedronTest, SmoothDoesNotCrash) {
  Deltahedron D = ico.GCtransform(2, 0);
  // Just verify it runs without crashing and preserves sizes
  D.smooth(0.5);
  EXPECT_EQ((int)D.points.size(), D.N);
  for(int u = 0; u < D.N; u++){
    EXPECT_FALSE(std::isnan(D.points[u][0]));
  }
}

// ================================================================
// ===== General GC(k,l) tests (l != 0) ==========================
// ================================================================

// ===== Node count tests for chiral transforms =====

TEST_F(DeltahedronTest, GC_1_1_NodeCount) {
  Deltahedron result = ico.GCtransform(1, 1);
  EXPECT_EQ(result.N, expected_gc_nodes(V0, E0, F0, 1, 1));  // 32
  EXPECT_EQ((int)result.points.size(), result.N);
}

TEST_F(DeltahedronTest, GC_2_1_NodeCount) {
  Deltahedron result = ico.GCtransform(2, 1);
  EXPECT_EQ(result.N, expected_gc_nodes(V0, E0, F0, 2, 1));  // 72
  EXPECT_EQ((int)result.points.size(), result.N);
}

// ===== Original vertices present at sqrt(T)*P_original =====
// Node IDs are NOT preserved through unfold/fold, so we search by position.

TEST_F(DeltahedronTest, GC_1_1_OriginalVerticesPresent) {
  Deltahedron result = ico.GCtransform(1, 1);
  double sqrtT = sqrt(3.0);

  for(int u = 0; u < ico.N; u++){
    coord3d expected = ico.points[u] * sqrtT;
    bool found = false;
    for(int v = 0; v < result.N; v++){
      if((result.points[v] - expected).norm() < 1e-10){
        found = true;
        break;
      }
    }
    EXPECT_TRUE(found) << "original vertex " << u
      << " not found at sqrt(3)*P in GC(1,1) output";
  }
}

TEST_F(DeltahedronTest, GC_2_1_OriginalVerticesPresent) {
  Deltahedron result = ico.GCtransform(2, 1);
  double sqrtT = sqrt(7.0);

  for(int u = 0; u < ico.N; u++){
    coord3d expected = ico.points[u] * sqrtT;
    bool found = false;
    for(int v = 0; v < result.N; v++){
      if((result.points[v] - expected).norm() < 1e-10){
        found = true;
        break;
      }
    }
    EXPECT_TRUE(found) << "original vertex " << u
      << " not found at sqrt(7)*P in GC(2,1) output";
  }
}

// ===== Surface preservation: all vertices lie on scaled parent triangles =====

TEST_F(DeltahedronTest, GC_1_1_SurfacePreservation) {
  Deltahedron result = ico.GCtransform(1, 1);
  double sqrtT = sqrt(3.0);

  for(int u = 0; u < result.N; u++){
    double min_dist = INFINITY;
    for(const auto& tri : ico.triangles){
      Tri3D T(ico.points[tri[0]]*sqrtT, ico.points[tri[1]]*sqrtT, ico.points[tri[2]]*sqrtT);
      double d = T.distance(result.points[u]);
      if(d < min_dist) min_dist = d;
    }
    EXPECT_NEAR(min_dist, 0.0, 1e-10)
      << "vertex " << u << " not on any scaled parent triangle plane";
  }
}

// ===== Topology consistency: matches Triangulation::GCtransform =====

TEST_F(DeltahedronTest, GC_1_1_TopologyConsistency) {
  Deltahedron result = ico.GCtransform(1, 1);
  Triangulation topo = ico.Triangulation::GCtransform(1, 1);

  EXPECT_EQ(result.N, topo.N);
  EXPECT_EQ(result.triangles.size(), topo.triangles.size());
}

TEST_F(DeltahedronTest, GC_2_1_TopologyConsistency) {
  Deltahedron result = ico.GCtransform(2, 1);
  Triangulation topo = ico.Triangulation::GCtransform(2, 1);

  EXPECT_EQ(result.N, topo.N);
  EXPECT_EQ(result.triangles.size(), topo.triangles.size());
}

// ===== No NaN or zero-norm coordinates =====

TEST_F(DeltahedronTest, GC_1_1_NoNaNPoints) {
  Deltahedron result = ico.GCtransform(1, 1);
  for(int u = 0; u < result.N; u++){
    for(int d = 0; d < 3; d++){
      EXPECT_FALSE(std::isnan(result.points[u][d]))
        << "NaN at vertex " << u << " coord " << d;
    }
    EXPECT_GT(result.points[u].norm(), 1e-10)
      << "zero-norm point at vertex " << u;
  }
}

TEST_F(DeltahedronTest, GC_2_1_NoNaNPoints) {
  Deltahedron result = ico.GCtransform(2, 1);
  for(int u = 0; u < result.N; u++){
    for(int d = 0; d < 3; d++){
      EXPECT_FALSE(std::isnan(result.points[u][d]))
        << "NaN at vertex " << u << " coord " << d;
    }
    EXPECT_GT(result.points[u].norm(), 1e-10)
      << "zero-norm point at vertex " << u;
  }
}

// ================================================================
// ===== C44 isomer tests =========================================
// ================================================================

// Load C44 geometry database from binary files once
static C44Geometries load_C44() {
  return C44Geometries::load(
    C44_DATA_DIR "/c044all.rspi.u8",
    C44_DATA_DIR "/c044all.geometry.f32");
}

static const C44Geometries& C44_db() {
  static C44Geometries db = load_C44();
  return db;
}

// Target triangle edge length: sqrt(3) * C-C bond length (1.45 Angstrom)
static const double target_L = sqrt(3.0) * 1.45;

// Reconstruct a Deltahedron from stored C44 isomer data:
//   spiral -> Triangulation -> dual_graph -> FullereneGraph -> Polyhedron -> dual
// Coordinates are computed from the reconstructed graph (via zero_order_geometry)
// rather than loaded from stored data, because the spiral round-trip does not
// preserve the vertex labeling from buckygen.
static Deltahedron make_C44_deltahedron(int idx) {
  const auto& db = C44_db();
  const uint8_t* pi = db.pentagon_indices(idx);

  // Build spiral code: 24 faces (N/2+2 for C44), 12 pentagons + 12 hexagons
  vector<int> spiral(24, 6);
  for(int i = 0; i < 12; i++)
    spiral[pi[i]] = 5;

  // Reconstruct triangulation and fullerene graph from spiral
  Triangulation T(spiral);
  FullereneGraph G = T.dual_graph();
  G.layout2d = G.tutte_layout();
  // scalerad * 1.5 = average cubic edge length.
  // Want cubic edge ~ target_L / sqrt(3) = 1.45 Angstrom.
  double scalerad = target_L / (1.5 * sqrt(3.0));
  vector<coord3d> fpoints = G.zero_order_geometry(scalerad);

  // Build fullerene polyhedron, take dual to get deltahedron
  Polyhedron P(G, fpoints, 6);
  Polyhedron dual = P.dual();
  return Deltahedron(dual);
}

// C44 dual: V=24, E=66, F=44
static constexpr int C44_V = 24, C44_E = 66, C44_F = 44;

// ===== Node count tests across all C44 isomers =====

TEST(C44DeltahedronTest, AllC44_GC_2_0_NodeCount) {
  int expected = expected_gc_nodes(C44_V, C44_E, C44_F, 2, 0);
  for(int idx = 0; idx < C44_db().n_isomers; idx++){
    SCOPED_TRACE("isomer " + to_string(idx));
    Deltahedron D = make_C44_deltahedron(idx);
    Deltahedron result = D.GCtransform(2, 0);
    EXPECT_EQ(result.N, expected);
    EXPECT_EQ((int)result.points.size(), result.N);
  }
}

TEST(C44DeltahedronTest, AllC44_GC_1_1_NodeCount) {
  int expected = expected_gc_nodes(C44_V, C44_E, C44_F, 1, 1);
  for(int idx = 0; idx < C44_db().n_isomers; idx++){
    SCOPED_TRACE("isomer " + to_string(idx));
    Deltahedron D = make_C44_deltahedron(idx);
    Deltahedron result = D.GCtransform(1, 1);
    EXPECT_EQ(result.N, expected);
    EXPECT_EQ((int)result.points.size(), result.N);
  }
}

// ===== No NaN or zero-norm coordinates =====

TEST(C44DeltahedronTest, AllC44_GC_2_0_NoNaN) {
  for(int idx = 0; idx < C44_db().n_isomers; idx++){
    SCOPED_TRACE("isomer " + to_string(idx));
    Deltahedron D = make_C44_deltahedron(idx);
    Deltahedron result = D.GCtransform(2, 0);
    for(int u = 0; u < result.N; u++){
      for(int d = 0; d < 3; d++)
        EXPECT_FALSE(std::isnan(result.points[u][d]));
      EXPECT_GT(result.points[u].norm(), 1e-10);
    }
  }
}

TEST(C44DeltahedronTest, AllC44_GC_1_1_NoNaN) {
  for(int idx = 0; idx < C44_db().n_isomers; idx++){
    SCOPED_TRACE("isomer " + to_string(idx));
    Deltahedron D = make_C44_deltahedron(idx);
    Deltahedron result = D.GCtransform(1, 1);
    for(int u = 0; u < result.N; u++){
      for(int d = 0; d < 3; d++)
        EXPECT_FALSE(std::isnan(result.points[u][d]));
      EXPECT_GT(result.points[u].norm(), 1e-10);
    }
  }
}

// ===== Topology consistency =====

TEST(C44DeltahedronTest, AllC44_GC_2_0_TopologyConsistency) {
  for(int idx = 0; idx < C44_db().n_isomers; idx++){
    SCOPED_TRACE("isomer " + to_string(idx));
    Deltahedron D = make_C44_deltahedron(idx);
    Deltahedron result = D.GCtransform(2, 0);
    Triangulation topo = D.halma_transform(1);  // halma_transform(k-1)
    EXPECT_EQ(result.N, topo.N);
    EXPECT_EQ(result.triangles.size(), topo.triangles.size());
  }
}

TEST(C44DeltahedronTest, AllC44_GC_1_1_TopologyConsistency) {
  for(int idx = 0; idx < C44_db().n_isomers; idx++){
    SCOPED_TRACE("isomer " + to_string(idx));
    Deltahedron D = make_C44_deltahedron(idx);
    Deltahedron result = D.GCtransform(1, 1);
    Triangulation topo = D.Triangulation::GCtransform(1, 1);
    EXPECT_EQ(result.N, topo.N);
    EXPECT_EQ(result.triangles.size(), topo.triangles.size());
  }
}

// ===== Surface preservation =====

TEST(C44DeltahedronTest, AllC44_GC_2_0_SurfacePreservation) {
  for(int idx = 0; idx < C44_db().n_isomers; idx++){
    SCOPED_TRACE("isomer " + to_string(idx));
    Deltahedron D = make_C44_deltahedron(idx);
    Deltahedron result = D.GCtransform(2, 0);

    for(int u = 0; u < result.N; u++){
      double min_dist = INFINITY;
      for(const auto& tri : D.triangles){
        Tri3D T(D.points[tri[0]]*2, D.points[tri[1]]*2, D.points[tri[2]]*2);
        double d = T.distance(result.points[u]);
        if(d < min_dist) min_dist = d;
      }
      EXPECT_NEAR(min_dist, 0.0, 1e-10)
        << "vertex " << u << " not on any scaled parent triangle";
    }
  }
}

TEST(C44DeltahedronTest, AllC44_GC_1_1_SurfacePreservation) {
  double sqrtT = sqrt(3.0);
  for(int idx = 0; idx < C44_db().n_isomers; idx++){
    SCOPED_TRACE("isomer " + to_string(idx));
    Deltahedron D = make_C44_deltahedron(idx);
    Deltahedron result = D.GCtransform(1, 1);

    for(int u = 0; u < result.N; u++){
      double min_dist = INFINITY;
      for(const auto& tri : D.triangles){
        Tri3D T(D.points[tri[0]]*sqrtT, D.points[tri[1]]*sqrtT, D.points[tri[2]]*sqrtT);
        double d = T.distance(result.points[u]);
        if(d < min_dist) min_dist = d;
      }
      EXPECT_NEAR(min_dist, 0.0, 1e-10)
        << "vertex " << u << " not on any scaled parent triangle";
    }
  }
}

// ===== Original vertices present =====

TEST(C44DeltahedronTest, AllC44_GC_2_0_OriginalVerticesPresent) {
  for(int idx = 0; idx < C44_db().n_isomers; idx++){
    SCOPED_TRACE("isomer " + to_string(idx));
    Deltahedron D = make_C44_deltahedron(idx);
    Deltahedron result = D.GCtransform(2, 0);

    // Halma path preserves node IDs: result.points[u] == 2 * D.points[u]
    for(int u = 0; u < D.N; u++){
      for(int d = 0; d < 3; d++){
        EXPECT_NEAR(result.points[u][d], 2.0 * D.points[u][d], 1e-10);
      }
    }
  }
}

TEST(C44DeltahedronTest, AllC44_GC_1_1_OriginalVerticesPresent) {
  double sqrtT = sqrt(3.0);
  for(int idx = 0; idx < C44_db().n_isomers; idx++){
    SCOPED_TRACE("isomer " + to_string(idx));
    Deltahedron D = make_C44_deltahedron(idx);
    Deltahedron result = D.GCtransform(1, 1);

    // Unfold/fold renumbers nodes, so search by position
    for(int u = 0; u < D.N; u++){
      coord3d expected = D.points[u] * sqrtT;
      bool found = false;
      for(int v = 0; v < result.N; v++){
        if((result.points[v] - expected).norm() < 1e-10){
          found = true;
          break;
        }
      }
      EXPECT_TRUE(found) << "original vertex " << u
        << " not found at sqrt(3)*P in GC(1,1) output";
    }
  }
}

// ================================================================
// ===== Larger GC transforms: GC(3,2) and GC(5,3) ===============
// ================================================================

// Helper: run all standard checks on GC(k,l) across all C44 isomers.
// Returns number of failures for diagnostics.
static void check_C44_GC(int k, int l, const char* label) {
  double scale = sqrt((double)(k*k + k*l + l*l));
  int expected_N = expected_gc_nodes(C44_V, C44_E, C44_F, k, l);

  for(int idx = 0; idx < C44_db().n_isomers; idx++){
    SCOPED_TRACE(string(label) + " isomer " + to_string(idx));
    Deltahedron D = make_C44_deltahedron(idx);
    Deltahedron result = D.GCtransform(k, l);

    // Node count
    EXPECT_EQ(result.N, expected_N);
    EXPECT_EQ((int)result.points.size(), result.N);

    // No NaN or zero-norm
    for(int u = 0; u < result.N; u++){
      for(int d = 0; d < 3; d++)
        EXPECT_FALSE(std::isnan(result.points[u][d]));
      EXPECT_GT(result.points[u].norm(), 1e-10);
    }

    // Topology consistency with Triangulation::GCtransform
    Triangulation topo = D.Triangulation::GCtransform(k, l);
    EXPECT_EQ(result.N, topo.N);
    EXPECT_EQ(result.triangles.size(), topo.triangles.size());

    // Surface preservation: every vertex lies on a scaled parent triangle
    for(int u = 0; u < result.N; u++){
      double min_dist = INFINITY;
      for(const auto& tri : D.triangles){
        Tri3D T(D.points[tri[0]]*scale, D.points[tri[1]]*scale, D.points[tri[2]]*scale);
        double d = T.distance(result.points[u]);
        if(d < min_dist) min_dist = d;
      }
      EXPECT_NEAR(min_dist, 0.0, 1e-10)
        << "vertex " << u << " not on any scaled parent triangle";
    }

    // Original vertices present (search by position)
    for(int u = 0; u < D.N; u++){
      coord3d expected = D.points[u] * scale;
      bool found = false;
      for(int v = 0; v < result.N; v++){
        if((result.points[v] - expected).norm() < 1e-10){
          found = true;
          break;
        }
      }
      EXPECT_TRUE(found) << "original vertex " << u
        << " not found at scale*P in " << label << " output";
    }
  }
}

// GC(3,2): T=19, scale=sqrt(19)
TEST(C44DeltahedronTest, AllC44_GC_3_2) {
  check_C44_GC(3, 2, "GC(3,2)");
}

// GC(5,3): T=49, scale=7
TEST(C44DeltahedronTest, AllC44_GC_5_3) {
  check_C44_GC(5, 3, "GC(5,3)");
}

// GC(10,10): T=300, scale=sqrt(300)=10*sqrt(3)
TEST(C44DeltahedronTest, AllC44_GC_10_10) {
  check_C44_GC(10, 10, "GC(10,10)");
}

// GC(10,11): T=331, scale=sqrt(331)
TEST(C44DeltahedronTest, AllC44_GC_10_11) {
  check_C44_GC(10, 11, "GC(10,11)");
}

// ================================================================
// ===== Optimization tests =======================================
// ================================================================

// Helper: compute edge length statistics (mean, variance, max deviation)
struct EdgeStats {
  double mean, variance, max_dev;
  static EdgeStats compute(const Deltahedron& D) {
    EdgeStats s;
    vector<edge_t> edges = D.undirected_edges();
    vector<double> lengths(edges.size());
    for(int i = 0; i < (int)edges.size(); i++)
      lengths[i] = coord3d::dist(D.points[edges[i].first],
                                  D.points[edges[i].second]);
    s.mean = 0;
    for(double l : lengths) s.mean += l;
    s.mean /= lengths.size();

    s.variance = 0;
    s.max_dev = 0;
    for(double l : lengths){
      double dev = l - s.mean;
      s.variance += dev * dev;
      s.max_dev = max(s.max_dev, fabs(dev));
    }
    s.variance /= lengths.size();
    return s;
  }
};

// Comprehensive optimization quality statistics.
struct OptStats {
  // Edge lengths
  double L_mean, L_min, L_max, L_cv;  // cv = coefficient of variation (std/mean)
  // Triangle angles (in degrees)
  double ang_mean, ang_min, ang_max, ang_std;
  // Convexity
  double h_min;
  // Convergence
  bool converged;
  int iters;

  static OptStats compute(const Deltahedron& D, bool conv) {
    OptStats s;
    s.converged = conv;
    s.iters = D.iterations_used;

    // Edge lengths
    vector<edge_t> edges = D.undirected_edges();
    double sum = 0, sum2 = 0;
    s.L_min = INFINITY; s.L_max = 0;
    for(const auto& e : edges){
      double l = coord3d::dist(D.points[e.first], D.points[e.second]);
      sum += l; sum2 += l*l;
      s.L_min = min(s.L_min, l);
      s.L_max = max(s.L_max, l);
    }
    int ne = (int)edges.size();
    s.L_mean = sum / ne;
    double L_var = sum2/ne - s.L_mean*s.L_mean;
    s.L_cv = (s.L_mean > 0) ? sqrt(max(0.0, L_var)) / s.L_mean : 0;

    // Triangle angles
    double asum = 0, asum2 = 0;
    int na = 0;
    s.ang_min = 180; s.ang_max = 0;
    for(const auto& tri : D.triangles){
      for(int c = 0; c < 3; c++){
        coord3d va = D.points[tri[(c+1)%3]] - D.points[tri[c]];
        coord3d vb = D.points[tri[(c+2)%3]] - D.points[tri[c]];
        double ang = coord3d::angle(va, vb) * 180.0 / M_PI;
        asum += ang; asum2 += ang*ang;
        s.ang_min = min(s.ang_min, ang);
        s.ang_max = max(s.ang_max, ang);
        na++;
      }
    }
    s.ang_mean = asum / na;
    s.ang_std = sqrt(max(0.0, asum2/na - s.ang_mean*s.ang_mean));

    // h_min must be set by caller (min_convexity_height defined later in file)
    s.h_min = 0;
    return s;
  }

  static void print_header() {
    printf("\n%4s  %5s  %8s %8s %8s  %7s %7s %7s  %8s %s\n",
           "iso", "iter", "L_cv", "L_min", "L_max",
           "ang_min", "ang_max", "ang_std",
           "h_min", "");
    printf("----  -----  -------- -------- --------  ------- ------- -------  -------- ----\n");
  }

  void print_row(int idx) const {
    printf("%4d  %5d  %8.5f %8.4f %8.4f  %7.2f %7.2f %7.3f  %+8.4f %s\n",
           idx, iters, L_cv, L_min, L_max,
           ang_min, ang_max, ang_std,
           h_min,
           converged ? "" : " [NC]");
  }

  static void print_summary(const vector<OptStats>& all) {
    double worst_cv = 0, worst_ang_min = 180, worst_ang_max = 0;
    double worst_ang_std = 0, worst_h = INFINITY;
    int n_nc = 0, max_iters = 0, sum_iters = 0;
    for(const auto& s : all){
      worst_cv = max(worst_cv, s.L_cv);
      worst_ang_min = min(worst_ang_min, s.ang_min);
      worst_ang_max = max(worst_ang_max, s.ang_max);
      worst_ang_std = max(worst_ang_std, s.ang_std);
      worst_h = min(worst_h, s.h_min);
      max_iters = max(max_iters, s.iters);
      sum_iters += s.iters;
      if(!s.converged) n_nc++;
    }
    printf("---- summary (%d isomers, %d not converged, %d total iters, max %d) ----\n",
           (int)all.size(), n_nc, sum_iters, max_iters);
    printf("  worst L_cv=%.5f  ang_range=[%.2f, %.2f]  worst ang_std=%.3f  worst h_min=%+.4f\n",
           worst_cv, worst_ang_min, worst_ang_max, worst_ang_std, worst_h);
  }
};

// Compute minimum signed height across all qualifying (deg<=6, all-low-deg-neighbors) vertices.
// Positive = fully convex, negative = worst concavity depth.
static double min_convexity_height(const Deltahedron& D) {
  double min_h = INFINITY;
  for(int v = 0; v < D.N; v++){
    int d = (int)D.neighbours[v].size();
    if(d > 6) continue;
    bool all_low = true;
    for(int i = 0; i < d; i++)
      if((int)D.neighbours[D.neighbours[v][i]].size() > 6){ all_low = false; break; }
    if(!all_low) continue;

    coord3d centroid(0,0,0);
    for(int i = 0; i < d; i++) centroid += D.points[D.neighbours[v][i]];
    centroid /= (double)d;
    coord3d n_fan(0,0,0);
    for(int i = 0; i < d; i++){
      coord3d e1 = D.points[D.neighbours[v][i]] - D.points[v];
      coord3d e2 = D.points[D.neighbours[v][(i+1)%d]] - D.points[v];
      n_fan += e1.cross(e2);
    }
    double n_len = n_fan.norm();
    if(n_len < 1e-15) continue;
    double h = (D.points[v] - centroid).dot(n_fan / n_len);
    if(h < min_h) min_h = h;
  }
  return min_h;
}

// Check that all deg<=6 vertices (with all deg<=6 neighbors) are convex.
static void check_convexity(const Deltahedron& D, double tol = 1e-6) {
  for(int v = 0; v < D.N; v++){
    int d = (int)D.neighbours[v].size();
    if(d > 6) continue;

    // Skip if any neighbor has deg > 6
    bool all_low = true;
    for(int i = 0; i < d; i++)
      if((int)D.neighbours[D.neighbours[v][i]].size() > 6){ all_low = false; break; }
    if(!all_low) continue;

    // Neighbor centroid
    coord3d centroid(0,0,0);
    for(int i = 0; i < d; i++)
      centroid += D.points[D.neighbours[v][i]];
    centroid /= (double)d;

    // Average outward normal from triangle fan
    coord3d n_fan(0,0,0);
    for(int i = 0; i < d; i++){
      int ni  = D.neighbours[v][i];
      int ni1 = D.neighbours[v][(i+1) % d];
      coord3d e1 = D.points[ni]  - D.points[v];
      coord3d e2 = D.points[ni1] - D.points[v];
      n_fan += e1.cross(e2);
    }
    double n_len = n_fan.norm();
    if(n_len < 1e-15) continue;
    coord3d n_hat = n_fan / n_len;

    double h = (D.points[v] - centroid).dot(n_hat);
    EXPECT_GE(h, -tol) << "vertex " << v << " (deg " << d
      << ") is concave: h=" << h;
  }
}

// ===== FD gradient check =====

TEST_F(DeltahedronTest, GradientCheck_Icosahedron) {
  // Perturbed icosahedron ensures all gradient terms are exercised
  vector<coord3d> noisy(ico.points);
  srand(123);
  for(int u = 0; u < ico.N; u++)
    for(int d = 0; d < 3; d++)
      noisy[u][d] += 0.1 * (2.0 * rand() / RAND_MAX - 1.0);

  double err = ico.gradient_check(noisy);
  printf("FD gradient check: max relative error = %.2e\n", err);
  EXPECT_LT(err, 1e-5) << "Analytic gradient must match FD within 1e-5";
}

TEST_F(DeltahedronTest, GradientCheck_LargePerturbation) {
  // Large perturbation creates concavities, exercising E_conv gradient
  vector<coord3d> noisy(ico.points);
  srand(789);
  for(int u = 0; u < ico.N; u++)
    for(int d = 0; d < 3; d++)
      noisy[u][d] += 0.3 * (2.0 * rand() / RAND_MAX - 1.0);

  double err = ico.gradient_check(noisy);
  printf("FD gradient check (large perturbation): max relative error = %.2e\n", err);
  EXPECT_LT(err, 1e-5) << "Analytic gradient must match FD within 1e-5";
}

// ===== Hessian check: analytical vs FD for optimize_patch's Hessian =====
// All three energy terms (E_bond, E_angle, E_conv) have exact analytical Hessians.

TEST_F(DeltahedronTest, HessianCheck_Icosahedron) {
  // All vertices free, small perturbation
  vector<coord3d> noisy(ico.points);
  srand(456);
  for(int u = 0; u < ico.N; u++)
    for(int d = 0; d < 3; d++)
      noisy[u][d] += 0.05 * (2.0 * rand() / RAND_MAX - 1.0);

  vector<bool> free_mask(ico.N, true);
  double err = ico.hessian_check(noisy, free_mask, {}, 0, 1e-5, true);
  printf("Hessian check (icosahedron, all free): max relative error = %.2e\n", err);
  EXPECT_LT(err, 5e-3) << "Analytical Hessian must match FD within 5e-3";
}

TEST_F(DeltahedronTest, HessianCheck_PartiallyFixed) {
  // Fix half the vertices (boundary), free the other half
  vector<coord3d> noisy(ico.points);
  srand(789);
  for(int u = 0; u < ico.N; u++)
    for(int d = 0; d < 3; d++)
      noisy[u][d] += 0.05 * (2.0 * rand() / RAND_MAX - 1.0);

  vector<bool> free_mask(ico.N, false);
  for(int i = 0; i < ico.N; i += 2) free_mask[i] = true;
  double err = ico.hessian_check(noisy, free_mask, {}, 0, 1e-5, true);
  printf("Hessian check (icosahedron, partial free): max relative error = %.2e\n", err);
  EXPECT_LT(err, 5e-3) << "Analytical Hessian must match FD within 5e-3";
}

TEST_F(DeltahedronTest, HessianCheck_WithConcavity) {
  // Large perturbation to create concavities, exercising E_conv Hessian
  vector<coord3d> noisy(ico.points);
  srand(321);
  for(int u = 0; u < ico.N; u++)
    for(int d = 0; d < 3; d++)
      noisy[u][d] += 0.3 * (2.0 * rand() / RAND_MAX - 1.0);

  vector<bool> free_mask(ico.N, true);
  // Use tighter FD step for curved landscape
  double err = ico.hessian_check(noisy, free_mask, {}, 0, 1e-6, true);
  printf("Hessian check (icosahedron, with concavity): max relative error = %.2e\n", err);
  EXPECT_LT(err, 5e-3) << "Analytical Hessian must match FD within 5e-3";
}

// ===== Icosahedron: already equilateral, should converge immediately =====

TEST_F(DeltahedronTest, Optimize_Icosahedron_StaysFixed) {
  Deltahedron D = ico;  // copy
  bool converged = (OptResult::CONVERGED == D.optimize(ico.points));
  EXPECT_TRUE(converged);

  // Should be very close to original (up to possible uniform scale)
  EdgeStats stats = EdgeStats::compute(D);
  EXPECT_LT(stats.variance, 1e-16)
    << "icosahedron edge variance should be near-zero after optimization";
}

// ===== Perturbed icosahedron: should recover equilateral =====

TEST_F(DeltahedronTest, Optimize_PerturbedIcosahedron) {
  Deltahedron D = ico;

  // Add noise to the coordinates
  vector<coord3d> noisy(ico.points);
  srand(42);
  for(int u = 0; u < D.N; u++){
    for(int d = 0; d < 3; d++)
      noisy[u][d] += 0.05 * (2.0 * rand() / RAND_MAX - 1.0);
  }

  EdgeStats before = EdgeStats::compute(Deltahedron(D, noisy));
  EXPECT_GT(before.variance, 1e-6) << "perturbation should create variance";

  bool converged = (OptResult::CONVERGED == D.optimize(noisy));
  EXPECT_TRUE(converged);

  EdgeStats after = EdgeStats::compute(D);
  EXPECT_LT(after.variance, before.variance * 0.01)
    << "optimization should reduce edge length variance by 100x";
  EXPECT_LT(after.variance, 1e-6)
    << "optimized edge length variance should be very small";
}

// Compute per-vertex signed height for all vertices.
// Non-qualifying vertices (deg>6 or neighbor with deg>6) get h=NaN.
static vector<double> vertex_convexity_heights(const Deltahedron& D) {
  vector<double> heights(D.N, NAN);
  for(int v = 0; v < D.N; v++){
    int d = (int)D.neighbours[v].size();
    if(d > 6) continue;
    bool all_low = true;
    for(int i = 0; i < d; i++)
      if((int)D.neighbours[D.neighbours[v][i]].size() > 6){ all_low = false; break; }
    if(!all_low) continue;

    coord3d centroid(0,0,0);
    for(int i = 0; i < d; i++) centroid += D.points[D.neighbours[v][i]];
    centroid /= (double)d;
    coord3d n_fan(0,0,0);
    for(int i = 0; i < d; i++){
      coord3d e1 = D.points[D.neighbours[v][i]] - D.points[v];
      coord3d e2 = D.points[D.neighbours[v][(i+1)%d]] - D.points[v];
      n_fan += e1.cross(e2);
    }
    double n_len = n_fan.norm();
    if(n_len < 1e-15) continue;
    heights[v] = (D.points[v] - centroid).dot(n_fan / n_len);
  }
  return heights;
}

// Dump deltahedron geometry + per-vertex convexity to a file for visualization.
static void dump_geometry(const Deltahedron& D, const string& filename) {
  FILE* f = fopen(filename.c_str(), "w");
  if(!f) return;
  vector<double> h = vertex_convexity_heights(D);
  // Header
  fprintf(f, "%d %d\n", D.N, (int)D.triangles.size());
  // Vertices: x y z degree h
  for(int v = 0; v < D.N; v++){
    fprintf(f, "%.12f %.12f %.12f %d %.12f\n",
            D.points[v][0], D.points[v][1], D.points[v][2],
            (int)D.neighbours[v].size(), std::isnan(h[v]) ? 999.0 : h[v]);
  }
  // Triangles: v0 v1 v2
  for(const auto& tri : D.triangles)
    fprintf(f, "%d %d %d\n", tri[0], tri[1], tri[2]);
  fclose(f);
}

// ===== C44: rank isomers by difficulty (initial concavity) =====
// Run with --gtest_filter='C44DeltahedronTest.RankByDifficulty'
// to identify the hardest isomers for the optimization test.
TEST(C44DeltahedronTest, DISABLED_RankByDifficulty) {
  vector<pair<double,int>> ranked;
  for(int idx = 0; idx < C44_db().n_isomers; idx++){
    Deltahedron D = make_C44_deltahedron(idx);
    double h = min_convexity_height(D);
    ranked.push_back({h, idx});
  }
  sort(ranked.begin(), ranked.end());
  printf("\nC44 isomers ranked by initial min convexity height (worst first):\n");
  for(auto& [h, idx] : ranked)
    printf("  isomer %2d: min_h = %+.6f\n", idx, h);
}

// Verify that round-tripping through spiral gives consistent vertex labeling.
// Dump geometry of most concave C44 isomers for visual inspection.
TEST(C44DeltahedronTest, DISABLED_DumpConcaveGeometries) {
  static const vector<int> worst = {65, 44, 57, 30, 49};
  for(int idx : worst){
    Deltahedron D = make_C44_deltahedron(idx);

    {
      Polyhedron P(D, D.points);
      string fname = "/tmp/c44_iso" + to_string(idx) + "_before.mol2";
      Polyhedron::to_file(P, fname);
      printf("Wrote %s (min_h=%.4f)\n", fname.c_str(), min_convexity_height(D));
    }

    D.optimize(D.points, target_L);
    {
      Polyhedron P(D, D.points);
      string fname = "/tmp/c44_iso" + to_string(idx) + "_after.mol2";
      Polyhedron::to_file(P, fname);
      printf("Wrote %s (min_h=%.4f)\n", fname.c_str(), min_convexity_height(D));
    }
  }
}

// ===== C44 duals: optimize the ~20 hardest isomers =====

TEST(C44DeltahedronTest, Optimize_AllC44) {
  // Representative subset (includes hardest isomers 65, 80 from RankByDifficulty)
  static const vector<int> subset = {0, 65, 80};

  OptStats::print_header();
  vector<OptStats> all_stats;

  for(int idx : subset){
    SCOPED_TRACE("isomer " + to_string(idx));
    Deltahedron D = make_C44_deltahedron(idx);
    bool converged = (OptResult::CONVERGED == D.optimize(D.points, target_L));
    OptStats s = OptStats::compute(D, converged);
    s.h_min = min_convexity_height(D);
    s.print_row(idx);
    all_stats.push_back(s);

    EXPECT_LT(s.L_cv, 0.03) << "C44 #" << idx << " L_cv=" << s.L_cv;
    EXPECT_GT(s.h_min, -0.02) << "C44 #" << idx << " h_min=" << s.h_min;
    EXPECT_GT(s.ang_min, 45.0) << "C44 #" << idx << " ang_min=" << s.ang_min;
    EXPECT_LT(s.ang_max, 75.0) << "C44 #" << idx << " ang_max=" << s.ang_max;

    for(int u = 0; u < D.N; u++)
      for(int d = 0; d < 3; d++)
        EXPECT_FALSE(std::isnan(D.points[u][d]));
  }
  OptStats::print_summary(all_stats);
}

// ================================================================
// ===== Larger fullerene dual optimization tests =================
// ================================================================

// Helper: build Deltahedron from IsomerDB entry using zero_order_geometry (unoptimized).
// This gives a rough starting guess -- the whole point is to test whether
// the deltahedron optimizer can handle a poor initial geometry.
static Deltahedron make_deltahedron_from_db(int N, const IsomerDB::Entry& entry) {
  FullereneGraph G = IsomerDB::makeIsomer(N, entry);
  if(G.layout2d.empty())
    G.layout2d = G.tutte_layout();
  // scalerad * 1.5 = average cubic edge length.
  // Want cubic edge ~ target_L / sqrt(3) = 1.45 Angstrom.
  double scalerad = target_L / (1.5 * sqrt(3.0));
  vector<coord3d> points = G.zero_order_geometry(scalerad);
  Polyhedron P(G, points, 6);
  Polyhedron dual = P.dual();
  return Deltahedron(dual);
}

// Parametrized test for optimizing fullerene duals of a given size.
// Tests all (or up to max_isomers) isomers for the given N.
// Only prints per-isomer rows for non-converged or problematic isomers (verbose=true prints all).
static void test_optimize_fullerene_duals(int N, int max_isomers = -1, bool IPR = true,
                                          long long max_work = 0, bool verbose = false) {
  IsomerDB db = IsomerDB::readPDB(N, IPR);
  ASSERT_GT(db.entries.size(), 0u) << "No isomers found for C" << N;

  int n_test = (max_isomers > 0) ? min(max_isomers, (int)db.entries.size())
                                 : (int)db.entries.size();

  printf("\nC%d%s (%d isomers, V_dual=%d):\n",
         N, IPR ? " IPR" : "", n_test, N/2+2);
  OptStats::print_header();
  vector<OptStats> all_stats;
  vector<int> failed_indices;

  for(int idx = 0; idx < n_test; idx++){
    SCOPED_TRACE("C" + to_string(N) + (IPR ? " IPR" : "") + " isomer " + to_string(idx));
    Deltahedron D = make_deltahedron_from_db(N, db.entries[idx]);

    int V_dual = N/2 + 2;
    EXPECT_EQ(D.N, V_dual) << "dual vertex count should be N/2+2";

    bool converged = (OptResult::CONVERGED == D.optimize(D.points, target_L, 1e-10, {}, max_work));
    OptStats s = OptStats::compute(D, converged);
    s.h_min = min_convexity_height(D);

    bool problematic = !converged || s.h_min < -0.01 || s.L_cv > 0.01;
    if(verbose || problematic) s.print_row(idx);
    if(!converged) failed_indices.push_back(idx);
    all_stats.push_back(s);

    EXPECT_LT(s.L_cv, 0.03) << "C" << N << " #" << idx << " L_cv=" << s.L_cv;
    EXPECT_GT(s.h_min, -0.02) << "C" << N << " #" << idx << " h_min=" << s.h_min;
    EXPECT_GT(s.ang_min, 45.0) << "C" << N << " #" << idx << " ang_min=" << s.ang_min;
    EXPECT_LT(s.ang_max, 75.0) << "C" << N << " #" << idx << " ang_max=" << s.ang_max;

    for(int u = 0; u < D.N; u++)
      for(int d = 0; d < 3; d++)
        EXPECT_FALSE(std::isnan(D.points[u][d]));
  }
  OptStats::print_summary(all_stats);

  if(!failed_indices.empty()){
    printf("Non-converged isomer indices (%d): {", (int)failed_indices.size());
    for(int i = 0; i < (int)failed_indices.size(); i++)
      printf("%s%d", i ? ", " : "", failed_indices[i]);
    printf("}\n");
  }
}

// ================================================================
// ===== Optimizer test suite =====================================
// ================================================================
//
// Three tiers:
//   1. Optimize_SmallFullerenes: All isomers C20-C50 (total ~833)
//   2. Optimize_C60_Hard: 50 hardest C60 isomers (from k_conv sweep)
//   3. Optimize_C100_Sample: 50 random IPR C100 isomers (scaling check)
//
// Together these run in ~1 minute and cover small/medium/large sizes.

// Tier 1: Small fullerene duals.
// All isomers for C20-C34 (20 total), sampled for C36-C50 (10 per size).
// Uses max_iter=1000 since small isomers converge quickly.
TEST(OptimizeTest, SmallFullerenes) {
  // Fullerenes exist for all even N >= 20 except N=22
  static const int sizes[] = {20, 24, 26, 28, 30, 32, 34, 36, 38, 40, 42, 44, 46, 48, 50};
  static const int sample_per_size = 0;  // 0 = test all isomers
  int total = 0, total_nc = 0, total_concave = 0;

  for(int N : sizes){
    IsomerDB db = IsomerDB::readPDB(N, false);
    int n_all = (int)db.entries.size();
    if(n_all == 0) continue;

    // Test all isomers (or sample if sample_per_size > 0)
    vector<int> indices;
    if(sample_per_size <= 0 || n_all <= sample_per_size){
      for(int i = 0; i < n_all; i++) indices.push_back(i);
    } else {
      // Deterministic sample
      srand(N * 1000 + 42);
      vector<int> pool(n_all);
      for(int i = 0; i < n_all; i++) pool[i] = i;
      for(int i = 0; i < sample_per_size; i++){
        int j = i + rand() % (n_all - i);
        swap(pool[i], pool[j]);
      }
      indices.assign(pool.begin(), pool.begin() + sample_per_size);
      sort(indices.begin(), indices.end());
    }

    int n_test = (int)indices.size();
    vector<OptStats> stats;
    vector<int> nc;
    for(int idx : indices){
      Deltahedron D = make_deltahedron_from_db(N, db.entries[idx]);
      bool converged = (OptResult::CONVERGED == D.optimize(D.points, target_L));
      OptStats s = OptStats::compute(D, converged);
      s.h_min = min_convexity_height(D);
      stats.push_back(s);
      if(!converged) nc.push_back(idx);
    }

    // Find worst stats
    double worst_Lcv = 0, worst_hmin = INFINITY;
    double worst_ang_min = 180, worst_ang_max = 0, worst_ang_std = 0;
    int n_concave = 0;
    for(auto& s : stats){
      worst_Lcv = max(worst_Lcv, s.L_cv);
      worst_hmin = min(worst_hmin, s.h_min);
      worst_ang_min = min(worst_ang_min, s.ang_min);
      worst_ang_max = max(worst_ang_max, s.ang_max);
      worst_ang_std = max(worst_ang_std, s.ang_std);
      if(s.h_min < -0.01) n_concave++;
    }

    printf("  C%-3d: %3d/%3d  %2d NC  L_cv=%.5f  h_min=%+.4f  ang=[%.2f,%.2f] std=%.3f  %d concave\n",
           N, n_test, n_all, (int)nc.size(), worst_Lcv, worst_hmin,
           worst_ang_min, worst_ang_max, worst_ang_std, n_concave);

    total += n_test;
    total_nc += (int)nc.size();
    total_concave += n_concave;

    // Geometry quality assertions: both triangle quality and convexity
    for(int i = 0; i < n_test; i++){
      EXPECT_LT(stats[i].L_cv, 0.03)
        << "C" << N << " isomer " << indices[i] << " L_cv=" << stats[i].L_cv;
      EXPECT_GT(stats[i].h_min, -0.02)
        << "C" << N << " isomer " << indices[i] << " h_min=" << stats[i].h_min;
      EXPECT_GT(stats[i].ang_min, 45.0)
        << "C" << N << " isomer " << indices[i] << " ang_min=" << stats[i].ang_min;
      EXPECT_LT(stats[i].ang_max, 75.0)
        << "C" << N << " isomer " << indices[i] << " ang_max=" << stats[i].ang_max;
    }
  }
  printf("  Total: %d isomers, %d NC, %d concave\n", total, total_nc, total_concave);
}

// Tier 2: 50 hardest C60 isomers.
// Selected from full C60 sweep: the 49 non-converged at k_conv=2 (the optimal
// k_conv value), plus isomer 388 (deepest concave local minimum).
TEST(OptimizeTest, C60_Hard) {
  static const vector<int> hard = {
    11, 118, 119, 146, 175, 248, 287, 296, 297, 298, 300, 301, 327,
    353, 388, 394, 395, 458, 503, 529, 536, 539, 541, 542, 545, 555,
    558, 568, 570, 618, 791, 798, 943, 986, 1044, 1093, 1099, 1136,
    1156, 1197, 1246, 1281, 1343, 1373, 1405, 1428, 1511, 1676, 1687,
    1783
  };

  IsomerDB db = IsomerDB::readPDB(60, false);
  ASSERT_GT((int)db.entries.size(), hard.back());

  printf("\nC60 hard subset (%d isomers, V_dual=32):\n", (int)hard.size());
  OptStats::print_header();
  vector<OptStats> all_stats;
  vector<int> nc;

  for(int idx : hard){
    Deltahedron D = make_deltahedron_from_db(60, db.entries[idx]);
    bool converged = (OptResult::CONVERGED == D.optimize(D.points, target_L));
    OptStats s = OptStats::compute(D, converged);
    s.h_min = min_convexity_height(D);

    bool problematic = !converged || s.h_min < -0.01 || s.L_cv > 0.01;
    if(problematic) s.print_row(idx);
    if(!converged) nc.push_back(idx);
    all_stats.push_back(s);
  }
  OptStats::print_summary(all_stats);

  int n_convex = 0, n_concave = 0, n_stuck = 0;
  for(auto& s : all_stats){
    if(s.L_cv > 0.001) n_stuck++;
    else if(s.h_min < -0.01) n_concave++;
    else n_convex++;
  }
  printf("NC=%d, convex=%d, concave=%d, stuck=%d\n",
         (int)nc.size(), n_convex, n_concave, n_stuck);

  // Quality assertions: all isomers should have acceptable geometry
  for(int i = 0; i < (int)all_stats.size(); i++){
    EXPECT_LT(all_stats[i].L_cv, 0.03)
      << "C60 #" << hard[i] << " L_cv=" << all_stats[i].L_cv;
    EXPECT_GT(all_stats[i].h_min, -0.02)
      << "C60 #" << hard[i] << " h_min=" << all_stats[i].h_min;
    EXPECT_GT(all_stats[i].ang_min, 45.0)
      << "C60 #" << hard[i] << " ang_min=" << all_stats[i].ang_min;
    EXPECT_LT(all_stats[i].ang_max, 75.0)
      << "C60 #" << hard[i] << " ang_max=" << all_stats[i].ang_max;
  }
}

// Tier 2b: IPR C60-Ih pentakis dodecahedron validation.
TEST(OptimizeTest, C60_IPR) {
  IsomerDB db = IsomerDB::readPDB(60, false);
  ASSERT_GT((int)db.entries.size(), 1811);

  Deltahedron D = make_deltahedron_from_db(60, db.entries[1811]);
  bool converged = (OptResult::CONVERGED == D.optimize(D.points, target_L));
  OptStats s = OptStats::compute(D, converged);
  s.h_min = min_convexity_height(D);
  s.print_row(1811);

  EXPECT_TRUE(converged) << "IPR C60 should converge";
  EXPECT_LT(s.L_cv, 1e-6) << "IPR C60 should have uniform edge lengths";
  EXPECT_LT(s.ang_std, 0.01) << "IPR C60 should have all angles = 60 deg";
  EXPECT_GT(s.h_min, 0.2) << "IPR C60 should be strongly convex";

  // Check pentakis dodecahedron geometry: r5/r6 ratio
  coord3d cm;
  for(int v = 0; v < D.N; v++) cm += D.points[v];
  cm /= (double)D.N;

  double r5_sum = 0, r6_sum = 0;
  int n5 = 0, n6 = 0;
  for(int v = 0; v < D.N; v++){
    double r = (D.points[v] - cm).norm();
    if((int)D.neighbours[v].size() == 5){ r5_sum += r; n5++; }
    else                                 { r6_sum += r; n6++; }
  }
  ASSERT_EQ(n5, 12);
  ASSERT_EQ(n6, 20);

  double ratio = (r5_sum/n5) / (r6_sum/n6);
  double exact_ratio = 1.169839420;
  EXPECT_NEAR(ratio, exact_ratio, 1e-4)
    << "r5/r6 = " << ratio << ", expected " << exact_ratio;

  printf("  IPR C60: r5/r6 = %.6f (exact: %.6f)\n", ratio, exact_ratio);
}

// Tier 3: 50 random IPR C100 isomers (deterministic seed).
TEST(OptimizeTest, C100_Sample) {
  IsomerDB db = IsomerDB::readPDB(100, true);  // 450 IPR isomers
  int n_total = (int)db.entries.size();
  ASSERT_GT(n_total, 0);

  // Deterministic sample of 50 from 450
  const int n_sample = min(50, n_total);
  vector<int> sample_indices;
  srand(42);
  if(n_sample >= n_total){
    for(int i = 0; i < n_total; i++) sample_indices.push_back(i);
  } else {
    // Fisher-Yates partial shuffle
    vector<int> pool(n_total);
    for(int i = 0; i < n_total; i++) pool[i] = i;
    for(int i = 0; i < n_sample; i++){
      int j = i + rand() % (n_total - i);
      swap(pool[i], pool[j]);
    }
    sample_indices.assign(pool.begin(), pool.begin() + n_sample);
    sort(sample_indices.begin(), sample_indices.end());
  }

  printf("\nC100 IPR sample (%d of %d isomers, V_dual=52):\n", n_sample, n_total);
  OptStats::print_header();
  vector<OptStats> all_stats;
  vector<int> nc;

  for(int idx : sample_indices){
    Deltahedron D = make_deltahedron_from_db(100, db.entries[idx]);
    bool converged = (OptResult::CONVERGED == D.optimize(D.points, target_L));
    OptStats s = OptStats::compute(D, converged);
    s.h_min = min_convexity_height(D);

    bool problematic = !converged || s.h_min < -0.01 || s.L_cv > 0.01;
    if(problematic) s.print_row(idx);
    if(!converged) nc.push_back(idx);
    all_stats.push_back(s);
  }
  OptStats::print_summary(all_stats);
  printf("NC=%d of %d\n", (int)nc.size(), n_sample);

  // Quality assertions
  for(int i = 0; i < (int)all_stats.size(); i++){
    EXPECT_LT(all_stats[i].L_cv, 0.03)
      << "C100 #" << sample_indices[i] << " L_cv=" << all_stats[i].L_cv;
    EXPECT_GT(all_stats[i].h_min, -0.02)
      << "C100 #" << sample_indices[i] << " h_min=" << all_stats[i].h_min;
    EXPECT_GT(all_stats[i].ang_min, 45.0)
      << "C100 #" << sample_indices[i] << " ang_min=" << all_stats[i].ang_min;
    EXPECT_LT(all_stats[i].ang_max, 75.0)
      << "C100 #" << sample_indices[i] << " ang_max=" << all_stats[i].ang_max;
  }
}

// ================================================================
// ===== Full sweeps (DISABLED — too slow for regular use) ========
// ================================================================


TEST(SweepTest, DISABLED_Optimize_AllC60) {
  test_optimize_fullerene_duals(60, -1, false);
}

TEST(SweepTest, DISABLED_Optimize_AllC70) {
  test_optimize_fullerene_duals(70, -1, false);
}
