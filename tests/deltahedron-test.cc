#include <gtest/gtest.h>
#include <cmath>
#include "fullerenes/deltahedron.hh"
#include "fullerenes/polyhedron.hh"
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
  ico.halma_transform(m, &face_grids);
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
  ico.halma_transform(m, &face_grids);
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

// Reconstruct a Deltahedron from stored C44 isomer data:
//   spiral -> Triangulation -> dual_graph -> Polyhedron -> dual -> Deltahedron
static Deltahedron make_C44_deltahedron(int idx) {
  const auto& db = C44_db();
  const uint8_t* pi = db.pentagon_indices(idx);
  const float* pts = db.points(idx);

  // Build spiral code: 24 faces (N/2+2 for C44), 12 pentagons + 12 hexagons
  vector<int> spiral(24, 6);
  for(int i = 0; i < 12; i++)
    spiral[pi[i]] = 5;

  // Reconstruct triangulation (fullerene dual) from spiral
  Triangulation T(spiral);

  // Get fullerene graph and attach stored optimized coordinates
  FullereneGraph G = T.dual_graph();
  vector<coord3d> fpoints(44);
  for(int v = 0; v < 44; v++)
    fpoints[v] = coord3d(pts[v*3], pts[v*3+1], pts[v*3+2]);

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
