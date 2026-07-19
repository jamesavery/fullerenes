// PLY mesh I/O round-trip: Polyhedron::{from,to}_ply and Deltahedron::{from,to}_ply
// preserve geometry and connectivity through write -> read, for both binary and
// ascii and for triangle and n-gon faces. Exercises std::tmpfile (no on-disk
// fixtures) and the outward-orientation normalisation in from_ply.

#include "fullerenes/polyhedron.hh"
#include "fullerenes/deltahedron.hh"

#include <gtest/gtest.h>
#include <algorithm>
#include <array>
#include <cmath>
#include <cstdio>
#include <vector>

namespace {

// Signed enclosed volume * 6 (fan triangulation): sign = orientation handedness,
// positive for outward-facing (CCW-on-outside) faces. Independent re-derivation
// of from_ply's internal normalisation, so the test can actually check it.
double signed_volume6(const Polyhedron& P) {
  double V = 0;
  for (const auto& f : P.faces())
    for (size_t i = 1; i + 1 < f.size(); i++)
      V += P.points[f[0]].dot(P.points[f[i]].cross(P.points[f[i+1]]));
  return V;
}

// Recovered triangle set as sorted vertex triples, for order-insensitive content
// comparison.
std::vector<std::array<int,3>> sorted_triangles(const std::vector<tri_t>& tris) {
  std::vector<std::array<int,3>> out;
  for (const auto& t : tris) {
    std::array<int,3> a{t[0], t[1], t[2]};
    std::sort(a.begin(), a.end());
    out.push_back(a);
  }
  std::sort(out.begin(), out.end());
  return out;
}

// Round-trip a mesh through Polyhedron::{to,from}_ply and assert the recovered
// polyhedron has the same vertex count, face count, and (within tolerance)
// vertex positions. ASSERT_EQ on N first so a vertex-count regression stops the
// test before the position loop reads out of bounds.
void check_polyhedron_roundtrip(const Polyhedron& P, bool binary) {
  FILE* tmp = std::tmpfile();
  ASSERT_NE(tmp, nullptr);
  ASSERT_TRUE(Polyhedron::to_ply(P, tmp, binary));
  std::rewind(tmp);
  Polyhedron Q = Polyhedron::from_ply(tmp);
  std::fclose(tmp);

  ASSERT_EQ(Q.N, P.N) << (binary ? "binary" : "ascii");
  EXPECT_EQ(Q.faces().size(), P.faces().size());
  double max_err = 0;
  for (int i = 0; i < P.N; i++)
    for (int d = 0; d < 3; d++)
      max_err = std::max(max_err, std::fabs(Q.points[i][d] - P.points[i][d]));
  EXPECT_LT(max_err, 1e-6) << (binary ? "binary" : "ascii");
}

} // namespace

// Build a Polyhedron straight from a hand-written PLY (covers from_ply on n-gon
// faces), then round-trip it both ways.
TEST(PlyIO, CubeRoundTrip) {
  const char* cube_ply =
    "ply\nformat ascii 1.0\n"
    "element vertex 8\n"
    "property float x\nproperty float y\nproperty float z\n"
    "element face 6\n"
    "property list uchar int vertex_indices\nend_header\n"
    "0 0 0\n1 0 0\n1 1 0\n0 1 0\n0 0 1\n1 0 1\n1 1 1\n0 1 1\n"
    "4 0 3 2 1\n4 4 5 6 7\n4 0 1 5 4\n4 1 2 6 5\n4 2 3 7 6\n4 3 0 4 7\n";
  FILE* tmp = std::tmpfile();
  ASSERT_NE(tmp, nullptr);
  std::fputs(cube_ply, tmp);
  std::rewind(tmp);
  Polyhedron cube = Polyhedron::from_ply(tmp);
  std::fclose(tmp);

  EXPECT_EQ(cube.N, 8);
  EXPECT_EQ(cube.faces().size(), 6u);          // n-gon (quad) faces preserved
  EXPECT_GT(cube.volume(), 0.9);               // unit cube, ~1.0
  check_polyhedron_roundtrip(cube, /*binary=*/true);
  check_polyhedron_roundtrip(cube, /*binary=*/false);
}

// from_ply normalises to outward: an inward-wound cube comes back with the same
// (positive) volume and outward faces.
TEST(PlyIO, NormalisesToOutward) {
  const char* inward_cube =
    "ply\nformat ascii 1.0\n"
    "element vertex 8\n"
    "property float x\nproperty float y\nproperty float z\n"
    "element face 6\n"
    "property list uchar int vertex_indices\nend_header\n"
    "0 0 0\n1 0 0\n1 1 0\n0 1 0\n0 0 1\n1 0 1\n1 1 1\n0 1 1\n"
    "4 1 2 3 0\n4 7 6 5 4\n4 4 5 1 0\n4 5 6 2 1\n4 6 7 3 2\n4 7 4 0 3\n";  // reversed winding
  FILE* tmp = std::tmpfile();
  ASSERT_NE(tmp, nullptr);
  std::fputs(inward_cube, tmp);
  std::rewind(tmp);
  Polyhedron cube = Polyhedron::from_ply(tmp);
  std::fclose(tmp);
  EXPECT_TRUE(cube.is_consistently_oriented());
  EXPECT_GT(cube.volume(), 0.9);               // unit cube, ~1.0 (volume() is fabs)
  // The real test of outward normalisation: the SIGNED volume must be positive.
  // (volume()/is_consistently_oriented() above cannot distinguish inward winding.)
  EXPECT_GT(signed_volume6(cube), 0.0);
}

// A tetrahedron is a triangulation: Deltahedron::{from,to}_ply round-trips it.
TEST(PlyIO, DeltahedronRoundTrip) {
  const char* tet_ply =
    "ply\nformat ascii 1.0\n"
    "element vertex 4\n"
    "property float x\nproperty float y\nproperty float z\n"
    "element face 4\n"
    "property list uchar int vertex_indices\nend_header\n"
    "0 0 0\n1 0 0\n0 1 0\n0 0 1\n"
    "3 0 2 1\n3 0 1 3\n3 1 2 3\n3 2 0 3\n";
  FILE* tmp = std::tmpfile();
  ASSERT_NE(tmp, nullptr);
  std::fputs(tet_ply, tmp);
  std::rewind(tmp);
  Deltahedron D = Deltahedron::from_ply(tmp);
  std::fclose(tmp);
  EXPECT_EQ(D.N, 4);
  EXPECT_EQ((int)D.triangles().size(), 4);

  FILE* tmp2 = std::tmpfile();
  ASSERT_NE(tmp2, nullptr);
  ASSERT_TRUE(Deltahedron::to_ply(D, tmp2, /*binary=*/true));
  std::rewind(tmp2);
  Deltahedron E = Deltahedron::from_ply(tmp2);
  std::fclose(tmp2);
  ASSERT_EQ(E.N, D.N);
  // Same triangle CONTENT, not just count -- a winding/order regression in the
  // write->read would otherwise pass on equal sizes alone.
  EXPECT_EQ(sorted_triangles(E.triangles()), sorted_triangles(D.triangles()));
}

// from_ply welds duplicate coincident vertices (a zero-length edge and its two
// zero-area faces -- the tomography-mesh defect). This "doubled apex" is a valid
// closed manifold (5 verts, 6 faces, consistently wound) where apex vertex 4 is
// an exact copy of vertex 3 at (0,0,1), joined by the zero-length edge (3,4);
// welding merges 4 into 3, drops the two degenerate bridge faces, and leaves the
// clean tetrahedron {0,1,2,3}.
TEST(PlyIO, WeldsCoincidentVertices) {
  const char* doubled_apex =
    "ply\nformat ascii 1.0\n"
    "element vertex 5\n"
    "property float x\nproperty float y\nproperty float z\n"
    "element face 6\n"
    "property list uchar int vertex_indices\nend_header\n"
    "1 0 0\n-0.5 0.866 0\n-0.5 -0.866 0\n0 0 1\n0 0 1\n"       // vertex 4 == vertex 3
    "3 0 2 1\n3 0 1 3\n3 1 2 4\n3 2 0 4\n3 3 1 4\n3 0 3 4\n";  // (3,4) is zero-length
  FILE* tmp = std::tmpfile();
  ASSERT_NE(tmp, nullptr);
  std::fputs(doubled_apex, tmp);
  std::rewind(tmp);
  Deltahedron D = Deltahedron::from_ply(tmp);
  std::fclose(tmp);

  // The coincident apex is gone: a clean 4-vertex, 4-face tetrahedron.
  EXPECT_EQ(D.N, 4);
  EXPECT_EQ((int)D.triangles().size(), 4);
  // No zero-length edge survives (the defect the weld exists to remove): every
  // triangle side has positive 3D length.
  double min_len = 1e300;
  for (const tri_t& t : D.triangles())
    for (int e = 0; e < 3; e++)
      min_len = std::min(min_len, (D.points[t[(e + 1) % 3]] - D.points[t[e]]).norm());
  EXPECT_GT(min_len, 1e-9);
}

// A clean mesh (no coincident vertices) is untouched by the weld -- the vertex
// count and faces are exactly as read.
TEST(PlyIO, WeldLeavesCleanMeshUnchanged) {
  const char* tet_ply =
    "ply\nformat ascii 1.0\n"
    "element vertex 4\n"
    "property float x\nproperty float y\nproperty float z\n"
    "element face 4\n"
    "property list uchar int vertex_indices\nend_header\n"
    "0 0 0\n1 0 0\n0 1 0\n0 0 1\n"
    "3 0 2 1\n3 0 1 3\n3 1 2 3\n3 2 0 3\n";
  FILE* tmp = std::tmpfile();
  ASSERT_NE(tmp, nullptr);
  std::fputs(tet_ply, tmp);
  std::rewind(tmp);
  Deltahedron D = Deltahedron::from_ply(tmp);
  std::fclose(tmp);
  EXPECT_EQ(D.N, 4);
  EXPECT_EQ((int)D.triangles().size(), 4);
}

// A non-triangulation PLY is rejected by Deltahedron::from_ply with mesh_io_error.
TEST(PlyIO, DeltahedronRejectsNonTriangulation) {
  const char* cube_ply =
    "ply\nformat ascii 1.0\n"
    "element vertex 8\n"
    "property float x\nproperty float y\nproperty float z\n"
    "element face 6\n"
    "property list uchar int vertex_indices\nend_header\n"
    "0 0 0\n1 0 0\n1 1 0\n0 1 0\n0 0 1\n1 0 1\n1 1 1\n0 1 1\n"
    "4 0 3 2 1\n4 4 5 6 7\n4 0 1 5 4\n4 1 2 6 5\n4 2 3 7 6\n4 3 0 4 7\n";
  FILE* tmp = std::tmpfile();
  ASSERT_NE(tmp, nullptr);
  std::fputs(cube_ply, tmp);
  std::rewind(tmp);
  EXPECT_THROW(Deltahedron::from_ply(tmp), mesh_io_error);
  std::fclose(tmp);
}

// A malformed header throws mesh_io_error rather than aborting.
TEST(PlyIO, MalformedThrows) {
  FILE* tmp = std::tmpfile();
  ASSERT_NE(tmp, nullptr);
  std::fputs("ply\nformat ascii 1.0\nelement vertex 3\nend_header\n", tmp);  // no x/y/z
  std::rewind(tmp);
  EXPECT_THROW(Polyhedron::from_ply(tmp), mesh_io_error);
  std::fclose(tmp);
}

// An out-of-range face index is rejected (was an out-of-bounds heap write).
TEST(PlyIO, OutOfRangeFaceIndexRejected) {
  const char* bad =
    "ply\nformat ascii 1.0\n"
    "element vertex 4\nproperty float x\nproperty float y\nproperty float z\n"
    "element face 1\nproperty list uchar int vertex_indices\nend_header\n"
    "0 0 0\n1 0 0\n0 1 0\n0 0 1\n"
    "3 0 1 99\n";                              // index 99 >= nv=4
  FILE* tmp = std::tmpfile();
  ASSERT_NE(tmp, nullptr);
  std::fputs(bad, tmp);
  std::rewind(tmp);
  EXPECT_THROW(Polyhedron::from_ply(tmp), mesh_io_error);
  std::fclose(tmp);
}

// A directed edge shared by two faces (non-manifold / inconsistent winding) is
// rejected rather than silently building a mis-wired graph.
TEST(PlyIO, NonManifoldRejected) {
  const char* bad =
    "ply\nformat ascii 1.0\n"
    "element vertex 4\nproperty float x\nproperty float y\nproperty float z\n"
    "element face 2\nproperty list uchar int vertex_indices\nend_header\n"
    "0 0 0\n1 0 0\n0 1 0\n0 0 1\n"
    "3 0 1 2\n3 0 1 3\n";                      // both traverse arc 0->1
  FILE* tmp = std::tmpfile();
  ASSERT_NE(tmp, nullptr);
  std::fputs(bad, tmp);
  std::rewind(tmp);
  EXPECT_THROW(Polyhedron::from_ply(tmp), mesh_io_error);
  std::fclose(tmp);
}
