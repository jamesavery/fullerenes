#include <gtest/gtest.h>
#include "fullerenes/polyhedron.hh"

using namespace std;

// The convex hull of a regular dodecahedron (C20) is itself.
TEST(ConvexHull, C20IsOriented) {
  Polyhedron P = Polyhedron::C20();
  Polyhedron CH = P.convex_hull();

  EXPECT_TRUE(CH.is_consistently_oriented());
  EXPECT_EQ(CH.N, P.N);
}

// Volume should be positive (outward normals convention).
TEST(ConvexHull, C20VolumePositive) {
  Polyhedron P = Polyhedron::C20();
  Polyhedron CH = P.convex_hull();

  EXPECT_GT(CH.volume(), 0);
}

// Convex hull of an icosahedron (dual of C20).
TEST(ConvexHull, IcosahedronIsOriented) {
  Polyhedron ico = Polyhedron::C20().dual();
  Polyhedron CH = ico.convex_hull();

  EXPECT_TRUE(CH.is_consistently_oriented());
  EXPECT_EQ(CH.N, ico.N);
  EXPECT_GT(CH.volume(), 0);
}

// Faces computed from the oriented graph should all be triangles
// (convex hull is a triangulation).
TEST(ConvexHull, FacesAreTriangles) {
  Polyhedron CH = Polyhedron::C20().convex_hull();

  auto fs = CH.faces();
  EXPECT_FALSE(fs.empty());
  for (const auto& f : fs)
    EXPECT_EQ(f.size(), 3u);
}

// Euler formula: V - E + F = 2 for a convex polyhedron.
TEST(ConvexHull, EulerFormula) {
  Polyhedron CH = Polyhedron::C20().convex_hull();

  int V = CH.N;
  int E = CH.count_edges();
  int F = CH.faces().size();
  EXPECT_EQ(V - E + F, 2);
}
