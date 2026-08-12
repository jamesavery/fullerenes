// GraphView::oriented_surface / is_consistently_oriented / component_surfaces.
//
// The predicate's whole content is Euler's relation faces == E - N + 2 - 2g, so
// every case here is a graph whose face count is known by hand, checked against
// a genus the test states out loud.  The corrupted cases are the ones that
// matter: they are perfectly consistent rotation systems -- OF THE WRONG
// SURFACE -- and the pre-2026-08-09 predicate accepted every one of them.

#include <gtest/gtest.h>
#include <stdexcept>
#include "fullerenes/graph.hh"
#include "fullerenes/fullerenegraph.hh"

using Code = OrientedSurface::Code;

// K4 embedded in the sphere: 4 vertices, 6 edges, 4 triangular faces.
static Graph k4_planar() {
  return Graph{{1,2,3},{0,3,2},{0,1,3},{0,2,1}};
}

// The 3x3 quotient of the triangular lattice: the flat torus triangulated by 9
// vertices, 27 edges and 18 triangles (chi = 0, genus 1).  Vertex (i,j) sits at
// i*a + j*b with a = (1,0), b = (1/2, sqrt(3)/2), and its six lattice
// neighbours are +a, +b, b-a, -a, -b, a-b at 0, 60, 120, 180, 240, 300 degrees
// -- that CCW list IS the rotation system.
static Graph torus_triangulation_3x3() {
  auto id = [](int i, int j){ return node_t((((i%3)+3)%3)*3 + (((j%3)+3)%3)); };
  const int di[6] = { 1, 0, -1, -1,  0,  1 };
  const int dj[6] = { 0, 1,  1,  0, -1, -1 };

  Graph T(size_t(9), 6);
  for(int i=0;i<3;i++)
    for(int j=0;j<3;j++)
      for(int k=0;k<6;k++)
        T.push_back(id(i,j), id(i+di[k], j+dj[k]));
  return T;
}

// Reverse the cyclic order at one degree-3 vertex, by transposing two of its
// three entries.  Nothing else changes: adjacency stays symmetric, every degree
// stays put, and the result is still a valid rotation system -- of another
// surface.
static void transpose_row(Graph& G, node_t u, int i, int j) {
  std::swap(G[u][i], G[u][j]);
}

// ---------------------------------------------------------------------------
// The core predicate
// ---------------------------------------------------------------------------

TEST(OrientedSurface, PlanarK4IsGenus0) {
  const OrientedSurface S = k4_planar().oriented_surface(0);
  EXPECT_EQ(S.code, Code::Ok);
  EXPECT_EQ(S.faces, 4);       // 4 - 6 + 4 = 2
  EXPECT_EQ(S.genus, 0);
  EXPECT_TRUE(k4_planar().is_consistently_oriented());
}

// Acceptance case 1: row(0) = [2,1,3] claiming genus 0 must FAIL, 2 faces
// against the 4 that genus 0 requires.  The old predicate returned true.
TEST(OrientedSurface, CorruptedK4IsATorusNotAPlane) {
  Graph G = k4_planar();
  transpose_row(G, 0, 0, 1);                 // {1,2,3} -> {2,1,3}

  const OrientedSurface S = G.oriented_surface(0);
  EXPECT_EQ(S.code, Code::GenusMismatch);
  EXPECT_EQ(S.faces, 2);                     // one 9-arc face and one 3-arc face
  EXPECT_EQ(S.genus, 1);                     // 4 - 6 + 2 = 0 = 2 - 2*1
  EXPECT_FALSE(G.is_consistently_oriented());

  // It is not garbage -- it is a consistently oriented TORUS, and says so.
  EXPECT_TRUE(G.is_consistently_oriented(1));
  EXPECT_EQ(G.oriented_surface(1).code, Code::Ok);
}

TEST(OrientedSurface, C20IsGenus0) {
  const FullereneGraph C20 = FullereneGraph::C20();
  const OrientedSurface S = C20.oriented_surface(0);
  EXPECT_EQ(S.code, Code::Ok);
  EXPECT_EQ(S.faces, 12);                    // 20 - 30 + 12 = 2
  EXPECT_EQ(S.genus, 0);
}

// Acceptance case 2: transposing two entries of one degree-3 row of C20 merges
// its three pentagons into one 15-arc face, 12 -> 10, chi = 0.
TEST(OrientedSurface, CorruptedC20IsATorusNotAFullerene) {
  Graph G(FullereneGraph::C20());
  transpose_row(G, 0, 0, 1);                 // {1,4,7} -> {4,1,7}

  const OrientedSurface S = G.oriented_surface(0);
  EXPECT_EQ(S.code, Code::GenusMismatch);
  EXPECT_EQ(S.faces, 10);
  EXPECT_EQ(S.genus, 1);
  EXPECT_FALSE(G.is_consistently_oriented());
  EXPECT_TRUE(G.is_consistently_oriented(1));
}

// Acceptance case 3: the genus is NOT hardcoded.  A genuine torus triangulation
// declaring genus 1 passes, and the same object declaring genus 0 fails.
TEST(OrientedSurface, TorusTriangulationDeclaringGenus1Passes) {
  const Graph T = torus_triangulation_3x3();
  ASSERT_EQ(T.count_edges(), 27u);

  const OrientedSurface S = T.oriented_surface(1);
  EXPECT_EQ(S.code, Code::Ok);
  EXPECT_EQ(S.faces, 18);                    // 9 - 27 + 18 = 0
  EXPECT_EQ(S.genus, 1);
  EXPECT_TRUE(T.is_consistently_oriented(1));

  EXPECT_FALSE(T.is_consistently_oriented(0));
  EXPECT_EQ(T.oriented_surface(0).code, Code::GenusMismatch);
}

// Degenerate input FAILS rather than passing vacuously: it was an edgeless view
// slipping through this predicate that used to let orient_triangulation SIGSEGV
// (that function is now gone -- see tests/volume-moments-test.cc, which makes
// the same case against the surface integrals that reached it).
TEST(OrientedSurface, DegenerateGraphsFail) {
  const Graph edgeless(size_t(8), 3);        // 8 vertices, no edges
  EXPECT_FALSE(edgeless.is_consistently_oriented());
  EXPECT_EQ(edgeless.oriented_surface().code, Code::Degenerate);

  const Graph empty;
  EXPECT_FALSE(empty.is_consistently_oriented());
  EXPECT_EQ(empty.oriented_surface().code, Code::Degenerate);
}

// The structural half of the test, still named separately from the topological
// half: an arc whose reverse is missing is not a rotation system at all.
TEST(OrientedSurface, AsymmetricAdjacencyIsItsOwnOutcome) {
  Graph G = k4_planar();
  G.erase_at(0, 0);                          // drop arc 0->1, keep 1->0

  const OrientedSurface S = G.oriented_surface(0);
  EXPECT_EQ(S.code, Code::AsymmetricAdjacency);
  EXPECT_FALSE(G.is_consistently_oriented());
}

// The SIMPLE precondition is enforced, not merely documented: with two arcs
// between the same pair, find(v,u) resolves both to one slot, phi stops being a
// permutation, and the walk would otherwise enter a cycle missing its start and
// spin forever.  It throws instead.
TEST(OrientedSurface, ParallelArcsThrowRatherThanHang) {
  const Graph doubled{{1,1},{0,0}};            // two vertices, two parallel arcs
  ASSERT_TRUE(doubled.adjacency_is_symmetric());   // so it does reach the walk
  EXPECT_THROW((void)doubled.oriented_surface(0), std::logic_error);
}

// ---------------------------------------------------------------------------
// The precondition, enforced
// ---------------------------------------------------------------------------

// require_oriented_surface is what the operations that need an embedding call
// instead of assert()ing (an assert can end a month-long unsupervised sweep, and
// vanishes entirely under -DNDEBUG).  What it must deliver over the assert is a
// diagnosis: WHO refused, and WHAT surface it was handed instead.
TEST(RequireOrientedSurface, NamesTheOperationAndTheSurfaceItGot) {
  const Graph K4 = k4_planar();
  EXPECT_NO_THROW(require_oriented_surface(K4, "test"));

  Graph G = k4_planar();
  transpose_row(G, 0, 0, 1);

  try {
    require_oriented_surface(G, "MyOperation");
    FAIL() << "a genus-1 rotation system claiming genus 0 must be refused";
  } catch(const unoriented_surface_error& e) {
    EXPECT_EQ(e.surface.code,  Code::GenusMismatch);
    EXPECT_EQ(e.surface.faces, 2);
    EXPECT_EQ(e.surface.genus, 1);
    EXPECT_EQ(e.expected_genus, 0);

    const string what(e.what());
    printf("  %s\n", e.what());
    EXPECT_NE(what.find("MyOperation"), string::npos);   // who refused
    EXPECT_NE(what.find("2 faces"),     string::npos);   // what it got
    EXPECT_NE(what.find("genus 1"),     string::npos);
  }

  // A torus is refused only against a claim it does not meet.
  EXPECT_NO_THROW(require_oriented_surface(G, "MyOperation", 1));

  // The two structural codes get their own sentences rather than a genus count.
  const Graph edgeless(size_t(4), 3);
  EXPECT_THROW(require_oriented_surface(edgeless, "MyOperation"), unoriented_surface_error);
}

// ---------------------------------------------------------------------------
// The multi-component wrapper
// ---------------------------------------------------------------------------

// Acceptance case 4: three planar trees.  Every tree has exactly one face
// (2E arcs in one orbit), so E - N + 2 - F = 0 for each.
TEST(ComponentSurfaces, PlanarForestOfThreePasses) {
  const Graph forest{
    {1}, {0,2}, {1},                        // path   0-1-2
    {4,5,6}, {3}, {3}, {3},                 // star   3-{4,5,6}
    {8}, {7}                                // edge   7-8
  };
  const vector<OrientedSurface> S = forest.component_surfaces();
  ASSERT_EQ(S.size(), 3u);
  for(const OrientedSurface& s: S){
    EXPECT_EQ(s.code, Code::Ok);
    EXPECT_EQ(s.faces, 1);
    EXPECT_EQ(s.genus, 0);
  }
}

// Acceptance case 5: corrupting ONE component's rotation fails exactly that
// component, and the wrapper says which.
TEST(ComponentSurfaces, OneCorruptedComponentFailsAlone) {
  Graph G(size_t(12), 3);                    // three disjoint copies of K4
  const Graph K4 = k4_planar();
  for(int c=0;c<3;c++)
    for(node_t u=0;u<4;u++)
      for(node_t v: K4[u]) G.push_back(4*c+u, 4*c+v);

  {
    const vector<OrientedSurface> S = G.component_surfaces();
    ASSERT_EQ(S.size(), 3u);
    for(const OrientedSurface& s: S){ EXPECT_EQ(s.code, Code::Ok); EXPECT_EQ(s.faces, 4); }
  }

  transpose_row(G, 4, 0, 1);                 // vertex 0 of the middle copy

  const vector<OrientedSurface> S = G.component_surfaces();
  ASSERT_EQ(S.size(), 3u);
  EXPECT_EQ(S[0].code, Code::Ok);            // connected_components() is index-ordered,
  EXPECT_EQ(S[2].code, Code::Ok);            // so component 1 is the middle copy
  EXPECT_EQ(S[1].code, Code::GenusMismatch);
  EXPECT_EQ(S[1].faces, 2);
  EXPECT_EQ(S[1].genus, 1);

  // ... and declaring that one to be a torus makes the whole forest agree.
  const vector<OrientedSurface> claimed = G.component_surfaces({0,1,0});
  for(const OrientedSurface& s: claimed) EXPECT_EQ(s.code, Code::Ok);
}

// The wrapper agrees with the core predicate on connected input -- the contract
// that lets the core predicate charge nothing for connectivity.
TEST(ComponentSurfaces, AgreesWithTheCorePredicateOnConnectedInput) {
  const Graph T = torus_triangulation_3x3();
  const vector<OrientedSurface> S = T.component_surfaces({1});
  ASSERT_EQ(S.size(), 1u);
  EXPECT_EQ(S[0].code, T.oriented_surface(1).code);
  EXPECT_EQ(S[0].faces, T.oriented_surface(1).faces);
  EXPECT_EQ(S[0].genus, T.oriented_surface(1).genus);
}
