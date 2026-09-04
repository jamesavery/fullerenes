// Regression coverage for the four ORIENTATION BOUNDARIES (layout2d.hh): the
// input formats that carry no rotation system by nature, and are therefore the
// only places in the library allowed to write one.  Every case below drives the
// boundary's own entry point rather than layout2d::planar_orient /
// orient_neighbours directly -- calling those from anywhere else is exactly what
// the boundary rule forbids, and a test that did it would not be testing the
// path any caller can reach.
//
//   1. Polyhedron::from_mol2           -- FacesPreservedThroughMol2
//   2. oriented_graph_from_adjacency   -- AdjacencyMatrixBoundary (via the
//                                         Fortran ABI's new_fullerene_graph_)
//   3. Polyhedron(points, tolerance)   -- PolyhedronFromCoordinates,
//                                         PolyhedronFromCoordinatesFaceStructure
//   4. set_layout2d_                   -- SetLayout2dFromGoodLayout

#include <gtest/gtest.h>
#include "fullerenes/fullerenegraph.hh"
#include "fullerenes/triangulation.hh"
#include "fullerenes/polyhedron.hh"
#include "fullerenes/layout2d.hh"

#include <cstdio>

using namespace std;

// The Fortran C ABI, whose C++ side is src/c++/graph_fortran.cc.  Boundaries 2
// and 4 have no other entry point -- oriented_graph_from_adjacency is
// file-static and set_layout2d_ is only ever called from Fortran.
extern "C" {
  FullereneGraph* new_fullerene_graph_(const int *nmax, const int *N, const int *adjacency);
  void            delete_fullerene_graph_(FullereneGraph**);
  void            set_layout2d_(PlanarGraph **g, const double *layout2d);
}

// Build an unoriented graph by extracting undirected edges and rebuilding neighbour lists.
static Graph to_unoriented(const Graph& G) {
  Graph nb(G.N, GRAPH_DMAX);
  for(node_t u = 0; u < G.N; u++)
    for(node_t v : G.nbrs(u))
      if(u < v) {
        nb.push_back(u, v);
        nb.push_back(v, u);
      }
  return Graph(nb);
}

// Check that two graphs have the same undirected edge set
static bool same_edges(const Graph& A, const Graph& B) {
  if(A.N != B.N) return false;
  set<edge_t> ea, eb;
  for(node_t u = 0; u < A.N; u++)
    for(node_t v : A.nbrs(u))
      if(u < v) ea.insert({u,v});
  for(node_t u = 0; u < B.N; u++)
    for(node_t v : B.nbrs(u))
      if(u < v) eb.insert({u,v});
  return ea == eb;
}

// Check that two oriented graphs have identical CW neighbour lists (up to cyclic rotation)
static bool same_orientation(const Graph& A, const Graph& B) {
  if(A.N != B.N) return false;
  for(node_t u = 0; u < A.N; u++) {
    auto na = A.nbrs(u);
    auto nb_list = B.nbrs(u);
    int deg = na.size();
    if((int)nb_list.size() != deg) return false;
    int offset = -1;
    for(int i = 0; i < deg; i++)
      if(nb_list[i] == na[0]) { offset = i; break; }
    if(offset < 0) return false;
    for(int i = 0; i < deg; i++)
      if(na[i] != nb_list[(offset + i) % deg]) return false;
  }
  return true;
}

// Build test graphs from known-good sources (no database needed).
static vector<pair<string, FullereneGraph>> make_test_graphs() {
  vector<pair<string, FullereneGraph>> graphs;

  // C20: dodecahedron
  graphs.push_back({"C20", FullereneGraph::C20()});

  // C20 dual = icosahedron (12 vertices, all degree 5)
  Triangulation ico(vector<int>(12, 5));

  // C80 via GC(2,0)
  {
    Triangulation gc = ico.GCtransform(2, 0);
    graphs.push_back({"C80-GC(2,0)", FullereneGraph(gc.dual_graph())});
  }

  // C60 via GC(1,1) (leapfrog = Buckminsterfullerene)
  {
    Triangulation gc = ico.GCtransform(1, 1);
    graphs.push_back({"C60-GC(1,1)", FullereneGraph(gc.dual_graph())});
  }

  // C180 via GC(3,0)
  {
    Triangulation gc = ico.GCtransform(3, 0);
    graphs.push_back({"C180-GC(3,0)", FullereneGraph(gc.dual_graph())});
  }

  // C140 via GC(2,1)
  {
    Triangulation gc = ico.GCtransform(2, 1);
    graphs.push_back({"C140-GC(2,1)", FullereneGraph(gc.dual_graph())});
  }

  return graphs;
}

class OrientationTest : public ::testing::Test {
protected:
  static vector<pair<string, FullereneGraph>> test_graphs;
  static void SetUpTestSuite() { test_graphs = make_test_graphs(); }
};

vector<pair<string, FullereneGraph>> OrientationTest::test_graphs;

// Test 1: All test graphs start out oriented
TEST_F(OrientationTest, TestGraphsAreOriented) {
  for(const auto& [name, G] : test_graphs) {
    SCOPED_TRACE(name + " (N=" + to_string(G.N) + ")");
    EXPECT_TRUE(G.is_consistently_oriented());
  }
}

// Test 2: to_unoriented produces a graph with same edges but not consistently oriented
TEST_F(OrientationTest, UnorientedProperties) {
  for(const auto& [name, G] : test_graphs) {
    SCOPED_TRACE(name);
    Graph unoriented = to_unoriented(G);
    EXPECT_TRUE(same_edges(G, unoriented));
    // For cubic graphs, to_unoriented builds neighbours in (u<v) order,
    // which almost certainly breaks the CW face-tracing property.
    // We don't require is_consistently_oriented() to be false — it could be
    // accidentally true for highly symmetric graphs.
  }
}

// Test 3: Tutte layout works on oriented graphs and is crossing-free
TEST_F(OrientationTest, TutteLayoutOriented) {
  for(const auto& [name, G] : test_graphs) {
    SCOPED_TRACE(name);
    PlanarGraph PG(G);
    vector<coord2d> layout = PG.tutte_layout();
    ASSERT_EQ(layout.size(), (size_t)G.N);
    EXPECT_TRUE(layout2d::layout_is_crossingfree(PG, layout))
      << "Tutte layout of oriented " << name << " should be crossing-free";
  }
}

// Test 4 (OrientNeighboursIdempotent) is GONE.  It applied
// layout2d::orient_neighbours to an already-oriented graph with that graph's own
// layout and checked nothing changed -- a no-op on input no boundary ever sees.
// Its subject was the internal repair calls, which no longer exist.

// Sizes and distributions of the faces of an oriented planar graph, sorted:
// what survives a change of embedding is the multiset, not the face list.
static vector<int> face_sizes(const PlanarGraphView& G, int face_max) {
  vector<int> sizes;
  for(const face_t& f: G.compute_faces_oriented(face_max)) sizes.push_back(f.size());
  sort(sizes.begin(), sizes.end());
  return sizes;
}

// -------------------------------------------------------------------------
// BOUNDARY 3: Polyhedron(points, tolerance) -- coordinates only.
//
// Bonds are inferred by distance and pushed in ascending index order, so what
// reaches the constructor is symmetric but has no rotation system at all.  The
// two cases below are what used to be tests 5 and 5b, moved off
// layout2d::planar_orient and onto the boundary that calls it -- and onto the
// input that boundary exists for.  Until the orientation predicate became a real
// test this path silently produced a genus-5 "polyhedron" (C20: 2 faces).
// -------------------------------------------------------------------------

// A regular icosahedron: 12 vertices of degree 5, edge 2, next-nearest distance
// 2*phi -- comfortably outside the 1.2 bond tolerance.
static vector<coord3d> icosahedron_points() {
  const double phi = (1 + sqrt(5.0))/2;
  vector<coord3d> pts;
  for(double s1: {-1.0, 1.0})
    for(double s2: {-1.0, 1.0}){
      pts.push_back(coord3d(0, s1, s2*phi));
      pts.push_back(coord3d(s1, s2*phi, 0));
      pts.push_back(coord3d(s2*phi, 0, s1));
    }
  return pts;
}

TEST_F(OrientationTest, PolyhedronFromCoordinates) {
  const Polyhedron C20 = Polyhedron::C20();
  vector<pair<string, vector<coord3d>>> clouds{
    {"dodecahedron", vector<coord3d>(C20.points.begin(), C20.points.end())},
    {"icosahedron",  icosahedron_points()},
  };

  for(const auto& [name, pts]: clouds){
    SCOPED_TRACE(name);
    const Polyhedron P(pts);                    // the boundary
    EXPECT_TRUE(P.is_consistently_oriented())
      << name << ": construction from coordinates must establish an orientation";
    EXPECT_GT(P.volume(), 0.0)
      << name << ": and the outward one -- the signed-volume half of the boundary";
  }
}

TEST_F(OrientationTest, PolyhedronFromCoordinatesFaceStructure) {
  const Polyhedron C20 = Polyhedron::C20();
  const Polyhedron P(vector<coord3d>(C20.points.begin(), C20.points.end()));
  ASSERT_TRUE(P.is_consistently_oriented());
  EXPECT_EQ(face_sizes(P, 6), vector<int>(12, 5)) << "the dodecahedron has 12 pentagons";

  const Polyhedron I(icosahedron_points());
  ASSERT_TRUE(I.is_consistently_oriented());
  EXPECT_EQ(face_sizes(I, 3), vector<int>(20, 3)) << "the icosahedron has 20 triangles";
}

// -------------------------------------------------------------------------
// BOUNDARY 4: set_layout2d_ -- an externally computed drawing, re-imported
// across the Fortran ABI.  Was test 6, which called layout2d::orient_neighbours
// directly.
//
// STILL RED, and deliberately so: the ASSERT below is on
// layout2d::layout_is_crossingfree, which is written in slope-intercept form and
// calls the exact collinearities of an icosahedrally symmetric Tutte layout
// (C60-GC(1,1), C180-GC(3,0)) intersections.  That defect is tracked separately
// (claude-projects/curvature-flow/refactor-debt.md,
// 2026-08-09-layout-is-crossingfree-slope-form) and is not this refactor's to
// fix; the assertion stays because the drawing being crossing-free is the
// precondition this case is about.
// -------------------------------------------------------------------------
TEST_F(OrientationTest, SetLayout2dFromGoodLayout) {
  for(const auto& [name, G] : test_graphs) {
    SCOPED_TRACE(name);

    PlanarGraph oriented(G);
    vector<coord2d> layout = oriented.tutte_layout();
    ASSERT_TRUE(layout2d::layout_is_crossingfree(oriented, layout));

    vector<double> layout_data(2*G.N);
    for(node_t u = 0; u < G.N; u++){
      layout_data[2*u]   = layout[u].first;
      layout_data[2*u+1] = layout[u].second;
    }

    PlanarGraph PG(to_unoriented(G)), *p = &PG;
    set_layout2d_(&p, layout_data.data());

    EXPECT_TRUE(PG.is_consistently_oriented())
      << "a crossing-free drawing must re-import as an orientation";
    EXPECT_TRUE(same_orientation(G, PG))
      << "and as the original one, since the drawing came from the original";
  }
}

// The other half of boundary 4, which is the reason it is the only boundary that
// rolls back: its input may already be oriented, so a bad drawing can make things
// WORSE.  A collinear drawing is one -- every neighbour of u lies at angle 0 or
// pi from it, so the CCW sort has nothing to sort by.
TEST_F(OrientationTest, SetLayout2dRefusesACorruptingDrawing) {
  const FullereneGraph& G = test_graphs[0].second;          // C20

  vector<double> layout_data(2*G.N);
  vector<coord2d> collinear(G.N);
  for(node_t u = 0; u < G.N; u++){
    collinear[u]       = coord2d(u, 0);
    layout_data[2*u]   = u;
    layout_data[2*u+1] = 0;
  }

  {  // the fixture is supposed to corrupt: check that it does, at the raw sort
    PlanarGraph probe(G);
    layout2d::orient_neighbours(probe, collinear);
    ASSERT_FALSE(probe.is_consistently_oriented())
      << "the collinear drawing was expected to destroy the orientation";
  }

  PlanarGraph PG(G), *p = &PG;
  ASSERT_TRUE(PG.is_consistently_oriented());
  set_layout2d_(&p, layout_data.data());

  EXPECT_TRUE(PG.is_consistently_oriented())
    << "the boundary must not install a drawing that breaks the invariant";
  EXPECT_TRUE(same_orientation(G, PG))
    << "the graph must keep exactly the rows it arrived with";
}

// -------------------------------------------------------------------------
// BOUNDARY 1: Polyhedron::from_mol2 -- a bare list of unordered bonds.  Was
// test 7 (faces preserved through a re-orientation), moved onto the format
// whose records genuinely carry no order.
// -------------------------------------------------------------------------
TEST_F(OrientationTest, FacesPreservedThroughMol2) {
  const Polyhedron P = Polyhedron::C20();
  ASSERT_TRUE(P.is_consistently_oriented());

  FILE* f = tmpfile();
  ASSERT_NE(f, nullptr);
  ASSERT_TRUE(Polyhedron::to_mol2(P, f));
  rewind(f);
  const Polyhedron Q = Polyhedron::from_mol2(f);
  fclose(f);

  EXPECT_TRUE(Q.is_consistently_oriented())
    << "from_mol2 must establish an orientation, the file having none";
  EXPECT_TRUE(same_edges(P, Q));
  EXPECT_EQ(face_sizes(Q, 6), face_sizes(P, 6))
    << "and the same faces: 12 pentagons, whichever way round they are wound";
}

// -------------------------------------------------------------------------
// BOUNDARY 2: oriented_graph_from_adjacency -- an adjacency MATRIX, which fixes
// the edge set and nothing else.  Was test 8, which open-coded the boundary's
// body; it now goes through the Fortran ABI entry point that is the only way to
// reach it.
// -------------------------------------------------------------------------
TEST_F(OrientationTest, AdjacencyMatrixBoundary) {
  for(const auto& [name, G] : test_graphs) {
    SCOPED_TRACE(name);

    const int N = G.N;
    vector<int> adj(size_t(N)*N, 0);
    for(node_t u = 0; u < N; u++)
      for(node_t v : G.nbrs(u))
        adj[size_t(u)*N + v] = 1;

    FullereneGraph *H = new_fullerene_graph_(&N, &N, adj.data());
    ASSERT_NE(H, nullptr) << name << ": the boundary refused a planar adjacency matrix";

    EXPECT_TRUE(H->is_consistently_oriented());
    EXPECT_TRUE(same_edges(G, *H));
    EXPECT_EQ(face_sizes(*H, 6), face_sizes(PlanarGraph(G), 6))
      << name << ": 12 pentagons and the rest hexagons";

    delete_fullerene_graph_(&H);
  }
}

// NOT TESTED, and not for want of trying: the boundary's REFUSAL path, i.e. an
// adjacency matrix that is not planar (K5) reaching new_fullerene_graph_ and
// coming back NULL.  It cannot be exercised from a test binary, because
// PlanarGraph::tutte_layout abort()s -- src/c++/layout.cc:244 and :253 -- before
// planar_orient can return false, taking the process with it.  That is one entry
// in the library-wide assert/abort inventory
// (claude-projects/curvature-flow/refactor-debt.md,
// 2026-08-09-library-assert-abort-inventory); when tutte_layout throws instead,
// this case is three lines.
