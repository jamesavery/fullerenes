#include <gtest/gtest.h>
#include "fullerenes/fullerenegraph.hh"
#include "fullerenes/triangulation.hh"
#include "fullerenes/layout2d.hh"

using namespace std;

// Build an unoriented graph by extracting undirected edges and rebuilding neighbour lists.
static Graph to_unoriented(const Graph& G) {
  neighbours_t nb(G.N);
  for(node_t u = 0; u < G.N; u++)
    for(node_t v : G.nbrs(u))
      if(u < v) {
        nb[u].push_back(v);
        nb[v].push_back(u);
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

// Test 4: orient_neighbours on (oriented graph + oriented layout) is idempotent
TEST_F(OrientationTest, OrientNeighboursIdempotent) {
  for(const auto& [name, G] : test_graphs) {
    SCOPED_TRACE(name);
    PlanarGraph PG(G);
    vector<coord2d> layout = PG.tutte_layout();

    // Apply orient_neighbours to an already-oriented graph
    PlanarGraph PG2(G);
    layout2d::orient_neighbours(PG2, layout);

    EXPECT_TRUE(PG2.is_consistently_oriented())
      << "orient_neighbours on oriented graph should stay oriented";
    EXPECT_TRUE(same_orientation(G, PG2))
      << "orient_neighbours on oriented graph should not change orientation";
  }
}

// Test 5: planar_orient produces consistent orientation from unoriented input
TEST_F(OrientationTest, PlanarOrientFromUnoriented) {
  for(const auto& [name, G] : test_graphs) {
    SCOPED_TRACE(name);
    Graph unoriented = to_unoriented(G);
    bool ok = layout2d::planar_orient(unoriented);
    ASSERT_TRUE(ok) << name << " should be planar";
    EXPECT_TRUE(unoriented.is_consistently_oriented())
      << "planar_orient should produce consistent orientation for " << name;
  }
}

// Test 5b: planar_orient produces correct face structure (12 pentagons, rest hexagons)
TEST_F(OrientationTest, PlanarOrientFaceStructure) {
  for(const auto& [name, G] : test_graphs) {
    SCOPED_TRACE(name);
    Graph unoriented = to_unoriented(G);
    layout2d::planar_orient(unoriented);
    ASSERT_TRUE(unoriented.is_consistently_oriented());

    PlanarGraph PG(unoriented);
    vector<face_t> faces = PG.compute_faces_oriented(6);
    int n5 = 0, n6 = 0;
    for(const auto& f : faces) {
      if(f.size() == 5) n5++;
      else if(f.size() == 6) n6++;
    }
    EXPECT_EQ(n5, 12) << name << ": should have 12 pentagons";
    EXPECT_EQ(n6, G.N/2 + 2 - 12) << name << ": wrong hexagon count";
  }
}

// Test 6: orient_neighbours from an oriented layout always produces correct orientation,
// regardless of the input graph's neighbour ordering
TEST_F(OrientationTest, OrientFromGoodLayout) {
  for(const auto& [name, G] : test_graphs) {
    SCOPED_TRACE(name);

    // First get a good layout from the oriented graph
    PlanarGraph oriented(G);
    vector<coord2d> layout = oriented.tutte_layout();
    ASSERT_TRUE(layout2d::layout_is_crossingfree(oriented, layout));

    // Now apply orient_neighbours to the unoriented graph with this good layout
    Graph unoriented = to_unoriented(G);
    PlanarGraph PG(unoriented);
    layout2d::orient_neighbours(PG, layout);

    EXPECT_TRUE(PG.is_consistently_oriented())
      << "orient_neighbours with good layout should always work";
    EXPECT_TRUE(same_orientation(G, PG))
      << "Should match original orientation when using good layout";
  }
}

// Test 7: Face structure is preserved through orient from good layout
TEST_F(OrientationTest, FacesPreservedFromGoodLayout) {
  for(const auto& [name, G] : test_graphs) {
    SCOPED_TRACE(name);

    PlanarGraph orig(G);
    vector<face_t> faces_orig = orig.compute_faces_oriented(6);

    // Get good layout, then orient unoriented graph
    vector<coord2d> layout = orig.tutte_layout();
    Graph unoriented = to_unoriented(G);
    PlanarGraph PG(unoriented);
    layout2d::orient_neighbours(PG, layout);
    vector<face_t> faces_recon = PG.compute_faces_oriented(6);

    ASSERT_EQ(faces_orig.size(), faces_recon.size())
      << name << ": face count mismatch";

    vector<int> sizes_orig, sizes_recon;
    for(const auto& f : faces_orig) sizes_orig.push_back(f.size());
    for(const auto& f : faces_recon) sizes_recon.push_back(f.size());
    sort(sizes_orig.begin(), sizes_orig.end());
    sort(sizes_recon.begin(), sizes_recon.end());
    EXPECT_EQ(sizes_orig, sizes_recon)
      << name << ": face size distributions should match";

    // For fullerenes: exactly 12 pentagons, rest hexagons
    int n5 = count(sizes_recon.begin(), sizes_recon.end(), 5);
    int n6 = count(sizes_recon.begin(), sizes_recon.end(), 6);
    EXPECT_EQ(n5, 12) << name << ": should have 12 pentagons";
    EXPECT_EQ(n6, G.N/2 + 2 - 12) << name << ": wrong hexagon count";
  }
}

// Test 8: Adjacency matrix -> planar_orient path (simulates Fortran interface)
TEST_F(OrientationTest, AdjacencyMatrixPlanarOrient) {
  for(const auto& [name, G] : test_graphs) {
    SCOPED_TRACE(name);

    // Build adjacency matrix
    int N = G.N;
    vector<int> adj(N * N, 0);
    for(node_t u = 0; u < N; u++)
      for(node_t v : G.nbrs(u))
        adj[u * N + v] = 1;

    // Reconstruct from adjacency matrix (loses orientation)
    neighbours_t nb(N);
    for(int i = 0; i < N; i++)
      for(int j = i + 1; j < N; j++)
        if(adj[i * N + j]) {
          nb[i].push_back(j);
          nb[j].push_back(i);
        }

    Graph recon(nb);
    bool ok = layout2d::planar_orient(recon);
    ASSERT_TRUE(ok) << "Should be planar";
    EXPECT_TRUE(recon.is_consistently_oriented());
    EXPECT_TRUE(same_edges(G, recon));

    // Verify face structure
    PlanarGraph PG(recon);
    vector<face_t> faces = PG.compute_faces_oriented(6);
    int n5 = 0, n6 = 0;
    for(const auto& f : faces) {
      if(f.size() == 5) n5++;
      else if(f.size() == 6) n6++;
    }
    EXPECT_EQ(n5, 12) << name << ": should have 12 pentagons";
    EXPECT_EQ(n6, G.N/2 + 2 - 12) << name << ": wrong hexagon count";
  }
}
