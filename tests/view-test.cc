// Tests for span-based view construction across the full class hierarchy.
//
// Validates that Graph, PlanarGraph, CubicGraph, FullereneGraph,
// Triangulation, Polyhedron, and Deltahedron all work correctly both
// with owned storage and as views over externally-owned memory.

#include <gtest/gtest.h>
#include <fullerenes/fullerenegraph.hh>
#include <fullerenes/triangulation.hh>
#include <fullerenes/polyhedron.hh>
#include <fullerenes/isomerdb.hh>

// ---------------------------------------------------------------------------
// Graph: owned vs view
// ---------------------------------------------------------------------------

TEST(GraphView, OwnedGraphBasics) {
    // Graph with owned storage
    Graph G = FullereneGraph::C20();
    EXPECT_EQ(G.N, 20);
    EXPECT_EQ(G.dmax, 3);
    EXPECT_TRUE(G.owns_memory());  // owns its data

    for (node_t u = 0; u < G.N; u++)
        EXPECT_EQ(G.degree(u), 3);
}

TEST(GraphView, ViewFromExternalMemory) {
    // Create owned graph, then create a view over its memory
    Graph owned = FullereneGraph::C20();

    Graph view(owned.N, owned.dmax,
               std::span<node_t>(owned.owned_neighbours),
               std::span<uint8_t>(owned.owned_deg));

    EXPECT_FALSE(view.owns_memory());  // view, not owned
    EXPECT_EQ(view.N, 20);
    EXPECT_EQ(view.dmax, 3);

    // View reads the same data
    for (node_t u = 0; u < view.N; u++) {
        EXPECT_EQ(view.degree(u), owned.degree(u));
        auto vn = view.nbrs(u);
        auto on = owned.nbrs(u);
        for (int i = 0; i < 3; i++)
            EXPECT_EQ(vn[i], on[i]);
    }
}

TEST(GraphView, ViewMutatesExternalMemory) {
    // Mutations through a view write to the external buffer
    std::vector<node_t> values(4 * 3, node_t(-1));
    std::vector<uint8_t> deg(4, 0);

    Graph view(4, 3,
               std::span<node_t>(values),
               std::span<uint8_t>(deg));

    view.push_back(0, 1);
    view.push_back(0, 2);
    view.push_back(0, 3);

    // Changes visible in external buffer
    EXPECT_EQ(deg[0], 3);
    EXPECT_EQ(values[0], 1);
    EXPECT_EQ(values[1], 2);
    EXPECT_EQ(values[2], 3);
}

TEST(GraphView, CopyOfViewIsOwned) {
    Graph owned = FullereneGraph::C20();
    Graph view(owned.N, owned.dmax,
               std::span<node_t>(owned.owned_neighbours),
               std::span<uint8_t>(owned.owned_deg));

    // Copy a view -> produces an owned copy
    Graph copy(view);
    EXPECT_FALSE(copy.owns_memory());  // copy of view is still a view

    // But constructing via base_t makes an owned copy
    Graph deep_copy(static_cast<const Graph::base_t&>(view));
    EXPECT_TRUE(deep_copy.owns_memory());  // owned copy
    EXPECT_EQ(deep_copy.degree(0), view.degree(0));
}

TEST(GraphView, MovePreservesView) {
    Graph owned = FullereneGraph::C20();
    Graph view(owned.N, owned.dmax,
               std::span<node_t>(owned.owned_neighbours),
               std::span<uint8_t>(owned.owned_deg));

    Graph moved(std::move(view));
    EXPECT_FALSE(moved.owns_memory());  // still a view
    EXPECT_EQ(moved.N, 20);
    EXPECT_EQ(moved.degree(0), 3);
}

// ---------------------------------------------------------------------------
// CubicGraph / FullereneGraph: view propagation through hierarchy
// ---------------------------------------------------------------------------

TEST(HierarchyView, CubicGraphFromViewGraph) {
    Graph owned = FullereneGraph::C20();
    // Already dmax=3, so CubicGraph won't restride — preserves view
    Graph view(owned.N, owned.dmax,
               std::span<node_t>(owned.owned_neighbours),
               std::span<uint8_t>(owned.owned_deg));

    CubicGraph cg(view);
    // CubicGraph copies the view (via PlanarGraph -> Graph copy ctor)
    // Since owned_neighbours of view is empty, copy is also a view
    EXPECT_FALSE(cg.owns_memory());
    EXPECT_EQ(cg.N, 20);
    EXPECT_EQ(cg.dmax, 3);
}

TEST(HierarchyView, FullereneGraphFromViewGraph) {
    Graph owned = FullereneGraph::C20();
    Graph view(owned.N, owned.dmax,
               std::span<node_t>(owned.owned_neighbours),
               std::span<uint8_t>(owned.owned_deg));

    FullereneGraph fg(view);
    EXPECT_FALSE(fg.owns_memory());
    EXPECT_EQ(fg.N, 20);
    EXPECT_TRUE(fg.is_consistently_oriented());
}

// ---------------------------------------------------------------------------
// Polyhedron: view of both adjacency and coordinates
// ---------------------------------------------------------------------------

TEST(PolyhedronView, OwnedPolyhedron) {
    Polyhedron P = Polyhedron::C20();
    EXPECT_EQ(P.N, 20);
    EXPECT_TRUE(P.points.owns_memory());  // owns coordinates

    auto fs = P.faces();
    EXPECT_EQ(int(fs.size()), 12);  // C20 has 12 pentagonal faces
    for (auto& f : fs) EXPECT_EQ(int(f.size()), 5);
}

TEST(PolyhedronView, ViewCoordinates) {
    Polyhedron owned = Polyhedron::C20();

    // Create a view that shares the coordinate memory
    std::span<coord3d> pts_span(owned.points.owned.data(), owned.points.owned.size());
    Polyhedron view(static_cast<const PlanarGraph&>(owned), pts_span, 6);

    EXPECT_FALSE(view.points.owns_memory());  // coordinate view
    EXPECT_EQ(view.N, 20);

    // Coordinates are the same
    for (node_t u = 0; u < view.N; u++) {
        EXPECT_DOUBLE_EQ(view.points[u][0], owned.points[u][0]);
        EXPECT_DOUBLE_EQ(view.points[u][1], owned.points[u][1]);
        EXPECT_DOUBLE_EQ(view.points[u][2], owned.points[u][2]);
    }

    // Faces computed from adjacency
    auto fs = view.faces();
    EXPECT_EQ(int(fs.size()), 12);
}

TEST(PolyhedronView, ViewMutatesCoordinates) {
    Polyhedron owned = Polyhedron::C20();
    std::span<coord3d> pts_span(owned.points.owned.data(), owned.points.owned.size());
    Polyhedron view(static_cast<const PlanarGraph&>(owned), pts_span, 6);

    // Mutating through the view changes the original
    view.points[0] = coord3d(99.0, 99.0, 99.0);
    EXPECT_DOUBLE_EQ(owned.points[0][0], 99.0);
}

// ---------------------------------------------------------------------------
// Batch slicing: the GPU use case
// ---------------------------------------------------------------------------

TEST(BatchSlicing, GraphBatch) {
    // Simulate GPU batch: one flat allocation, multiple graph views
    const int B = 3;  // batch size
    const int N = 20;
    const int dmax = 3;

    // Build B fullerene graphs
    FullereneGraph fg = FullereneGraph::C20();

    // Allocate batch storage
    std::vector<node_t> all_values(B * N * dmax);
    std::vector<uint8_t> all_deg(B * N);

    // Fill batch by copying the same graph B times
    for (int b = 0; b < B; b++) {
        std::copy(fg.neighbours.begin(), fg.neighbours.end(),
                  all_values.begin() + b * N * dmax);
        std::copy(fg.deg.begin(), fg.deg.end(),
                  all_deg.begin() + b * N);
    }

    // Create views into the batch
    for (int b = 0; b < B; b++) {
        std::span<node_t> vals(all_values.data() + b * N * dmax, N * dmax);
        std::span<uint8_t> degs(all_deg.data() + b * N, N);

        Graph view(N, dmax, vals, degs);
        EXPECT_FALSE(view.owns_memory());
        EXPECT_EQ(view.N, N);

        // Run graph algorithms on the view
        EXPECT_TRUE(view.is_connected());
        EXPECT_TRUE(view.is_consistently_oriented());
        EXPECT_EQ(int(view.count_edges()), 30);  // 3*20/2

        // Construct fullerene graph from the view
        FullereneGraph fg_view(view);
        EXPECT_FALSE(fg_view.owns_memory());

        // Spiral extraction works on view
        vector<int> spiral;
        jumplist_t jumps;
        EXPECT_TRUE(fg_view.get_rspi_from_fg(spiral, jumps));
    }
}

TEST(BatchSlicing, PolyhedronBatch) {
    // Simulate GPU batch with both adjacency and coordinates
    const int B = 3;
    const int N = 20;
    const int dmax = 3;

    Polyhedron P0 = Polyhedron::C20();

    // Allocate batch storage
    std::vector<node_t> all_values(B * N * dmax);
    std::vector<uint8_t> all_deg(B * N);
    std::vector<coord3d> all_points(B * N);

    for (int b = 0; b < B; b++) {
        std::copy(P0.neighbours.begin(), P0.neighbours.end(),
                  all_values.begin() + b * N * dmax);
        std::copy(P0.deg.begin(), P0.deg.end(),
                  all_deg.begin() + b * N);
        std::copy(P0.points.begin(), P0.points.end(),
                  all_points.begin() + b * N);
    }

    // Create per-graph views
    for (int b = 0; b < B; b++) {
        std::span<node_t> vals(all_values.data() + b * N * dmax, N * dmax);
        std::span<uint8_t> degs(all_deg.data() + b * N, N);
        std::span<coord3d> pts(all_points.data() + b * N, N);

        // Build a Graph view, then wrap as PlanarGraph, then Polyhedron
        Graph g_view(N, dmax, vals, degs);
        Polyhedron poly(PlanarGraph(g_view), pts, 6);

        EXPECT_FALSE(poly.owns_memory());
        EXPECT_FALSE(poly.points.owns_memory());
        EXPECT_EQ(poly.N, N);

        // Geometry works on the view
        EXPECT_GT(poly.volume(), 0);
        EXPECT_GT(poly.surface_area(), 0);

        auto fs = poly.faces();
        EXPECT_EQ(int(fs.size()), 12);  // 12 pentagonal faces
    }
}

TEST(BatchSlicing, MutationThroughViewWritesBack) {
    const int N = 20;
    const int dmax = 3;

    Polyhedron P0 = Polyhedron::C20();

    // Flat batch storage (single element for simplicity)
    std::vector<node_t> values(P0.neighbours.begin(), P0.neighbours.end());
    std::vector<uint8_t> deg(P0.deg.begin(), P0.deg.end());
    std::vector<coord3d> pts(P0.points.begin(), P0.points.end());

    // Create view
    Graph g_view(N, dmax,
                 std::span<node_t>(values),
                 std::span<uint8_t>(deg));
    Polyhedron poly(PlanarGraph(g_view), std::span<coord3d>(pts), 6);

    // Move vertex 0
    poly.points[0] = coord3d(1.0, 2.0, 3.0);

    // Verify it wrote through to the flat buffer
    EXPECT_DOUBLE_EQ(pts[0][0], 1.0);
    EXPECT_DOUBLE_EQ(pts[0][1], 2.0);
    EXPECT_DOUBLE_EQ(pts[0][2], 3.0);
}

// ---------------------------------------------------------------------------
// Triangulation: view from Graph
// ---------------------------------------------------------------------------

TEST(TriangulationView, FromGraphView) {
    FullereneGraph fg = FullereneGraph::C20();
    Triangulation T_owned(fg.dual_graph());

    // Create a Graph view of the triangulation's adjacency
    Graph view(T_owned.N, T_owned.dmax,
               std::span<node_t>(T_owned.owned_neighbours),
               std::span<uint8_t>(T_owned.owned_deg));

    Triangulation T_view(view);
    EXPECT_FALSE(T_view.owns_memory());
    EXPECT_EQ(T_view.N, T_owned.N);

    // triangles() works on the view
    auto tris = T_view.triangles();
    EXPECT_EQ(int(tris.size()), 2 * T_view.N - 4);
}

// ---------------------------------------------------------------------------
// Fill constructor
// ---------------------------------------------------------------------------

TEST(GraphFill, FillConstructor) {
    // Graph(N, initial_row) fills all rows with the same values
    Graph G(10, std::vector<node_t>(3, -1));
    EXPECT_EQ(G.N, 10);
    EXPECT_EQ(G.dmax, 3);
    for (node_t u = 0; u < G.N; u++) {
        EXPECT_EQ(G.degree(u), 3);
        for (int i = 0; i < 3; i++)
            EXPECT_EQ(G[u][i], -1);
    }
}

// ---------------------------------------------------------------------------
// Edge case: empty graph
// ---------------------------------------------------------------------------

TEST(GraphView, EmptyGraph) {
    Graph G;
    EXPECT_EQ(G.N, 0);
    EXPECT_FALSE(G.owns_memory());
    EXPECT_TRUE(G.neighbours.empty());
}

TEST(GraphView, EmptyView) {
    std::vector<node_t> values;
    std::vector<uint8_t> deg;
    Graph view(0, 3,
               std::span<node_t>(values),
               std::span<uint8_t>(deg));
    EXPECT_EQ(view.N, 0);
    EXPECT_FALSE(view.owns_memory());
}
