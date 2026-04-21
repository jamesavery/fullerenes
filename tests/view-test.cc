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
#include <fullerenes/batch/batchable.hh>

// ---------------------------------------------------------------------------
// Phase 1: intrinsic batchability contract
// ---------------------------------------------------------------------------
static_assert(batch::batchable_view<GraphView>);
static_assert(batch::batchable_view<PlanarGraphView>);
static_assert(batch::batchable_view<CubicGraphView>);
static_assert(batch::batchable_view<FullereneGraphView>);
static_assert(batch::batchable_view<TriangulationView>);
static_assert(batch::batchable_view<FullereneDualView>);
static_assert(batch::batchable_view<PolyhedronView<double>>);
static_assert(batch::batchable_view<PolyhedronView<float>>);
static_assert(batch::batchable_view<DeltahedronView<double>>);

static_assert(GraphView::n_fields         == 3);
static_assert(CubicGraphView::n_fields    == 3);
static_assert(TriangulationView::n_fields == 3);
static_assert(PolyhedronView<double>::n_fields == 4);
static_assert(DeltahedronView<double>::n_fields == 4);

// Graph-like views must share layout with each other.
static_assert([]{
    return batch::layout_compatible<GraphView, CubicGraphView>(60, 3)
        && batch::layout_compatible<GraphView, TriangulationView>(32, 6)
        && batch::layout_compatible<CubicGraphView, FullereneGraphView>(60, 3)
        && batch::layout_compatible<TriangulationView, FullereneDualView>(32, 6);
}());

// Geometry views share their graph fields with plain graph views.
static_assert([]{
    return batch::layout_compatible<PolyhedronView<double>, PlanarGraphView>(60, 10);
}());

// ---------------------------------------------------------------------------
// Phase 1 contract: runtime checks of to_tuple() / size factors
// ---------------------------------------------------------------------------
TEST(BatchContract, GraphSizeFactors) {
    auto f = GraphView::get_size_factors(60, 10);
    EXPECT_EQ(f[0], 10u);  // neighbours = dmax per vertex
    EXPECT_EQ(f[1], 1u);   // deg = 1 per vertex
    EXPECT_EQ(f[2], 10u);  // twin = dmax per vertex
}

TEST(BatchContract, PolyhedronSizeFactors) {
    auto f = PolyhedronView<double>::get_size_factors(60, 10);
    ASSERT_EQ(f.size(), 4u);
    EXPECT_EQ(f[0], 10u);
    EXPECT_EQ(f[1], 1u);
    EXPECT_EQ(f[2], 10u);
    EXPECT_EQ(f[3], 1u);  // points = 1 per vertex
}

TEST(BatchContract, ToTupleAliasesGraphFields) {
    Graph G = FullereneGraph::C20();
    auto t = G.to_tuple();
    static_assert(std::tuple_size_v<decltype(t)> == 3);
    EXPECT_EQ(std::get<0>(t).data(), G.neighbours.data());
    EXPECT_EQ(std::get<1>(t).data(), G.deg.data());
    EXPECT_EQ(std::get<0>(t).size(), (size_t)(G.N * G.dmax));
}

TEST(BatchContract, PolyhedronToTupleIncludesPoints) {
    Polyhedron P = Polyhedron::C20();
    auto t = P.to_tuple();
    static_assert(std::tuple_size_v<decltype(t)> == 4);
    EXPECT_EQ(std::get<3>(t).data(), P.points.data());
    EXPECT_EQ(std::get<3>(t).size(), (size_t)P.N);
}

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
    // CubicGraph deep-copies from the view (via Owned<PlanarGraphView>)
    EXPECT_TRUE(cg.owns_memory());
    EXPECT_EQ(cg.N, 20);
    EXPECT_EQ(cg.dmax, 3);
}

TEST(HierarchyView, FullereneGraphFromViewGraph) {
    Graph owned = FullereneGraph::C20();
    Graph view(owned.N, owned.dmax,
               std::span<node_t>(owned.owned_neighbours),
               std::span<uint8_t>(owned.owned_deg));

    FullereneGraph fg(view);
    EXPECT_TRUE(fg.owns_memory());
    EXPECT_EQ(fg.N, 20);
    EXPECT_TRUE(fg.is_consistently_oriented());
}

// ---------------------------------------------------------------------------
// Polyhedron: view of both adjacency and coordinates
// ---------------------------------------------------------------------------

TEST(PolyhedronView, OwnedPolyhedron) {
    Polyhedron P = Polyhedron::C20();
    EXPECT_EQ(P.N, 20);
    EXPECT_TRUE(P.owns_memory());  // owns coordinates
    auto fs = P.faces();
    EXPECT_EQ(int(fs.size()), 12);  // C20 has 12 pentagonal faces
    for (auto& f : fs) EXPECT_EQ(int(f.size()), 5);
}

TEST(PolyhedronView, ViewCoordinates) {
    Polyhedron owned = Polyhedron::C20();

    // Create a view that shares the coordinate memory
    std::span<coord3d> pts_span(owned.owned_points.data(), owned.owned_points.size());
    Polyhedron view(static_cast<const PlanarGraphView&>(owned), pts_span, 6);

    // Polyhedron always deep-copies now (both adjacency and points)
    EXPECT_TRUE(view.owns_memory());
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
    // With Owned<PolyhedronView>, points is always owned (deep copy).
    // Test that writing through the span modifies the owned storage.
    owned.points[0] = coord3d(99.0, 99.0, 99.0);
    EXPECT_DOUBLE_EQ(owned.owned_points[0][0], 99.0);
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

        // Construct fullerene graph from the view (deep copies)
        FullereneGraph fg_view(view);
        EXPECT_TRUE(fg_view.owns_memory());

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

        // Build a Graph view, then construct Polyhedron (deep copies)
        Graph g_view(N, dmax, vals, degs);
        Polyhedron poly(static_cast<const PlanarGraphView&>(PlanarGraph(g_view)), std::vector<coord3d>(pts.begin(), pts.end()), 6);

        EXPECT_TRUE(poly.owns_memory());
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

    // Polyhedron now always deep-copies, so writes go to owned storage.
    // Verify that writing through the points span updates owned_points.
    Graph g_view(N, dmax,
                 std::span<node_t>(values),
                 std::span<uint8_t>(deg));
    Polyhedron poly(static_cast<const PlanarGraphView&>(PlanarGraph(g_view)),
                    std::vector<coord3d>(pts.begin(), pts.end()), 6);

    poly.points[0] = coord3d(1.0, 2.0, 3.0);

    EXPECT_DOUBLE_EQ(poly.owned_points[0][0], 1.0);
    EXPECT_DOUBLE_EQ(poly.owned_points[0][1], 2.0);
    EXPECT_DOUBLE_EQ(poly.owned_points[0][2], 3.0);
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
    EXPECT_TRUE(T_view.owns_memory());  // deep-copies from view
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

// ---------------------------------------------------------------------------
// Phase 4: Batch<V> / BatchView<V>
// ---------------------------------------------------------------------------
#include <fullerenes/batch/batch.hh>

TEST(Batch, CubicGraphBatchPushBackAndIndex) {
    Graph G20 = FullereneGraph::C20();
    ASSERT_EQ(G20.N, 20);
    CubicGraph CG(G20);

    batch::Batch<CubicGraphView> B(CG.N, /*capacity=*/4, CG.dmax);
    EXPECT_EQ(B.size(), 0);
    EXPECT_EQ(B.capacity(), 4);
    EXPECT_EQ(B.N(), 20);
    EXPECT_EQ(B.dmax(), 3);

    B.push_back(CG);
    B.push_back(CG);
    EXPECT_EQ(B.size(), 2);

    auto v = B.view();
    EXPECT_EQ(v.size(), 2);

    for (std::size_t i = 0; i < 2u; ++i) {
        CubicGraphView entry = v[i];
        EXPECT_EQ(entry.N, 20);
        EXPECT_EQ(entry.dmax, 3);
        ASSERT_EQ(entry.neighbours.size(), (size_t)entry.N * entry.dmax);
        // Same adjacency as the source.
        for (int u = 0; u < entry.N; ++u) {
            ASSERT_EQ(entry.deg[u], CG.deg[u]);
            for (int j = 0; j < entry.deg[u]; ++j)
                EXPECT_EQ(entry.neighbours[u * entry.dmax + j],
                          CG.neighbours [u * CG.dmax     + j]);
        }
    }
}

TEST(Batch, PolyhedronBatchPreservesCoordinates) {
    Polyhedron P = Polyhedron::C20();
    ASSERT_EQ(P.N, 20);

    batch::Batch<PolyhedronView<double>> B(P.N, 2, P.dmax);
    B.push_back(P);
    ASSERT_EQ(B.size(), 1);

    PolyhedronView<double> entry = B[0];
    ASSERT_EQ(entry.N, 20);
    ASSERT_EQ(entry.points.size(), 20u);
    for (int u = 0; u < entry.N; ++u)
        for (int k = 0; k < 3; ++k)
            EXPECT_DOUBLE_EQ(entry.points[u][k], P.points[u][k]);
}

TEST(Batch, BatchViewSliceAndIterate) {
    Graph G20 = FullereneGraph::C20();
    CubicGraph CG(G20);

    batch::Batch<CubicGraphView> B(CG.N, 5, CG.dmax);
    for (int i = 0; i < 5; ++i) B.push_back(CG);

    auto all = B.view();
    ASSERT_EQ(all.size(), 5);

    auto mid = all.slice(1, 3);   // entries [1,4)
    EXPECT_EQ(mid.size(), 3);
    EXPECT_EQ(mid.N(),    CG.N);

    int seen = 0;
    for (auto v : mid) {
        EXPECT_EQ(v.N, CG.N);
        ++seen;
    }
    EXPECT_EQ(seen, 3);
}

// ---------------------------------------------------------------------------
// Phase 5: BatchState and LayoutScratch<T>
// ---------------------------------------------------------------------------
#include <fullerenes/batch/batch_state.hh>
#include <fullerenes/batch/layout_scratch.hh>

TEST(BatchState, PushBackAndView) {
    batch::BatchState S(4);
    EXPECT_EQ(S.size(), 0);
    EXPECT_EQ(S.capacity(), 4);

    int i0 = S.push_back(/*id=*/1001, StatusFlag{StatusEnum::EMPTY}, /*iter=*/0, /*vi=*/0);
    int i1 = S.push_back(/*id=*/1002, StatusFlag{StatusEnum::CUBIC_INITIALIZED}, /*iter=*/7, /*vi=*/1);
    EXPECT_EQ(i0, 0);
    EXPECT_EQ(i1, 1);
    EXPECT_EQ(S.size(), 2);

    auto v = S.view();
    ASSERT_EQ(v.size(), 2);
    EXPECT_EQ(v.id[0], 1001u);
    EXPECT_EQ(v.id[1], 1002u);
    EXPECT_EQ(v.iteration[1], 7);
    EXPECT_EQ(v.valid_index[1], 1);
}

TEST(BatchState, SliceAndMutateThroughView) {
    batch::BatchState S(5);
    for (int i = 0; i < 5; ++i) S.push_back(uint64_t(100 + i), StatusFlag{}, i, i);

    auto mid = S.slice(1, 3);
    EXPECT_EQ(mid.size(), 3);
    EXPECT_EQ(mid.id[0], 101u);
    EXPECT_EQ(mid.id[2], 103u);

    // Write-through: mid.iteration is a span into S's storage.
    for (auto& it : mid.iteration) it = 99;
    EXPECT_EQ(S.iteration()[0], 0);
    EXPECT_EQ(S.iteration()[1], 99);
    EXPECT_EQ(S.iteration()[3], 99);
    EXPECT_EQ(S.iteration()[4], 4);
}

TEST(LayoutScratch, SizingAndPerEntryAccess) {
    constexpr int Nf = 12;
    batch::LayoutScratch<double> L(Nf, /*capacity=*/3);
    L.resize(3);
    EXPECT_EQ(L.size(), 3);
    EXPECT_EQ(L.per_entry(), Nf);

    // Write a sentinel into entry 1, verify entries 0/2 are untouched.
    auto e1 = L[1];
    ASSERT_EQ(int(e1.size()), Nf);
    for (int k = 0; k < Nf; ++k) e1[k] = coord2<double>(double(k), double(-k));

    auto e0 = L[0];
    auto e2 = L[2];
    for (int k = 0; k < Nf; ++k) {
        EXPECT_EQ(L[1][k].first,  double(k));
        EXPECT_EQ(L[1][k].second, double(-k));
    }
    // untouched entries default-initialized to (0,0) by BatchAlloc
    EXPECT_EQ(e0[0].first,  0.0);
    EXPECT_EQ(e2[Nf-1].first, 0.0);
}

TEST(LayoutScratch, SliceSharesStorage) {
    constexpr int Nf = 5;
    batch::LayoutScratch<float> L(Nf, 4);
    L.resize(4);

    auto v = L.view();
    for (int i = 0; i < 4; ++i)
        for (int k = 0; k < Nf; ++k)
            v[i][k] = coord2<float>(float(i), float(k));

    auto sub = L.slice(1, 2);          // entries [1,3)
    EXPECT_EQ(sub.size(), 2);
    EXPECT_EQ(sub[0][0].first, 1.0f);
    EXPECT_EQ(sub[1][Nf-1].second, float(Nf - 1));
}

TEST(Batch, WritesThroughEntryViewUpdateStorage) {
    Graph G20 = FullereneGraph::C20();
    CubicGraph CG(G20);

    batch::Batch<CubicGraphView> B(CG.N, 2, CG.dmax);
    B.push_back(CG);
    B.push_back(CG);

    // Mutate entry 0 through the view, then re-read through the view.
    {
        auto v = B.view();
        CubicGraphView e0 = v[0];
        e0.neighbours[0] = 42;  // clobber one entry
    }
    {
        auto v = B.view();
        CubicGraphView e0 = v[0];
        EXPECT_EQ(e0.neighbours[0], 42);
        CubicGraphView e1 = v[1];
        EXPECT_NE(e1.neighbours[0], 42);  // untouched
    }
}

