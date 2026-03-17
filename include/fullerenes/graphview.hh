#pragma once

// View hierarchy for the graph class system.
//
// All view types are trivially copyable (TC) -- just spans + scalars.
// This enables zero-copy interop with SYCL kernels, numpy arrays, and
// batch GPU pipelines.
//
// The view hierarchy provides full polymorphism:
//   GraphView -> PlanarGraphView -> CubicGraphView -> FullereneGraphView
//   GraphView -> PlanarGraphView -> TriangulationView -> FullereneDualView
//
// Geometry views (PolyhedronView, DeltahedronView) add a span<coord3d>
// for vertex coordinates.
//
// Algorithms are compiled once on the view types. Owned<View> (owned.hh)
// adds vector storage for any view type.

#include "fullerenes/dense_graph.hh"
#include "fullerenes/geometry.hh"
#include "fullerenes/auxiliary.hh"
#include "fullerenes/matrix.hh"

// --- Forward declarations for owned types ---
// Type aliases (Graph = Owned<GraphView>, etc.) are deferred until
// the old class hierarchy is fully replaced. During migration, the
// old classes inherit from these views.

// ---------------------------------------------------------------------------
// GraphView: general undirected/directed graph with adjacency lists.
// ---------------------------------------------------------------------------
struct GraphView : Spanify::RSRAdjacencyView<node_t> {
    using RSRAdjacencyView::RSRAdjacencyView;
    static constexpr uint8_t default_dmax = 10;
    GraphView() = default;

    // --- Edge operations (mutate through spans) ---
    bool insert_edge(const arc_t& e, const node_t suc_uv=-1, const node_t suc_vu=-1);
    bool remove_edge(const edge_t& e);
    bool edge_exists(const edge_t& e) const;
    void flip_all_orientations();

    // --- Vertex-pair navigation (O(degree) per call) ---
    // Bring arc-index overloads from base into scope (otherwise hidden
    // by the vertex-pair overloads declared here).
    using RSRAdjacencyView::next;
    using RSRAdjacencyView::prev;
    using RSRAdjacencyView::next_on_face;

    int  arc_ix(node_t u, node_t v) const;
    node_t next(node_t u, node_t v) const;
    node_t prev(node_t u, node_t v) const;
    node_t next_on_face(node_t u, node_t v) const;
    node_t prev_on_face(node_t u, node_t v) const;

    // --- Topology queries ---
    bool is_consistently_oriented() const;
    bool adjacency_is_symmetric() const;
    bool has_separating_triangles() const;

    // --- Connectivity and shortest paths ---
    bool is_connected(const set<node_t>& subgraph = set<node_t>()) const;
    void single_source_shortest_paths(node_t source, int *distances, size_t max_depth = INT_MAX) const;
    matrix<int> all_pairs_shortest_paths(const vector<node_t>& V,
                                         const unsigned int max_depth = INT_MAX) const;
    matrix<int> all_pairs_shortest_paths(const unsigned int max_depth = INT_MAX) const;
    vector<vector<node_t>> connected_components() const;

    vector<node_t> shortest_cycle(node_t s, const int max_depth) const;
    vector<node_t> shortest_cycle(const vector<node_t>& prefix, const int max_depth) const;
    vector<int> multiple_source_shortest_paths(const vector<node_t>& sources, const unsigned int max_depth=INT_MAX) const;

    // --- Hamiltonian paths ---
    int hamiltonian_count() const;
    int hamiltonian_count(node_t current_node, vector<bool>& used_edges, vector<bool>& used_nodes, vector<node_t>& path, const vector<int>& distances) const;

    // --- Geometry helpers ---
    coord2d centre2d(const vector<coord2d>& layout) const;
    coord3d centre3d(std::span<const coord3d> layout) const;

    // --- Degree and edge queries ---
    int max_degree() const;
    vector<edge_t>  undirected_edges() const;
    vector<arc_t> directed_edges() const;
    size_t count_edges() const;
};

// ---------------------------------------------------------------------------
// PlanarGraphView: planar graph with oriented embedding.
// ---------------------------------------------------------------------------
struct PlanarGraphView : GraphView {
    using GraphView::GraphView;
    static constexpr uint8_t default_dmax = 6;

    // Algorithm methods will be migrated here from PlanarGraph when
    // PlanarGraph switches from inheriting Graph to inheriting PlanarGraphView.
    // For now, PlanarGraph : Graph keeps the method definitions.
};

// ---------------------------------------------------------------------------
// CubicGraphView: 3-regular planar graph.
// ---------------------------------------------------------------------------
struct CubicGraphView : PlanarGraphView {
    using PlanarGraphView::PlanarGraphView;
    static constexpr uint8_t default_dmax = 3;
};

// ---------------------------------------------------------------------------
// FullereneGraphView: fullerene graph (3-regular, 12 pentagons, rest hex).
// ---------------------------------------------------------------------------
struct FullereneGraphView : CubicGraphView {
    using CubicGraphView::CubicGraphView;
    static constexpr uint8_t default_dmax = 3;
};

// ---------------------------------------------------------------------------
// TriangulationView: planar triangulation (max degree 6).
// ---------------------------------------------------------------------------
struct TriangulationView : PlanarGraphView {
    using PlanarGraphView::PlanarGraphView;
    static constexpr uint8_t default_dmax = 6;
};

// ---------------------------------------------------------------------------
// FullereneDualView: dual of a fullerene (triangulation with 12 degree-5
// vertices, rest degree 6).
// ---------------------------------------------------------------------------
struct FullereneDualView : TriangulationView {
    using TriangulationView::TriangulationView;
    static constexpr uint8_t default_dmax = 6;
};

// ---------------------------------------------------------------------------
// PolyhedronView: planar graph with 3D vertex coordinates.
// ---------------------------------------------------------------------------
struct PolyhedronView : PlanarGraphView {
    std::span<coord3d> points;
    static constexpr uint8_t default_dmax = 10;

    PolyhedronView() = default;

    // Construct from adjacency view + coordinate span.
    PolyhedronView(const PlanarGraphView& g, std::span<coord3d> pts)
        : PlanarGraphView(g), points(pts) {}

    // Full view constructor (adjacency + coordinates).
    PolyhedronView(node_t N, int dmax,
                   std::span<node_t> neighbours, std::span<uint8_t> deg,
                   std::span<coord3d> pts,
                   std::span<uint8_t> twin = {})
        : PlanarGraphView(N, dmax, neighbours, deg, twin), points(pts) {}
};

// ---------------------------------------------------------------------------
// DeltahedronView: triangulation with 3D vertex coordinates (equilateral
// triangle embedding).
// ---------------------------------------------------------------------------
struct DeltahedronView : TriangulationView {
    std::span<coord3d> points;
    static constexpr uint8_t default_dmax = 6;

    DeltahedronView() = default;

    // Construct from adjacency view + coordinate span.
    DeltahedronView(const TriangulationView& t, std::span<coord3d> pts)
        : TriangulationView(t), points(pts) {}

    // Full view constructor (adjacency + coordinates).
    DeltahedronView(node_t N, int dmax,
                    std::span<node_t> neighbours, std::span<uint8_t> deg,
                    std::span<coord3d> pts,
                    std::span<uint8_t> twin = {})
        : TriangulationView(N, dmax, neighbours, deg, twin), points(pts) {}
};

// ---------------------------------------------------------------------------
// TC verification: all view types must be trivially copyable.
// ---------------------------------------------------------------------------
static_assert(std::is_trivially_copyable_v<GraphView>,
    "GraphView must be trivially copyable");
static_assert(std::is_trivially_copyable_v<PlanarGraphView>,
    "PlanarGraphView must be trivially copyable");
static_assert(std::is_trivially_copyable_v<CubicGraphView>,
    "CubicGraphView must be trivially copyable");
static_assert(std::is_trivially_copyable_v<FullereneGraphView>,
    "FullereneGraphView must be trivially copyable");
static_assert(std::is_trivially_copyable_v<TriangulationView>,
    "TriangulationView must be trivially copyable");
static_assert(std::is_trivially_copyable_v<FullereneDualView>,
    "FullereneDualView must be trivially copyable");
static_assert(std::is_trivially_copyable_v<PolyhedronView>,
    "PolyhedronView must be trivially copyable");
static_assert(std::is_trivially_copyable_v<DeltahedronView>,
    "DeltahedronView must be trivially copyable");

// Verify correct hierarchy relationships
static_assert(std::is_base_of_v<GraphView, PlanarGraphView>);
static_assert(std::is_base_of_v<PlanarGraphView, CubicGraphView>);
static_assert(std::is_base_of_v<CubicGraphView, FullereneGraphView>);
static_assert(std::is_base_of_v<PlanarGraphView, TriangulationView>);
static_assert(std::is_base_of_v<TriangulationView, FullereneDualView>);
static_assert(std::is_base_of_v<PlanarGraphView, PolyhedronView>);
static_assert(std::is_base_of_v<TriangulationView, DeltahedronView>);
