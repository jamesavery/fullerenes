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
#include "fullerenes/spiral.hh"

// --- Forward declarations for owned types ---
// Type aliases (Graph = Owned<GraphView>, etc.) are deferred until
// the old class hierarchy is fully replaced. During migration, the
// old classes inherit from these views.
struct Graph;

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
// Forward declarations for owned return types.
class PlanarGraph;

struct PlanarGraphView : GraphView {
    using GraphView::GraphView;
    static constexpr uint8_t default_dmax = 6;

    typedef spiral_nomenclature::construction_scheme_t construction_scheme_t;

    // --- Topology queries ---
    bool is_a_fullerene(bool verbose=true) const;
    bool is_cubic() const;
    bool is_triangulation() const;
    bool is_cut_vertex(node_t v) const;

    // --- Face computation ---
    vector<face_t> compute_faces_oriented(int Fmax=INT_MAX) const;
    vector<face_t> compute_faces(unsigned int Fmax=INT_MAX) const;
    face_t get_face_oriented(const arc_t& e, int Fmax=INT_MAX) const;
    arc_t get_face_representation(arc_t e, int Fmax=INT_MAX) const;
    vector<arc_t> compute_face_representations(int Fmax=INT_MAX) const;

    int face_size(node_t u, node_t v) const {
        int d = 1;
        node_t u0 = u;
        while (v != u0) {
            node_t w = v;
            v = next_on_face(u, v);
            u = w;
            d++;
        }
        return d;
    }

    // --- Dual and derived graphs (return owned types) ---
    PlanarGraph dual_graph(unsigned int Fmax=INT_MAX) const;
    PlanarGraph leapfrog_dual() const;
    PlanarGraph enveloping_triangulation(construction_scheme_t& scheme) const;
    PlanarGraph enveloping_triangulation(const construction_scheme_t& scheme) const;

    // --- Combinatorics ---
    size_t count_perfect_matchings() const;
    vector<node_t> vertex_numbers(vector<vector<node_t>>& perms, const vector<node_t>& loc) const;

    // --- Triangulation helpers ---
    vector<tri_t> triangulation(int face_max=INT_MAX) const;
    vector<tri_t> triangulation(const vector<face_t>& faces) const;
    vector<tri_t> centroid_triangulation(const vector<face_t>& faces) const;
    vector<tri_t>& orient_triangulation(vector<tri_t>& tris) const;

    // --- Layout ---
    vector<coord2d> tutte_layout(node_t s=0, node_t t=-1, node_t r=-1, unsigned int face_max=INT_MAX) const;
    vector<coord2d> tutte_layout(const face_t& outer_face) const;
    vector<coord2d> tutte_layout_iterative(const face_t& outer_face, const vector<coord2d>& outer_coords) const;
    vector<coord2d> tutte_layout_direct(const face_t& outer_face, const vector<coord2d>& outer_coords) const;

    // --- Geometry ---
    vector<coord3d> zero_order_geometry(double scalerad=4) const;
};

// ---------------------------------------------------------------------------
// CubicGraphView: 3-regular planar graph.
// ---------------------------------------------------------------------------
// Forward declarations for return types.
class Triangulation;

struct CubicGraphView : PlanarGraphView {
    using PlanarGraphView::PlanarGraphView;
    static constexpr uint8_t default_dmax = 3;

    // --- Spiral methods ---
    bool get_spiral_from_cg(node_t f1, node_t f2, node_t f3,
                            vector<int>& spiral, jumplist_t& jumps,
                            bool general=true) const;
    bool get_spiral_from_cg(vector<int>& spiral, jumplist_t& jumps,
                            bool canonical=true, bool general=true,
                            bool pentagon_start=true) const;

    vector<node_t> vertex_numbers(const Triangulation& T,
                                  const vector<vector<node_t>>& perm,
                                  const vector<node_t>& loc) const;
};

// ---------------------------------------------------------------------------
// FullereneGraphView: fullerene graph (3-regular, 12 pentagons, rest hex).
// ---------------------------------------------------------------------------
// Forward declaration for FullereneGraph return types.
class FullereneGraph;

struct FullereneGraphView : CubicGraphView {
    using CubicGraphView::CubicGraphView;
    static constexpr uint8_t default_dmax = 3;

    // --- Fullerene-specific methods ---
    FullereneGraph halma_fullerene(int n, bool do_layout=false) const;
    FullereneGraph leapfrog_fullerene(bool do_layout=false) const;
    FullereneGraph GCtransform(unsigned k=1, unsigned l=0) const;

    bool get_rspi_from_fg(node_t f1, node_t f2, node_t f3,
                          vector<int>& rspi, jumplist_t& jumps,
                          bool general=true) const;
    bool get_rspi_from_fg(vector<int>& rspi, jumplist_t& jumps,
                          bool general=true, bool pentagon_start=true) const;

    matrix<int> pentagon_distance_mtx() const;
    vector<coord3d> zero_order_geometry(double scalerad=4) const;
    vector<coord3d> optimized_geometry(std::span<const coord3d> initial_geometry,
                                       int opt_method=3, double ftol=1e-12) const;
};

// ---------------------------------------------------------------------------
// TriangulationView: planar triangulation (max degree 6).
// ---------------------------------------------------------------------------
struct TriangulationView : PlanarGraphView {
    using PlanarGraphView::PlanarGraphView;
    static constexpr uint8_t default_dmax = 6;

    // --- Nested types ---
    struct simple_geodesic {
        Eisenstein g;
        int axis;
        simple_geodesic(int a, int b=0, int axis=0) : g(a,b), axis(axis) {}
    };

    struct geodesic {
        vector<Eisenstein> g;
        double d;
        int axis;
    };

    // --- Triangulation methods ---
    PlanarGraph dual_graph() const;
    vector<face_t> cubic_faces() const;
    unordered_map<arc_t,arc_t> arc_translation() const;

    size_t max_degree() const {
        size_t max_deg = 0;
        for (node_t u=0; u<N; u++) max_deg = std::max(max_deg, (size_t)degree(u));
        return max_deg;
    }

    vector<uint8_t> n_degrees() const {
        vector<uint8_t> nd(max_degree(),0);
        for (node_t u=0; u<N; u++) nd[degree(u)-1]++;
        return nd;
    }

    PlanarGraph inverse_leapfrog_dual() const;
    pair<node_t,node_t> adjacent_tris(const arc_t& e) const;
    vector<tri_t> compute_faces_oriented() const;

    Triangulation GCtransform(unsigned k=1, unsigned l=0) const;
    Triangulation halma_transform(int m, vector<map<edge_t,node_t>>* face_grids = nullptr) const;

    // --- Spiral methods ---
    bool get_spiral_implementation(node_t f1, node_t f2, node_t f3,
                                   vector<int>& v, jumplist_t& j,
                                   vector<node_t>& permutation, bool general=true,
                                   const vector<int>& S0=vector<int>(),
                                   const jumplist_t& J0=jumplist_t()) const;
    bool get_spiral(node_t f1, node_t f2, node_t f3,
                    vector<int>& v, jumplist_t& j,
                    vector<node_t>& permutation, bool general=true) const;
    bool get_spiral(vector<int>& v, jumplist_t& j,
                    vector<vector<node_t>>& permutations,
                    bool only_rarest_special=true, bool general=true, bool CW_only=false) const;
    bool get_spiral(vector<int>& v, jumplist_t& j,
                    bool rarest_start=true, bool general=true, bool CW_only=false) const;
    general_spiral get_general_spiral(bool rarest_start=true, bool CW_only=false) const;
    void get_all_spirals(vector<vector<int>>& spirals, vector<jumplist_t>& jumps,
                         vector<vector<node_t>>& permutations,
                         bool only_special=false, bool general=false) const;

    void symmetry_information(int N_generators, Graph& coxeter_diagram, vector<int>& coxeter_labels) const;
    vector<node_t> vertex_numbers(vector<vector<node_t>>& perms, const vector<node_t>& loc) const;

    vector<tri_t> triangles() const { return compute_faces_oriented(); }
    int n_triangles() const { return 2*N - 4; }

    // --- Geodesic methods ---
    matrix<int> pentagon_distance_mtx() const;
    matrix<int> simple_square_surface_distances(vector<node_t> only_nodes={}, bool calculate_self_geodesics=false) const;
    matrix<double> surface_distances(vector<node_t> only_nodes={}, bool calculate_self_geodesics=false) const;
    matrix<geodesic> surface_geodesics(vector<node_t> only_nodes={}, bool calculate_self_geodesics=false) const;
    matrix<simple_geodesic> simple_geodesics(vector<node_t> only_nodes={}, bool calculate_self_geodesics=false) const;

    node_t end_of_the_line(node_t u0, int i, int a, int b) const;
    vector<vector<node_t>> quads_of_the_line(node_t u0, int i, int a, int b) const;

    Triangulation sort_nodes() const;
};

// ---------------------------------------------------------------------------
// FullereneDualView: dual of a fullerene (triangulation with 12 degree-5
// vertices, rest degree 6).
// ---------------------------------------------------------------------------
struct FullereneDualView : TriangulationView {
    using TriangulationView::TriangulationView;
    static constexpr uint8_t default_dmax = 6;

    bool get_rspi(node_t f1, node_t f2, node_t f3,
                  vector<int>& r, jumplist_t& j, bool general=true) const;
    bool get_rspi(vector<int>& r, jumplist_t& j,
                  bool general=true, bool pentagon_start=true) const;
    general_spiral get_rspi(bool rarest_start=true) const;

    spiral_nomenclature name(bool rarest_start=true) const;
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
