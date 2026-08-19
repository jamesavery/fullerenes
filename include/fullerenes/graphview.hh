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

#include <stdexcept>
#include <array>
#include <cmath>

#include "fullerenes/dense_graph.hh"
#include "fullerenes/geometry.hh"
#include "fullerenes/auxiliary.hh"
#include "fullerenes/matrix.hh"
#include "fullerenes/spiral.hh"
#include "fullerenes/permutation.hh"

// --- Forward declarations for owned types ---
// Type aliases (Graph = Owned<GraphView>, etc.) are deferred until
// the old class hierarchy is fully replaced. During migration, the
// old classes inherit from these views.
struct Graph;

// ---------------------------------------------------------------------------
// The orientable surface a rotation system embeds a graph in.
//
// A rotation system -- the cyclic order of each vertex's neighbour row -- always
// embeds its graph in SOME orientable surface (Edmonds; Heffter-Ringel-Youngs).
// With adjacency symmetric, the arc-successor map
//
//     phi(u->v) = v -> next(v, u)          ("the arc after u->v around its face")
//
// is a PERMUTATION of the arc set, so its orbits always partition the arcs and
// always ARE the faces of an embedding -- for a valid rotation system, a
// corrupted one, or a randomly shuffled one alike.  A bare rotation system
// therefore has nothing to violate, and "is it consistently oriented?" is a
// question only once the caller says WHICH surface the object claims to be.
// Euler then reads the realised genus off the face count: for one connected
// component
//
//     faces == E - N + 2 - 2*genus.
//
// So transposing two entries of one C20-dual row does not produce an
// inconsistent orientation -- it produces a perfectly consistent orientation OF
// A TORUS, 10 faces where 12 were claimed.  `faces` and `genus` below are what
// the rotation system actually realises; `code` compares that with the claim.
// ---------------------------------------------------------------------------
struct OrientedSurface {
    enum class Code {
        Ok,                    // the rotation system realises the surface claimed
        Degenerate,            // nothing to embed: N == 0 or E == 0
        AsymmetricAdjacency,   // some arc u->v has no reverse v->u, so phi is not a
                               // permutation and the graph has no faces at all
        GenusMismatch,         // consistently oriented, but of a different surface
    };
    // Only Ok and GenusMismatch leave the numbers meaningful; the two structural
    // codes describe an object with no faces, and a default-constructed
    // OrientedSurface is the Degenerate failure, not a blank success.
    Code code  = Code::Degenerate;
    int  faces = 0;    // the phi-orbit count
    int  genus = -1;   // (E - N + 2 - faces)/2 -- the surface actually realised
};

// ---------------------------------------------------------------------------
// The @pre "consistently oriented" of an operation that walks faces, violated by
// its caller -- the exception channel, not a modeled outcome (style-failures.md):
// such an operation has no answer to give, and every one of them is deep enough
// that only the outermost per-isomer harness can decide what to do.
//
// It carries the OrientedSurface rather than only a sentence, so a handler can
// ASK what it was handed -- and so the sentence can say it.  ("Not oriented" on
// its own names no offender and no surface, which is why the asserts it replaced
// were worth nothing to whoever read the crash.)
// ---------------------------------------------------------------------------
struct unoriented_surface_error : public std::logic_error {
    OrientedSurface surface;          // the surface the rotation system realises
    int             expected_genus;   // the one the operation needed
    unoriented_surface_error(const std::string& what, const OrientedSurface& s, int g)
        : std::logic_error(what), surface(s), expected_genus(g) {}
};

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

    // Inline so device code can call them: these five are the vocabulary for
    // walking a rotation system.
    //
    // Total -- -1 for "not an arc", never an exception.  A caller for which a
    // missing arc is a contract violation raises that itself.
    int arc_ix(node_t u, node_t v) const { return find(u, v); }

    // Successor to v in the oriented neighbours of u.
    node_t next(node_t u, node_t v) const {
        const std::span<const node_t> nu = (*this)[u];
        const int j = arc_ix(u, v);
        return j >= 0 ? nu[(j + 1) % nu.size()] : -1;
    }

    // Predecessor to v in the oriented neighbours of u.
    node_t prev(node_t u, node_t v) const {
        const std::span<const node_t> nu = (*this)[u];
        const int j = arc_ix(u, v);
        return j >= 0 ? nu[(j - 1 + nu.size()) % nu.size()] : -1;
    }

    // Successor / predecessor to v in the face containing the arc u->v.
    node_t next_on_face(node_t u, node_t v) const { return prev(v, u); }
    node_t prev_on_face(node_t u, node_t v) const { return next(v, u); }

    // --- Topology queries ---

    // Does this rotation system embed the graph in the surface the caller
    // claims -- non-degenerate, symmetric, and of genus `genus`?  At the
    // default genus 0 this is "consumers can work with this planar graph
    // safely", so an edgeless or empty graph FAILS rather than passing
    // vacuously.
    //
    // TWO PRECONDITIONS, DELIBERATELY NOT CHECKED:
    //  - CONNECTED.  Euler's relation below is per component; paying
    //    connected_components() here would put an O(N) sweep on a hot path
    //    whose entire point is to be TOLD the topology.  Use
    //    component_surfaces() when connectivity is unknown.
    //  - SIMPLE (no parallel arcs, no self-loops).  Face tracing goes through
    //    find(v,u), which resolves only the FIRST occurrence of a neighbour in
    //    a row, so parallel arcs all resolve to one slot and phi stops being a
    //    permutation.  Not silent: the walk THROWS std::logic_error the moment
    //    it steps onto an already-visited arc, which under a permutation cannot
    //    happen.  Whether the Graph hierarchy should admit parallel arcs and
    //    loops at all is an open question (the Delta Complex type exists for
    //    that): see claude-projects/curvature-flow/refactor-debt.md, entry
    //    2026-08-09-graph-hierarchy-parallel-arcs-and-loops.
    //
    // Undecidable here, and not attempted: global handedness.  The mirror
    // rotation system is equally valid, so CW vs CCW needs coordinates
    // (orient_polyhedron_neighbours' signed-volume test resolves it on
    // PolyhedronView).  For a 3-connected planar graph nothing else is left --
    // by Steinitz and Whitney it has exactly two genus-0 rotation systems,
    // sigma and its mirror -- so at genus 0 this predicate is as tight as a
    // combinatorial test can be.
    //
    // @anchor graph-is-consistently-oriented
    // @pre    connected: is_connected()
    // @post   result == (oriented_surface(genus).code == OrientedSurface::Code::Ok)
    // @time   O(E_dmax)
    bool is_consistently_oriented(int genus = 0) const;

    // is_consistently_oriented(genus) with the surface attached: the faces the
    // rotation system actually has and the genus they imply, so a caller can
    // say WHY it failed and by how much.  Same preconditions.
    //
    // @anchor graph-oriented-surface
    // @pre    connected: is_connected()
    // @error  Degenerate          when N == 0 || count_edges() == 0
    // @error  AsymmetricAdjacency when !adjacency_is_symmetric()
    // @error  GenusMismatch       when result.genus != genus
    // @post on Ok: result.genus == genus &&
    //                  result.faces == int(count_edges()) - N + 2 - 2*genus
    // @throws std::logic_error when the graph is not simple (see above)
    // @time   O(E_dmax)
    OrientedSurface oriented_surface(int genus = 0) const;

    // One OrientedSurface per connected component, in connected_components()
    // order -- the wrapper for input whose connectivity is unknown or
    // deliberately plural, so the caller sees WHICH component is wrong rather
    // than only that one is.  `genus` is the claim per component and defaults
    // to all-zero; entries past its end are read as 0.
    //
    // Costs ONE connected_components() call, which yields the count and the
    // partition together, and no subgraph is ever built: per-component N_i and
    // E_i come off the partition, and since a phi-orbit never leaves its
    // component a single arc walk attributes every face to one.  On connected
    // input it therefore agrees with oriented_surface(genus[0]) exactly -- which
    // is the contract that lets the core predicate charge nothing for
    // connectivity.  One exception to the per-component resolution: a missing
    // reverse arc breaks the arc set rather than one component's topology, so
    // AsymmetricAdjacency is reported on EVERY component.
    //
    // @anchor graph-component-surfaces
    // @post   result.size() == connected_components().size()
    // @post   connected: implies(is_connected(), result.size() == 1)
    // @time   O(E_dmax)
    vector<OrientedSurface> component_surfaces(const vector<int>& genus = {}) const;

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

    // --- Hamiltonian cycles ---
    // Number of undirected Hamiltonian cycles, each counted once (the
    // convention of IsomerDB's ncycham).  Babic's backtracking with the
    // exclusion rule of the legacy hamilton.f HamiltonCyc, generalised to any
    // degree.  Definition staged in claude-projects/unfortran/src/hamiltonian.cc
    // until promotion.
    int64_t hamiltonian_count() const;

    // --- Geometry helpers ---
    coord2d centre2d(const vector<coord2d>& layout) const;
    coord3d centre3d(std::span<const coord3d> layout) const;

    // --- Degree and edge queries ---
    int max_degree() const;
    vector<edge_t>  undirected_edges() const;
    vector<arc_t> directed_edges() const;
    size_t count_edges() const;

    // Sort the vertex indices {0, ..., N-1} by `less` (a stable sort) and
    // return the resulting permutation pi with pi[u_old] = u_new (i.e.
    // pi maps old labels to new labels).  *this is unchanged.  `less` is
    // any callable invocable as `bool less(node_t a, node_t b)`; lambdas
    // typically capture `this` to access vertex data such as degrees.
    // (Materialise the sorted graph via apply_permutation(pi) on an owner.)
    template<typename Less>
    Permutation argsort(Less less) const {
        std::vector<int> idx(N);
        for (int i = 0; i < N; ++i) idx[i] = i;
        std::stable_sort(idx.begin(), idx.end(), less);  // idx[u_new] = u_old
        Permutation pi(N);
        for (int u_new = 0; u_new < N; ++u_new) pi[idx[u_new]] = u_new;
        return pi;
    }
};

// Enforce the "consistently oriented" precondition of `operation`: return if G's
// rotation system embeds it in the genus-`genus` surface, throw otherwise.  One
// function, so every operation that refuses an unoriented graph refuses it in
// the same sentence -- and so that sentence reports the surface it GOT
// (faces and implied genus) beside the one it wanted, which is the whole reason
// these are throws now and not asserts.
//
// @anchor graph-require-oriented-surface
// @pre    connected: G.is_connected()   (inherited from oriented_surface)
// @throws unoriented_surface_error when G.oriented_surface(genus).code != Ok
// @time   O(E_dmax)
void require_oriented_surface(const GraphView& G, const char* operation, int genus = 0);

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

    // The size of the face containing the arc u->v; 0 if the walk runs off a
    // missing arc or does not close within max_face steps, so a non-oriented
    // rotation system cannot spin forever.
    static constexpr int max_face = 64;
    int face_size(node_t u, node_t v) const {
        int d = 1;
        const node_t u0 = u;
        while (v != u0) {
            const node_t w = v;
            v = next_on_face(u, v);
            if (v < 0) return 0;
            u = w;
            if (++d > max_face) return 0;
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
    //
    // Two ways to cut polygon faces into triangles, both of which PRESERVE the
    // faces' orientation arc for arc: a triangle comes out as itself, and every
    // interior arc a face gains appears once in each direction.  So both return a
    // consistently oriented triangle surface whenever `faces` is one -- which,
    // faces() being compute_faces_oriented(), it always is.  That is a property of
    // the construction, not of a repair pass: the orientation-repair call these
    // two used to end with (orient_triangulation, deleted 2026-08-09) could only
    // ever confirm what the loop had already produced, and did it by abort()ing
    // the process when it disagreed.
    //
    // @pre     oriented: `faces` are the oriented faces of *this
    // @post    the result is a closed consistently oriented triangle surface
    // @throws  std::logic_error on a face with fewer than 3 vertices

    // Fan: face (v0,...,vk) -> the k-1 triangles (v0,vj,vj+1).  Adds no vertices,
    // so the result indexes *this's vertices alone.
    vector<tri_t> triangulation(int face_max=INT_MAX) const;
    vector<tri_t> triangulation(const vector<face_t>& faces) const;

    // Centroid: face i of more than 3 vertices -> one triangle per edge, meeting
    // at a NEW vertex N+i, its centroid; a triangle face stays one triangle.  The
    // caller supplies the coordinates of those N+i (PolyhedronView::
    // centroid_surface does), and allocates one per face whether it is used or
    // not, so face i's centroid is always vertex N+i.  Unlike the fan it does not
    // depend on which vertex the face list starts at.
    vector<tri_t> centroid_triangulation(const vector<face_t>& faces) const;

    // --- Layout ---
    vector<coord2d> tutte_layout(node_t s=0, node_t t=-1, node_t r=-1, unsigned int face_max=INT_MAX) const;
    vector<coord2d> tutte_layout(const face_t& outer_face) const;
    vector<coord2d> tutte_layout_iterative(const face_t& outer_face, const vector<coord2d>& outer_coords) const;
    vector<coord2d> tutte_layout_direct(const face_t& outer_face, const vector<coord2d>& outer_coords) const;

    // --- Geometry ---
    // Tutte layout -> sphere shells.  One specific construction, available to
    // callers that want it by name (e.g. to compare against the
    // eisenstein-paint methods); subtypes inherit it.
    // @pre  the rotation system is consistently oriented (otherwise
    //           tutte_layout falls back to shortest_cycle for the outer face)
    // @post returns N points on a sphere of average bond length scalerad*1.5,
    //           centred on their centroid
    vector<coord3d> tutte_sphere_geometry(double scalerad=4) const;

    // The best initial geometry for this type.  Subtypes point it at the
    // tightest construction their structure allows -- FullereneGraphView at
    // eisenstein_paint_geometry -- so which method runs, and what the argument
    // means, follow the type.  Callers wanting a specific method name it.
    vector<coord3d> zero_order_geometry(double scalerad=4) const {
        return tutte_sphere_geometry(scalerad);
    }
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
    using InitialGeometry = vector<coord3d> (FullereneGraphView::*)(double) const;
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

    vector<coord3d> zero_order_geometry(double bond_length=1.44) const { return eisenstein_paint_geometry(bond_length); }

    // Warm-start 3D geometry from the cubic-metric Eisenstein paint
    // (eisenstein_paint::cubic_geometry): AlexandrovIDTCubic realizes the
    // cubic polyhedral metric's 20..60 pentagon-incident vertices
    // exactly as the unique convex embedding realizing Alexandrov's theorem. 
    // The exact integer Eisenstein paint interpolates the rest over the realized
    // dual polytope's flat cells and projects onto the coarse embedding's faces.  
    // Throws std::runtime_error if the paint pipeline fails (must never happen).
    vector<coord3d> eisenstein_paint_geometry(double bond_length=1.44) const;
    vector<coord3d> tutte_sphere_geometry(double scalerad=4) const;

    // Force-field geometry optimization (Wu / extended Wu, harmonic) --
    // the production optimizer for cubic fullerene graphs.  opt_method
    // keeps the legacy variant numbering: 1/2 Wu, 3/4 extended Wu, 5/6
    // softened extended Wu; even variants add a Coulomb origin repulsion
    // phase.  Default 3 = hard-harmonic ExtWu.  Implemented by
    // wu::forcefield + wu::optimize (wu_forcefield.hh); replaces the
    // legacy Fortran SA_OptFF path, against which it is cross-validated
    // pointwise and end-to-end by the parity tests in
    // claude-projects/unfortran/tests (which link the Fortran archives
    // directly -- libfullerenes itself is Fortran-free).
    vector<coord3d> optimized_geometry(std::span<const coord3d> initial_geometry,
                                       int opt_method=3, double ftol=1e-12) const;

    vector<coord3d> optimized_geometry(
        const InitialGeometry f_initial = &FullereneGraphView::eisenstein_paint_geometry, 
        double initial_scalar = 1.44,
        int opt_method=3, double ftol=1e-12) const {
            vector<coord3d> initial_geometry = (this->*f_initial)(1.44);          
            return optimized_geometry(initial_geometry, opt_method, ftol);
        }
};

// ---------------------------------------------------------------------------
// TriangulationView: planar triangulation (max degree 6).
// ---------------------------------------------------------------------------
struct TriangulationView : PlanarGraphView {
    using PlanarGraphView::PlanarGraphView;
    static constexpr uint8_t default_dmax = 6;

    // --- Nested types ---
    // A simple geodesic u -> v: the Eisenstein displacement (a, b) (with
    // a, b >= 0) walked from u, starting along the `axis`-th out-arc at u
    // (the a-axis). Reaches v = end_of_the_line(u, axis, a, b); squared
    // length is g.norm2() = a^2 + a*b + b^2.
    struct simple_geodesic {
        Eisenstein g;
        int axis;
        simple_geodesic(int a=0, int b=0, int axis=0) : g(a,b), axis(axis) {}
    };

    // A general geodesic u -> v: a sequence of simple geodesics
    // u -> K_1 -> K_2 -> ... -> v, broken at intermediate cones.
    struct geodesic {
        vector<simple_geodesic> segments;
        geodesic() = default;
        geodesic(int) {}   // matrix<geodesic>(m, n) zero-init compatibility
    };

    // One iteration of the simple_geodesics search from a traced source: the
    // (axis,a,b) ray probed vertex v at squared distance d2; H_before is the
    // running min H(U,V) *before* this probe; improved == (d2 <= H_before) ==
    // "the witness geodesic was (re)assigned". Aggregate; emitted in the exact
    // search order (axis outer, then a, then b) so a viz can replay the min.
    struct geodesic_step {
        int    axis, a, b;
        node_t v;
        int    d2, H_before;
        bool   improved;
    };

    // The Thurston star unfolding from a source cone: the surface cut along
    // the shortest simple geodesics source -> every other cone (a star tree)
    // and developed into the Eisenstein plane.  Boundary-only and exact in
    // Z[w]: every vertex of the development is a lattice point.
    //
    // The payload is the STANDARD labelled outline (Unfolding::outline's
    // type, CW like it), so the fold/unfold pipeline consumes it directly:
    // entries alternate (source copy, source), (leaf position, cone), so the
    // even entries are the k-1 source copies, the odd entries the k-1 other
    // cones, and the i-th squared cut length is
    // (outline[2i+1].first - outline[2i].first).norm2().  The only fields
    // beyond the outline are the outcome code and the globally-shortest
    // flag, which is a property of the input's geodesics and cannot be read
    // off the picture.  Construction, proofs and sweep evidence:
    // claude-projects/visualization/STAR-UNFOLDING.md.
    struct star_unfolding {
        enum class Code { Ok, SourceNotCone, NoSimpleGeodesic, CollinearCuts };
        Code   code = Code::Ok;
        string metadata;                    // names the offender when code != Ok

        node_t source = -1;
        vector<pair<Eisenstein, node_t>> outline;   // CW 2(k-1)-gon, see above
        bool cuts_globally_shortest = true;         // false: some cut is only
                                                    // simple-shortest (degenerate)

        static star_unfolding error(Code c, string m) {
            star_unfolding r;
            r.code     = c;
            r.metadata = std::move(m);
            return r;
        }
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
    matrix<double> surface_distances(vector<node_t> only_nodes={}, bool calculate_self_geodesics=false,
                                     matrix<geodesic>* geodesics_out=nullptr) const;
    matrix<geodesic> surface_geodesics(vector<node_t> only_nodes={}, bool calculate_self_geodesics=false) const;
    // trace_u (optional): if >= 0, every (axis,a,b) ray probed from source
    // trace_u is appended to *trace_out in search order, and *M_out receives
    // its search radius M[U] = max graph distance from trace_u. Additive: with
    // trace_u<0 the run is byte-identical to the untraced search.
    matrix<simple_geodesic> simple_geodesics(vector<node_t> only_nodes={}, bool calculate_self_geodesics=false,
                                             node_t trace_u=-1, vector<geodesic_step>* trace_out=nullptr,
                                             int* M_out=nullptr) const;

    // Concatenate simple geodesics along a node-index path:
    //   [U, K_1, ..., V] -> { simple(U, K_1), simple(K_1, K_2), ..., simple(K_n, V) }.
    // Returns an empty geodesic for paths of length <= 1 (U == V or unreachable).
    static geodesic compose_simple_geodesics(const vector<int>& path,
                                             const matrix<simple_geodesic>& simple);

    // Thurston star unfolding from `source` (any vertex labelling).
    // @anchor tri-star-unfold
    // @pre  degree(source) != 6                        (else Code::SourceNotCone)
    // @post on Ok: outline is the closed CW 2(k-1)-gon of the cut development
    //       (k = cone count), alternating source copies and leaf cones, each
    //       leaf appearing exactly once; the walk closes exactly, and the
    //       (a,b) lattice shoelace equals -(2N - 4) (CW traversal covering
    //       every face once) -- both asserted.
    // @error Code::SourceNotCone    when degree(source) == 6
    // @error Code::NoSimpleGeodesic when some cone has no simple geodesic
    //       from source
    // @error Code::CollinearCuts    when two cut directions coincide at source
    // @throws std::runtime_error on closure/area failure -- a deep invariant
    //       that cannot occur for valid input (swept clean over every isomer
    //       and every source, C20-C100).
    // @ref  claude-projects/visualization/STAR-UNFOLDING.md
    star_unfolding star_unfold(node_t source) const;
    // Sweep form: shares precomputed simple_geodesics / surface_distances
    // over `cones` (the degree != 6 vertices) across many sources.
    star_unfolding star_unfold(node_t source, const vector<node_t>& cones,
                               const matrix<simple_geodesic>& simple,
                               const matrix<double>& dist) const;

    node_t end_of_the_line(node_t u0, int i, int a, int b) const;
    // coords_out (optional): per quad vertex, its integer Eisenstein lattice coord
    // (origin u0 at (0,0), axis i along (1,0)); filled in lockstep with the vertex
    // quads. Lets callers lay the unfolded strip in the plane / lift it to 3D.
    vector<vector<node_t>> quads_of_the_line(node_t u0, int i, int a, int b,
                                             vector<vector<Eisenstein>>* coords_out=nullptr) const;

    Triangulation sort_nodes() const;

    // Return the permutation that sorts cone vertices (degree != 6) before
    // flat (degree == 6) vertices, preserving original order within each
    // group.  pi[u_old] = u_new.  *this is unchanged; apply via
    // `apply_permutation(pi)` on a copy to materialise the sorted graph.
    Permutation sort_flat_last() const;
};

// ---------------------------------------------------------------------------
// FullereneDualView: dual of a fullerene (triangulation with 12 degree-5
// vertices, rest degree 6).
// ---------------------------------------------------------------------------

// The pentagon and hexagon neighbour indices of a fullerene (Raghavachari;
// Fowler & Manolopoulos, An Atlas of Fullerenes; legacy spiral.f
// DualAnalyze): P[k] = number of pentagons with k pentagon neighbours
// (k = 0..5, sum 12), H[k] = number of hexagons with k hexagon neighbours
// (k = 0..6, sum N/2 - 10).
struct NeighbourIndices {
    std::array<uint8_t,6> P{};
    std::array<uint8_t,7> H{};
};
// Np = number of pentagon-pentagon adjacencies = sum_k k P[k] / 2
// (legacy util.f IPentInd).
inline int pentagon_adjacencies(const NeighbourIndices& ni) {
    int m = 0; for (int k = 0; k < 6; k++) m += k * ni.P[k]; return m / 2;
}
// sigma_h = standard deviation of the hexagon-neighbour count over the
// hexagons with k >= kmin neighbours (legacy util.f HexInd, kmin = 0; an
// older program version summed only k >= 3).  0 when there are none.
inline double hexagon_index(const NeighbourIndices& ni, int kmin = 0) {
    long n = 0, hk = 0, hk2 = 0;
    for (int k = kmin; k < 7; k++) { n += ni.H[k]; hk += k * ni.H[k]; hk2 += k * k * ni.H[k]; }
    return n ? std::sqrt(std::fabs(double(hk2)/n - std::pow(double(hk)/n, 2))) : 0.0;
}

struct FullereneDualView : TriangulationView {
    using TriangulationView::TriangulationView;
    static constexpr uint8_t default_dmax = 6;

    // The pentagon/hexagon neighbour indices of this dual.
    NeighbourIndices neighbour_indices() const;

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
// Forward declaration for return types.
class Polyhedron;

// ---------------------------------------------------------------------------
// The moments of a mass distribution: total mass, centre of mass, and the
// CENTRAL second moment about that centre.  One shape for both mass models
// (MassModel below), so the two inertia branches are the same sentence with one
// word changed.
//
// `mass` is whatever the weighting makes it, and is the ONLY field whose
// meaning the two branches state differently: for volume_moments() -- the
// enclosed SOLID of unit density, whose mass IS its volume -- it is the enclosed
// volume, and for the point-mass sibling vertex_moments() it is the total point
// mass, i.e. the point COUNT.  (Until 2026-08-09 the type was VolumeMoments and
// the field `volume`, which told the truth for the first and not the second.)
// ---------------------------------------------------------------------------
struct MassMoments {
    // Why the fields do not describe a body.  Only Ok leaves them meaningful:
    // every other code leaves the zeroed sentinel below, and centroid {0,0,0} is
    // a perfectly plausible-looking point, so a caller must read `code` (or test
    // `mass > 0`) rather than eyeball the numbers.
    enum class Code {
        Ok,                        // the fields are the moments of a body
        Degenerate,                // no mass to take moments of: < 4 points, no faces or
                                   // zero enclosed volume (volume_moments); no points at
                                   // all (vertex_moments)
        NotPositiveSemidefinite,   // covariance has a negative or non-finite minimum
                                   // eigenvalue -- it is not a second moment, so the
                                   // input was not the mass distribution it claimed
    };
    Code     code = Code::Degenerate;
    double   mass = 0;           // the distribution's total mass: the enclosed VOLUME for
                                 // volume_moments() (the surface's orientation sign is
                                 // normalised away), the point COUNT for vertex_moments().
                                 // > 0 iff code == Ok
    coord3d  centroid{0,0,0};    // the distribution's centre of mass
    matrix3d covariance{};       // ∫ (x-centroid)(x-centroid)^T dm, about the CENTROID
                                 // (the sum over the points for vertex_moments())
};

// Exact moments of the solid enclosed by an oriented triangle surface, by
// signed-tetrahedron sums over the faces: nothing is sampled and nothing is
// vertex-weighted, so a densely triangulated patch does not move the answer.
// Per tetrahedron (a,b,c) with apex o, exactly:
//     V            = det[a b c]/6
//     ∫ x    dV    = det (a+b+c)/24
//     ∫ x x^T dV   = det M S M^T,  M = [a b c] as COLUMNS,
//                    S = diag(1/60) + offdiag(1/120)  (reference tetrahedron)
// S = (J+I)/120 with J all-ones, so M S M^T = (ss^T + aa^T + bb^T + cc^T)/120
// with s = a+b+c: three outer products instead of a 3x3 triple product.
//
// This is the KERNEL: the arithmetic over a point set and an oriented triangle
// list, with no view type in the signature -- PolyhedronView<double>::
// volume_moments() (which triangulates its polygon faces first) is one caller,
// and any triangulated view (DeltahedronView, which is not a PolyhedronView, so
// it reaches this through its own triangles()) or caller-supplied triangle list
// is another.
//
// CALLER'S OBLIGATION, stated here rather than in a formal slot because the
// library carries no predicate for it: `tris` must be a CLOSED, CONSISTENTLY
// ORIENTED triangle surface over `points` -- every directed edge appearing
// exactly once.  A surface that is not returns Code::NotPositiveSemidefinite in
// the cases that are detectable and is otherwise meaningless; the winding
// DIRECTION does not matter, only its consistency (the orientation sign is
// normalised away).  The two failure modes below are the ones the kernel
// detects: too little input to enclose anything (< 4 points, no triangles, or a
// signed volume that is not positive), and a second moment that is not one.
//
// @anchor solid-volume-moments
// @post on Ok:         result.mass > 0 && result.covariance.eigenvalues()[2] >= 0 &&
//                          (result.covariance - result.covariance.transpose()).norm() == 0
// @post on Degenerate: result.mass == 0 && result.covariance.norm() == 0
// @post on NotPositiveSemidefinite:
//                      result.mass == 0 && result.covariance.norm() == 0
MassMoments volume_moments(std::span<const coord3d> points, const vector<tri_t>& tris);

// The moments of the point-mass distribution "unit mass at every point" -- the
// sibling of volume_moments() for MassModel::Atoms, in the same triple so the
// two mass models compose identically.  `mass` carries the total mass, which
// for unit masses is points.size() -- a COUNT, not a volume, which is why the
// shared type is MassMoments; `covariance` is Σ (x-mean)(x-mean)^T about the
// point mean.  Needs no topology at all.
//
// @anchor point-mass-moments
// @post on Ok:         result.mass == (double)points.size() &&
//                          result.covariance.eigenvalues()[2] >= 0 &&
//                          (result.covariance - result.covariance.transpose()).norm() == 0
// @post on Degenerate: result.mass == 0 && result.covariance.norm() == 0
// @post on NotPositiveSemidefinite:
//                      result.mass == 0 && result.covariance.norm() == 0
MassMoments vertex_moments(std::span<const coord3d> points);

// ---------------------------------------------------------------------------
// Which mass distribution an inertia tensor / principal frame describes.
// The two are NOT interchangeable -- see PolyhedronView::inertia_matrix.
// ---------------------------------------------------------------------------
enum class MassModel {
    Solid,   // uniform density over the enclosed solid (needs a closed oriented surface)
    Atoms,   // uniform mass at the vertices (the molecular convention; needs no topology)
};

// ---------------------------------------------------------------------------
// Why a principal frame could not be built.  Every failure returns the identity
// matrix, which is a perfectly legal rotation, so a caller that does not look at
// the code cannot tell a frame from a fallback.
// ---------------------------------------------------------------------------
struct InertialFrame {
    enum class Code {
        Ok,
        DegenerateMass,            // the distribution has no body: zero mass / zero tensor
        NotPositiveSemidefinite,   // the second moment is not a second moment (see MassMoments)
        NonFiniteTensor,           // NaN/inf in the tensor -- non-finite coordinates upstream
        NotUnitary,                // the eigenvector matrix is not orthogonal to 1e-2
    };
    Code     code = Code::Ok;
    matrix3d axes = matrix3d::unit_matrix();   // row i is the i-th principal axis
};

template<typename T = double>
struct PolyhedronView : PlanarGraphView {
    std::span<coord3<T>> points;
    static constexpr uint8_t default_dmax = 10;

    PolyhedronView() = default;

    // Construct from adjacency view + coordinate span.
    PolyhedronView(const PlanarGraphView& g, std::span<coord3<T>> pts)
        : PlanarGraphView(g), points(pts) {}

    // Full view constructor (adjacency + coordinates).
    PolyhedronView(node_t N, int dmax,
                   std::span<node_t> neighbours, std::span<uint8_t> deg,
                   std::span<coord3<T>> pts,
                   std::span<uint8_t> twin = {})
        : PlanarGraphView(N, dmax, neighbours, deg, twin), points(pts) {}

    // -----------------------------------------------------------------------
    // Batchability contract (extends graph fields with `points`).
    // Canonical tuple order: {neighbours, deg, twin, points}.
    //   points : 1 coord per vertex  (N total)
    // -----------------------------------------------------------------------
    static constexpr std::size_t n_fields = 4;

    auto to_tuple() {
        return std::forward_as_tuple(this->neighbours, this->deg, this->twin, points);
    }
    auto to_tuple() const {
        return std::forward_as_tuple(this->neighbours, this->deg, this->twin, points);
    }

    static constexpr std::array<std::size_t, n_fields>
    get_size_factors(int /*N*/, int dmax) {
        return { (std::size_t)dmax, (std::size_t)1, (std::size_t)dmax, (std::size_t)1 };
    }

    // --- Geometry queries ---
    //
    // The three surface integrals all integrate centroid_surface(), and all
    // return a bare double because the Fortran C ABI (get_surface_area_ /
    // get_volume_, src/fortran/volume.f) and the pybind generator both classify
    // this API by return type -- a coded struct here would drop .surface_area and
    // .volume from the Python bindings.  So the outcome hangs off the surface
    // instead, exactly as principal_axes() hangs its outcome off
    // principal_frame(): on ANY centroid_surface() code other than Ok there is no
    // surface to integrate and all three return 0.0.  Zero area and zero volume
    // are legal-looking answers, so a caller that must tell them apart from a
    // real measurement calls centroid_surface() and reads its code.
    double surface_area() const;
    double volume() const { return volume_divergence(); }
    double volume_tetra() const;
    double volume_divergence() const;
    double diameter() const;
    pair<coord3<T>,coord3<T>> bounding_box() const;
    coord3<T> width_height_depth() const;

    // The closed triangle surface every surface integral on this class runs
    // over: the polygon faces split at their centroids (centroid_triangulation),
    // with those centroids appended to the vertex list.  One object rather than
    // two calls because the halves index each other -- `tris` indexes `points` --
    // and because the point of using it FOUR times (surface_area, volume_tetra,
    // volume_divergence, volume_moments) is that they integrate the SAME surface.
    struct CentroidSurface {
        // Why an empty surface is not an empty answer.  Every failure leaves the
        // two lists EMPTY, and an empty triangle list integrates to a perfectly
        // plausible 0 -- for area, for volume, for every moment -- so a caller
        // that eyeballs the number cannot tell "this body has no volume" from
        // "this object was never a body".  Read `code`.
        enum class Code {
            Ok,             // points and tris are a closed, consistently oriented surface
            NoFaces,        // nothing to integrate: no vertices or no edges, hence no faces
            NotOrientable,  // it has faces, but they are not a sphere's: the rotation
                            // system realises some other surface, or is not one at all
                            // (see OrientedSurface, which is what decides this)
        };
        Code            code = Code::NoFaces;
        vector<coord3d> points;   // the polyhedron's vertices, then one centroid per face
        vector<tri_t>   tris;     // closed and consistently oriented iff code == Ok
    };

    // @anchor polyhedron-centroid-surface
    // @error  NoFaces       when oriented_surface().code == Degenerate
    // @error  NotOrientable when oriented_surface().code is neither Ok nor Degenerate
    // @post on Ok:            !tris.empty() && points.size() == N + faces().size()
    // @post on NoFaces:       points.empty() && tris.empty()
    // @post on NotOrientable: points.empty() && tris.empty()
    // @time   O(E_dmax)
    CentroidSurface centroid_surface() const;

    // Moments of the enclosed solid: centroid_surface() -- the same surface
    // surface_area(), volume_tetra() and volume_divergence() integrate over --
    // through the volume_moments() kernel above.  So `.mass` (which for this
    // weighting IS the enclosed volume) and volume_tetra() are the same integral,
    // agreeing exactly in exact arithmetic; in floating
    // point they part company where the apex choice earns its keep, since
    // volume_tetra() sums tetrahedra from the ORIGIN rather than from the vertex
    // mean.  Measured on a 3x5x7 box (2026-08-08, tests/volume-moments-test.cc):
    // both exact at offsets up to 1e6, but at offset 1e8 the origin-apex sum has
    // relative error 9.5e5 while this one is still exact.  On a cage shifted by
    // |d| = 18 the two agree to 3e-15 relative.
    //
    // The one member of the four that CAN say so returns MassMoments::Degenerate
    // when centroid_surface() has no surface to hand it, forwarding that code
    // rather than integrating an empty triangle list to a zero that means the
    // same thing silently.
    // @pre closed: is_consistently_oriented() && !faces().empty()
    // @post on a centroid_surface() code other than Ok: result.code == Degenerate
    MassMoments volume_moments() const;

    // The CENTRAL inertia tensor I = tr(M) Id - M of the second moment M of one
    // of two mass distributions (MassModel) -- taken about that distribution's
    // own centre, so it is invariant under translating the polyhedron (before
    // 2026-08-07 the vertex sum was about the ORIGIN, and an off-centre cage got
    // a frame that depended on where it sat).
    //   MassModel::Solid (DEFAULT): the enclosed solid of uniform density, about
    //     the volume centroid -- exact per triangle, and independent of how
    //     finely the surface is triangulated to the extent that the faces are
    //     planar (a real cage face is not: fan-triangulating instead of splitting
    //     at the centroid moves the C60 volume by 3.3e-5 relative).  This is the
    //     answer the type name promises, and the one a principal frame should be
    //     built on.
    //   MassModel::Atoms: uniform mass at the VERTICES, about the vertex mean --
    //     the molecular convention (mass sits at the atoms), right for a fullerene
    //     cage, and an explicit opt-in because a densely triangulated region
    //     silently drags the axis.
    // The two are NOT interchangeable, so choose deliberately: measured over 40
    // C60/C70 cages plus a squashed icosahedral C80 the two tensors agree to 5.9%
    // (Frobenius, unit-trace) but the FRAMES do not -- an axis both models resolve
    // turns by up to 3.40 deg, 3.70 deg at a marginal eigenvalue gap, and 6.27 deg
    // raw inside a near-degenerate eigenplane (where no direction is determined at
    // all, so that last number measures nothing).  A symmetric cage is the
    // exception: the squashed icosahedral C80 agrees to 0.0006 deg.  Numbers from
    // tests/volume-moments-test.cc, which prints them on every run.  The repo's
    // molecular call sites therefore pass MassModel::Atoms explicitly rather than
    // inherit the default; of the eight, five are live (fullerene_polyhedron()
    // and programs/{fullerene-isomers-polyhedra,fullerene-polyhedron,
    // polyhedron-optimize,halma-polyhedron}.cc) and three are not built: density.cc
    // is commented out of programs/CMakeLists.txt, and gaudi.cc / gs-ex.cc appear
    // in no build file and no longer compile against the current API
    // (gaudi.cc:498 pushes onto points, now a std::span; gs-ex.cc:309 calls
    // dual(10,true) against a nullary dual()).  So the pinning is VERIFIED for
    // five call sites and asserted for three.
    // NOTE: the SYCL path (src/sycl/geometry-properties.cc) has its own
    // vertex-weighted, origin-based inertia_matrix/principal_axes and does NOT
    // follow this default.
    //
    // @anchor polyhedron-inertia-matrix
    // @pre closed: implies(m == MassModel::Solid,
    //                      is_consistently_oriented() && !faces().empty())
    // @post symmetric: (result - result.transpose()).norm() == 0
    matrix3d inertia_matrix(MassModel m = MassModel::Solid) const;

    // The principal frame of inertia_matrix(m): row i of the result is the unit
    // eigenvector of the i-th eigenvalue, ordered by |lambda| ascending (see
    // matrix3d::eigensystem), which for a positive-semidefinite tensor is
    // ascending lambda -- so row 0 is the smallest moment, the LONGEST axis.
    //
    // Every failure returns the identity, which is a legal rotation and
    // therefore indistinguishable from a real frame: a caller that needs to know
    // calls principal_frame(m) instead and reads its code.
    //
    // @anchor polyhedron-principal-axes
    // @pre closed: implies(m == MassModel::Solid,
    //                      is_consistently_oriented() && !faces().empty())
    // @post orthogonal: (result*result.transpose() - matrix3d::unit_matrix()).norm() < 1e-2
    matrix3d principal_axes(MassModel m = MassModel::Solid) const;

    // principal_axes(m) with the outcome attached -- the same matrix, plus the
    // reason when it is the identity fallback rather than a frame.
    //
    // @anchor polyhedron-principal-frame
    // @pre closed: implies(m == MassModel::Solid,
    //                      is_consistently_oriented() && !faces().empty())
    // @post on Ok: (result.axes*result.axes.transpose() -
    //                   matrix3d::unit_matrix()).norm() < 1e-2
    // @post on DegenerateMass:          result.axes == matrix3d::unit_matrix()
    // @post on NotPositiveSemidefinite: result.axes == matrix3d::unit_matrix()
    // @post on NonFiniteTensor:         result.axes == matrix3d::unit_matrix()
    // @post on NotUnitary:              result.axes == matrix3d::unit_matrix()
    InertialFrame principal_frame(MassModel m = MassModel::Solid) const;
    bool is_invalid() const;

    Polyhedron incremental_convex_hull() const;
    Polyhedron dual() const;
    Polyhedron leapfrog_dual() const;

    bool optimize(int opt_method=3, double ftol=1e-10);
    bool optimize_other(bool optimize_angles=true, map<edge_t,double> zero_values_dist={});

    vector<face_t> faces(int face_max=INT_MAX) const { return compute_faces(face_max); }

    void scale(const coord3<T>& s) { for(node_t u=0;u<N;u++) points[u] *= s; }
    void move(const coord3<T>& d)  { for(node_t u=0;u<N;u++) points[u] += d; }
    // Apply a linear map to every vertex in place: points <- M points.
    void transform(const matrix3d& M) { for(node_t u=0;u<N;u++) points[u] = M * points[u]; }
    void move_to_origin() { move(-mean(std::span<const coord3<T>>(points.data(), N))); }
    // Rotate into the principal frame of inertia_matrix(m), in place.  Inherits
    // that method's precondition: MassModel::Solid integrates the enclosed solid,
    // so it needs a closed oriented surface, where MassModel::Atoms needs no
    // topology at all.  An unbuildable frame silently applies the identity --
    // call principal_frame(m) first if that must be distinguished.
    // @pre closed: implies(m == MassModel::Solid,
    //                      is_consistently_oriented() && !faces().empty())
    void align_with_axes(MassModel m = MassModel::Solid) { transform(principal_axes(m)); }

    vector<coord2d> polar_angles() const {
        vector<coord2d> angles(N);
        for(node_t u=0;u<N;u++) angles[u] = coord3d(points[u][0],points[u][1],points[u][2]).polar_angle();
        return angles;
    }

    static vector<coord3<T>> polar_mapping(const vector<coord2d>& angles) {
        vector<coord3<T>> surface(angles.size());
        for(size_t u=0;u<surface.size();u++){
            const double theta = angles[u].first, phi = angles[u].second;
            surface[u] = coord3<T>((T)(cos(theta)*sin(phi)), (T)(sin(theta)*sin(phi)), (T)cos(phi));
        }
        return surface;
    }
};

// ---------------------------------------------------------------------------
// Explicit specialization declarations for PolyhedronView<double>.
// Prevent "specialization after instantiation" errors: the compiler sees these
// declarations before any implicit instantiation of PolyhedronView<double>.
// Definitions are in src/c++/polyhedron.cc and polyhedron-optimize.cc.
// ---------------------------------------------------------------------------
template<> double PolyhedronView<double>::diameter() const;
template<> double PolyhedronView<double>::surface_area() const;
template<> double PolyhedronView<double>::volume_divergence() const;
template<> double PolyhedronView<double>::volume_tetra() const;
template<> pair<coord3d,coord3d> PolyhedronView<double>::bounding_box() const;
template<> coord3d PolyhedronView<double>::width_height_depth() const;
template<> PolyhedronView<double>::CentroidSurface PolyhedronView<double>::centroid_surface() const;
template<> MassMoments PolyhedronView<double>::volume_moments() const;
template<> matrix3d PolyhedronView<double>::inertia_matrix(MassModel m) const;
template<> matrix3d PolyhedronView<double>::principal_axes(MassModel m) const;
template<> InertialFrame PolyhedronView<double>::principal_frame(MassModel m) const;
template<> bool PolyhedronView<double>::is_invalid() const;
template<> Polyhedron PolyhedronView<double>::incremental_convex_hull() const;
template<> Polyhedron PolyhedronView<double>::dual() const;
template<> Polyhedron PolyhedronView<double>::leapfrog_dual() const;
template<> bool PolyhedronView<double>::optimize(int opt_method, double ftol);
template<> bool PolyhedronView<double>::optimize_other(bool optimize_angles, map<edge_t,double> zero_values_dist);

// ---------------------------------------------------------------------------
// DeltahedronView: triangulation with 3D vertex coordinates (equilateral
// triangle embedding).
// ---------------------------------------------------------------------------
// Forward declaration for return types.
class Deltahedron;

// --- Curvature-cone detection (DeltahedronView<double>::find_cones) ---
//
// A closed genus-0 triangulation carries total Gaussian curvature 4*pi
// (Gauss-Bonnet); a fullerene shell concentrates it in 12 cones of pi/3 each.
// find_cones locates those cones as surface points with the angle-defect mass
// each integrates -- the inverse of placing pentagons on a sphere. See
// src/c++/deltahedron-cones.cc for the algorithm.

// A curvature cone located on the triangulation surface.
struct SurfaceCone {
  int                  face;                 // triangle index (into triangles())
  std::array<double,3> bary;                 // barycentric coords within `face`
  coord3d              position;             // 3D position = bary . triangle corners
  double               integrated_K;         // signed angle defect summed over the cone's
                                             // R-disk, EXCLUDING vertices already claimed by
                                             // an earlier cone (so cones partition the mass);
                                             // ~pi/3 ideal, slightly less if disks overlap
  bool                 meanshift_converged;
  int                  meanshift_iters;
};

struct ConeFinderParams {
  // THE default number of Taubin passes for cone detection.  Every localizer that
  // seeds from find_cones (the sub-project's signed / bayesian / joint tracks) forwards
  // its own taubin_iters here, so this constant is the single source of that default.
  static constexpr int default_taubin_iters = 10;

  double R_pent;                  // disk + exclusion radius (a RADIUS), mesh units (required)
  double threshold_frac  = 0.1;   // accept a cone iff disk-K >= this * pi/3
  int    max_centres     = 12;    // upper bound (fullerene: exactly 12 pentagons)
  int    taubin_iters    = default_taubin_iters;  // Taubin passes on the internal working copy
  int    meanshift_iters = 200;   // max mean-shift iterations per centre
  explicit ConeFinderParams(double R) : R_pent(R) {}
};

struct DelaunayTriangulation;   // forward decl for DeltahedronView::intrinsic_metric()

template<typename T = double>
struct DeltahedronView : TriangulationView {
    std::span<coord3<T>> points;
    static constexpr uint8_t default_dmax = 6;

    DeltahedronView() = default;

    // Construct from adjacency view + coordinate span.
    DeltahedronView(const TriangulationView& t, std::span<coord3<T>> pts)
        : TriangulationView(t), points(pts) {}

    // Full view constructor (adjacency + coordinates).
    DeltahedronView(node_t N, int dmax,
                    std::span<node_t> neighbours, std::span<uint8_t> deg,
                    std::span<coord3<T>> pts,
                    std::span<uint8_t> twin = {})
        : TriangulationView(N, dmax, neighbours, deg, twin), points(pts) {}

    // -----------------------------------------------------------------------
    // Batchability contract (extends graph fields with `points`).
    // Canonical tuple order: {neighbours, deg, twin, points}.
    // -----------------------------------------------------------------------
    static constexpr std::size_t n_fields = 4;

    auto to_tuple() {
        return std::forward_as_tuple(this->neighbours, this->deg, this->twin, points);
    }
    auto to_tuple() const {
        return std::forward_as_tuple(this->neighbours, this->deg, this->twin, points);
    }

    static constexpr std::array<std::size_t, n_fields>
    get_size_factors(int /*N*/, int dmax) {
        return { (std::size_t)dmax, (std::size_t)1, (std::size_t)dmax, (std::size_t)1 };
    }

    // --- Quality metrics ---
    double max_angle_relerr() const;
    int count_concave() const;

    // --- Geometry operations ---
    vector<face_t> compute_dual_faces() const;
    void smooth(double q);
    // Taubin lambda|mu low-pass: `iters` passes of a forward step smooth(lambda)
    // then a shrink-correcting step smooth(mu) (mu < -lambda). Modifies points.
    void taubin_smooth(int iters, double lambda=0.5, double mu=-0.53);
    // Discrete Gaussian curvature K(v) = 2*pi - sum of incident triangle angles
    // (the angle defect), per vertex, on the current embedding. Size N. The
    // tris overload reuses an already-computed triangulation (triangles() is O(E)).
    std::vector<double> angle_defects() const;
    std::vector<double> angle_defects(const std::vector<tri_t>& tris) const;
    // Locate up to params.max_centres curvature cones (e.g. the 12 fullerene
    // pentagons of a tomogram shell), detected on a Taubin-smoothed copy (this
    // mesh is untouched) and reported on the original geometry, in a deterministic
    // order. Fewer than max_centres is under-segmentation, a legitimate outcome.
    // If smoothed_points is non-null it receives the smoothed copy. See
    // deltahedron-cones.cc.
    // @pre     radius:    params.R_pent > 0 (and finite)
    // @pre     fraction:  params.threshold_frac in [0,1]
    // @pre     finite:    every vertex coordinate is finite (no NaN/inf)
    // @post    bounded:   result.size() <= params.max_centres
    // @post    gated:     every returned cone has integrated_K >= threshold_frac * pi/3
    // @post    on_mesh:   every cone's face is a valid triangle index, position lies on it
    // @throws  std::invalid_argument when a @pre is violated (bad params or non-finite geometry)
    std::vector<SurfaceCone> find_cones(const ConeFinderParams& params,
                                        std::vector<coord3d>* smoothed_points = nullptr) const;
    // The intrinsic surface metric: the iDT whose edge lengths are this mesh's
    // Euclidean edge lengths. The bridge from a 3D embedding to the purely-intrinsic
    // curvature flow (curvature-flow.hh); the only step that reads 3D coordinates.
    DelaunayTriangulation intrinsic_metric() const;
    Deltahedron GCtransform(unsigned k, unsigned l) const;
    Deltahedron halma_transform(int m) const;

    // optimize() and optimize_patch() stay on Deltahedron (use optimizer state fields).
    int reflect_concave(std::span<coord3<T>> pts, double threshold=0,
                        const vector<bool>& fixed={},
                        vector<bool>* reflected_mask=nullptr) const;
    int reflect_all_concave(std::span<coord3<T>> pts, double threshold=0,
                            const vector<bool>& fixed={},
                            vector<bool>* reflected_mask=nullptr) const;
    // Project concave vertices onto the convex hull of the point set.
    int project_onto_convex_hull(std::span<coord3<T>> pts) const;
    double gradient_check(std::span<const coord3<T>> geometry, double target_L=0, double eps=1e-6) const;
    double hessian_check(std::span<const coord3<T>> geometry, const vector<bool>& free_mask,
                         const vector<bool>& interior_mask={}, double target_L=0,
                         double eps=1e-5, bool verbose=false) const;

    // --- AET (equilateral-triangle target) energy stack ---
    // The energy/gradient/HVP that optimize() and optimize_patch() are
    // built on (E_bond + E_angle + E_curv [+ E_flat] [+ E_conv]), exposed
    // as thin wrappers over the internal term implementations -- the
    // single source of truth -- for external optimizer frameworks
    // (incubation seam for claude-projects/optimize).  `edges` is the
    // cached undirected_edges() list, precomputed once by the caller so
    // repeated evaluations keep the internal callers' cost profile.
    // aet_energy_gradient: returns E; if grad is non-null it is zeroed
    // and filled.  aet_hv_product: Hv = (grad^2 E) v (Hv zeroed here).
    double aet_energy_gradient(const vector<edge_t>& edges,
                               std::span<const coord3<T>> x, vector<coord3<T>>* grad,
                               double target_L, double k_bond, double k_angle,
                               double k_curv, double k_flat,
                               double k_conv = 0, double sigma_conv = 0,
                               const vector<bool>& conv_mask = {}) const;
    void aet_hv_product(const vector<edge_t>& edges,
                        std::span<const coord3<T>> x, const vector<coord3<T>>& v,
                        vector<coord3<T>>& Hv, double target_L,
                        double k_bond, double k_angle, double k_curv, double k_flat,
                        const vector<bool>& fixed = {},
                        double k_conv = 0, double sigma_conv = 0) const;
    // Signed convexity height h per vertex (h > 0 convex, h < 0 concave;
    // fixed or degree > 6 vertices get h = 1).  Part of the AET seam: the
    // convex-constrained acceptance gate consumes it.
    void aet_h_values(std::span<const coord3<T>> x, vector<double>& h,
                      const vector<bool>& fixed = {}) const;
};

// ---------------------------------------------------------------------------
// Explicit specialization declarations for DeltahedronView<double>.
// Definitions live in src/c++/deltahedron.cc.
// ---------------------------------------------------------------------------
template<> double DeltahedronView<double>::max_angle_relerr() const;
template<> int    DeltahedronView<double>::count_concave() const;
template<> vector<face_t> DeltahedronView<double>::compute_dual_faces() const;
template<> void   DeltahedronView<double>::smooth(double q);
template<> void   DeltahedronView<double>::taubin_smooth(int iters, double lambda, double mu);
template<> std::vector<double> DeltahedronView<double>::angle_defects() const;
template<> std::vector<double> DeltahedronView<double>::angle_defects(const std::vector<tri_t>& tris) const;
template<> std::vector<SurfaceCone> DeltahedronView<double>::find_cones(
    const ConeFinderParams& params, std::vector<coord3d>* smoothed_points) const;
template<> DelaunayTriangulation DeltahedronView<double>::intrinsic_metric() const;
template<> Deltahedron DeltahedronView<double>::GCtransform(unsigned k, unsigned l) const;
template<> Deltahedron DeltahedronView<double>::halma_transform(int m) const;
template<> int DeltahedronView<double>::reflect_concave(std::span<coord3d> pts, double threshold,
                                                        const vector<bool>& fixed,
                                                        vector<bool>* reflected_mask) const;
template<> int DeltahedronView<double>::reflect_all_concave(std::span<coord3d> pts, double threshold,
                                                            const vector<bool>& fixed,
                                                            vector<bool>* reflected_mask) const;
template<> int DeltahedronView<double>::project_onto_convex_hull(std::span<coord3d> pts) const;
template<> double DeltahedronView<double>::gradient_check(std::span<const coord3d> geometry,
                                                          double target_L, double eps) const;
template<> double DeltahedronView<double>::hessian_check(std::span<const coord3d> geometry,
                                                         const vector<bool>& free_mask,
                                                         const vector<bool>& interior_mask,
                                                         double target_L, double eps,
                                                         bool verbose) const;
template<> double DeltahedronView<double>::aet_energy_gradient(
    const vector<edge_t>& edges, std::span<const coord3d> x, vector<coord3d>* grad,
    double target_L, double k_bond, double k_angle, double k_curv, double k_flat,
    double k_conv, double sigma_conv, const vector<bool>& conv_mask) const;
template<> void DeltahedronView<double>::aet_hv_product(
    const vector<edge_t>& edges, std::span<const coord3d> x, const vector<coord3d>& v,
    vector<coord3d>& Hv, double target_L, double k_bond, double k_angle,
    double k_curv, double k_flat, const vector<bool>& fixed,
    double k_conv, double sigma_conv) const;
template<> void DeltahedronView<double>::aet_h_values(
    std::span<const coord3d> x, vector<double>& h,
    const vector<bool>& fixed) const;

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
static_assert(std::is_trivially_copyable_v<PolyhedronView<double>>,
    "PolyhedronView<double> must be trivially copyable");
static_assert(std::is_trivially_copyable_v<PolyhedronView<float>>,
    "PolyhedronView<float> must be trivially copyable");
static_assert(std::is_trivially_copyable_v<DeltahedronView<double>>,
    "DeltahedronView<double> must be trivially copyable");
static_assert(std::is_trivially_copyable_v<DeltahedronView<float>>,
    "DeltahedronView<float> must be trivially copyable");


// Verify correct hierarchy relationships
static_assert(std::is_base_of_v<GraphView, PlanarGraphView>);
static_assert(std::is_base_of_v<PlanarGraphView, CubicGraphView>);
static_assert(std::is_base_of_v<CubicGraphView, FullereneGraphView>);
static_assert(std::is_base_of_v<PlanarGraphView, TriangulationView>);
static_assert(std::is_base_of_v<TriangulationView, FullereneDualView>);
static_assert(std::is_base_of_v<PlanarGraphView, PolyhedronView<double>>);
static_assert(std::is_base_of_v<PlanarGraphView, PolyhedronView<float>>);
static_assert(std::is_base_of_v<TriangulationView, DeltahedronView<double>>);
static_assert(std::is_base_of_v<TriangulationView, DeltahedronView<float>>);
