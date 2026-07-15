#ifndef DELAUNAY_HH
#define DELAUNAY_HH

#include "triangulation.hh"
#include "geometry.hh"

#include <array>
#include <cstdio>       // FILE* for the .idt serialization (to_ascii / from_ascii)
#include <functional>


// Diamond: the local geometry around an edge in a metrized triangulation.
//
//      B
//     / \          upper triangle: sides (e, a, b)
//    a   b         a = side adjacent to u, b = side adjacent to v
//   /     \  .
//  u---e---v       e = diagonal being tested/flipped
//   \     /
//    c   d         lower triangle: sides (e, c, d)
//     \ /          c = side adjacent to u, d = side adjacent to v
//      D
//
// All geometric predicates (Delaunay, convexity, flipped length) depend only
// on these five edge lengths.  No vertex IDs or topology needed.
struct Diamond {
  double e, a, b, c, d;

  bool   is_delaunay()    const;  // cot(angle_B) + cot(angle_D) >= 0
  bool   is_convex()      const;  // angle sum < pi at both u and v
  double flipped_length() const;  // length of BD (the other diagonal)

  // Cocircular ("tight") test: cot(angle_B) + cot(angle_D) == 0 exactly,
  // i.e. the four points u, v, B, D are concyclic on the surface.  In this
  // case both triangulations of the diamond are equally-valid Delaunay
  // refinements of the same cell.  Uses exact integer arithmetic on
  // length-squared (valid when all five lengths square to non-negative
  // integers, e.g. equilateral triangulations and their flips).
  // See CANONICAL-TESSELATION.md for the derivation.
  bool is_cocircular() const;

  // Floating-point cocircular test for general (non-equilateral) metrics,
  // where length-squared is not integer so the exact predicate above does not
  // apply: tight iff |cot(angle_B) + cot(angle_D)| < tol. Scale-invariant
  // (cotangents are dimensionless), so tol is a pure angle threshold.
  bool is_cocircular(double tol) const;
};

// ============================================================================
// Canonical Delaunay tesselation
//
// On a piecewise-flat surface, the iDT triangulation is not unique: when
// k >= 4 cone points are concyclic, the cocircular k-gon admits multiple
// equally-valid Delaunay triangulations, and which one the algorithm emits
// depends on the input vertex labelling.  The Bobenko-Springborn (2007)
// theorem guarantees that the underlying TESSELATION (cell decomposition
// where every cocircular polygon is left intact) IS unique.
//
// CanonicalTesselation is exactly that invariant: a sorted multi-set of
// CCW-oriented polygons, each in user-supplied vertex labels with
// integer length-squared annotation per edge.  Two iDT outputs that
// differ only by trivial flips inside cocircular cells (and/or by a
// label-equivariant relabelling, once mapped through the same external
// label-map) compare equal here, even though their raw triangulations
// differ.
// ============================================================================

struct CanonicalTesselation {
  // Polygon = CCW boundary of one Delaunay cell.  Each entry is
  // (vertex_label, length²-of-outgoing-edge-to-next-vertex).  Polygons
  // are normalized to lex-min cyclic rotation.
  using Polygon = std::vector<std::pair<int, long long>>;

  // Sorted multi-set of cells (a deterministic linear order is enough for
  // equality / ordering).
  std::vector<Polygon> cells;

  bool operator==(const CanonicalTesselation& o) const { return cells == o.cells; }
  bool operator!=(const CanonicalTesselation& o) const { return cells != o.cells; }
  bool operator< (const CanonicalTesselation& o) const { return cells <  o.cells; }

  // Counts (for quick queries; empty iDT yields zeros).
  int n_cells() const { return (int)cells.size(); }

  // Stable 64-bit fingerprint for hash maps and at-a-glance comparison.
  size_t fingerprint() const;

  // Compact human-readable form; one cell per line.
  std::string to_string() const;
};

// How a geodesic disk measures distance on the metric DCEL:
//   Edge   -- Dijkstra along mesh edges (a graph upper bound on the geodesic).
//   Unfold -- fast marching: the wavefront crosses unfolded triangle interiors, so
//             paths cut across faces (the true intrinsic geodesic, to first order).
enum class DiskMetric { Edge, Unfold };

// One cell of a geodesic-disk partition: the seeding vertex and the vertices it
// claims (nearest source, within R), each with its geodesic distance. See
// DelaunayTriangulation::geodesic_disks.
struct GeodesicDisk {
  int                                source;
  std::vector<std::pair<int,double>> members;   // (vertex, geodesic distance <= R)
};

// ============================================================================
// DCEL-based intrinsic Delaunay triangulation (delta-complex).
//
// Half-edge (DCEL) representation that correctly handles multi-edges and
// self-loops.  Every edge is identified by a half-edge index, so flip_edge(h)
// is unambiguous even when multiple edges connect the same vertex pair.
//
// Twin convention: half-edges 2k and 2k+1 are always twins.
// Face orientation: he_next traverses each face CCW.
// Vertex circulation: cw(h) rotates CW around origin(h).
// ============================================================================

struct DelaunayTriangulation {
  // --- Counts ---
  int nv = 0;   // live vertices
  int nh = 0;   // allocated half-edges (including dead slots)
  int nf = 0;   // allocated faces (including dead slots)

  // --- Half-edge topology (indexed 0..nh-1) ---
  vector<int>    he_next;    // next half-edge CCW in same face
  vector<int>    he_origin;  // origin vertex (-1 = dead)
  vector<int>    he_face;    // face to the left

  // --- Metric (indexed 0..nh-1) ---
  vector<double> he_length;  // edge length (same for h and twin(h))
  vector<double> he_angle;   // angle at origin(h) in face(h)

  // --- Per-vertex (indexed 0..nv-1) ---
  vector<int>    v_out;          // one outgoing half-edge (-1 = dead vertex)
  vector<double> v_cone_angle;   // cone angle = total vertex angle (deg*pi/3 when equilateral)
  vector<int>    v_orig_degree;  // degree in original triangulation

  // --- Diagnostics ---
  // If > 0, remove_flat_vertices prints live-vertex progress to stderr every
  // this-many removed vertices (and a header/footer). 0 = silent. Output is
  // flushed each line so it is visible even on a slow run.
  int verbose_removal = 0;

  // --- Per-face (indexed 0..nf-1) ---
  vector<int> f_he;        // one boundary half-edge (-1 = dead face)

  // --- Free lists ---
  vector<int> free_edges;  // recycled edge slots (half-edge id / 2)
  vector<int> free_faces;  // recycled face slots

  // --- Clean accessors ---
  int  twin(int h)  const { return h ^ 1; }
  int  edge(int h)  const { return h >> 1; }
  int  prev(int h)  const { return he_next[he_next[h]]; }  // only for triangulations
  int  dest(int h)  const { return he_origin[twin(h)]; }
  bool alive(int h) const { return he_origin[h] >= 0; }
  // Bigon edge: both half-edges of h bound the same face.  Arises in
  // Δ-complexes around low-degree cone vertices (an i–j edge of an "iji"
  // face); dihedral quantities across such an edge are undefined.
  bool is_bigon(int h) const { return he_face[h] == he_face[twin(h)]; }

  // CW rotation around origin(h): next outgoing half-edge clockwise.
  int cw(int h) const { return he_next[twin(h)]; }

  // CCW rotation around origin(h): next outgoing half-edge counterclockwise.
  int ccw(int h) const { return twin(prev(h)); }

  // The three half-edges of (triangular) face f, in he_next order starting
  // from its representative.  Pre: f is live (f_he[f] >= 0).
  std::array<int,3> face_halfedges(int f) const {
    const int h = f_he[f];
    return {h, he_next[h], prev(h)};
  }
  // The three corner vertices of face f, CCW (origins of face_halfedges(f)).
  std::array<int,3> face_vertices(int f) const {
    const auto h = face_halfedges(f);
    return {he_origin[h[0]], he_origin[h[1]], he_origin[h[2]]};
  }

  int vertex_degree(int v) const;  // count outgoing half-edges from v

  // Range over the outgoing half-edges of v (the cw ring from v_out[v]); empty if v
  // has no incident live edge. The canonical vertex traversal.
  struct IncidentHalfEdges {
    const DelaunayTriangulation& D; int v;
    struct iterator {
      const DelaunayTriangulation* D; int start, h; bool done;
      int       operator*()  const { return h; }
      iterator& operator++()       { h = D->cw(h); done = (h == start); return *this; }
      bool      operator!=(const iterator&) const { return !done; }
    };
    iterator begin() const { int h0 = D.v_out[v]; return {&D, h0, h0, h0 < 0}; }
    iterator end()   const { return {}; }
  };
  IncidentHalfEdges incident(int v) const { return {*this, v}; }

  // Range over one (even) half-edge per live edge. The canonical edge
  // traversal, skipping dead slots.
  struct LiveEdges {
    const DelaunayTriangulation& D;
    struct iterator {
      const DelaunayTriangulation* D; int h;
      void advance() { while (h < D->nh && !D->alive(h)) h += 2; }
      int       operator*()  const { return h; }
      iterator& operator++()       { h += 2; advance(); return *this; }
      bool      operator!=(const iterator& o) const { return h != o.h; }
    };
    iterator begin() const { iterator it{&D, 0}; it.advance(); return it; }
    iterator end()   const { return {&D, D.nh}; }
  };
  LiveEdges edges() const { return {*this}; }

  // --- Construction ---

  // A metric on the edges: length of the undirected edge {u,v}.
  using EdgeLengthFn = std::function<double(node_t u, node_t v)>;

  // Build the iDT topology from the combinatorial triangulation T with the
  // EQUILATERAL metric (all edge lengths 1, all angles pi/3, cone angle
  // deg*pi/3). This is the special case from_intrinsic_metric(T, 1).
  static DelaunayTriangulation from_triangulation(const Triangulation& T);

  // Build the iDT topology from T with a PRESCRIBED intrinsic metric: edge
  // lengths come from `length(u,v)`, angles and cone angles are then derived
  // from those lengths (recompute_all_angles + vertex_angle_sum). Generalises
  // from_triangulation to non-equilateral surfaces.
  static DelaunayTriangulation from_intrinsic_metric(const Triangulation& T,
                                                     const EdgeLengthFn& length);

  // --- Geometry ---
  Diamond diamond(int h) const;
  void recompute_face_angles(int f);
  // Recompute he_angle for every face, then refresh the v_cone_angle cache.
  // THE entry point for restoring angle coherence after a metric (length) change.
  void recompute_all_angles();
  // Recompute he_angle for every face only (no cone-cache refresh) -- for hot
  // callers that read just he_angle; pair with recompute_cone_angles when needed.
  void recompute_all_face_angles();
  // Refresh only the v_cone_angle cache from the current he_angle. Cone angle is
  // flip-invariant, so flips do not need this; a length change does.
  void recompute_cone_angles();

  // Total angle at vertex v = sum of incident corner angles. The Gaussian
  // curvature / angle defect is 2*pi - vertex_angle_sum(v); a flat vertex
  // has angle sum 2*pi. Requires angles current (call recompute_all_angles
  // after changing any length).
  double vertex_angle_sum(int v) const;

  // A vertex is flat (zero cone curvature, hence removable) iff its cone
  // angle is 2*pi to within flat_tol. Metric-based replacement for the
  // equilateral-only "degree == 6" test. flat_tol must lie above the input
  // metric's curvature noise floor and below the smallest real cone curvature:
  // 1e-6 suits exact (equilateral) metrics; a numerically solved metric (e.g.
  // a CEPS conformal solve, whose layout-cut seams leave ~5e-4 spurious
  // curvature at genuinely flat vertices) needs ~1e-2.
  bool is_flat(int v, double flat_tol = 1e-6) const;

  // Discrete Gaussian curvature K(v) = 2*pi - cone angle, per live vertex (the angle
  // defect). Reads the cached v_cone_angle, so requires it current.
  std::vector<double> curvature() const;

  // --- Delaunay operations ---
  bool is_delaunay_edge(int h) const;
  // Flip the diagonal of the diamond around edge h.  Accepts any edge with a
  // strictly convex diamond and a positive flipped length, including the
  // B == D case, which produces a self-loop edge at B (strictly Delaunay
  // when the flipped edge was non-Delaunay; paper lem:selfloop-delaunay).
  // @post on true:  the DCEL invariants hold with the new diagonal B-D in
  //                 place of h; only the two diamond faces are rewired; the
  //                 surface metric is unchanged (paper thm:metric-invariance).
  // @post on false: the complex is untouched.
  bool flip_edge(int h);
  // Global Delaunay restore, seeding with every live edge.
  // @post is_delaunay()
  // @throws std::runtime_error when the 200*nv flip budget is exhausted
  //         (impossible by the energy argument, paper thm:lawson-converge;
  //         fail-loud backstop).
  int  lawson_sweep();
  // Seeded core: restore Delaunay starting from a frontier of edges (one half-edge id
  // per edge), propagating to each flip's rim. For a hot caller that changed only a
  // local region of an otherwise-Delaunay mesh; reaches the same iDT the no-arg global
  // sweep would. `in_stack` (sized nh/2) is the caller-owned stack-membership scratch;
  // it is left all-false on normal return, so a hot loop reuses one buffer without
  // re-zeroing. The no-arg lawson_sweep() seeds with all live edges and its own scratch.
  // @pre  seeded: every non-Delaunay edge is in seed_edges or arises on the
  //       rim of a seed cascade (paper lem:seeded-sweep).
  // @post is_delaunay()
  // @throws std::runtime_error when the 200*nv flip budget is exhausted.
  int  lawson_sweep(const vector<int>& seed_edges, vector<bool>& in_stack);
  int  count_non_delaunay() const;
  int  flip_to_delaunay();
  bool is_delaunay() const;

  // --- Vertex removal ---
  void remove_flat_vertex(int v);
  // Remove every flat vertex (is_flat(v, flat_tol)), leaving the cones. See
  // is_flat for choosing flat_tol (1e-6 exact, ~1e-2 for a CEPS metric).
  //
  // Observer hook: if `on_pop` is non-empty, it is invoked as on_pop(v) at
  // EVERY work-list pop of a live flat vertex v -- after the live/flat guard
  // passes and BEFORE that vertex's self-loop cleanup, i.e. at the moment the
  // mesh is globally Delaunay per the driver invariant. It fires in every drain
  // path (the initial drain and each restructure round). The observer MUST NOT
  // mutate the mesh (it is an observation hook only); this is taken on trust,
  // not checked. An empty on_pop (the default) is a zero-overhead no-op.
  //
  // @throws std::runtime_error on an un-flippable self-loop at a flat vertex
  //         that survives the cocircular tie-break.  An EXACTLY TIED loop
  //         diamond (end angle = pi; realizable, paper rem:ties-real) is
  //         resolved by flipping the tied side's closing spoke, whose
  //         diamond is exactly cocircular at EVERY fan size (the
  //         cocircular-fan lemma, paper thm:cocircular-fan +
  //         thm:tie-break), so on exact sphere input of minimum degree 3
  //         this throw is provably dead (paper thm:main); it guards the
  //         floating-point boundary strip and out-of-scope inputs.
  // @throws std::runtime_error on Lawson budget exhaustion (propagated from
  //         lawson_sweep).
  // @throws std::runtime_error on a stuck reduction: a live flat vertex
  //         survives the fixed-point rounds (same scope note as above).
  void remove_flat_vertices(double flat_tol = 1e-6,
                            const std::function<void(int)>& on_pop = {});

  // Renumber the live vertices to 0..n_live-1 (dropping removed ones) and
  // shrink nv, rewriting he_origin and the per-vertex arrays. Needed after a
  // removal that leaves live vertices scattered (the metric-based path, which
  // unlike compute(T) does not sort flats last), because AlexandrovSolver
  // sizes its system by nv and assumes every index 0..nv-1 is a live cone.
  // Returns new_to_old: the original index of each surviving vertex.
  std::vector<int> compact_vertices();

  // --- Edge/face allocation ---
  int  alloc_edge();         // returns half-edge id of first half-edge
  int  alloc_face();         // returns face id
  void dealloc_edge(int h);  // mark edge as dead, add to free list
  void dealloc_face(int f);  // mark face as dead, add to free list

  // Allocate an edge and set its endpoints and length.
  // Returns the half-edge h with origin(h) = u, origin(twin h) = v,
  // length = L on both sides.  Faces remain unassigned.
  int  alloc_directed_edge(int u, int v, double L);

  // Allocate a vertex with the given cone angle and original degree, growing the
  // three per-vertex arrays (v_out, v_cone_angle, v_orig_degree) in lockstep.
  // Returns its id; v_out is left at -1 (the caller wires its edges).
  int  alloc_vertex(double cone_angle, int orig_degree);

  // Wire three half-edges into a CCW triangle face and compute its angles
  // from the stored edge lengths.  Returns the new face id.
  // Preconditions: h0, h1, h2 already have their origin and length set;
  // their endpoints form a triangle with origin(h0)=u, dest(h0)=v=origin(h1),
  // dest(h1)=w=origin(h2), dest(h2)=u.
  int  wire_triangle(int h0, int h1, int h2);

  // Split the face left of h0 (= a->b, with he_next giving b->c, c->a) into three
  // triangles fanning to a new flat vertex P, at spoke lengths {P->a, P->b, P->c}.
  // The three boundary edges are kept, so the metric is preserved exactly and P has
  // cone angle 2*pi. Refreshes v_cone_angle at P and the three corners. Returns P.
  // The 1->3 sibling of the (internal) edge bisection.
  int  split_face(int h0, std::array<double,3> spokes);

  // --- Full algorithm ---
  // Computes the unique intrinsic Delaunay triangulation of the input
  // surface (Bobenko-Springborn 2007).  The output is generally a
  // delta-complex and may contain multi-edges, self-loops, and bigons
  // around cone vertices (deg-2 cones) -- all valid iDT features.
  // On sphere input with minimum vertex degree 3 (any cubic-polyhedron
  // dual, fullerene duals included) success is unconditional
  // (paper thm:main).
  // @post result.nv == number of deg != 6 vertices of T, all live
  // @post result.is_delaunay()
  static DelaunayTriangulation compute(const Triangulation& T);

  // As compute(T), but for a PRESCRIBED intrinsic metric `length(u,v)`.
  // Flat vertices are identified by cone angle (is_flat), not degree, so it
  // works on any surface whose flat vertices need not be degree-6. Removes
  // the flat vertices, leaving the cone vertices in a delta-complex iDT
  // ready for AlexandrovSolver. Equivalent to compute(T) when length == 1.
  // flat_tol sets the cone/flat threshold (see is_flat): 1e-6 for an exact
  // metric, ~1e-2 for a numerically solved one.
  // If new_to_old is non-null it receives compact_vertices()' map: the
  // original T-label of each surviving cone vertex 0..nv-1.
  static DelaunayTriangulation compute(const Triangulation& T,
                                       const EdgeLengthFn& length,
                                       double flat_tol = 1e-6,
                                       std::vector<int>* new_to_old = nullptr);

  // --- Surface metric (intrinsic; promoted from the delta-complex project) ---
  // Per-cone-pair geodesic distances and geodesics on the metric delta-complex,
  // by priority-bounded BFS over triangle strips unfolded into the Eisenstein
  // plane (exact-integer for the simple distances; floating point only in the
  // APSP step).  The iDT IS the cone graph (flat vertices already removed), so
  // every vertex 0..nv-1 is a cone: these take no node subset and no source
  // mesh, and the BFS bound is derived from the iDT's own weighted graph
  // diameter.  Validated on every fullerene dual C20-C160 (211,203,353 isomers,
  // 0 failures).  Mirror of TriangulationView's surface-metric methods
  // (graphview.hh).
  //
  // Default scope: cone-to-cone.  Self-geodesics (closed loops based at one
  // cone) are excluded by an API pin in the walk and not computed; pass
  // calculate_self_geodesics = true to compute the diagonal too -- the BFS
  // runs in self-mode for each cone, picking up seed-edge self-loops at seed
  // setup and any u_start apex placements the unmodified BFS validity gate
  // accepts.  Wrap-around closures rely on multi-seed coverage of the full
  // cone angle.  Mirrors TriangulationView's calculate_self_geodesics flag.
  // See claude-projects/delta-complex/DELTA-COMPLEX-SURFACE-METRIC.md
  // §"Self-geodesics" for the full account.

  // A simple geodesic u -> v, exact in Z[w]: the displacement g (u at origin,
  // |g|^2 = squared length) in the unfolding seeded at half-edge `axis` with
  // its far endpoint placed at B_seed (the split-prime representative).  The
  // triple (axis, B_seed, g) determines the developed face strip uniquely.
  struct simple_geodesic { Eisenstein g; int axis; Eisenstein B_seed; };

  // A general geodesic u -> K_1 -> ... -> v: simple geodesics broken at
  // intermediate cones.  Mirror of TriangulationView::geodesic.
  struct geodesic {
    std::vector<simple_geodesic> segments;
    geodesic() = default;
    geodesic(int) {}   // matrix<geodesic>(m, n) zero-init compatibility
  };

  // Exact squared simple-geodesic distances; LLONG_MAX where no simple geodesic
  // exists.  Diagonal: 0 in cone-to-cone mode; the closed-geodesic squared length
  // (or the search-radius sentinel) in self mode.
  // @throws std::logic_error on a deep invariant (degenerate face or
  // non-Loeschian edge), which cannot occur on a valid iDT.
  matrix<long long>       simple_square_surface_distances(bool calculate_self_geodesics = false) const;
  // The realizing simple geodesic per ordered cone pair (in u's unfolding frame).
  matrix<simple_geodesic> simple_geodesics(bool calculate_self_geodesics = false) const;
  // Squared surface distances (APSP-smoothed across intermediate cones).  If
  // geodesics_out != nullptr, also fills it with the composed per-pair geodesics.
  matrix<double>          surface_distances(bool calculate_self_geodesics = false,
                                            matrix<geodesic>* geodesics_out = nullptr) const;
  // The composed surface geodesic for every cone pair.
  matrix<geodesic>        surface_geodesics(bool calculate_self_geodesics = false) const;

  // Partition the surface into geodesic R-disks: one multi-source front seeded at all
  // `sources` at distance 0, so each vertex is claimed by its NEAREST source -- a
  // geodesic Voronoi cell clipped at R. Disjoint by construction (no vertex in two
  // disks). Edge metric relaxes along edges; Unfold also runs the triangle wavefront
  // across each incident face within a cell (cross-cell seams stay on the edge metric).
  std::vector<GeodesicDisk> geodesic_disks(const std::vector<int>& sources, double R,
                                           DiskMetric metric) const;

  // Concatenate the simple geodesics along a cone path [u, K_1, ..., v] into one
  // multi-segment geodesic; empty for paths of length <= 1.
  static geodesic compose_simple_geodesics(const std::vector<int>& path,
                                           const matrix<simple_geodesic>& simple);

  // --- Weighted graph algorithms on the 1-skeleton ---
  // Single-source weighted shortest paths on the DCEL with he_length as edge
  // weights.  Returns dist[v] for v in 0..nv-1; +infinity for unreachable v.
  // Standard Dijkstra with a binary heap; O((nv + alive_edges) log nv).
  std::vector<double> single_source_shortest_paths(int src) const;

  // An upper bound on the weighted graph diameter, via a textbook double-sweep
  // of single_source_shortest_paths (Dijkstra from any vertex u_0 -> find the
  // farthest vertex u_1 -> Dijkstra from u_1).  For metric graphs the
  // double-sweep result d(u_1, u_2) is a lower bound on the true diameter D
  // satisfying d(u_1, u_2) >= D/2, so the returned value is 2 * d(u_1, u_2)
  // and is guaranteed to be an upper bound on D.  Used to size BFS priority
  // cutoffs in the surface-metric routines.  O((nv + alive_edges) log nv).
  double diameter_upper_bound() const;

  // Smallest degree among live (non-removed) vertices, or INT_MAX if
  // none.  A value below 3 is one (but not the only) non-simplicial
  // signature -- use is_simplicial() for a complete check.
  int min_live_degree() const;

  // True iff the iDT's 1-skeleton is a simple graph: no self-loops,
  // no multi-edges.  Equivalently: the map h |-> (origin(h), dest(h))
  // is injective on live half-edges.  Self-loops fail because both
  // twins encode the arc (v,v); multi-edges fail because two non-twin
  // half-edges encode the same arc.
  // Non-simplicial outputs are valid iDT delta-complexes (Bobenko-
  // Springborn 2007), not algorithm failures.
  // Cost: O(E log E).
  bool is_simplicial() const;

  // True iff the DCEL is structurally well-formed: every live half-edge
  // belongs to exactly one face cycle traversed via he_next, and every
  // such cycle has length 3.  This is the DCEL counterpart of
  // Graph::is_consistently_oriented (which walks faces via next(v,u)).
  // A well-formed iDT can still be non-simplicial (multi-edges, self-
  // loops); a not-well-formed iDT signals a bug in the algorithm.
  // Cost: O(E).
  bool is_well_formed() const;

  // Bisect all multi-edges by inserting midpoint vertices.
  // Multi-edges (multiple geodesics between the same cone-point pair) can't be
  // realized as distinct straight edges in R³.  Bisecting each with a midpoint
  // makes them geometrically distinct.  Returns the number of vertices added.
  int bisect_multi_edges();

  // --- Tesselation invariant ---
  // True iff edge h is cocircular ("tight"): the diamond's four cone points
  // share a common circumcircle, so flipping h yields an equally-valid
  // Delaunay triangulation.  Both half-edges of an edge return the same value.
  bool is_cocircular_edge(int h) const;

  // Per-half-edge cocircular mask: tight[h] == tight[twin(h)]; dead half-edges
  // are false.  O(num_edges) integer-arithmetic predicates.
  std::vector<bool> cocircular_edges() const;

  // Float cocircular mask for general metrics (Diamond::is_cocircular(tol)).
  // Use as the `tight` argument to canonical_tesselation when the metric is
  // not equilateral (the integer predicate above is then invalid).
  std::vector<bool> cocircular_edges(double tol) const;

  // Canonical Delaunay tesselation in `vertex_labels` coordinates.
  // `vertex_labels[k]` is the external label assigned to D's live vertex k;
  // typically the vertex's position in the original Triangulation it came
  // from (as recovered through `Triangulation::sort_flat_last()`).
  // Cocircular interior edges are merged so each cell becomes one polygon;
  // polygons are normalized (lex-min rotation), the multi-set is sorted.
  CanonicalTesselation canonical_tesselation(
      const std::vector<int>& vertex_labels) const;

  // General tesselation extraction with a caller-supplied "interior" mask.
  // `tight[h]` true ⇒ edge h is interior to a cell (it is collapsed away);
  // `tight[h]` false ⇒ edge h is on a cell boundary.  Both half-edges of an
  // edge must agree (`tight[h] == tight[h^1]`).
  //
  // The cell-walk is identical to canonical_tesselation; only the source of
  // the interior-edge mask changes.  This lets the Alexandrov solver pass
  // a Bobenko-Izmestiev "inessential" mask (edges with θ_e = π in the
  // weighted-Delaunay tesselation at κ=0) to obtain the polytope's
  // 2-skeleton T̄ (the polygonal flat 2-faces of the convex polytope),
  // distinct from but compatible with the cocircular tesselation.
  //
  // Both notions agree exactly on flat-2-face diagonals whose four cone
  // points are also cocircular in the surface metric; they may differ
  // when the flat 2-face is non-cyclic, or when 4 surface-cocircular
  // points are spread across multiple polytope faces.  Cross-validation
  // between the two is a useful sanity check on the GCP geometry at κ=0.
  CanonicalTesselation canonical_tesselation(
      const std::vector<int>& vertex_labels,
      const std::vector<bool>& tight) const;

  // --- Validation ---
  bool check_consistency() const;

  // --- Serialization (.idt) ---
  // Text round-trip of the intrinsic Delaunay delta-complex (the ".idt" format, magic
  // "iDT-DCEL 1"). ALIVE-only and half-edge based, so it is faithful to the multi-edges,
  // self-loops and self-adjacent faces of a non-simplicial iDT that a vertex-pair adjacency
  // cannot express -- and it stays compact even after remove_flat_vertices leaves the arrays
  // full of dead slots (compact_vertices shrinks only nv). Only the topology + edge lengths are
  // stored, at full double precision; faces, v_out and angles are rebuilt on read. Format spec
  // and field layout are in delaunay.cc.
  //
  // @pre  D is compacted: every vertex index in [0, nv) is live (throws otherwise).
  // @post returns ferror(file) == 0.
  static bool to_ascii(const DelaunayTriangulation& D, FILE* file);

  // @post the result passes check_consistency(); a malformed stream is rejected, never returned
  //       as a silently-wrong triangulation.
  // @throws std::runtime_error naming the fault on malformed input.
  static DelaunayTriangulation from_ascii(FILE* file);
};


#endif
