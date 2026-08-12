#ifndef DELAUNAY_HH
#define DELAUNAY_HH

#include "triangulation.hh"
#include "geometry.hh"
#include "delaunay_view.hh"        // DelaunayView, Diamond, dcel_capacities
#include "delaunay_transport.hh"   // DelaunayPointTrackerView, TrackerTransport

#include <array>
#include <cstdio>       // FILE* for the .idt serialization (to_ascii / from_ascii)
#include <functional>

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
// Point tracker (optional; flip-tape transport).
// See claude-projects/delaunay-fillin/DESIGN-cubic-exact-paint.md.  When
// active, flip_edge and the flat-vertex removal transport the tracked
// points through each operation's isometric planar development, so a point
// recorded anywhere on the surface keeps a valid (face, barycentric)
// location as the triangulation changes.  compact_vertices renumbers only
// vertices, so tracked state is untouched by it; the tracker is NOT
// serialized by to_ascii/from_ascii.  Inactive by default: every hook is a
// single branch -- zero cost and zero behavior change for non-tracking use.
//
// Barycentrics are anchored to HALF-EDGE SLOTS, never vertex ids: slot i of
// face f is the origin of the i-th half-edge of the cycle starting at
// f_he[f].  Delta-complex faces may repeat corner labels (self-loop and
// bigon faces arise transiently); slots are always unambiguous.
//
// Failure contract: a transport-detected failure (the clamp band, a
// degenerate development, an exhausted arena) trips the DCEL's Status latch
// and the owner converts it to a throw at its boundary (throw_on_status).
// A REMOVAL's transport failure is commit-or-nothing: the plan hooks run
// before any mutation and splice_fan early-outs on the latch, so complex and
// tracker are unchanged.  A FLIP's is not: the view's flip_edge has no
// status check between the plan hook and the first write, so the diagonal is
// rewired, the commit hook then refuses, and the complex is left poisoned --
// as it is by a throw from the DCEL surgery itself (splice_fan deep
// invariants), per the file's deep-invariant convention.
// bisect_multi_edges and split_face are NOT transport-hooked and throw
// when tracking is active.  Copies of a tracking complex snapshot the
// tracker: transport applies only to the copy you mutate.
//
// Seeding through remove_flat_vertices requires the removed vertices to be
// flat to roughly CLAMP_TOL precision (their star development must close);
// the loose flat_tol regime (~1e-2, CEPS metrics) is incompatible with
// tracking -- the seed placement throws rather than proceeding.
// ============================================================================

// The host OWNER of the tracker state: the vectors behind
// DelaunayPointTrackerView's spans plus the transport's carry scratch (the
// DelaunayStorage owner/view pattern -- storage here, algorithm on the view,
// rebind() the @inv bound establisher).  The transported state, its
// invariants and the clamp policy live on the view
// (delaunay_transport.hh); this type is storage and growth only.
struct DelaunayPointTracker {
  bool active = false;

  std::vector<int>    owned_label, owned_face, owned_next;
  std::vector<int>    owned_head, owned_tail;      // the by-face partition
  std::vector<double> owned_bary;
  std::vector<int>    owned_carry_idx, owned_carry_cand;
  std::vector<double> owned_carry_x, owned_carry_y, owned_carry_b;

  DelaunayPointTrackerView view;    // the transported state
  TransportCarry           carry;   // the plan-to-commit scratch

  // @inv bound: each span is exactly the span over its vector.  The counts
  // (view.n, the clamp accounting) survive the rebind, as the DCEL's own
  // free-list counts do.
  void rebind() {
    view.label = owned_label; view.face = owned_face; view.next = owned_next;
    view.head  = owned_head;  view.tail = owned_tail; view.bary = owned_bary;
    carry.cidx  = owned_carry_idx; carry.rcand = owned_carry_cand;
    carry.cx    = owned_carry_x;   carry.cy    = owned_carry_y;
    carry.rb    = owned_carry_b;
  }

  // Monotone growth to `points_cap` tracked points over `buckets_cap` face
  // slots; the carry follows the points' capacity (transport_carry_cap).
  // New slots carry the empty representation (-1 / 0).  Ends in rebind().
  void ensure_capacity(long points_cap, long buckets_cap) {
    if ((long)owned_label.size() < points_cap) {
      owned_label.resize(points_cap, -1);
      owned_face .resize(points_cap, -1);
      owned_next .resize(points_cap, -1);
      owned_bary .resize(3 * points_cap, 0.0);
    }
    if ((long)owned_head.size() < buckets_cap) {
      owned_head.resize(buckets_cap, -1);
      owned_tail.resize(buckets_cap, -1);
    }
    const long C = transport_carry_cap((long)owned_label.size());
    if ((long)owned_carry_idx.size() < C) {
      owned_carry_idx .resize(C, -1);
      owned_carry_cand.resize(C, -1);
      owned_carry_x   .resize(C, 0.0);
      owned_carry_y   .resize(C, 0.0);
      owned_carry_b   .resize(3 * C, 0.0);
    }
    rebind();
  }

  // Headroom for a removal pass: every removed vertex seeds exactly ONE
  // tracked point, ON TOP of whatever the caller already registered through
  // track_point.  The owner grows (it holds vectors), so a tracked removal
  // never has to refuse for want of seed slots -- which is what the
  // enable-time sizing alone would do once the caller registers more points
  // than the pass leaves survivors.
  void reserve_for_removals(long n_removals, long buckets_cap) {
    ensure_capacity((long)view.n + n_removals, buckets_cap);
  }
};

// ============================================================================
// DelaunayStorage: everything DelaunayTriangulation OWNS -- the nine vectors
// backing the view's spans, the free lists, the tracker, the diagnostics
// flag.  A plain struct so the compiler-generated memberwise copy/move is
// complete BY CONSTRUCTION: a member added here is copied, moved, and (via
// owned_tuple's arity assert) repointed, with no hand-written list to forget
// it in.
// ============================================================================
struct DelaunayStorage {
  std::vector<int>    owned_he_next, owned_he_origin, owned_he_face;
  std::vector<double> owned_he_length, owned_he_angle;
  std::vector<int>    owned_v_out;
  std::vector<double> owned_v_cone_angle;
  std::vector<int>    owned_v_orig_degree;
  std::vector<int>    owned_f_he;

  // Backing storage of the view's bounded free lists (the SpanStacks carry
  // the counts; these carry the slots, capacity-managed by the ensure_*
  // helpers: one entry per possible edge / face id).
  std::vector<int> owned_free_edges, owned_free_faces;

  DelaunayPointTracker tracker;   // value member: copies/moves with the complex

  // If > 0, remove_flat_vertices prints live-vertex progress to stderr every
  // this-many removed vertices (and a header/footer). 0 = silent. Output is
  // flushed each line so it is visible even on a slow run.
  int verbose_removal = 0;

  // The nine backing vectors in the view's canonical field order
  // (DelaunayView::to_tuple) -- the zip repoint() folds over.
  auto owned_tuple() {
    return std::forward_as_tuple(owned_he_next, owned_he_origin, owned_he_face,
                                 owned_he_length, owned_he_angle,
                                 owned_v_out, owned_v_cone_angle,
                                 owned_v_orig_degree, owned_f_he);
  }
};

// ============================================================================
// HostDelaunayWorkspace: a workspace that owns its storage -- a byte slab
// sized by the layout, with the span workspace carved into it.  IS-A
// DelaunayWorkspace, so it passes anywhere the view machinery wants one.
// Host-only (allocates); device code carves DelaunayWorkspace from its own
// arena via DelaunayWorkspace::Layout.
// ============================================================================
struct HostDelaunayWorkspace : DelaunayWorkspace {
  std::vector<std::byte> slab;
  explicit HostDelaunayWorkspace(DelaunayWorkspace::Layout l) : slab(l.bytes()) {
    static_cast<DelaunayWorkspace&>(*this) = l.make(std::span<std::byte>(slab));
  }
};

// ============================================================================
// DCEL-based intrinsic Delaunay triangulation (delta-complex).
//
// Half-edge (DCEL) representation that correctly handles multi-edges and
// self-loops.  Every edge is identified by a half-edge index, so flip_edge(h)
// is unambiguous even when multiple edges connect the same vertex pair.
//
// The counts, SoA arrays, navigation, and intrinsic-geometry methods live on
// the inherited DelaunayView (delaunay_view.hh); the vectors behind the spans
// live on the inherited DelaunayStorage.  The aliasing contract, stated once:
//
//   @inv bound:  each view span is exactly the span over its storage vector
//                (same data pointer, same size) -- established by repoint(),
//                which folds the two field tuples together and rebinds the
//                tracker's own spans.
//   @inv sized:  he_*.size() == nh, v_*.size() >= nv, f_he.size() == nf.
//
// Resizer inventory (everything that may reallocate, each ending in
// repoint()): the monotone ensure_* helpers (post-build growth, and the
// construction paths, which allocate through them before the view-level
// build), from_ascii (which establishes the arrays wholesale), and
// compact_vertices (the one shrinking operation).  The tracker's storage
// has its own monotone grower (DelaunayPointTracker::ensure_capacity, ending
// in its own rebind), driven by enable_point_tracking and track_point.
// Everything else mutates in place through the spans.
//
// Why not Owned<View> (owned.hh): that template is specialized to the RSR
// adjacency layout (neighbours/deg/twin + optional points); a tuple-driven
// generalization could subsume this type and SpanVector, and should, the day
// a third owner-with-repoint appears -- this is the second.
//
// MOVED-FROM STATE: a moved-from complex is the valid EMPTY complex
// (nv = nh = nf = 0, spans re-bound to its emptied vectors -- never
// dangling).  Consistency predicates pass it vacuously; it is
// indistinguishable from a genuinely empty complex by design.
// ============================================================================

struct DelaunayTriangulation : DelaunayView, DelaunayStorage {
  // Backward-compatible spellings of the hoisted tracker types.
  using TrackedPoint = DelaunayTrackedPoint;
  using PointTracker = DelaunayPointTracker;

  // Re-bind the inherited spans to the owned vectors (the @inv bound
  // establisher).  Must run after any owned_* reallocation and in every
  // copy/move; folds the view's field tuple onto the storage's, so the
  // field list exists in exactly two places (the two tuples), whose arity
  // equality is asserted below the class.
  void repoint() {
    auto v = DelaunayView::to_tuple();
    auto s = owned_tuple();
    [&]<std::size_t... I>(std::index_sequence<I...>) {
      ((std::get<I>(v) = std::get<I>(s)), ...);
    }(std::make_index_sequence<n_fields>{});
    nv_cap = (int)owned_v_out.size();
    nh_cap = (int)owned_he_next.size();
    nf_cap = (int)owned_f_he.size();
    free_edges.rebind(owned_free_edges);   // counts survive the rebind
    free_faces.rebind(owned_free_faces);
    tracker.rebind();                      // the tracker's own @inv bound
  }

  // --- Growth (owner-only; the view never resizes) ---
  // Monotone: ensure the array family holds at least `need` slots; never
  // shrinks.  New slots carry the DEAD-SLOT REPRESENTATION the predicates
  // read (-1 = dead half-edge / dead vertex / unassigned face; metric 0),
  // and every call ends in repoint().
  void ensure_halfedges(int need) {
    if (need <= (int)owned_he_next.size()) return;
    owned_he_next.resize(need, -1);
    owned_he_origin.resize(need, -1);
    owned_he_face.resize(need, -1);
    owned_he_length.resize(need, 0);
    owned_he_angle.resize(need, 0);
    owned_free_edges.resize(need / 2);   // DcelCapacities::free_edges_cap = nh/2
    repoint();
  }
  void ensure_faces(int need) {
    if (need <= (int)owned_f_he.size()) return;
    owned_f_he.resize(need, -1);
    owned_free_faces.resize(need);       // DcelCapacities::free_faces_cap = nf
    repoint();
  }
  void ensure_vertices(int need) {
    if (need <= (int)owned_v_out.size()) return;
    owned_v_out.resize(need, -1);
    owned_v_cone_angle.resize(need, 0.0);
    owned_v_orig_degree.resize(need, 0);
    repoint();
  }

  // Growth-shadowing allocators: the view's allocators are capacity-guarded
  // and never grow; the owner pre-ensures capacity when the free lists are
  // exhausted, then delegates.  During reduction the free lists always cover
  // demand, so these grow only under post-build insertion (split_face,
  // bisect_multi_edges, from_ascii).
  int alloc_edge() {
    if (free_edges.empty() && nh + 2 > nh_cap) ensure_halfedges(nh + 2);
    return DelaunayView::alloc_edge();
  }
  int alloc_face() {
    if (free_faces.empty() && nf + 1 > nf_cap) ensure_faces(nf + 1);
    return DelaunayView::alloc_face();
  }
  int alloc_directed_edge(int u, int v, double L) {
    if (free_edges.empty() && nh + 2 > nh_cap) ensure_halfedges(nh + 2);
    return DelaunayView::alloc_directed_edge(u, v, L);
  }
  int wire_triangle(int h0, int h1, int h2) {
    if (free_faces.empty() && nf + 1 > nf_cap) ensure_faces(nf + 1);
    return DelaunayView::wire_triangle(h0, h1, h2);
  }

  // --- Rule of 5 ---
  // A defaulted copy would leave the new object's spans aliasing the
  // SOURCE's vectors, so every copy/move re-binds via repoint().  The
  // storage base carries all owned state memberwise (see its banner), so
  // these are one-liners; copy assignment is copy-and-swap for the strong
  // guarantee (a throwing vector copy must not leave half-assigned storage
  // under live spans).
  DelaunayTriangulation() = default;
  DelaunayTriangulation(const DelaunayTriangulation& o)
      : DelaunayView(o), DelaunayStorage(o) { repoint(); }
  DelaunayTriangulation(DelaunayTriangulation&& o) noexcept
      : DelaunayView(o), DelaunayStorage(std::move(o)) {
    repoint();
    o.nv = o.nh = o.nf = 0;   // leave the source as the empty complex
    o.free_edges.clear();     // counts over emptied storage
    o.free_faces.clear();
    o.status = Status::Ok;    // the empty complex is valid
    o.status_site = nullptr;
    o.status_witness = -1;
    o.repoint();
    o.tracker.active = false; // and the empty complex tracks nothing
    o.tracker.view.reset();
  }
  DelaunayTriangulation& operator=(const DelaunayTriangulation& o) {
    if (this != &o) {
      DelaunayTriangulation tmp(o);
      *this = std::move(tmp);
    }
    return *this;
  }
  DelaunayTriangulation& operator=(DelaunayTriangulation&& o) noexcept {
    if (this != &o) {
      DelaunayView::operator=(o);
      DelaunayStorage::operator=(std::move(o));
      repoint();
      o.nv = o.nh = o.nf = 0;
      o.free_edges.clear();
      o.free_faces.clear();
      o.status = Status::Ok;
      o.status_site = nullptr;
      o.status_witness = -1;
      o.repoint();
      o.tracker.active = false;
      o.tracker.view.reset();
    }
    return *this;
  }

  // Activate transport (sizes the per-face index).  Call before the
  // operations whose points you want carried (e.g. before
  // remove_flat_vertices; the compute(..., track_removed) overload does this
  // and seeds every removed flat vertex with label = its input-label).
  // One tracking session per complex.  Sizes the tracker's storage at
  // tracked_points_cap(nv) points over the face capacity (removal seeds one
  // point per removed flat vertex; track_point grows past it).
  // @throws std::runtime_error when the tracker already carries points
  //         (stale state would be silently transported).
  void enable_point_tracking();

  // Register a point at (face, barycentric); returns its index in
  // tracker.view.
  // @pre tracker.active; face live; b finite, >= -CLAMP_TOL componentwise,
  //      summing to 1 within 1e-9 (clamped + renormalized on registration)
  // @throws std::runtime_error on any violated precondition.
  int track_point(int label, int face, double b0, double b1, double b2);

  // (Navigation -- twin/edge/prev/dest/alive/is_bigon/cw/ccw, face_halfedges/
  // face_vertices, vertex_degree, and the incident()/edges() ranges -- is
  // inherited from DelaunayView.)

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

  // (Intrinsic geometry -- diamond, recompute_face_angles/recompute_all_angles/
  // recompute_all_face_angles/recompute_cone_angles, vertex_angle_sum, and
  // is_flat -- is inherited from DelaunayView.)

  // Discrete Gaussian curvature K(v) = 2*pi - cone angle, per live vertex (the angle
  // defect). Reads the cached v_cone_angle, so requires it current.
  std::vector<double> curvature() const;

  // --- Delaunay operations ---
  // (is_delaunay_edge / is_delaunay / count_non_delaunay and the
  // workspace-based sweeps are inherited from DelaunayView.  The wrappers
  // below keep the historical allocation-free-caller-facing signatures:
  // they materialize a workspace, select the transport policy by
  // tracker.active, run the view machinery, and convert a non-Ok status to
  // the documented throws.)

  // Flip the diagonal of the diamond around edge h, transporting tracked
  // points when the tracker is active.  @post as DelaunayView::flip_edge.
  bool flip_edge(int h);
  // Global Delaunay restore, seeding with every live edge.
  // @post is_delaunay()
  // @throws std::runtime_error when the 200*nv flip budget is exhausted
  //         (impossible by the energy argument, paper thm:lawson-converge;
  //         fail-loud backstop).
  int  lawson_sweep();
  int  flip_to_delaunay();

  // --- Vertex removal ---
  void remove_flat_vertex(int v);
  // Remove every flat vertex (is_flat(v, flat_tol)), leaving the cones. See
  // is_flat for choosing flat_tol (1e-6 for an equilateral-derived metric,
  // ~1e-2 for a CEPS metric).
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

  // The exact-regime removal driver (paper sec:exactness): derives the Lsq
  // carry from the current he_length values and instantiates
  // ExactIntegerMetric, so every decision -- flips, ears, flatness (the
  // integer cone excess eps_v == 0), and the tie-break side order -- is a
  // function of integer arithmetic alone.  The equilateral compute(T)
  // calls this; the banded remove_flat_vertices above remains the
  // general-metric driver.
  // Throws as remove_flat_vertices, plus InvariantViolated on a corrupt
  // Lsq carry (lattice development or flip placement).
  // @throws std::runtime_error when the current metric is not an exact one:
  //         a live he_length that does not square to a positive in-envelope
  //         integer, or a live vertex whose float curvature disagrees with
  //         its integer cone excess (e.g. split_face / bisect-inserted
  //         vertices) -- loud, never a silently wrong reduction.
  void remove_flat_vertices_exact(const std::function<void(int)>& on_pop = {});

  // Renumber the live vertices to 0..n_live-1 (dropping removed ones) and
  // shrink nv, rewriting he_origin and the per-vertex arrays. Needed after a
  // removal that leaves live vertices scattered (the metric-based path, which
  // unlike compute(T) does not sort flats last), because AlexandrovSolver
  // sizes its system by nv and assumes every index 0..nv-1 is a live cone.
  // Returns new_to_old: the original index of each surviving vertex.
  std::vector<int> compact_vertices();

  // --- Edge/face allocation ---
  // (alloc_edge / alloc_face / alloc_directed_edge / wire_triangle are the
  // growth-shadowing owner versions defined inline above; dealloc_edge /
  // dealloc_face are inherited from DelaunayView.)

  // Allocate a vertex with the given cone angle and original degree, growing the
  // three per-vertex arrays (v_out, v_cone_angle, v_orig_degree) in lockstep.
  // Returns its id; v_out is left at -1 (the caller wires its edges).
  int  alloc_vertex(double cone_angle, int orig_degree);

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
  // (paper thm:main).  Runs the EXACT integer metric (paper sec:exactness):
  // every decision through the reduction -- the tie-break side order
  // included -- is a function of integer arithmetic alone, and every
  // output he_length is the square root of an integer.  The Lsq carry is
  // transient: post-compute owner surgery (flip_edge / lawson_sweep /
  // bisect_multi_edges / split_face) runs the banded predicates on the
  // sqrt-of-integer lengths.
  // @post result.nv == number of deg != 6 vertices of T, all live
  // @post result.is_delaunay()
  static DelaunayTriangulation compute(const Triangulation& T);

  // compute(T) with the general-metric (banded float) predicates applied
  // to the same equilateral input.  The frozen cross-reference for the
  // parallel-primitives mirrors and the exact-vs-banded corpus verdict;
  // production callers use compute(T).
  static DelaunayTriangulation compute_banded(const Triangulation& T);

  // As compute(T), but for a PRESCRIBED intrinsic metric `length(u,v)`,
  // running the BANDED float regime.  Flat vertices are identified by cone
  // angle (is_flat), not degree, so it works on any surface whose flat
  // vertices need not be degree-6. Removes the flat vertices, leaving the
  // cone vertices in a delta-complex iDT ready for AlexandrovSolver.  On
  // unit lengths it agrees with compute(T) wherever the tolerance slivers
  // are empty (the corpus-verified case); the regimes differ, so agreement
  // is an empirical verdict, not an identity.
  // flat_tol sets the cone/flat threshold (see is_flat): 1e-6 for an
  // equilateral-derived metric, ~1e-2 for a numerically solved one.
  // If new_to_old is non-null it receives compact_vertices()' map: the
  // original T-label of each surviving cone vertex 0..nv-1.
  // With track_removed, point tracking is enabled before the reduction and
  // every removed flat vertex is seeded as a tracked point (label = its
  // T-label; removal never renumbers, so labels are the input's), ending
  // located in the result's cells -- see the tracker banner above.
  static DelaunayTriangulation compute(const Triangulation& T,
                                       const EdgeLengthFn& length,
                                       double flat_tol = 1e-6,
                                       std::vector<int>* new_to_old = nullptr,
                                       bool track_removed = false);

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
  //
  // ORDER CONTRACT: the result is 1:1 with `sources` and in the same order --
  // result[i].source == sources[i] for every i (each source seeds its own disk at the
  // minimal distance 0, so it always owns itself). Callers may index the two in lockstep.
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

  // (min_live_degree is inherited from DelaunayView.  A value below 3 is
  // one (but not the only) non-simplicial signature -- use is_simplicial()
  // for a complete check.)

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
  // Exact-integer predicate, valid on the Eisenstein-lattice domain
  // (equilateral metrics and their reductions); a diamond whose faces fit
  // no lattice answers FALSE -- use the tol overload for general metrics.
  bool is_cocircular_edge(int h) const;

  // Per-half-edge cocircular mask: tight[h] == tight[twin(h)]; dead half-edges
  // are false.  O(num_edges) integer-arithmetic predicates (lattice domain,
  // as is_cocircular_edge).
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

  // (check_consistency -- the nine numbered class invariants -- is
  // inherited from DelaunayView.)

  // Convert a non-Ok view status to the documented throws (the owner error
  // boundary; message = status name + operation).  No-op on Ok.
  void throw_on_status(const char* op) const;

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

// ---------------------------------------------------------------------------
// Contract sanity: the storage tuple's arity must match the view's field
// count (repoint() folds them index-by-index), and the inheritance edges of
// the owner pattern must hold.
// ---------------------------------------------------------------------------
static_assert(std::tuple_size_v<decltype(std::declval<DelaunayStorage>().owned_tuple())>
                  == DelaunayView::n_fields,
              "DelaunayStorage::owned_tuple arity must equal DelaunayView::n_fields");
static_assert(std::is_base_of_v<DelaunayView, DelaunayTriangulation>);
static_assert(std::is_base_of_v<DelaunayStorage, DelaunayTriangulation>);

#endif
