#pragma once

// DelaunayView -- trivially-copyable span view of the delta-complex DCEL.
//
// The owning type is DelaunayTriangulation (delaunay.hh), which inherits this
// view and points the spans into its own vectors (the graphview.hh migration
// pattern).  The view's OWN code allocates nothing, throws nothing, and does
// no I/O.  Since the exact-metric landing it includes eisenstein.hh (via
// delaunay_geometry.hh), which is NOT yet device-clean -- it brings <vector>,
// throwing host conveniences, and a global using-directive, though every
// primitive the view actually calls is scalar, header-inline, and
// exception-free.  Splitting the device-clean Eisenstein core is queued
// refactor debt (parallel-primitives refactor-debt.md,
// 2026-07-26-eisenstein-device-core); until it lands, a GPU translation unit
// cannot include this header alone.
//
// Twin convention:      twin(h) = h^1  (half-edges 2k, 2k+1 are twins).
// Face orientation:     he_next traverses each face CCW.
// Vertex circulation:   cw(h) rotates CW around origin(h).
//
// CONST DOES NOT PROTECT THE ARRAYS: the members are span<T>, so element
// access through a const view (or a const owner) yields mutable references.
// The vectors' old const-propagation is gone; a read-only method mutating
// the complex is no longer a compile error.  Const-ness of a view is
// documentation, stated here because the compiler no longer states it.
//
// Scope: the SoA arrays, navigation, the intrinsic geometry, the Diamond
// predicates, the capacity formulas, check_consistency, the canonical field
// order, the view-level build (build_topology / from_triangulation over a
// caller-supplied arc-to-halfedge map) -- and the whole mutation machinery
// (allocation over the bounded free lists, flips, sweeps, vertex removal)
// over DelaunayWorkspace, with the Status latch as the run-path error
// channel, the Transport policy carrying the point tracker's hooks, and the
// Metric policy selecting the predicate regime (exact integer vs banded
// float).  Owner-level (delaunay.hh): storage, growth, the tracked
// TrackerTransport policy, serialization, and every allocating convenience.
//
// (See the device caveat in the banner above: eisenstein.hh, reached via
// delaunay_geometry.hh, is the one non-device-clean include.)

#include <algorithm>
#include <array>
#include <cmath>
#include <concepts>
#include <cstddef>
#include <cstdint>
#include <numbers>
#include <optional>
#include <span>
#include <stdexcept>
#include <string>
#include <tuple>
#include <type_traits>
#include <utility>

#include "fullerenes/span_vector.hh"   // Spanify::SpanStack, Spanify::BitSpan

// The DCEL's element type is int, bound as the lib's 32-bit index width
// (PROMOTION-DESIGN.md Q4): machine-check the premise where the arrays live.
static_assert(std::is_same_v<int, std::int32_t>,
    "DelaunayView assumes int is the lib's 32-bit index type");

#include "fullerenes/delaunay_geometry.hh"  // Diamond, DiamondSq, FanPolygon, delaunay_detail
#include "fullerenes/eisenstein.hh"         // Eisenstein, wedge, apex_factor,
                                            // place_third_eis_total (exact metric)

// ============================================================================
// dcel_capacities(nv0): the workspace sizes a bounded DCEL over a closed
// genus-0 triangulation with nv0 vertices is built to (nh = 2E = 6*nv0 - 12,
// nf = 2*nv0 - 4; the free lists hold at most one entry per edge / face).
// Reduction never exceeds the build sizes (flips rewire in place; removal is
// net-decreasing), so these are also the lifetime capacities of the whole
// compute() pipeline.  General delta-complexes (arbitrary cone count / genus,
// or post-build vertex insertion) do NOT satisfy these identities and must
// size explicitly.
// ============================================================================
struct DcelCapacities {
  long nh_cap;          // half-edge slots   = 6*nv0 - 12
  long nf_cap;          // face slots        = 2*nv0 - 4
  long free_edges_cap;  // recycled edge ids = nh_cap / 2  (dead until stage 3)
  long free_faces_cap;  // recycled face ids = nf_cap      (dead until stage 3)
};

// @pre nv0 >= 3 (the smallest closed triangulation); smaller inputs -- e.g.
// an unfilled batch slot -- yield all-zero capacities rather than the
// negative values the Euler identity would produce (which a size_t cast
// would turn into near-SIZE_MAX allocations).
constexpr DcelCapacities dcel_capacities(int nv0) {
  if (nv0 < 3) return { 0, 0, 0, 0 };
  const long nh = 6L * nv0 - 12;
  const long nf = 2L * nv0 - 4;
  return { nh, nf, nh / 2, nf };
}

// Pin the Euler identity on the icosahedron (nv0 = 12: E = 30, F = 20) and
// the degenerate guard.
static_assert(dcel_capacities(12).nh_cap == 60 && dcel_capacities(12).nf_cap == 20
              && dcel_capacities(12).free_edges_cap == 30
              && dcel_capacities(12).free_faces_cap == 20,
              "dcel_capacities must satisfy the Euler identity");
static_assert(dcel_capacities(0).nh_cap == 0 && dcel_capacities(2).nh_cap == 0,
              "dcel_capacities must guard nv0 < 3");

// The build sizes of the DCEL over ANY closed triangulation with degree sum
// deg_sum = sum_v deg(v) = 2E: nh = 2E half-edge slots, nf = 2E/3 faces
// (3F = 2E, every face a triangle).  dcel_capacities is this specialized to
// the connected genus-0 class, where Euler fixes 2E = 6*nv0 - 12.
constexpr DcelCapacities dcel_build_capacities(long deg_sum) {
  return { deg_sum, deg_sum / 3, deg_sum / 2, deg_sum / 3 };
}
static_assert(dcel_build_capacities(6L * 12 - 12).nh_cap == dcel_capacities(12).nh_cap
              && dcel_build_capacities(6L * 12 - 12).nf_cap == dcel_capacities(12).nf_cap,
              "genus-0 specialization: degree sum 6*nv0-12 is dcel_capacities");

// ============================================================================
// DelaunayWorkspace: the bounded scratch the mutation machinery runs over.
// Pure aggregate of spans and counts; carve from one arena via Factory, or
// host-own via DelaunayWorkspaceOwner (delaunay.hh).
// ============================================================================
struct DelaunayWorkspace {
  int k_max = 0;

  // Fan scratch (extract_fan / ear_clip_fan) -- capacity k_max.  fan.P is
  // the exact lattice development (develop_fan_lattice; exact metric only).
  FanPolygon fan;
  std::span<int> poly;        // [k_max]     ear-clip working polygon; also the
                              //             tie-break's sector-arc scratch
  std::span<int> rpoly;       // [k_max]     residual-polygon scratch
  FanTriangulation tri;       // diagonals/triangles, both capacity k_max

  // Driver scratch (remove_flat_vertices).  Every collector is a SpanStack,
  // so capacity guarding and the calibration peak() are uniform.
  Spanify::SpanStack<int> work;            // cap nv0
  Spanify::BitSpan        on_queue;        // nv0 bits
  Spanify::SpanStack<int> ring;            // cap 2*k_max (both sides of cleanup)
  Spanify::SpanStack<int> tie_side[2];     // cap 2*k_max each (tie-break side fans)
  Spanify::SpanStack<int> new_faces;       // cap k_max: the surgery's face list
  Spanify::SpanStack<int> lawson_stack;    // cap nh0/2
  Spanify::BitSpan        sweep_in_stack;  // nh0/2 bits

  // ------------------------------------------------------------------------
  // Layout: the workspace's shape -- three sizes and the one function that
  // lays the fields over a cursor.  Carves one workspace from a caller byte
  // arena
  // (device-legal, no allocation).  k_max bounds the star degree seen during
  // reduction; on a delta-complex deg(v) can reach the outgoing-arc count,
  // so k_max = nh0 is the always-sufficient (and, with the O(k) layout,
  // cheap) choice; k_max = 0 suffices for sweep-only use.
  //
  // ONE layout list: lay() carves the fields in order over a Carve cursor.
  // A counting cursor (null base) yields the exact byte total; a real arena
  // yields the workspace, bounds-checked -- bytes() and make() cannot drift.
  // ------------------------------------------------------------------------
  struct Layout {
    int  nv0   = 0;   // vertex count of the input triangulation
    int  k_max = 0;   // star-degree capacity (see banner)
    long nh_explicit = -1;   // override for complexes grown past the build
                             // size (bisect/split); -1 = the Euler formula

    long nh0() const {
      return nh_explicit >= 0 ? nh_explicit : dcel_capacities(nv0).nh_cap;
    }

    // Bump cursor over one arena; null base = counting mode.
    struct Carve {
      std::byte*  base = nullptr;
      std::size_t off  = 0, cap = 0;
      bool        ok   = true;
      template <class T> std::span<T> take(long n) {
        off = (off + alignof(T) - 1) / alignof(T) * alignof(T);
        std::byte* q = base ? base + off : nullptr;
        off += (std::size_t)n * sizeof(T);
        ok = ok && (!base || off <= cap);
        return {reinterpret_cast<T*>(q), (std::size_t)n};
      }
    };

    // THE field layout (the only list).
    DelaunayWorkspace lay(Carve& c) const {
      const long lawson_cap     = nh0() / 2;
      const long on_queue_words = Spanify::BitSpan::words_for(nv0);
      const long sweep_words    = Spanify::BitSpan::words_for(lawson_cap);

      DelaunayWorkspace ws;
      ws.k_max = k_max;
      ws.fan.nb        = c.take<int>(k_max);
      ws.fan.spoke_he  = c.take<int>(k_max);
      ws.fan.inner_rim = c.take<int>(k_max);
      ws.fan.spokes    = c.take<double>(k_max);
      ws.fan.rims      = c.take<double>(k_max);
      ws.fan.cum       = c.take<double>(k_max + 1);
      ws.fan.P         = c.take<Eisenstein>(k_max);
      ws.poly          = c.take<int>(k_max);
      ws.rpoly         = c.take<int>(k_max);
      ws.tri.diagonals = c.take<FanTriangulation::Diagonal>(k_max);
      ws.tri.triangles = c.take<FanTriangulation::Triangle>(k_max);
      ws.work          = Spanify::SpanStack<int>(c.take<int>(nv0));
      ws.on_queue      = Spanify::BitSpan(c.take<std::uint64_t>(on_queue_words));
      ws.ring          = Spanify::SpanStack<int>(c.take<int>(2L * k_max));
      ws.tie_side[0]   = Spanify::SpanStack<int>(c.take<int>(2L * k_max));
      ws.tie_side[1]   = Spanify::SpanStack<int>(c.take<int>(2L * k_max));
      ws.new_faces     = Spanify::SpanStack<int>(c.take<int>(k_max));
      ws.lawson_stack  = Spanify::SpanStack<int>(c.take<int>(lawson_cap));
      ws.sweep_in_stack = Spanify::BitSpan(c.take<std::uint64_t>(sweep_words));
      return ws;
    }

    std::size_t bytes() const {
      Carve c;
      lay(c);
      return c.off;
    }

    // On an undersized arena the carve fails loudly: the returned workspace
    // has k_max = 0, so the first fan extraction trips CapacityExceeded.
    DelaunayWorkspace make(std::span<std::byte> arena) const {
      Carve c{arena.data(), 0, arena.size()};
      DelaunayWorkspace ws = lay(c);
      return c.ok ? ws : DelaunayWorkspace{};
    }
  };
};

// ============================================================================
// Transport policy: the point-tracker's hook points on the two
// topology-changing operations.  NoTransport (the default, and the only
// policy device code sees) is a set of empty inline no-ops -- the
// instantiated bodies are the untracked operations exactly.  The owner's
// tracked path supplies a host-side policy (TrackerTransport, delaunay.cc)
// reproducing the flip-tape transport through the same surgery bodies.
// Plan hooks run BEFORE any mutation and may throw (host policies only);
// commit hooks run after the surgery and must not fail.
// ============================================================================
struct NoTransport {
  static constexpr bool tracking() { return false; }
  void plan_flip(int /*h*/, int /*fh*/, int /*ft*/) {}
  void commit_flip(int /*fh*/, int /*ft*/) {}
  void plan_star(int /*v*/, const FanPolygon&) {}
  void plan_charts(const FanPolygon&, const FanTriangulation&) {}
  void commit_removal(std::span<const int> /*new_faces*/, int /*seed_label*/) {}
};

// Observer policy for remove_flat_vertices: on_pop fires at every live-flat
// work-list pop (the documented observation hook), on_removed after each
// successful removal (the diagnostics hook).  Must not mutate the mesh.
struct NoRemovalObserver {
  void on_pop(int /*v*/) {}
  void on_removed(int /*v*/) {}
};

// ============================================================================
// Metric policy: compile-time selection of the predicate regime, threaded
// beside Transport through the flip/sweep/removal machinery (paper
// sec:exactness).  Each policy carries the COMPLETE predicate-and-length
// vocabulary the machinery consults -- flatness, the three edge
// predicates, the flipped length, the ear acceptance, the star
// preparation, the tie-break side order, and the edge-length write -- so
// the algorithms contain no regime branches at all; the two regimes are
// two small objects a proof can cite separately.  Definitions follow
// DelaunayView (their operations read it).
//
// The pipeline entry point is the selector -- compute(T) (equilateral
// input) instantiates ExactIntegerMetric, compute(T, lengths)
// BandedFloatMetric -- so the exact predicates are not callable by
// convention on a metric they do not apply to; the exact owner driver
// additionally derives and VERIFIES its carry (remove_flat_vertices_exact).
// ============================================================================
struct BandedFloatMetric;
struct ExactIntegerMetric;

// The input surface the view-level build reads: a triangulation presented
// as a ROTATION SYSTEM -- row u is the CCW cyclic order of u's neighbours,
// find(v,u) locates the reverse arc, and arcid(u,i) is the dense arc index
// the build's arc-to-halfedge map is keyed by.  Satisfied by
// TriangulationView and by a Triangulation owner; a concept rather than a
// concrete parameter type so this header stays free of the graphview.hh
// include set (the device-TU concern in the banner above).  The semantic
// side -- rings CCW-consistent and closed -- is what
// GraphView::is_consistently_oriented() checks; the build's guards catch
// the violations it cannot assume (unpaired arcs, non-triangle corners),
// not all of them.
template<typename TriView>
concept oriented_triangulation = requires(const TriView& T) {
  T.N;  T.dmax;  T[0];
  { T.find(0, 0) }  -> std::convertible_to<int>;
  { T.arcid(0, 0) } -> std::convertible_to<std::size_t>;
  { T.degree(0) }   -> std::convertible_to<int>;
};

// The dense arc-index space of T: arcid(u,i) ranges over [0, arc_space(T)).
template<oriented_triangulation TriView>
constexpr long arc_space(const TriView& T) { return (long)T.N * T.dmax; }

// Degree sum 2E = sum_v deg(v): the DCEL's half-edge count for T.
template<oriented_triangulation TriView>
constexpr long degree_sum(const TriView& T) {
  long s = 0;
  for (int v = 0; v < (int)T.N; v++) s += T.degree(v);
  return s;
}

// ============================================================================
// DelaunayView: span SoA view of the half-edge DCEL.
// ============================================================================
struct DelaunayView {
  // The run-path error channel: a sticky first-failure latch (the view never
  // throws).  First failure wins; every mutating operation below early-outs
  // once status != Ok, so one bad operation cannot cascade into corrupted
  // state or clobber the original diagnostic.  The owner converts a non-Ok
  // status to its documented throws at its boundary.
  enum class Status : int { Ok, CapacityExceeded, BudgetExceeded, InvariantViolated };
  Status      status         = Status::Ok;
  const char* status_site    = nullptr;   // string literal at the trip site
  int         status_witness = -1;        // offending h / v / size

  // Latch the first failure together with the obligation it violated; later
  // trips never clobber the original diagnostic (first-failure-wins is
  // structural, not by-inspection).  Returns false so guard sites can
  // `return trip(...)`.  A non-Ok status is TERMINAL for the complex: no
  // reset exists; recover by re-assigning from a clean complex (which copies
  // the clean status).
  bool trip(Status s, const char* site, int witness = -1) {
    if (status == Status::Ok) {
      status = s;
      status_site = site;
      status_witness = witness;
    }
    return false;
  }

  // --- Counts (dead slots included in nh/nf) ---
  int nv = 0;   // live vertices
  int nh = 0;   // allocated half-edges (including dead slots)
  int nf = 0;   // allocated faces (including dead slots)

  // --- Capacities (fixed at bind time; the view never resizes) ---
  int nv_cap = 0, nh_cap = 0, nf_cap = 0;

  // --- Half-edge topology (indexed 0..nh-1) ---
  std::span<int>    he_next;    // next half-edge CCW in same face
  std::span<int>    he_origin;  // origin vertex (-1 = dead)
  std::span<int>    he_face;    // face to the left

  // --- Half-edge metric (indexed 0..nh-1) ---
  std::span<double> he_length;  // edge length (same for h and twin(h))
  std::span<double> he_angle;   // angle at origin(h) in face(h)

  // --- Per-vertex (indexed 0..nv-1) ---
  std::span<int>    v_out;          // one outgoing half-edge (-1 = dead vertex)
  std::span<double> v_cone_angle;   // cone angle = total vertex angle
  std::span<int>    v_orig_degree;  // degree in the original triangulation

  // --- Per-face (indexed 0..nf-1) ---
  std::span<int>    f_he;       // one boundary half-edge (-1 = dead face)

  // --- Free lists (bounded LIFO == the historical vectors' LIFO use;
  //     capacity one entry per edge / per face -- DcelCapacities'
  //     free_*_cap -- so overflow is unreachable in contract, loud if
  //     reached) ---
  Spanify::SpanStack<int> free_edges;  // recycled edge slots (half-edge id / 2)
  Spanify::SpanStack<int> free_faces;  // recycled face slots

  // -------------------------------------------------------------------------
  // Navigation.
  // -------------------------------------------------------------------------
  int  twin(int h)  const { return h ^ 1; }
  int  edge(int h)  const { return h >> 1; }
  int  prev(int h)  const { return he_next[he_next[h]]; }  // only for triangulations
  int  dest(int h)  const { return he_origin[twin(h)]; }
  bool alive(int h) const { return he_origin[h] >= 0; }

  // Bigon edge: both half-edges of h bound the same face.  Arises in
  // delta-complexes around low-degree cone vertices (an i-j edge of an "iji"
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
  // The cycle slot (0, 1, 2) of half-edge h within its face's 3-cycle.
  // he_origin[h] is corner SLOT cycle_slot(h) of any chart built by the
  // cell-metric construction (eisenstein_paint's cell_metric reads its
  // corners off this same cycle) -- the label-free identification that
  // stays sound when corner LABELS repeat on delta-complex cells.
  // Pre: h is live and on a well-formed 3-cycle.
  int cycle_slot(int h) const {
    const auto cyc = face_halfedges(he_face[h]);
    for (int k = 0; k < 3; k++)
      if (cyc[k] == h) return k;
    throw std::logic_error("DelaunayView::cycle_slot: half-edge "
                           + std::to_string(h) + " not on face "
                           + std::to_string(he_face[h]) + "'s 3-cycle");
  }

  // Count outgoing half-edges from v; 0 for a dead vertex.
  int vertex_degree(int v) const {
    int deg = 0;
    for ([[maybe_unused]] int h : incident(v)) deg++;   // empty when v_out[v] < 0
    return deg;
  }

  // Range over the outgoing half-edges of v (the cw ring from v_out[v]);
  // empty if v has no incident live edge.  The canonical vertex traversal.
  // SELF-CONTAINED VALUE: the range and its iterator carry the one span the
  // walk needs (cw(h) = he_next[h^1]), so a range taken from a temporary
  // view -- e.g. a by-value batch entry -- never dangles.
  struct IncidentHalfEdges {
    std::span<const int> he_next; int h0;
    struct iterator {
      std::span<const int> he_next; int start, h; bool done;
      int       operator*()  const { return h; }
      iterator& operator++()       { h = he_next[h ^ 1]; done = (h == start); return *this; }
      bool      operator!=(const iterator&) const { return !done; }
    };
    iterator begin() const { return {he_next, h0, h0, h0 < 0}; }
    iterator end()   const { return {}; }
  };
  IncidentHalfEdges incident(int v) const { return {he_next, v_out[v]}; }

  // Range over one (even) half-edge per live edge.  The canonical edge
  // traversal, skipping dead slots.  Self-contained value, like incident():
  // alive(h) = he_origin[h] >= 0 is the only fact the walk reads.
  struct LiveEdges {
    std::span<const int> he_origin; int nh;
    struct iterator {
      std::span<const int> he_origin; int h, nh;
      void advance() { while (h < nh && he_origin[h] < 0) h += 2; }
      int       operator*()  const { return h; }
      iterator& operator++()       { h += 2; advance(); return *this; }
      bool      operator!=(const iterator& o) const { return h != o.h; }
    };
    iterator begin() const { iterator it{he_origin, 0, nh}; it.advance(); return it; }
    iterator end()   const { return {he_origin, nh, nh}; }
  };
  LiveEdges edges() const { return {he_origin, nh}; }

  // -------------------------------------------------------------------------
  // Intrinsic geometry.
  // -------------------------------------------------------------------------

  // The five half-edges whose lengths form the diamond of h, with the same
  // field roles as Diamond.  h: u->v; face left of h has third vertex
  // B = dest(next(h)); the twin face has D = dest(next(twin(h))).
  struct DiamondArcs { int e, a, b, c, d; };
  DiamondArcs diamond_arcs(int h) const {
    int t = twin(h);
    int h_vB = he_next[h];          // v->B  (side b, adjacent to v)
    int h_Bu = he_next[h_vB];       // B->u  (side a, adjacent to u)
    int h_uD = he_next[t];          // u->D  (side c, adjacent to u)
    int h_Dv = he_next[h_uD];       // D->v  (side d, adjacent to v)
    return {h, h_Bu, h_vB, h_uD, h_Dv};
  }

  // The five-length local geometry around edge h.
  Diamond diamond(int h) const {
    auto A = diamond_arcs(h);
    return {he_length[A.e], he_length[A.a], he_length[A.b],
            he_length[A.c], he_length[A.d]};
  }

  // The same diamond in exact squared-length coordinates, read from the
  // metric's Lsq carry (never from rounded doubles).
  DiamondSq diamond_sq(int h, std::span<const long long> Lsq) const {
    auto A = diamond_arcs(h);
    return {Lsq[A.e], Lsq[A.a], Lsq[A.b], Lsq[A.c], Lsq[A.d]};
  }

  // Recompute the three corner angles of face f from its edge lengths.
  void recompute_face_angles(int f) {
    if (f_he[f] < 0) return;
    // h_i: u_i -> u_{i+1} with length L_i.  Angle at origin(h_i) is the
    // corner between sides L_i (outgoing) and L_{i-1} (incoming), opposite
    // to L_{i+1}.
    const auto h = face_halfedges(f);
    double L[3] = { he_length[h[0]], he_length[h[1]], he_length[h[2]] };
    for (int i = 0; i < 3; i++)
      he_angle[h[i]] = delaunay_detail::triangle_angle(L[i], L[(i + 2) % 3], L[(i + 1) % 3]);
  }

  // Recompute every corner angle (he_angle) from the current edge lengths,
  // then refresh the per-vertex cone-angle cache (v_cone_angle) that
  // curvature() / is_flat() / remove_flat_vertices() read.  Both are derived
  // from he_length, so this is THE entry point that re-establishes angle
  // coherence after any change to the metric.  Delaunay flips do NOT need
  // it: the cone angle is flip-invariant (the quad's interior angle at each
  // diamond vertex is independent of which diagonal is present -- paper
  // lem:coneflip), so flip_edge keeps v_cone_angle correct on its own.
  void recompute_all_angles() {
    recompute_all_face_angles();
    recompute_cone_angles();
  }

  // Recompute he_angle for every face, WITHOUT refreshing the v_cone_angle
  // cache.  For a hot caller (a line search) that reads only he_angle per
  // trial; pair with recompute_cone_angles once the kept metric needs the
  // cone cache.
  void recompute_all_face_angles() {
    for (int f = 0; f < nf; f++)
      recompute_face_angles(f);
  }

  // Refresh the cone-angle cache at one vertex.  @pre he_angle current.
  void recompute_cone_angle(int v) {
    v_cone_angle[v] = vertex_angle_sum(v);
  }

  // The equilateral metric: every edge length 1, every corner angle pi/3,
  // cone angle Theta_v = deg(v) * pi/3 (so kappa_v = (6 - deg v) * pi/3).
  // Writes exactly the nh live-slot prefix and the nv vertices.
  void set_equilateral_metric() {
    constexpr double pi = std::numbers::pi_v<double>;
    for (int h = 0; h < nh; h++) he_length[h] = 1.0;
    for (int h = 0; h < nh; h++) he_angle[h]  = pi / 3.0;
    for (int v = 0; v < nv; v++)
      v_cone_angle[v] = v_orig_degree[v] * pi / 3.0;   // (deg*pi)/3: the
                                  // reference association, kept to the ULP
  }

  // Refresh the cone-angle cache v_cone_angle[v] = vertex_angle_sum(v) for
  // every live vertex.  @pre he_angle current.  O(nh).
  void recompute_cone_angles() {
    for (int v = 0; v < nv; v++)
      if (v_out[v] >= 0)
        recompute_cone_angle(v);
  }

  // Total angle at vertex v = sum of incident corner angles (the paper's
  // Theta_v).  @pre angles current (recompute_all_angles after any length
  // change).
  double vertex_angle_sum(int v) const {
    // he_angle[h] is the angle at origin(h) in face(h): one corner per
    // outgoing half-edge.  incident(v) is empty when v_out[v] < 0.
    double sum = 0.0;
    for (int h : incident(v)) sum += he_angle[h];
    return sum;
  }

  // Discrete Gaussian curvature (angle defect) at v: kappa_v = 2*pi - Theta_v.
  // @pre v_cone_angle current (recompute_all_angles after any length change).
  double curvature(int v) const {
    return delaunay_detail::two_pi - v_cone_angle[v];
  }

  // Total curvature Sigma_v kappa_v over the live vertices.  Gauss-Bonnet:
  // equal to 4*pi on any closed genus-0 complex -- the cheapest global check
  // of the cone-angle cache, exposed for tests and validators.
  double total_curvature() const {
    double sum = 0.0;
    for (int v = 0; v < nv; v++)
      if (v_out[v] >= 0) sum += curvature(v);
    return sum;
  }

  // A vertex is flat (zero cone curvature, hence removable) iff
  // |kappa_v| < flat_tol.  flat_tol must separate the metric's curvature
  // noise floor from the smallest real cone curvature (pi/3 ~ 1.047): 1e-6
  // suits exact (equilateral) metrics, where any tol in (0, pi/3) agrees
  // exactly with the degree-6 test; a numerically solved metric (e.g. a
  // CEPS conformal solve, whose layout-cut seams leave ~5e-4 spurious
  // curvature at genuinely flat vertices) needs ~1e-2.
  // @pre v_cone_angle current.
  bool is_flat(int v, double flat_tol = 1e-6) const {
    return std::abs(curvature(v)) < flat_tol;
  }

  // Integer cone excess at v: eps_v = 6 - deg0(v), the curvature in units
  // of pi/3 (kappa_v = eps_v * pi/3 on the equilateral metric): +1 per
  // pentagon, 0 per hexagon, -1 per heptagon.  Exact through the WHOLE
  // reduction: cone angle is flip-invariant (paper lem:coneflip) and only
  // eps = 0 vertices are ever removed, so v_orig_degree stays the
  // authority.  @pre equilateral regime (v_orig_degree meaningful).
  int cone_excess(int v) const { return 6 - v_orig_degree[v]; }

  // Integer Gauss-Bonnet: Sigma eps_v = 12 on any closed genus-0 complex
  // ("a fullerene has twelve pentagons") -- the exact-regime form of
  // total_curvature, an integer identity rather than a 4*pi float check.
  int total_cone_excess() const {
    int sum = 0;
    for (int v = 0; v < nv; v++)
      if (v_out[v] >= 0) sum += cone_excess(v);
    return sum;
  }

  // (Metric-dispatched flatness lives on the policies: ExactIntegerMetric
  // decides by the integer cone excess -- no tolerance at all, provably
  // identical to the banded test for any tol in (0, pi/3) -- and
  // BandedFloatMetric by the curvature band it carries.)

  // -------------------------------------------------------------------------
  // Allocation over the free lists (capacity-guarded; the view never grows --
  // an owner pre-ensures capacity for genuine growth, and during reduction
  // the free lists always cover demand, so the guards are dead in contract).
  // -------------------------------------------------------------------------
  int alloc_edge() {
    if (status != Status::Ok) return -1;
    int eid;
    if (!free_edges.empty()) {
      eid = free_edges.back();
      free_edges.pop_back();
    } else {
      eid = nh / 2;
      if (nh + 2 > nh_cap) {
        trip(Status::CapacityExceeded, "alloc_edge: half-edge capacity", nh);
        return -1;
      }
      nh += 2;
      he_next[2*eid] = -1;   he_next[2*eid+1] = -1;
      he_origin[2*eid] = -1; he_origin[2*eid+1] = -1;
      he_face[2*eid] = -1;   he_face[2*eid+1] = -1;
      he_length[2*eid] = 0;  he_length[2*eid+1] = 0;
      he_angle[2*eid] = 0;   he_angle[2*eid+1] = 0;
    }
    he_origin[2*eid] = -1;
    he_origin[2*eid+1] = -1;
    return 2 * eid;
  }

  int alloc_face() {
    if (status != Status::Ok) return -1;
    int fid;
    if (!free_faces.empty()) {
      fid = free_faces.back();
      free_faces.pop_back();
    } else {
      if (nf + 1 > nf_cap) {
        trip(Status::CapacityExceeded, "alloc_face: face capacity", nf);
        return -1;
      }
      fid = nf++;
      f_he[fid] = -1;
    }
    f_he[fid] = -1;
    return fid;
  }

  void dealloc_edge(int h) {
    if (status != Status::Ok) return;
    int eid = edge(h);
    he_origin[2*eid] = -1;
    he_origin[2*eid+1] = -1;
    he_next[2*eid] = -1;
    he_next[2*eid+1] = -1;
    he_face[2*eid] = -1;
    he_face[2*eid+1] = -1;
    he_length[2*eid] = 0;
    he_length[2*eid+1] = 0;
    if (!free_edges.push_back(eid))
      trip(Status::CapacityExceeded, "dealloc_edge: free-list overflow (double dealloc?)", eid);
  }

  void dealloc_face(int f) {
    if (status != Status::Ok) return;
    f_he[f] = -1;
    if (!free_faces.push_back(f))
      trip(Status::CapacityExceeded, "dealloc_face: free-list overflow (double dealloc?)", f);
  }

  // Allocate an edge and set its endpoints and length.  Returns the
  // half-edge h with origin(h) = u, origin(twin h) = v, length L on both
  // sides; faces remain unassigned.
  int alloc_directed_edge(int u, int v, double L) {
    if (status != Status::Ok) return -1;
    int h = alloc_edge();
    if (h < 0) return h;
    he_origin[h]     = u;
    he_origin[twin(h)] = v;
    he_length[h] = he_length[twin(h)] = L;
    return h;
  }

  // Wire three half-edges into a CCW triangle face and compute its angles
  // from the stored edge lengths.  Returns the new face id.
  // @pre h0, h1, h2 have origin and length set; their endpoints chain.
  int wire_triangle(int h0, int h1, int h2) {
    if (status != Status::Ok) return -1;
    // The historical local_arc.at() threw on a missing arc; a dead arc here
    // would otherwise write he_next[-1] (paper lem:wire-triangle scope).
    if (h0 < 0 || h1 < 0 || h2 < 0) {
      trip(Status::InvariantViolated, "wire_triangle: dead arc handle", std::min({h0, h1, h2}));
      return -1;
    }
    int fid = alloc_face();
    if (fid < 0) return fid;
    he_next[h0] = h1; he_next[h1] = h2; he_next[h2] = h0;
    he_face[h0] = he_face[h1] = he_face[h2] = fid;
    f_he[fid] = h0;
    recompute_face_angles(fid);
    return fid;
  }

  // -------------------------------------------------------------------------
  // Construction (the view-level build; PROMOTION-DESIGN.md 1.4).
  // -------------------------------------------------------------------------

  // Reset to the empty complex at capacity: every slot carries the dead
  // representation (-1 topology, 0 metric), counts zero, free lists empty.
  // The single establisher of the fresh state over a caller arena whose
  // previous contents are arbitrary.  (The free lists' calibration peaks
  // are telemetry and deliberately survive -- outside the state contract.)
  void reset_to_dead_slots() {
    for (int h = 0; h < nh_cap; h++) {
      he_next[h] = -1; he_origin[h] = -1; he_face[h] = -1;
      he_length[h] = 0.0; he_angle[h] = 0.0;
    }
    for (int v = 0; v < nv_cap; v++) {
      v_out[v] = -1; v_cone_angle[v] = 0.0; v_orig_degree[v] = 0;
    }
    for (int f = 0; f < nf_cap; f++) f_he[f] = -1;
    nv = 0; nh = 0; nf = 0;
    free_edges.clear();
    free_faces.clear();
  }

  // Build phase 1: number the undirected edges 0..E-1, recording each
  // arc's half-edge slot -- 2k for u->v, 2k+1 for v->u -- in the
  // arc-to-halfedge map, the reverse arc located by T.find(v,u).  Ids are
  // assigned in the canonical order (u ascending, v in u's ring order,
  // counting only u<v arcs), so the numbering -- and every id derived from
  // it -- is bit-identical to the historical owner-level build (frozen as
  // the reference_build_topology oracle in tests/delaunay-test.cc).
  // Sets nh.
  // @post paired: the ids written form couples {2k, 2k+1} on mutually
  //       reverse arc slots; -1 marks an arc still unassigned (what
  //       phase 2 tests).
  template<oriented_triangulation TriView>
  void assign_halfedge_ids(const TriView& T, std::span<int> he_of_arc)
  {
    const long n_arcs = arc_space(T);
    for (long a = 0; a < n_arcs; a++) he_of_arc[a] = -1;
    int eid = 0;
    for (int u = 0; u < (int)T.N; u++) {
      auto row = T[u];
      const int deg = (int)row.size();
      for (int i = 0; i < deg; i++) {
        const int v = row[i];
        if (u < v) {
          const int jr = T.find(v, u);   // position of u in v's row
          if (jr < 0) {
            trip(Status::InvariantViolated,
                 "assign_halfedge_ids: u<v arc without reverse (not an "
                 "oriented triangulation)", u);
            return;
          }
          // Backstop: dead under build_topology's degree-sum prologue gate
          // plus the simple-input @pre; a malformed ring can still reach it.
          if (2 * eid + 2 > nh_cap) {
            trip(Status::CapacityExceeded,
                 "assign_halfedge_ids: half-edge capacity", 2 * eid);
            return;
          }
          he_of_arc[T.arcid(u, i)]  = 2 * eid;
          he_of_arc[T.arcid(v, jr)] = 2 * eid + 1;
          eid++;
        }
      }
    }
    nh = 2 * eid;
  }

  // Build phase 2: origins (origin(arc u->*) = u).  The id lookup also
  // catches the pairing failure phase 1 cannot see: an arc u->v with u > v
  // whose reverse is absent was never assigned an id.  (The historical
  // build wrote he_origin[-1] here.)
  // @post total: every arc slot of T carries an id.  Together with
  //       phase 1's `paired`, this is what lets phase 3 dereference w->u
  //       unguarded.
  template<oriented_triangulation TriView>
  void set_origins(const TriView& T, std::span<int> he_of_arc)
  {
    for (int u = 0; u < (int)T.N; u++) {
      auto row = T[u];
      const int deg = (int)row.size();
      for (int i = 0; i < deg; i++) {
        const int h = he_of_arc[T.arcid(u, i)];
        if (h < 0) {
          trip(Status::InvariantViolated,
               "set_origins: arc with larger origin unpaired (not an "
               "oriented triangulation)", u);
          return;
        }
        he_origin[h] = u;
      }
    }
  }

  // Build phase 3: next pointers, face assignments, and the per-vertex
  // epilogue (v_out, v_orig_degree).  For vertex u with CCW neighbours
  // [..., v, w, ...]: arc u->v lies in face (u, v, w) with w the next
  // neighbour after v, and next(u->v) = v->w.  The corner needs arc v->w
  // to exist (guarded: its absence means a face is not a triangle); w->u
  // needs no guard by phase 1 `paired` + phase 2 `total`: u->w having an
  // id forces its reverse w->u to have one.
  //
  // The CANONICAL CORNER of a face is the one whose apex is its smallest
  // vertex; by the simple-input @pre exactly one of a face's three corners
  // qualifies, so each face is created once and fids run in (smallest
  // vertex, ring position) order -- the order the reference oracle freezes.
  template<oriented_triangulation TriView>
  void wire_faces(const TriView& T, std::span<int> he_of_arc)
  {
    for (int u = 0; u < (int)T.N; u++) {
      auto row = T[u];
      const int deg = (int)row.size();
      for (int j = 0; j < deg; j++) {
        const int v = row[j], w = row[(j + 1) % deg];
        const int jw = T.find(v, w);
        if (jw < 0) {
          trip(Status::InvariantViolated,
               "wire_faces: corner (u,v,w) missing arc v->w (input face "
               "is not a triangle)", v);
          return;
        }
        const int h_uv = he_of_arc[T.arcid(u, j)];
        const int h_vw = he_of_arc[T.arcid(v, jw)];
        he_next[h_uv] = h_vw;

        if (u < v && u < w) {
          // Backstop: dead under build_topology's prologue gate + @pre.
          if (nf + 1 > nf_cap) {
            trip(Status::CapacityExceeded, "wire_faces: face capacity", nf);
            return;
          }
          const int fid = nf++;
          const int h_wu = he_of_arc[T.arcid(w, T.find(w, u))];
          he_face[h_uv] = fid;
          he_face[h_vw] = fid;
          he_face[h_wu] = fid;
          f_he[fid] = h_uv;
        }
      }
      if (deg > 0) v_out[u] = he_of_arc[T.arcid(u, 0)];
      v_orig_degree[u] = deg;   // original degree: topological, metric-free
    }
  }

  // Build the metric-independent half-edge topology from a combinatorial
  // triangulation: the composition reset_to_dead_slots ->
  // assign_halfedge_ids -> set_origins -> wire_faces.  The metric arrays
  // (he_length, he_angle, v_cone_angle) are zeroed but left for the caller
  // to fill -- the shared core of from_triangulation (equilateral) and the
  // owner's from_intrinsic_metric (prescribed lengths).
  //
  // he_of_arc: the arc-to-halfedge map, dense over [0, arc_space(T)),
  // caller-supplied.  The prologue gates ALL capacities against the
  // input's own sizes -- nv_cap >= T.N, then dcel_build_capacities(
  // degree_sum(T)); a connected genus-0 input meets dcel_capacities(T.N)
  // with equality -- so every capacity failure refuses BEFORE touching the
  // arena, and the in-phase capacity trips are backstops, dead under this
  // prologue plus the @pre below.
  //
  // @pre each ring lists distinct neighbours (simple input; multi-edge
  //      delta-complexes are built by from_ascii, never by this build).
  //      The guards detect unpaired arcs and non-triangle corners, not
  //      every malformed adjacency: an ASYMMETRIC duplicate (u twice in
  //      v's ring, v once in u's) evades them.
  // @post on Ok: check_consistency's @inv 1-6 (paper sec:dcel -- twin
  //      closure, triangular faces, origin chaining, face agreement,
  //      v_out / f_he witnesses) hold on the live prefix.  @inv 7-9 are
  //      metric facts, so check_consistency() in full is the @post of
  //      from_triangulation and of the owner's from_intrinsic_metric,
  //      never of the bare topology build.
  //
  // Starting from reset_to_dead_slots(), an Ok build is a deterministic
  // function of (T, capacities) alone, even on a reused batch arena.
  // Failures latch the Status channel and are terminal for the binding:
  // capacity refusals leave the arena untouched; an invariant trip leaves
  // a partial build behind the latch, which gates every consumer.
  template<oriented_triangulation TriView>
  void build_topology(const TriView& T, std::span<int> he_of_arc)
  {
    if (status != Status::Ok) return;
    if ((long)he_of_arc.size() < arc_space(T)) {
      trip(Status::CapacityExceeded, "build_topology: he_of_arc workspace",
           (int)he_of_arc.size());
      return;
    }
    if (nv_cap < (int)T.N) {
      trip(Status::CapacityExceeded, "build_topology: vertex capacity",
           nv_cap);
      return;
    }
    const DcelCapacities need = dcel_build_capacities(degree_sum(T));
    if (nh_cap < need.nh_cap) {
      trip(Status::CapacityExceeded, "build_topology: half-edge capacity",
           nh_cap);
      return;
    }
    if (nf_cap < need.nf_cap) {
      trip(Status::CapacityExceeded, "build_topology: face capacity",
           nf_cap);
      return;
    }

    reset_to_dead_slots();
    nv = T.N;

    assign_halfedge_ids(T, he_of_arc);
    if (status != Status::Ok) return;
    set_origins(T, he_of_arc);
    if (status != Status::Ok) return;
    wire_faces(T, he_of_arc);
  }

  // The equilateral build: topology + unit metric (every edge length 1,
  // every corner pi/3, cone angle deg*pi/3).  One call per isomer over a
  // caller arena.  Static: the view's named constructor, mirroring the
  // owner's DelaunayTriangulation::from_triangulation.
  // @post on Ok: check_consistency() in full (topology @inv 1-6 from
  //       build_topology, metric @inv 7-9 from the equilateral fill).
  template<oriented_triangulation TriView>
  static void from_triangulation(DelaunayView& D, const TriView& T,
                                 std::span<int> he_of_arc)
  {
    D.build_topology(T, he_of_arc);
    if (D.status != Status::Ok) return;
    D.set_equilateral_metric();
  }

  // -------------------------------------------------------------------------
  // Delaunay operations.
  // -------------------------------------------------------------------------
  bool is_delaunay_edge(int h) const { return diamond(h).is_delaunay(); }

  bool is_delaunay() const {
    for (int h : edges())
      if (!is_delaunay_edge(h)) return false;
    return true;
  }

  int count_non_delaunay() const {
    int count = 0;
    for (int h : edges())
      if (!is_delaunay_edge(h)) count++;
    return count;
  }

  // (Exact-regime Delaunay verdicts on a finished complex go through
  // Diamond::squared() -- the carry is transient, so no ExactIntegerMetric
  // survives compute(); see the corpus tests for the shape.)

  // Smallest degree among live vertices, or INT32_MAX if none.
  int min_live_degree() const {
    int m = INT32_MAX;
    for (int v = 0; v < nv; v++) {
      if (v_out[v] < 0) continue;
      int d = vertex_degree(v);
      if (d < m) m = d;
    }
    return m;
  }

  // Flip the diagonal of the diamond around edge h.  Accepts any edge with a
  // strictly convex diamond and a positive flipped length, including the
  // B == D case, which produces a self-loop edge at B (strictly Delaunay
  // when the flipped edge was non-Delaunay; paper lem:selfloop-delaunay).
  // @post on true:  the DCEL invariants hold with the new diagonal B-D in
  //                 place of h; only the two diamond faces are rewired; the
  //                 surface metric is unchanged (paper thm:metric-invariance).
  // @post on false: the complex is untouched.
  // @error InvariantViolated (exact regime, via m.flipped): a convex
  //        diamond without a lattice diagonal, or a diagonal past the
  //        exact envelope -- both mean the Lsq carry is corrupt.
  // The transport policy's plan hook runs before the first write (it may
  // throw with everything untouched); commit after the rewire.  The metric
  // policy decides convexity and the new diagonal's length in both
  // coordinates (Length), and owns the paired he_length / Lsq write.
  template <class Metric = BandedFloatMetric, class Transport = NoTransport>
  bool flip_edge(int h, Metric&& m = Metric{}, Transport&& tr = Transport{}) {
    if (status != Status::Ok) return false;

    int t = twin(h);
    int h1 = he_next[h], h2 = he_next[h1];
    int h4 = he_next[t], h5 = he_next[h4];
    int u = he_origin[h],  v = he_origin[t];
    int B = he_origin[h2], D = he_origin[h5];

    if (!m.convex(*this, h)) return false;
    auto f = m.flipped(*this, h);
    if (!f) return false;

    // Rewire the diagonal: h becomes B->D, t becomes D->B.  Reuse the two
    // face slots fh, ft by rewiring in place (avoids dealloc/realloc).
    int fh = he_face[h], ft = he_face[t];
    tr.plan_flip(h, fh, ft);

    he_origin[h] = B;
    he_origin[t] = D;
    m.set_edge_length(*this, h, *f);

    // Face left of h: (B, D, v) via half-edges h -> h5 -> h1 (origins B, D, v).
    he_next[h]  = h5;  he_next[h5] = h1;  he_next[h1] = h;
    he_face[h]  = he_face[h5] = he_face[h1] = fh;
    f_he[fh] = h;

    // Face left of t: (D, B, u) via half-edges t -> h2 -> h4 (origins D, B, u).
    he_next[t]  = h2;  he_next[h2] = h4;  he_next[h4] = t;
    he_face[t]  = he_face[h2] = he_face[h4] = ft;
    f_he[ft] = t;

    recompute_face_angles(fh);
    recompute_face_angles(ft);

    // u and v lost their incident diagonal; find a new outgoing half-edge.
    if (v_out[u] == h) v_out[u] = h4;
    if (v_out[v] == t) v_out[v] = h1;
    // B == D case: h and t are now self-loops at B.  he_origin[h5] is never
    // reassigned by the flip and equals D == B; h5 is wired into face fh
    // above -- a valid outgoing half-edge from B.
    if (B == D) v_out[B] = h5;

    tr.commit_flip(fh, ft);
    return true;
  }

  // Shared drain of the seeded Lawson sweep.  The Bobenko-Springborn
  // discrete Dirichlet energy strictly decreases on every flip (paper
  // thm:lawson-converge), so the sweep terminates; the 200*nv budget is a
  // fail-loud backstop, not part of the algorithm.
  template <class Metric = BandedFloatMetric, class Transport = NoTransport>
  int lawson_sweep_drain(Spanify::SpanStack<int>& S, Spanify::BitSpan& in_stack,
                         Metric&& m = Metric{}, Transport&& tr = Transport{}) {
    if (status != Status::Ok) return 0;
    int flips = 0;
    int budget = 200 * nv;
    while (!S.empty()) {
      if (budget <= 0) {
        trip(Status::BudgetExceeded,
             "lawson_sweep: 200*nv flip budget (paper thm:lawson-converge precludes this)", nv);
        return flips;
      }

      int h = S.back(); S.pop_back();
      in_stack.clear(edge(h));

      if (!alive(h) || m.delaunay(*this, h)) continue;

      // Record rim edges before flipping (they'll be checked next).
      int h1 = he_next[h], h2 = he_next[h1];
      int h4 = he_next[twin(h)], h5 = he_next[h4];

      if (!flip_edge(h, m, tr)) {
        trip(Status::InvariantViolated,
             "lawson_sweep: flip_edge rejected a non-Delaunay edge (paper lem:ndimpliesconvex)", h);
        return flips;
      }
      flips++; budget--;

      for (int rim : {h1, h2, h4, h5}) {
        int eid = edge(rim);
        if (!in_stack.test(eid)) {
          if (!S.push_back(rim & ~1)) {
            trip(Status::CapacityExceeded, "lawson_sweep: rim push", rim);
            return flips;
          }
          in_stack.set(eid);
        }
      }
    }
    return flips;
  }

  // The seeded restore around a surgery rim (paper cor:local-restore): the
  // frontier is every edge incident to a live rim vertex, pushed in ring
  // order with in_stack dedup -- the same seed set, in the same order, the
  // historical collect-then-sweep produced, without materializing it.
  // @pre  as the seeded lawson_sweep below.
  template <class Metric = BandedFloatMetric, class Transport = NoTransport>
  int lawson_sweep_around(std::span<const int> rim, Spanify::SpanStack<int>& S,
                          Spanify::BitSpan& in_stack,
                          Metric&& m = Metric{}, Transport&& tr = Transport{}) {
    if (status != Status::Ok) return 0;
    S.clear();
    for (int w : rim) {
      if (w < 0 || w >= nv || v_out[w] < 0) continue;
      for (int h : incident(w)) {
        if (!alive(h)) continue;
        int eid = edge(h);
        if (!in_stack.test_and_set(eid)) {
          if (!S.push_back(eid << 1)) {
            trip(Status::CapacityExceeded, "lawson_sweep_around: seed push", eid);
            return 0;
          }
        }
      }
    }
    return lawson_sweep_drain(S, in_stack, m, tr);
  }

  // Seeded Lawson: restore Delaunay from a frontier of edges, propagating to
  // each flip's rim.  in_stack is the caller-owned stack-membership scratch
  // (sized nh/2 bits); @pre all-false on entry, @post all-false on normal
  // (non-status-tripped) return.
  // @pre  seeded: every non-Delaunay edge is in seed_edges or arises on the
  //       rim of a seed cascade (paper lem:seeded-sweep).
  // @post is_delaunay() on Ok.
  template <class Metric = BandedFloatMetric, class Transport = NoTransport>
  int lawson_sweep(std::span<const int> seed_edges, Spanify::SpanStack<int>& S,
                   Spanify::BitSpan& in_stack,
                   Metric&& m = Metric{}, Transport&& tr = Transport{}) {
    if (status != Status::Ok) return 0;
    S.clear();
    for (int h : seed_edges)
      if (alive(h)) {
        int eid = edge(h);
        if (!in_stack.test_and_set(eid)) {
          if (!S.push_back(eid << 1)) {
            trip(Status::CapacityExceeded, "lawson_sweep: seed push", eid);
            return 0;
          }
        }
      }
    return lawson_sweep_drain(S, in_stack, m, tr);
  }

  // Global Delaunay restore: seed with every live edge, in edge order.
  // @post is_delaunay() on Ok (paper thm:lawson-converge).
  template <class Metric = BandedFloatMetric, class Transport = NoTransport>
  int flip_to_delaunay(Spanify::SpanStack<int>& S, Spanify::BitSpan& in_stack,
                       Metric&& m = Metric{}, Transport&& tr = Transport{}) {
    if (status != Status::Ok) return 0;
    S.clear();
    for (int h : edges()) {
      int eid = edge(h);
      if (!in_stack.test_and_set(eid)) {
        if (!S.push_back(eid << 1)) {
          trip(Status::CapacityExceeded, "flip_to_delaunay: seed push", eid);
          return 0;
        }
      }
    }
    return lawson_sweep_drain(S, in_stack, m, tr);
  }

  // -------------------------------------------------------------------------
  // Vertex removal machinery.
  // -------------------------------------------------------------------------

  // Extract the fan polygon around v into ws.fan.  The ONE chokepoint for
  // fan-degree overflow: every downstream fan buffer is sized to the same
  // k_max, so once this check passes the rest is safe by construction.
  void extract_fan(int v, DelaunayWorkspace& ws) {
    if (status != Status::Ok) return;
    FanPolygon& fan = ws.fan;

    fan.k = vertex_degree(v);
    if (fan.k > (int)fan.nb.size()) {
      trip(Status::CapacityExceeded, "extract_fan: star degree exceeds k_max", fan.k);
      return;
    }

    fan.spoke_he[0] = v_out[v];
    for (int i = 1; i < fan.k; i++)
      fan.spoke_he[i] = ccw(fan.spoke_he[i-1]);

    for (int i = 0; i < fan.k; i++) {
      fan.nb[i]        = dest(fan.spoke_he[i]);
      fan.spokes[i]    = he_length[fan.spoke_he[i]];
      fan.inner_rim[i] = he_next[fan.spoke_he[i]];
      fan.rims[i]      = he_length[fan.inner_rim[i]];
    }

    fan.cum[0] = 0;
    for (int i = 0; i < fan.k; i++)
      fan.cum[i+1] = fan.cum[i]
                   + delaunay_detail::triangle_angle(fan.spokes[i],
                                                     fan.spokes[(i+1) % fan.k],
                                                     fan.rims[i]);
  }

  // Exact isometric development of the extracted fan onto Z[w] (exact
  // regime): apex v at the origin, ws.fan.P[i] the lattice position of
  // nb[i].  Face i is placed over base origin -> P[i] by the lattice apex
  // construction (place_third_eis_total, chirality +1 = the face's CCW
  // orientation) -- the same third-vertex placement the paint/atlas layer
  // uses.
  //
  // The anchor P[0] is found by scan: distinct sector-0 representatives
  // of the same norm need not be lattice-equivalent (split primes), and
  // an anchor in an orbit incompatible with the star's chart MAY fail
  // with some apex off the lattice -- so the scan tries every
  // Sector0Reps representative in the deterministic b-order and takes
  // the first that develops fully.  Several orbits can complete (e.g.
  // the all-norm-7 star develops from both conjugate orbits, paper
  // lem:orbit-anchor); ANY complete development is a valid isometric
  // development: chirality-fixed placement makes any two completing
  // developments pure ROTATIONS of each other, and every downstream
  // predicate (wedge, dot2, norm) is rotation-invariant.  Because the apex is
  // EXACTLY flat (eps = 0), a complete development must CLOSE -- replaying
  // face k-1 must reproduce P[0] -- so a complete unclosed development, or
  // no completing representative at all, means the Lsq carry is corrupt ->
  // InvariantViolated.
  // @pre extract_fan(v, ws) done; Lsq in the inductive envelope (so the
  // int narrowing inside is exact).  (Defined after ExactIntegerMetric.)
  inline void develop_fan_lattice(int v, DelaunayWorkspace& ws,
                                  const ExactIntegerMetric& m);

  // The tie-break's side order, exact (paper thm:tie-break): the theorem's
  // side is the one subtending at most pi at the flat vertex.  The two
  // sides partition the cone angle 2*pi exactly, so theta_0 <= theta_1 iff
  // theta_0 <= pi -- decided by developing side 0's sector onto the
  // lattice (the same anchored development as the fan, without closure)
  // and testing the boundary wedge.  Returns 0 or 1; trips on a corrupt
  // carry.  (Defined after ExactIntegerMetric.)
  inline int lattice_first_tie_side(int h_loop, DelaunayWorkspace& ws,
                                    const ExactIntegerMetric& m);

  // Ear-clip the fan polygon into triangles.  By Meisters' theorem a simple
  // polygon with k >= 4 always has an ear; a pass with no ear is a theorem
  // counterexample -> InvariantViolated, never a silent partial result.
  // The ear acceptance test is the metric's (m.ear): exact wedge signs and
  // lattice norms over the development ws.fan.P (@pre m.prepare_star done),
  // or the banded float test over the polar coordinates.
  template <class Metric = BandedFloatMetric>
  void ear_clip_fan(DelaunayWorkspace& ws, Metric&& mt = Metric{}) {
    if (status != Status::Ok) return;
    const FanPolygon& fan = ws.fan;
    FanTriangulation& tri = ws.tri;
    int k = fan.k;

    tri.n_diagonals = 0;
    tri.n_triangles = 0;

    int poly_n = k;
    for (int i = 0; i < k; i++) ws.poly[i] = i;

    while (poly_n > 3) {
      int n = poly_n;
      bool clipped = false;

      for (int j = 0; j < n; j++) {
        int pp = ws.poly[(j - 1 + n) % n], pi = ws.poly[j], pn = ws.poly[(j + 1) % n];
        Length len = mt.ear(*this, fan, pp, pi, pn);
        if (status != Status::Ok) return;
        if (len.len <= 0) continue;
        tri.diagonals[tri.n_diagonals++] = {pp, pi, pn, len};
        for (int m = j; m < poly_n - 1; m++) ws.poly[m] = ws.poly[m + 1];
        poly_n--;
        clipped = true;
        break;
      }

      if (!clipped) {
        trip(Status::InvariantViolated,
             "ear_clip_fan: no acceptable ear (paper thm:two-ears violated)", poly_n);
        return;
      }
    }

    // Compose diagonals into the triangle list, appending the base triangle.
    int rpoly_n = k;
    for (int i = 0; i < k; i++) ws.rpoly[i] = i;
    for (int di = 0; di < tri.n_diagonals; di++) {
      const FanTriangulation::Diagonal& ear = tri.diagonals[di];
      tri.triangles[tri.n_triangles++] = {ear.from, ear.ear, ear.to};
      int fpos = -1;
      for (int m = 0; m < rpoly_n; m++) if (ws.rpoly[m] == ear.ear) { fpos = m; break; }
      if (fpos < 0) {
        trip(Status::InvariantViolated,
             "ear_clip_fan: clipped ear absent from residual polygon", ear.ear);
        return;
      }
      for (int m = fpos; m < rpoly_n - 1; m++) ws.rpoly[m] = ws.rpoly[m + 1];
      rpoly_n--;
    }
    if (rpoly_n != 3) {
      trip(Status::InvariantViolated, "ear_clip_fan: residual polygon size != 3", rpoly_n);
      return;
    }
    tri.triangles[tri.n_triangles++] = {ws.rpoly[0], ws.rpoly[1], ws.rpoly[2]};
  }

  // Replace v's star with the ear triangulation.  Arcs are keyed by polygon
  // POSITIONS (so repeated vertex ids in a delta-complex fan stay
  // unambiguous), resolved by formula -- arc (i, i+1) is inner_rim[i] -- or
  // by scanning the <= k-3 recorded diagonals (O(k) per triangle, no k x k
  // matrix and no per-call reset).  The new faces are pushed onto
  // ws.new_faces (wire_triangle pins f_he to its first argument, so face
  // slot i of a new face is fan index t.vi -- the anchoring the point
  // transport uses).
  template <class Metric = BandedFloatMetric>
  void splice_fan(int v, DelaunayWorkspace& ws, Metric&& m = Metric{}) {
    if (status != Status::Ok) return;
    const FanPolygon& fan = ws.fan;
    const FanTriangulation& tri = ws.tri;
    int k = fan.k;
    // diag_he[di] = the allocated half-edge of diagonal di, from -> to.
    // tri.diagonals has capacity k_max >= k, so reuse ws.poly as the
    // parallel handle array (ear_clip_fan is done with it).
    std::span<int> diag_he = ws.poly;
    auto arc = [&](int i, int j) -> int {
      if (j == (i + 1) % k) return fan.inner_rim[i];
      for (int di = 0; di < tri.n_diagonals; di++) {
        const FanTriangulation::Diagonal& d = tri.diagonals[di];
        if (d.from == i && d.to == j) return diag_he[di];
        if (d.to == i && d.from == j) return twin(diag_he[di]);
      }
      return -1;   // caught by wire_triangle's dead-arc guard
    };

    for (int i = 0; i < k; i++) {
      dealloc_face(he_face[fan.spoke_he[i]]);
      dealloc_edge(fan.spoke_he[i]);
    }
    v_out[v] = -1;
    if (status != Status::Ok) return;

    for (int di = 0; di < tri.n_diagonals; di++) {
      const FanTriangulation::Diagonal& d = tri.diagonals[di];
      int h_d = alloc_directed_edge(fan.nb[d.from], fan.nb[d.to], d.len.len);
      if (status != Status::Ok) return;
      m.set_edge_length(*this, h_d, d.len);   // the metric's paired write
      diag_he[di] = h_d;
    }

    for (int ti = 0; ti < tri.n_triangles; ti++) {
      const FanTriangulation::Triangle& t = tri.triangles[ti];
      int f = wire_triangle(arc(t.v0, t.v1), arc(t.v1, t.v2), arc(t.v2, t.v0));
      if (status != Status::Ok) return;
      if (!ws.new_faces.push_back(f)) {
        trip(Status::CapacityExceeded, "splice_fan: new-face list", f);
        return;
      }
    }

    // Restore each rim vertex's outgoing pointer from a known-live local
    // edge: inner_rim[i] originates at nb[i], is never deallocated (only
    // spokes are), and was just re-wired into an ear triangle.
    for (int i = 0; i < k; i++) v_out[fan.nb[i]] = fan.inner_rim[i];
  }

  // Remove one flat vertex (deg < 3: not removable here -- the driver's
  // restructure rounds handle it).  Degree 3 is the k = 3 instance of the
  // general path: ear_clip_fan's scan loop is empty and the base triangle
  // is the merged face, so splice_fan performs the direct three-face merge
  // (no diagonal is allocated).  Transport plan hooks run before the
  // surgery.  Under the exact metric the star is developed onto the
  // lattice first and the ear machinery runs on wedge signs and norms.
  template <class Metric = BandedFloatMetric, class Transport = NoTransport>
  void remove_flat_vertex(int v, DelaunayWorkspace& ws,
                          Metric&& m = Metric{}, Transport&& tr = Transport{}) {
    if (status != Status::Ok) return;
    if (vertex_degree(v) < 3) return;

    constexpr bool tracking = std::remove_reference_t<Transport>::tracking();
    extract_fan(v, ws);
    if (status != Status::Ok) return;
    if constexpr (tracking) tr.plan_star(v, ws.fan);

    m.prepare_star(*this, ws, v);   // exact: lattice development; banded: no-op
    if (status != Status::Ok) return;
    ear_clip_fan(ws, m);
    if (status != Status::Ok) return;
    if constexpr (tracking) tr.plan_charts(ws.fan, ws.tri);

    ws.new_faces.clear();
    splice_fan(v, ws, m);
    if (status != Status::Ok) return;

    if constexpr (tracking) tr.commit_removal(ws.new_faces.live(), v);
  }

  // The t = 1 cocircular tie-break: an exactly tied self-loop diamond (end
  // angle = pi) is resolved by flipping a side-fan spoke whose diamond is
  // exactly cocircular (energy-neutral, Delaunay-preserving), keeping the
  // flip only if it convexifies the loop (cocircular flips are involutive).
  // Smaller-theta side first (the theorem's side, paper thm:tie-break) --
  // the side ORDER, like every other decision, is the metric's: the exact
  // regime decides it by the lattice sector test (lattice_first_tie_side),
  // so the exact pipeline's output is a function of integer arithmetic
  // alone; the banded regime compares the accumulated float angles.
  template <class Metric = BandedFloatMetric, class Transport = NoTransport>
  bool tie_break_self_loop(int v, int h_loop, DelaunayWorkspace& ws,
                           Metric&& m = Metric{}, Transport&& tr = Transport{}) {
    if (status != Status::Ok) return false;
    const int slots[2] = {h_loop, h_loop ^ 1};
    double theta[2] = {0.0, 0.0};
    for (int i = 0; i < 2; i++) {
      ws.tie_side[i].clear();
      const int start = cw(slots[i] ^ 1), stop = cw(slots[i]);
      long guard = 0;
      for (int g = start; g != stop; g = cw(g)) {
        theta[i] += he_angle[g];
        if (dest(g) != v) {
          if (!ws.tie_side[i].push_back(g))
            return trip(Status::CapacityExceeded, "tie_break_self_loop: side-fan capacity",
                        (int)ws.tie_side[i].size());
        }
        if (++guard > nh)
          return trip(Status::InvariantViolated, "tie_break_self_loop: ring walk overran", h_loop);
      }
    }
    const int first = m.first_tie_side(*this, ws, h_loop, theta);
    if (status != Status::Ok) return false;
    for (int pass = 0; pass < 2; pass++) {
      const int i = (pass == 0) ? first : 1 - first;
      const std::span<const int> side = ws.tie_side[i].live();
      for (int s = 0; s < (int)side.size(); s++) {
        const int g = side[s];
        if (!m.cocircular(*this, g)) continue;
        if (!flip_edge(g, m, tr)) continue; // inscribed quad: convex, but stay guarded
        if (m.convex(*this, h_loop)) return true;
        flip_edge(g, m, tr);                // undo: not the theorem's spoke
      }
    }
    return false;
  }

  // Flip away all self-loops at flat vertex v (they must be cleared before
  // removal, else splice_fan would wire a live edge to the dead vertex).
  // By the half-angle lemma every self-loop at a flat vertex in a Delaunay
  // state has diamond ends <= pi; strict convexity fails only on an EXACT
  // tie, which the tie-break resolves.  The pass budget converts unforeseen
  // cocircular pathology into a loud status instead of a hang; a surviving
  // self-loop is InvariantViolated (provably dead on exact min-degree-3
  // sphere input, paper thm:main).
  template <class Metric = BandedFloatMetric, class Transport = NoTransport>
  void flip_away_self_loops(int v, DelaunayWorkspace& ws,
                            Metric&& m = Metric{}, Transport&& tr = Transport{}) {
    if (status != Status::Ok) return;
    if (v_out[v] < 0) return;
    int budget = 8 * (vertex_degree(v) + 4);
    bool flipped_any = true;
    while (flipped_any) {
      if (--budget < 0) {
        trip(Status::BudgetExceeded,
             "flip_away_self_loops: pass budget (cocircular pathology guard)", v);
        return;
      }
      flipped_any = false;
      for (int h : incident(v)) {
        if (dest(h) != v) continue;
        if (flip_edge(h, m, tr) || flip_edge(h ^ 1, m, tr) ||
            tie_break_self_loop(v, h, ws, m, tr)) {
          flipped_any = true;
          break;                  // the ring changed: restart the scan
        }
        if (status != Status::Ok) return;
      }
    }
    for (int h : incident(v))
      if (dest(h) == v) {
        trip(Status::InvariantViolated,
             "flip_away_self_loops: self-loop survives the tie-break (paper thm:main scope)", v);
        return;
      }
  }

  // Remove every flat vertex, leaving the cones: the work-list driver with
  // per-removal seeded Delaunay restore, restructure rounds, and the
  // stuck-reduction fail-loud guard.  Flatness is the metric's: the exact
  // regime removes eps_v == 0 vertices (integer, no tolerance), the banded
  // regime |kappa_v| < flat_tol (see is_flat for choosing the band).
  // Returns the status (Ok iff fully reduced).
  //
  // @ref  alg:idt (the driver), cor:local-restore (per-removal restore),
  //       thm:full-scan-term (termination), cor:stuck-diagnosis (loud stop).
  // @inv  delaunay: at every work-list pop the complex is globally Delaunay
  //       -- established by the initial flip_to_delaunay, re-established by
  //       each removal's seeded restore; the hypothesis of
  //       cor:local-restore.
  // @variant n_live_flat strictly decreases on each successful removal, and
  //       a restructure round that removes nothing ends the fixed point.
  // @post on Ok: no live v is flat under the metric, and is_delaunay().
  //
  // Flatness is isometric-invariant, so a removal only changes the
  // REMOVABILITY of star neighbours -- each vertex is touched O(deg) times
  // per round, O(N) total.  The two effects the local re-push provably
  // cannot reach, hence the outer restructure rounds: (a) a flat made
  // removable by a seeded-sweep flip OUTSIDE the rim; (b) a deg <= 2 flat
  // that only a global restructure lifts to deg >= 3 (paper lem:flat-deg3).
  template <class Metric = BandedFloatMetric, class Transport = NoTransport,
            class Observer = NoRemovalObserver>
  Status remove_flat_vertices(DelaunayWorkspace& ws,
                              Metric&& m = Metric{},
                              Transport&& tr = Transport{},
                              Observer&& obs = Observer{}) {
    if (status != Status::Ok) return status;

    // Reused across the per-removal seeded sweeps: each sweep leaves it
    // all-false on normal return, so the zero-fill happens once.
    ws.sweep_in_stack.clear_all();
    ws.work.clear();
    ws.on_queue.clear_all();

    // Admit v to the work list iff it is a live flat not already queued --
    // THE admission predicate, shared by the local re-push and the global
    // re-seed.  Returns whether v was newly admitted.
    auto enqueue = [&](int v) -> bool {
      if (status != Status::Ok) return false;
      if (v < 0 || v >= nv || v_out[v] < 0 || !m.is_flat(*this, v)) return false;
      if (ws.on_queue.test_and_set(v)) return false;
      if (!ws.work.push_back(v))
        return trip(Status::CapacityExceeded, "remove_flat_vertices: work list", v);
      return true;
    };
    auto push = [&](int v) { enqueue(v); };

    // Attempt one removal; on success restore Delaunay locally and re-push
    // the affected rim flats.  Flatness is isometric-invariant, so a removal
    // only changes the REMOVABILITY of star neighbours.
    auto try_remove = [&](int v) -> bool {
      if (status != Status::Ok) return false;
      if (v < 0 || v >= nv || v_out[v] < 0 || !m.is_flat(*this, v)) return false;
      obs.on_pop(v);
      // Collect the star rim on BOTH sides of the self-loop cleanup (its
      // flips perturb the pre-flip ring; the removal re-triangulates the
      // post-flip rim) -- a superset of the edges whose Delaunay status the
      // surgery can change.
      ws.ring.clear();
      for (int h : incident(v)) {
        if (!ws.ring.push_back(dest(h)))
          return trip(Status::CapacityExceeded, "remove_flat_vertices: ring", (int)ws.ring.size());
      }
      flip_away_self_loops(v, ws, m, tr);
      if (status != Status::Ok) return false;
      for (int h : incident(v)) {
        if (!ws.ring.push_back(dest(h)))
          return trip(Status::CapacityExceeded, "remove_flat_vertices: ring", (int)ws.ring.size());
      }

      remove_flat_vertex(v, ws, m, tr);
      if (status != Status::Ok) return false;
      const bool removed_v = (v_out[v] < 0);
      if (removed_v)
        while (nv > 0 && v_out[nv - 1] < 0) nv--;

      // Restore Delaunay from the rim even when the removal failed
      // (deg <= 2): the cleanup may have flipped, and every pop must see a
      // globally-Delaunay mesh for the seeded restore to be exact
      // (paper cor:local-restore).
      lawson_sweep_around(ws.ring.live(), ws.lawson_stack, ws.sweep_in_stack, m, tr);
      if (status != Status::Ok) return false;
      if (!removed_v) return false;

      for (int w : ws.ring.live()) push(w);
      obs.on_removed(v);
      return true;
    };

    auto seed_all_flats = [&]() -> bool {
      bool any = false;
      for (int v = 0; v < nv; v++) {
        if (status != Status::Ok) return any;
        any |= enqueue(v);
      }
      return any;
    };

    auto drain = [&]() -> bool {
      bool any = false;
      while (!ws.work.empty()) {
        int v = ws.work.back(); ws.work.pop_back(); ws.on_queue.clear(v);
        if (try_remove(v)) any = true;
        if (status != Status::Ok) return any;
      }
      return any;
    };

    // Establish the iDT up front so every per-removal seeded restore starts
    // from a globally-Delaunay mesh.
    flip_to_delaunay(ws.lawson_stack, ws.sweep_in_stack, m, tr);
    if (status != Status::Ok) return status;

    seed_all_flats();
    if (status != Status::Ok) return status;
    drain();
    if (status != Status::Ok) return status;

    // Restructure rounds: a global sweep can unblock deg <= 2 stragglers
    // and expose flats the local re-push missed.
    while (true) {
      flip_to_delaunay(ws.lawson_stack, ws.sweep_in_stack, m, tr);
      if (status != Status::Ok) return status;
      if (!seed_all_flats()) break;   // no flats remain -> done
      if (status != Status::Ok) return status;
      if (!drain()) break;            // flats remain but none removable -> done
      if (status != Status::Ok) return status;
    }

    // Fail loud on a stuck reduction: returning a partial reduction silently
    // would hand callers a non-cone iDT.
    for (int v = 0; v < nv; v++)
      if (v_out[v] >= 0 && m.is_flat(*this, v)) {
        trip(Status::InvariantViolated,
             "remove_flat_vertices: stuck with live flat vertex (paper cor:stuck-diagnosis)", v);
        return status;
      }

    // Final Lawson: the output is globally Delaunay regardless of removal
    // order (paper thm:lawson-converge; a B == D self-loop-creating flip's
    // new loop is strictly Delaunay by paper lem:selfloop-delaunay).
    flip_to_delaunay(ws.lawson_stack, ws.sweep_in_stack, m, tr);
    return status;
  }

  // Structural well-formedness over caller-owned visited bits (>= nh):
  // every live half-edge in exactly one he_next cycle, every cycle length 3.
  bool is_well_formed(Spanify::BitSpan& visited) const {
    visited.clear_all();
    for (int h0 = 0; h0 < nh; h0++) {
      if (!alive(h0) || visited.test(h0)) continue;
      int h = h0;
      for (int step = 0; step < 3; step++) {
        if (!alive(h) || visited.test(h)) return false;
        visited.set(h);
        h = he_next[h];
      }
      if (h != h0) return false;
    }
    for (int h = 0; h < nh; h++)
      if (alive(h) && !visited.test(h)) return false;
    return true;
  }

  // -------------------------------------------------------------------------
  // Structural + metric validation: the DCEL class invariant, executable.
  // True iff all nine facts hold on every live element:
  //   @inv 1  twin closure: twin(h) < nh for live h
  //   @inv 2  triangular faces: he_next^3(h) == h
  //   @inv 3  origin chaining: dest(h) == origin(next(h))
  //   @inv 4  face agreement: face(next(h)) == face(h)
  //   @inv 5  v_out witnesses: v_out[v] live and originating at v
  //   @inv 6  f_he witnesses: f_he[f] live and lying in face f
  //   @inv 7  positive lengths on live half-edges
  //   @inv 8  twin length equality: length(h) == length(twin h)
  //   @inv 9  triangle inequality on every live face (1e-9 relative slack)
  // Allocation-free; O(nh).
  // -------------------------------------------------------------------------
  bool check_consistency() const {
    for (int h = 0; h < nh; h++)                                    // @inv 1
      if (alive(h) && twin(h) >= nh) return false;

    for (int h = 0; h < nh; h++)                                    // @inv 2
      if (alive(h) && he_next[he_next[he_next[h]]] != h) return false;

    for (int h = 0; h < nh; h++)                                    // @inv 3
      if (alive(h) && dest(h) != he_origin[he_next[h]]) return false;

    for (int h = 0; h < nh; h++)                                    // @inv 4
      if (alive(h) && he_face[he_next[h]] != he_face[h]) return false;

    for (int v = 0; v < nv; v++)                                    // @inv 5
      if (v_out[v] >= 0 && (!alive(v_out[v]) || he_origin[v_out[v]] != v))
        return false;

    for (int f = 0; f < nf; f++)                                    // @inv 6
      if (f_he[f] >= 0 && (!alive(f_he[f]) || he_face[f_he[f]] != f))
        return false;

    for (int h = 0; h < nh; h++)                                    // @inv 7
      if (alive(h) && he_length[h] <= 0) return false;

    for (int h = 0; h < nh; h += 2)                                 // @inv 8
      if (alive(h) && he_length[h] != he_length[twin(h)]) return false;

    for (int h = 0; h < nh; h++) {                                  // @inv 9
      if (!alive(h)) continue;
      int h1 = he_next[h], h2 = he_next[h1];
      double L0 = he_length[h], L1 = he_length[h1], L2 = he_length[h2];
      double eps = 1e-9 * (L0 + L1 + L2);
      if (L0 > L1 + L2 + eps) return false;
      if (L1 > L0 + L2 + eps) return false;
      if (L2 > L0 + L1 + eps) return false;
    }

    return true;
  }

  // -------------------------------------------------------------------------
  // Canonical field order.
  //
  // to_tuple() lists the 9 SoA arrays in THE canonical order.  Deliberately
  // OUTSIDE the tuple: the counts nv/nh/nf, the capacities, the Status
  // latch, and the free-list stacks -- allocation bookkeeping, not surface
  // state (state comparators handle them separately and by name):
  //   { he_next, he_origin, he_face, he_length, he_angle,
  //     v_out, v_cone_angle, v_orig_degree, f_he }
  // This order is the contract the owner's repoint() fold and any field-wise
  // state comparator build on.  Allocation SIZING is dcel_capacities' job --
  // this view deliberately does not model the parent's batch::batchable_view
  // concept (its per-vertex multipliers cannot express the DCEL's affine
  // sizes); a batch adapter derives its tables from dcel_capacities.
  // -------------------------------------------------------------------------
  static constexpr std::size_t n_fields = 9;

  auto to_tuple() {
    return std::forward_as_tuple(he_next, he_origin, he_face,
                                 he_length, he_angle,
                                 v_out, v_cone_angle, v_orig_degree, f_he);
  }
  auto to_tuple() const {
    return std::forward_as_tuple(he_next, he_origin, he_face,
                                 he_length, he_angle,
                                 v_out, v_cone_angle, v_orig_degree, f_he);
  }
};

static_assert(std::is_trivially_copyable_v<DelaunayView>,
              "DelaunayView must be trivially copyable");
static_assert(std::tuple_size_v<decltype(std::declval<DelaunayView>().to_tuple())>
                  == DelaunayView::n_fields,
              "to_tuple arity must equal n_fields");

// ============================================================================
// The two metric policies (declared before DelaunayView; banner there).
// Each is the complete vocabulary of one regime -- flatness, the three edge
// predicates, the flipped length, the ear acceptance, the star preparation,
// the tie-break side order, and the paired edge-length write.
// ============================================================================

// General metrics: the float Diamond predicates with the named tolerance
// bands (delaunay_detail), absorbing FP noise.  Carries the flatness band.
struct BandedFloatMetric {
  double flat_tol = delaunay_detail::default_flat_tol;

  bool is_flat(const DelaunayView& V, int v) const { return V.is_flat(v, flat_tol); }
  bool delaunay(const DelaunayView& V, int h) const { return V.diamond(h).is_delaunay(); }
  bool convex(const DelaunayView& V, int h) const { return V.diamond(h).is_convex(); }
  bool cocircular(const DelaunayView& V, int h) const {
    return V.diamond(h).is_cocircular(delaunay_detail::tie_cocircular_tol);
  }
  std::optional<Length> flipped(DelaunayView& V, int h) const {
    double f = V.diamond(h).flipped_length();
    if (!std::isfinite(f) || f <= 0) return std::nullopt;
    return Length{f, 0};
  }
  void set_edge_length(DelaunayView& V, int h, Length l) const {
    V.he_length[h] = V.he_length[V.twin(h)] = l.len;
  }
  void prepare_star(DelaunayView&, DelaunayWorkspace&, int /*v*/) const {}
  Length ear(DelaunayView&, const FanPolygon& fan, int pp, int pi, int pn) const {
    return {delaunay_detail::ear_length_if_acceptable(fan, pp, pi, pn), 0};
  }
  int first_tie_side(DelaunayView&, DelaunayWorkspace&, int /*h_loop*/,
                     const double theta[2]) const {
    return (theta[0] <= theta[1]) ? 0 : 1;
  }
};

// Equilateral metrics and everything flips and flat-vertex removals produce
// from them: every face stays an Eisenstein-lattice triangle, so every
// decision is an integer sign (DiamondSq's tau form), flip and ear
// diagonals are lattice norms carried in Lsq, flatness is the integer cone
// excess eps_v == 0, and the tie-break side order is the lattice sector
// test -- the whole reduction is a function of integer arithmetic alone.
// Corruption of the carry trips InvariantViolated through the view's latch;
// nothing here throws.
// @inv envelope: every live Lsq entry is in [1, exact_lsq_max] -- holds at
//       entry (the owner driver derives and verifies the carry) and is
//       maintained inductively by the two guarded write sites below.
struct ExactIntegerMetric {
  std::span<long long> Lsq;   // squared length per half-edge (the carry)

  bool is_flat(const DelaunayView& V, int v) const { return V.cone_excess(v) == 0; }
  bool delaunay(const DelaunayView& V, int h) const {
    return V.diamond_sq(h, Lsq).is_delaunay();
  }
  bool convex(const DelaunayView& V, int h) const {
    return V.diamond_sq(h, Lsq).is_convex();
  }
  bool cocircular(const DelaunayView& V, int h) const {
    return V.diamond_sq(h, Lsq).is_cocircular();
  }
  std::optional<Length> flipped(DelaunayView& V, int h) const {
    auto fsq = V.diamond_sq(h, Lsq).flipped_length_sq();
    if (!fsq) {
      V.trip(DelaunayView::Status::InvariantViolated,
             "flip_edge: convex diamond without a lattice diagonal (Lsq corrupt or out of envelope)", h);
      return std::nullopt;
    }
    if (*fsq <= 0) return std::nullopt;       // geometrically degenerate diagonal
    if (*fsq > exact_lsq_max) {
      V.trip(DelaunayView::Status::InvariantViolated,
             "flip_edge: flipped Lsq exceeds the exact envelope", h);
      return std::nullopt;
    }
    return Length{std::sqrt((double)*fsq), *fsq};
  }
  void set_edge_length(DelaunayView& V, int h, Length l) const {
    const int t = V.twin(h);
    V.he_length[h] = V.he_length[t] = l.len;
    Lsq[h] = Lsq[t] = l.sq;
  }
  void prepare_star(DelaunayView& V, DelaunayWorkspace& ws, int v) const {
    V.develop_fan_lattice(v, ws, *this);
  }
  Length ear(DelaunayView& V, const FanPolygon& fan, int pp, int pi, int pn) const {
    const long long dsq = delaunay_detail::ear_diag_sq_if_acceptable(fan.P, pp, pi, pn);
    if (dsq == 0) return {0, 0};
    if (dsq > exact_lsq_max) {
      V.trip(DelaunayView::Status::InvariantViolated,
             "ear_clip_fan: ear Lsq exceeds the exact envelope", pi);
      return {0, 0};
    }
    return {std::sqrt((double)dsq), dsq};
  }
  int first_tie_side(DelaunayView& V, DelaunayWorkspace& ws, int h_loop,
                     const double* /*theta*/) const {
    return V.lattice_first_tie_side(h_loop, ws, *this);
  }
};

// ---------------------------------------------------------------------------
// Exact lattice developments (bodies here: they read ExactIntegerMetric).
// ---------------------------------------------------------------------------

inline void DelaunayView::develop_fan_lattice(int v, DelaunayWorkspace& ws,
                                              const ExactIntegerMetric& m) {
  if (status != Status::Ok) return;
  const FanPolygon& fan = ws.fan;
  const int k = fan.k;
  auto Ls = [&](int i) { return (int)m.Lsq[fan.spoke_he[i]]; };
  auto Lr = [&](int i) { return (int)m.Lsq[fan.inner_rim[i]]; };

  auto attempt = [&](Eisenstein P0) -> Development {
    fan.P[0] = P0;
    for (int i = 0; i < k; i++) {
      auto C = place_third_eis_total({0, 0}, fan.P[i],
                                     Ls(i), Ls((i + 1) % k), Lr(i), +1);
      if (!C) return Development::OffLattice;
      if (i + 1 < k) fan.P[i + 1] = *C;
      else if (!(*C == fan.P[0])) return Development::Unclosed;
    }
    return Development::Closed;
  };

  for (Eisenstein z : Sector0Reps(Ls(0))) {
    switch (attempt(z)) {
      case Development::Closed:     return;
      case Development::OffLattice: continue;   // orbit fails here: next representative
      case Development::Unclosed:
        trip(Status::InvariantViolated,
             "develop_fan_lattice: complete development fails to close at a flat apex", v);
        return;
    }
  }
  trip(Status::InvariantViolated,
       "develop_fan_lattice: no closing lattice development for the star", v);
}

inline int DelaunayView::lattice_first_tie_side(int h_loop, DelaunayWorkspace& ws,
                                                const ExactIntegerMetric& m) {
  // Side 0's sector arcs, directly in CCW order: the outgoing arcs from
  // h_loop round to twin(h_loop), inclusive (the CW walk the tie-break's
  // theta accumulation runs is exactly this walk reversed).  The face
  // between CCW-consecutive arcs A_t, A_{t+1} is face(A_t), with spoke Lsq
  // at the arcs and rim Lsq at he_next(A_t) -- the same per-face shape
  // develop_fan_lattice places, and the same anchor-independence argument:
  // any two completing developments are pure rotations of each other, and
  // wedge / dot2 are rotation-invariant, so the sector verdict does not
  // depend on the representative found.
  int n = 0;
  for (int g = h_loop; ; g = ccw(g)) {
    if (n >= (int)ws.poly.size()) {
      trip(Status::CapacityExceeded, "lattice_first_tie_side: sector-arc scratch", n);
      return 0;
    }
    ws.poly[n++] = g;
    if (g == twin(h_loop)) break;
  }
  auto Ls = [&](int t) { return (int)m.Lsq[ws.poly[t]]; };
  auto Lr = [&](int t) { return (int)m.Lsq[he_next[ws.poly[t]]]; };

  for (Eisenstein z : Sector0Reps(Ls(0))) {
    Eisenstein P0 = z, Pt = z;
    bool on_lattice = true;
    for (int t = 0; t + 1 < n && on_lattice; t++) {
      auto C = place_third_eis_total({0, 0}, Pt, Ls(t), Ls(t + 1), Lr(t), +1);
      if (!C) on_lattice = false;
      else    Pt = *C;
    }
    if (on_lattice)
      return delaunay_detail::lattice_sector_at_most_pi(P0, Pt) ? 0 : 1;
  }
  trip(Status::InvariantViolated,
       "lattice_first_tie_side: no lattice development for the tie sector", h_loop);
  return 0;
}
