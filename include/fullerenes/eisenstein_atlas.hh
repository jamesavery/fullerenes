#pragma once

// =====================================================================
// Tier 3.5 -- the exact cross-cell atlas over the eisenstein_paint cells.
//
// The paint tables chart each live iDT 2-cell into its OWN Eisenstein
// frame and enumerate the lattice points it claims.  This module glues
// those charts into an atlas so that exact integer/rational geometry
// can move BETWEEN cells:
//
//   * per-HALF-EDGE chart transitions: he_trans[h] is the unique
//     orientation-preserving lattice isometry (unit rotation + Z[w]
//     translation -- the crystallographic restriction makes it exact)
//     from face(h)'s frame to face(twin h)'s frame, matched on the
//     crossed edge's two corner SLOTS in each face cycle (position in
//     the cycle, not cone identity, so repeated-corner cells cannot
//     confuse the match).  Keyed by the crossing half-edge, a
//     delta-complex cell pair sharing several edges is unambiguous by
//     construction (a face-pair key could not distinguish them) --
//     a defensive property: no fullerene-dual witness is known
//     (checked C20-C132, raw and realized iDTs), and on all of those
//     the face-pair collapse computes identical transitions;
//   * anchored T_sorted edges: two ADJACENT lattice points claimed by
//     ONE cell.  Both lie in the closed CONVEX cell, so the unit
//     segment between them is inside it, and the cell development --
//     an isometry -- makes them the true adjacent copies.  Anchor
//     selection iterates entries in SCAN order, so every anchor -- and
//     hence every resolved chain -- is deterministic and
//     platform-stable;
//   * a multi-source BFS over the T_sorted edge graph (edges adjacent
//     iff they share a face) that routes every face to an anchored
//     edge by midpoint hops.  Each hop segment is a triangle midline:
//     it stays inside ONE flat, cone-free face, so its exact trace
//     through the cell complex is always well-defined;
//   * exact straight-segment tracing (trace_segment) and per-face
//     chart resolution (resolve_face lazily on the host; resolve_all
//     precomputes every face, the read-only form a device consumer
//     needs), on top of which any barycentric sample of any T_sorted
//     face resolves to an exact (host cell, rational position) pair
//     (locate_sample).
//
// Every atlas table is FLAT: dense arc-indexed arrays (over
// arc_space(T), keyed by u*dmax + arc_ix(u,v) -- the vocabulary the
// DCEL build uses), plain vectors of plain value structs; CellAtlas
// owns them and AtlasView is the span mirror (itself trivially
// copyable; the element structs carry Eisenstein and are not -- an
// int32 flattening is the device port's job, tracked there).  Cell
// charts are read straight from the ParamTablesView
// (claim/corner_ids/frame_points) -- the atlas holds NO per-cell
// index of its own.
//
// All positions are Eisenstein integers or EisensteinRational; every
// predicate is exact integer arithmetic (__int128 for the
// cross-multiplied comparisons, internal to the .cc).  This
// implementation runs host-only.  A device port must satisfy itself
// that its target supports __int128: the arithmetic has been measured
// conformant on one NVIDIA GPU under one AdaptiveCpp SSCP flow
// (2026-07-29), which says nothing about other vendors or about
// DPC++/SPIR-V -- re-measure per target rather than assuming.  No
// floating point enters any decision.  Width
// contract: Eisenstein components are int; denominators and sample
// weights must fit int after common-denominator scaling (checked
// narrows throw, never truncate); the subsequent int products are
// covered by the chart envelope -- cell coordinates are bounded by
// the cell development (|components| <~ 26000, eisenstein.hh) times
// a checked-int factor, far inside int width.
//
// Failure contract: every failure here is a DEEP invariant violation
// (malformed pipeline state, a chart that does not verify), thrown
// uniformly as std::logic_error.  The pipeline's own parametrize gate
// (PaintError) is the modeled-failure boundary; past it, the atlas
// either works or the input violated its contract.
//
// LIFETIME/STABILITY: the view members alias the parametrization's
// and owner's buffers -- valid while those live and do not grow;
// treat them as read-only.
// =====================================================================

#include "fullerenes/eisenstein.hh"        // Eisenstein, EisensteinRational, LatticeIsometry
#include "fullerenes/eisenstein_paint.hh"  // SurfaceParametrization
#include "fullerenes/eisenstein_paint_tables.hh"  // ParamTablesView, LatticePoint
#include "fullerenes/triangulation.hh"
#include "fullerenes/delaunay.hh"

#include <array>
#include <span>
#include <vector>

namespace eisenstein_paint {

// A T_sorted face developed into one chart: its three corners (CCW, unit
// triangle -- verified exactly on resolution) plus a located interior/edge
// point that seeds sample traces.
struct ResolvedFace {
  bool ok = false;
  int  cell = -1;                    // host cell of `anchor`
  std::array<Eisenstein, 3> P{};     // face corners (CCW) in `cell`'s frame
  EisensteinRational anchor;         // located midpoint of one face edge (inside `cell`)
};

// An anchored T_sorted edge: u, v claimed at adjacent positions by `cell`.
struct AnchorEdge {
  int cell; int u, v; Eisenstein pu, pv;
  bool operator==(const AnchorEdge&) const = default;
};

// The trivially-copyable span mirror of the atlas (the read surface;
// resolve_face's lazy cache stays on the owning CellAtlas).
struct AtlasView {
  ParamTablesView   V;   // the charted tables (claim / corners / frames)
  DelaunayView      D;   // the charted cone iDT
  TriangulationView T;   // T_sorted

  std::span<const LatticeIsometry>   he_trans;        // [D.nh] face(h) -> face(twin h)
  std::span<const tri_t>             tface;           // T_sorted faces, CCW
  std::span<const int32_t>           arc_face;        // [N*dmax] arcid -> face left of arc
  std::span<const int32_t>           arc_edge;        // [N*dmax] arcid -> undirected edge id
  std::span<const std::array<int,2>> tedge;           // undirected edges (u < v)
  std::span<const std::array<int,2>> edge_faces;      // per edge, its two faces
  std::span<const int32_t>           anchor_of_edge;  // edge id -> anchors index, or -1
  std::span<const AnchorEdge>        anchors;
  std::span<const int32_t>           bfs_parent_edge; // edge id -> parent edge (-1 root)
  std::span<const int32_t>           bfs_via_face;    // edge id -> face shared with parent
  std::span<const int32_t>           bfs_depth;       // edge id -> chain length to its anchor
  std::span<const ResolvedFace>      resolved;        // complete after resolve_all

  // Dense arc lookups (the DCEL build's arc vocabulary).  A (u, v)
  // that is not a T arc yields -1 -- misses stay loud at the caller
  // (arc_ix reports them as -1; folding that into the index would
  // read a neighbouring row's slot and look like a valid answer).
  // Named arc_of, not arcid: the graph's own arcid(u, i) takes a
  // neighbour SLOT, and an identically-spelled (u, v) overload would
  // invite silent misuse.
  int arc_of(int u, int v) const {
    const int ix = T.arc_ix(u, v);
    return ix < 0 ? -1 : u * T.dmax + ix;
  }
  int face_of_arc(int u, int v) const {
    const int a = arc_of(u, v);
    return a < 0 ? -1 : arc_face[a];
  }
  int edge_of(int u, int v) const {
    const int a = arc_of(u, v);
    return a < 0 ? -1 : arc_edge[a];
  }
};

static_assert(std::is_trivially_copyable_v<AtlasView>);

// The atlas: owning flat vectors + the lazy resolution cache.  Build
// with build_atlas; read through view() (or the owning members -- same
// data).  P must outlive the atlas and not be modified.
struct CellAtlas {
  ParamTablesView   V;
  DelaunayView      D;
  TriangulationView T;

  std::vector<LatticeIsometry>   he_trans;
  std::vector<tri_t>             tface;
  std::vector<int32_t>           arc_face;
  std::vector<int32_t>           arc_edge;
  std::vector<std::array<int,2>> tedge;
  std::vector<std::array<int,2>> edge_faces;
  std::vector<int32_t>           anchor_of_edge;
  std::vector<AnchorEdge>        anchors;
  std::vector<int32_t>           bfs_parent_edge;
  std::vector<int32_t>           bfs_via_face;
  std::vector<int32_t>           bfs_depth;

  // Per-face resolution cache (resolve_face lazily, resolve_all in
  // full); sized once at build, so the view's span stays valid.
  std::vector<ResolvedFace> resolved;

  AtlasView view() const {   // members assigned BY NAME: a reordered
    AtlasView v;             // or added table cannot silently mis-wire
    v.V = V; v.D = D; v.T = T;
    v.he_trans        = he_trans;
    v.tface           = tface;
    v.arc_face        = arc_face;
    v.arc_edge        = arc_edge;
    v.tedge           = tedge;
    v.edge_faces      = edge_faces;
    v.anchor_of_edge  = anchor_of_edge;
    v.anchors         = anchors;
    v.bfs_parent_edge = bfs_parent_edge;
    v.bfs_via_face    = bfs_via_face;
    v.bfs_depth       = bfs_depth;
    v.resolved        = resolved;
    return v;
  }
};

// Build the atlas over P's charts.
//
// P may be parametrized on the raw dual_idt or on a realized
// polytope's iDT; a raw chart whose development folds is refused
// loudly by the anchoring gate (vertices developed adjacent that are
// not a mesh edge).  The atlas is INTRINSIC-ONLY and
// curvature-sign-agnostic.  Deterministic: the output is a pure
// function of P (no hash-iteration order anywhere; anchor selection
// is in scan order).
//
//   @throws std::logic_error on any deep invariant violation: a live
//           cell without a chart, a chart whose corner order departs
//           from its face cycle (the corner-k = k-th-cycle-origin
//           theorem is gated per cell), a folded development
//           (non-edge lattice adjacency), an uncovered edge graph.
CellAtlas build_atlas(const SurfaceParametrization& P);

// =====================================================================
// The T-free (intrinsic) atlas -- charts and transitions over a pure
// flat cone metric, where no dual triangulation T is known.
// =====================================================================

// Every admissible development of every cell (CellDevelopments: corners,
// frames, lattice scans) plus THE SELECTION: accepted[f] names the
// development that develops consistently with its neighbours, and
// he_trans the transitions that selection makes exact.  Owning, unlike
// CellAtlas: there is no SurfaceParametrization to alias -- these charts
// ARE the parametrization the intrinsic path has.
//
// The lattice points of the accepted charts are numbered globally by
// node_first (the exclusive prefix over their entry counts), so
// node_at(f, p) is exactly the union-find node intrinsic_dual glues on.
// The claimed VERTEX IDS are absent by construction: deriving them is
// intrinsic_dual's whole job.
struct IntrinsicAtlas {
  DelaunayView         D;           // the charted cone iDT
  CellDevelopments     dev;         // every admissible development, per cell
  std::vector<int32_t> accepted;    // [D.nf] chosen development, -1 unplaced
  std::vector<int32_t> node_first;  // [D.nf+1] exclusive prefix of chart entries
  std::vector<LatticeIsometry> he_trans;  // [D.nh] face(h) -> face(twin h)

  bool placed(int f) const {
    return f >= 0 && f < (int)accepted.size() && accepted[f] >= 0;
  }
  // @pre placed(f).
  DevelopmentView chart(int f) const { return dev.development(f, accepted[f]); }
  std::array<Eisenstein, 3> frame_points(int f) const { return chart(f).frame_points(); }
  std::array<int, 3> corner_ids(int f) const {
    const CellCorners c = dev.corners(f);
    return { c.c0, c.c1, c.c2 };
  }
  int n_nodes() const { return node_first.empty() ? 0 : node_first.back(); }

  // The global node id of lattice point p in cell f's chart, or -1 when
  // f does not claim p (dead, unplaced, or p outside the scan).  Pure
  // row arithmetic -- ChartView::claim's index, with no vertex id to
  // read.
  int node_at(int f, Eisenstein p) const {
    if (!placed(f)) return -1;
    const DevelopmentView c = chart(f);
    if (p.second < c.scan.b_min || p.second > c.scan.b_max) return -1;
    const ScanRow r = c.rows[p.second - c.scan.b_min];
    if (p.first < r.a_left || p.first > r.a_right) return -1;
    return node_first[f] + r.entry_off + (p.first - r.a_left);
  }
};

// Build the atlas over D WITHOUT a known dual triangulation T (pure flat cone
// metric): the cell developments are the SAME exact enumeration build_atlas
// charts through (cell_developments), but the choice among them is made
// INTRINSICALLY -- by demanding a unit-rotation transition to an already-placed
// neighbour across their shared iDT edge -- where the parametrization would use
// T's frame walkers.  No claims, no T: only the accepted charts and their
// transitions, which is all the intrinsic dual reconstruction (no T, no 3D)
// consumes.  Curvature-sign-agnostic like build_atlas.
//
// Seed a cell with one of its developments (a global reflection/rotation gauge
// the canonical spiral name is invariant under), then BFS: each new cell takes
// the development whose shared-edge transition to its placed neighbour is a unit
// rotation.  A cell reachable through several edges is placed by the first
// neighbour that resolves it; one with no consistent development (a folded cell
// on a noisy metric) is left unplaced rather than aborting.
//
//   @pre    D a post-flip simplicial iDT with integer edge length-squares.
//   @post   consistently reachable cells have accepted >= 0; he_trans holds the
//           transition across every half-edge between two placed cells.
//   @throws PaintError(EMBED) when a live cell's metric leaves the integrality
//           band or admits no lattice development at all (cell_developments'
//           float-to-exact boundary); std::logic_error if D has no placeable
//           cell, or if two placed cells' shared edge does not develop by a
//           lattice isometry (deep invariants).
IntrinsicAtlas build_intrinsic_atlas(const DelaunayView& D);

// Reconstruct the ORIENTED dual triangulation from an intrinsic atlas: take each
// placed cell's lattice points, unite the copies shared across every edge (via
// the chart transitions), and read off the CCW unit-triangle faces as an
// oriented graph.  The cones become its degree-5 (pentagon) vertices and the
// interior lattice points its degree-6 (hexagon) vertices -- the fullerene dual
// -- ready for spiral_nomenclature.  Curvature-sign-agnostic: vertex degrees are
// read off, not assumed, so fulleroid duals reconstruct unchanged.  This is the
// intrinsic, combinatorial-only (no T, no 3D) analogue of run()'s dual.
//
//   @pre    A from build_intrinsic_atlas with every cell placed (a clean iDT).
//   @post   an oriented Triangulation isomorphic to the dual of A.D's surface.
//   @throws std::logic_error on an orientation inconsistency or a shared lattice
//           point absent from a cell scan (deep invariants).
Triangulation intrinsic_dual(const IntrinsicAtlas& A);

// Trace the straight segment a -> b (both in `cell`'s frame; a inside the
// CLOSED cell) through the cell complex.  Returns b's host cell, b
// re-expressed in that cell's frame, and the carried frame isometry.
// Terminates: the crossing parameter is non-decreasing, and strictly
// increasing whenever the open segment meets no cone -- which every
// internal caller guarantees (anchored segments and face midlines);
// the in-binary crossing cap is a deep-invariant guard, not a bound
// any valid trace approaches.
//   @throws std::logic_error when stuck (no exit edge) or on cap trip.
struct CellTrace { int cell; EisensteinRational pos; LatticeIsometry carried; };
CellTrace trace_segment(const CellAtlas& A, int cell,
                        EisensteinRational a, EisensteinRational b);

// Resolve T_sorted face fi into a chart (lazy, cached): its three corners
// CCW as a unit triangle in some cell's frame, plus a located seed point.
//   @throws std::logic_error on any failed hop or failed verification.
const ResolvedFace& resolve_face(CellAtlas& A, int fi);

// Precompute every face's resolution (the device-facing read-only form:
// after this, `resolved` is a complete table and no call mutates A).
void resolve_all(CellAtlas& A);

// Exact location of the barycentric sample (a*P0 + b*P1 + c*P2) / den of
// T_sorted face fi: the host cell and the exact rational position in its
// frame.  @pre a, b, c >= 0 and a + b + c == den (checked: the sample
// must lie in the CLOSED face -- the anchor-to-sample trace stays in
// flat cone-free territory only there).
struct CellPoint { int cell; EisensteinRational pos; };
CellPoint locate_sample(CellAtlas& A, int fi, long a, long b, long c, long den);

}  // namespace eisenstein_paint
