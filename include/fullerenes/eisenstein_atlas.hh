#pragma once

// =====================================================================
// Tier 3.5 -- the exact cross-cell atlas over the eisenstein_paint cells.
//
// The Tier-3 primitives embed each live iDT 2-cell into its OWN
// Eisenstein-lattice frame (a chart) and enumerate the lattice points it
// claims.  This module glues those charts into an atlas so that exact
// integer/rational geometry can move BETWEEN cells:
//
//   * per-adjacency chart transitions: the unique orientation-preserving
//     lattice isometry (unit rotation + Z[w] translation -- the
//     crystallographic restriction makes it exact) matching the shared
//     iDT edge's two cone corners across the two frames;
//   * anchored T_sorted edges: two ADJACENT lattice points claimed by ONE
//     cell.  Both lie in the closed CONVEX cell, so the unit segment
//     between them is inside it, and the cell development -- an isometry
//     -- makes them the true adjacent copies;
//   * a multi-source BFS over the T_sorted edge graph (edges adjacent iff
//     they share a face) that routes every face to an anchored edge by
//     midpoint hops.  Each hop segment is a triangle midline: it stays
//     inside ONE flat, cone-free face, so its exact trace through the
//     cell complex is always well-defined -- no development ever wraps a
//     cone;
//   * exact straight-segment tracing (trace_segment) and lazy per-face
//     chart resolution (resolve_face), on top of which any barycentric
//     sample of any T_sorted face resolves to an exact
//     (host cell, rational Eisenstein position) pair (locate_sample /
//     locate_vertex).
//
// All positions are Eisenstein integers or EisensteinRational; every
// predicate is exact integer arithmetic (__int128 for the
// cross-multiplied comparisons, internal to the .cc).  No floating point
// enters any decision.  Width contract: Eisenstein components are int
// (the eisenstein.hh sector-geometry banner); denominators and sample
// weights must fit int after common-denominator scaling -- checked
// narrows throw instead of silently truncating.
//
// Failure contract: every failure here is a DEEP invariant violation
// (malformed pipeline state, a chart that does not verify), thrown
// uniformly as std::logic_error.  There is no named-outcome enum: the
// pipeline's own prepare gate (StageError) is the modeled-failure
// boundary; past it, the atlas either works or the input violated its
// contract.
// =====================================================================

#include "fullerenes/eisenstein.hh"        // Eisenstein, EisensteinRational, LatticeIsometry
#include "fullerenes/eisenstein_paint.hh"  // SurfaceParametrization + the cell primitives
#include "fullerenes/triangulation.hh"
#include "fullerenes/delaunay.hh"

#include <array>
#include <unordered_map>
#include <vector>

namespace eisenstein_paint {

// One live iDT 2-cell as the atlas sees it: its chart (corner cone ids +
// corner positions in its own Eisenstein frame) and its lattice claims.
// Derived from the SurfaceParametrization's Cell/LatticeMap -- the claim
// map is the hash-indexed view of the LatticeMap that the exact tracing
// needs for O(1) position lookup.
struct AtlasCell {
  bool ok = false;
  std::array<int, 3>        corners{ -1, -1, -1 };  // corner cone ids (CCW)
  std::array<Eisenstein, 3> P{};                    // corner positions in the cell's own frame
  std::unordered_map<Eisenstein, int> claim;        // lattice pos -> T_sorted vertex id
};

// A T_sorted face developed into one chart: its three corners (CCW, unit
// triangle -- verified exactly on resolution) plus a located interior/edge
// point that seeds sample traces.
struct ResolvedFace {
  bool ok = false;
  int  cell = -1;                    // host cell of `anchor`
  std::array<Eisenstein, 3> P{};     // face corners (CCW) in `cell`'s frame
  EisensteinRational anchor;         // located midpoint of one face edge (inside `cell`)
};

// The atlas.  An OPEN struct: the fields past `cells` are the internal
// combinatorics (T_sorted face/edge tables, transitions, anchored edges,
// BFS routing, the lazy resolved-face cache), public for diagnostic
// tooling but subject to change with the implementation -- treat them as
// internal-facing.  (Per-vertex occurrences live on the
// SurfaceParametrization; locate_vertex is a Layer-I function.)
struct CellAtlas {
  const DelaunayTriangulation* D = nullptr;  // the charted cone iDT
  TriangulationView            T;            // T_sorted (non-owning view)
  std::vector<AtlasCell>       cells;        // indexed by D face id

  // Directed pair key shared by `trans`, `arc_face`, `edge_id` (and by
  // callers that query them): both ids packed into one 64-bit value.
  static constexpr long long arc_key(int a, int b) {
    return ((long long)a << 32) ^ (unsigned)b;
  }
  // Undirected variant: the arc key of the sorted pair (edge_id's key).
  static constexpr long long edge_key(int u, int v) {
    return u < v ? arc_key(u, v) : arc_key(v, u);
  }

  // --- internal-facing fields (see the struct banner) -----------------
  std::unordered_map<long long, LatticeIsometry> trans;  // arc_key(f_from, f_to) -> transition
  // T_sorted combinatorics
  std::vector<tri_t>                 tface;       // T_sorted faces, CCW
  std::unordered_map<long long, int> arc_face;    // arc_key(u, v) -> face id (arc u->v)
  std::vector<std::array<int, 2>>    tedge;       // undirected edges (u < v)
  std::unordered_map<long long, int> edge_id;     // arc_key(u, v), u < v -> edge id
  std::vector<std::array<int, 2>>    edge_faces;  // per edge, its two faces
  // anchored edges + BFS routing
  struct AnchorEdge { int cell; int u, v; Eisenstein pu, pv; };  // u, v claimed by `cell`
  std::vector<int>        anchor_of_edge;         // edge id -> anchors index, or -1
  std::vector<AnchorEdge> anchors;
  std::vector<int>        bfs_parent_edge;        // edge id -> parent edge (-1 root)
  std::vector<int>        bfs_via_face;           // edge id -> face shared with parent
  // lazy per-face resolution cache (resolve_face)
  std::vector<ResolvedFace> resolved;

  // The chart transition mapping f_from's frame onto f_to's.
  //   @pre    f_from and f_to are edge-adjacent live cells of D.
  //   @throws std::logic_error when no transition exists (deep invariant).
  LatticeIsometry transition(int f_from, int f_to) const;
};

// Build the atlas over P's charts.
//
// P may be parametrized on the raw dual_idt (recommended for pure
// surface work -- no Alexandrov solve / 3D embedding required) or on a
// realized polytope's iDT.  The atlas is INTRINSIC-ONLY: everything
// here lives in the flat cone metric.  It is also curvature-sign-
// agnostic: nothing assumes cone degrees < 6 or any bound on the cell
// diameter, so fulleroid-style inputs (cone angles > 2pi) pass through
// unchanged.
//
// No chart is recomputed: the parametrization's cells and lattice maps
// are indexed into hash form and the cross-cell structure (transitions,
// anchored edges, BFS routing) is built on top.  P (and the D it points
// to) must outlive the atlas.
//
//   @throws std::logic_error on any deep invariant violation (dead cell
//           embedded, unclaimed lattice adjacency, uncovered edge graph, ...).
CellAtlas build_atlas(const SurfaceParametrization& P);

// Trace the straight segment a -> b (both in `cell`'s frame; a inside the
// CLOSED cell) through the cell complex.  Returns b's host cell, b
// rewritten into that cell's frame, and the composed transition applied
// along the way -- the "carry companion": apply `carried` to any other
// point expressed in the starting chart to move it into the destination
// chart alongside b.
//
// Exact: every predicate is integer (__int128 for the cross-multiplied
// comparisons); the crossing parameter is monotone, so the walk provably
// terminates within its in-binary cap (4 * #cells + 64 crossings).
//
//   @throws std::logic_error when no exit edge exists or the cap trips
//           (both deep invariants for a inside the closed cell).
struct CellTrace { int cell; EisensteinRational pos; LatticeIsometry carried; };
CellTrace trace_segment(const CellAtlas& A, int cell,
                        EisensteinRational a, EisensteinRational b);

// Resolve face f (T_sorted face id): develop its three corners into one
// chart and locate the midpoint of one of its edges.  Lazy and cached in
// A.resolved.  Routing: the BFS parent chain back to an anchored edge,
// executed forward as midpoint hops (each an exact trace_segment).  On
// arrival the chart is verified exactly -- a unit CCW triangle -- so any
// orientation-convention drift throws instead of corrupting downstream
// samples.
//
//   @throws std::logic_error on any failed hop or failed verification.
const ResolvedFace& resolve_face(CellAtlas& A, int f);

// Locate the barycentric sample (a*u + b*v + c*w) / den of face f, where
// (u, v, w) are the face's corners in T_sorted CCW order: resolves the
// face, then traces from its anchor to the sample.  Returns the sample's
// host cell + exact position in that cell's frame.
//   @pre    a, b, c >= 0, den > 0, and all four fit int (checked).
//   @throws std::logic_error on a violated width bound or any failed
//           resolve/trace (deep invariants).
struct CellPoint { int cell; EisensteinRational pos; };
CellPoint locate_sample(CellAtlas& A, int f, long a, long b, long c, long den);

// (locate_vertex is a Layer-I function on SurfaceParametrization -- a
// vertex is located by its first recorded occurrence, no atlas needed.)

}  // namespace eisenstein_paint
