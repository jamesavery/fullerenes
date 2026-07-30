#pragma once

// =====================================================================
// Eisenstein-paint, Layer II: from a cone-point Alexandrov embedding to
// a full fullerene embedding.  This is the only layer that touches 3D.
//
// Layer I (eisenstein_paint.hh) charts the surface intrinsically: every
// vertex gets exact integer barycentric coordinates w.r.t. its cell's
// corners.  Layer II supplies the corners' positions in R^3 and
// evaluates the charts:
//
//   realize_dual(S)    -- Bobenko-Izmestiev Alexandrov solve of the
//                         dual cone metric.  The returned D is the
//                         post-flip iDT at kappa = 0: every 2-cell is a
//                         FLAT face of the polytope (extrinsically
//                         Delaunay), T-bar(0) is simple, and the
//                         reconstruction is validated convex/embedded.
//                         That is precisely the precondition under
//                         which evaluating the integer barycentric
//                         charts places every vertex EXACTLY on the
//                         polytope surface.  (T(0) itself may carry
//                         flat-diagonal multi-edges on exactly-
//                         symmetric isomers -- the charts handle them;
//                         simpliciality is NOT required here.)
//   evaluate_sorted(A, ...) -- integer barycentric combination of the
//                         three corner anchors per cell, gcd-reduced on
//                         edge points so adjacent cells produce
//                         bit-identical output across shared edges
//                         (idempotent paint); result in sorted labels.
//   evaluate(A, ...)   -- the same, back-permuted to original labels.
//
// Cubic-metric geometry (cubic_geometry): the carbon-atom geometry
// EXACTLY on the cubic Alexandrov polytope (the convex realization of
// the CUBIC polyhedral metric -- flat unit pentagons/hexagons; see
// AlexandrovIDTCubic, delaunay_alexandrov.hh).  Construction (the
// flip-tape transport of DESIGN-cubic-exact-paint.md):
//
//   c1. realize_cubic(T): AlexandrovIDTCubic with point tracking.  The
//       kis complex's flat vertices -- hexagon centers, pentagon
//       centers, hexagon-only cubic vertices -- are tracked as
//       (cell, barycentric) locations through the flat-vertex removal
//       and every homotopy flip, ending in the kappa = 0 iDT whose
//       cells are flat pieces of the polytope.
//   c2. Pentagon-incident cubic vertices are the cones: exact solver
//       positions.
//   c3. Every tracked kis vertex evaluates barycentrically against its
//       cell's three cone positions -- exactly ON the polytope surface.
//
// No dual-metric solve, no anchor substitution, no centroid step (the
// pre-2026-07-24 heuristics): the only error is the Alexandrov
// solver's own convergence/reconstruction tolerance.  A cubic bond is
// an intrinsic unit segment, so its 3D chord is <= 1 (equality iff it
// does not cross a polytope crease); CubicGeometry reports bond stats
// as a self-diagnostic -- any bond above 1 + tolerance indicates a bug.
//
// Failure handling: composable functions throw PaintError; the facades
// dual_geometry / cubic_geometry never throw and return a Status.
// =====================================================================

#include "fullerenes/eisenstein_paint.hh"
#include "fullerenes/deltahedron.hh"
#include "fullerenes/geometry.hh"

#include <vector>

namespace eisenstein_paint {

// =====================================================================
// The realized dual polytope.
// =====================================================================

struct DualPolytope {
    // The intrinsic metric remains the source of truth: post-flip iDT
    // at kappa = 0 (flat cells, validated).  Charted by Layer I's
    // parametrize; the 3D positions below are only ever combined with
    // the intrinsic charts, never fed back into them.
    DelaunayTriangulation D;
    std::vector<coord3d>  cone_pos;    // one per cone, sorted labels

    // The solve's converged radii r_v (B-I: r_v = |p_v| from the apex), one
    // per D vertex in sorted labels.  Exposed because the polytope's
    // 2-skeleton is unreachable without them: AlexandrovSolver::
    // inessential_edges(D, r) and polytope_tesselation(D, r, labels) both
    // need r, and T-bar -- not the flip-dependent triangulation -- is the
    // unique, comparable form of a realized cone metric.
    std::vector<double>   r;

    // The polytope as a library Deltahedron (cone iDT 1-skeleton +
    // positions; 12 vertices with degrees up to 11 for fullerene duals
    // -- not a fullerene dual itself, no deg-5/6 assumption).  Oriented
    // neighbour rings are built directly from D's DCEL cycles.
    // @throws std::logic_error when T(0) is non-simplicial (flat-
    //         diagonal multi-edges: the 1-skeleton is not a simple
    //         graph; use AlexandrovSolver::polytope_tesselation for the
    //         polygonal T-bar(0) instead).
    Deltahedron deltahedron() const;
};

// Alexandrov realization of S's cone metric (iDT + B-I solve).
// Throws PaintError: IDT_COMPUTE / NON_SIMPLICIAL / ALEXANDROV per the
// solver's validation verdict.
DualPolytope realize_dual(const SortedDual& S);

// =====================================================================
// Chart evaluation (the paint).
// =====================================================================

// Interpolate anchor positions over F's lattice via integer barycentric,
// writing into pos3d[vertex_id] for every non-cone entry.  Idempotent
// on shared on-edge vertices (one barycentric weight is 0 -> reduce by
// gcd of the other two before combining, eliminating the cell-wedge
// denominator from the FP arithmetic so both adjacent cells produce
// bit-identical 3D output across the shared edge).
//
// Preconditions:
//   - corners/frame belong to one charted cell; entries are its lattice
//     points (scanline-major)
//   - anchors.size() == n_cones, finite values
//   - pos3d.size() >= the charted complex's vertex count
// The (frame, corners, entries) form is the mathematical core; the
// (V, f) form projects a cell out of the tables (@pre 0 <= f < V.nf,
// V.cell_live(f)) and reads n_cones from the view.
void interpolate_cell(CellFrame frame, CellCorners corners,
                      std::span<const LatticePoint> entries,
                      const std::vector<coord3d>& anchors,
                      int n_cones,
                      std::vector<coord3d>& pos3d,
                      int cell_id_for_diag = -1);
void interpolate_cell(const ParamTablesView& V, int f,
                      const std::vector<coord3d>& anchors,
                      std::vector<coord3d>& pos3d);

// Evaluate every chart of A against `anchors` (position of cone c at
// anchors[c], c < A.n_cones), in SORTED (T_sorted) labels: the result
// has A.T.N entries, cones c < n_cones at anchors[c] verbatim, every
// other vertex by barycentric interpolation in its cell's chart
// (on-edge vertices idempotently from both adjacent cells).  CAUTION:
// pair the result with the SORTED graph (A's T / SortedDual::T) --
// indexing it by original labels is plausible-looking wrong output;
// use evaluate below when original labels are wanted.  The
// complex A was charted on must have flat cells w.r.t. the surface the
// anchors live on -- i.e. A must be parametrize(P.D, S) for a realized
// polytope P whose cone positions are the anchors.  Throws
// PaintError(INTERPOLATE).
std::vector<coord3d> evaluate_sorted(const SurfaceParametrization& A,
                                     const std::vector<coord3d>& anchors);

// evaluate = back-permutation o evaluate_sorted: the same positions
// re-indexed to ORIGINAL labels via A's own stored permutation
// (= SortedDual::perm at parametrize time; perm[u_orig] = u_sorted).
std::vector<coord3d> evaluate(const SurfaceParametrization& A,
                              const std::vector<coord3d>& anchors);

// =====================================================================
// The realized cubic polytope.
// =====================================================================

struct CubicPolytope {
    // The post-flip cubic iDT at kappa = 0: every cell a flat piece of
    // the polytope, and the POINT TRACKER carrying every removed kis
    // vertex as a (cell, barycentric) location on the surface.  Tracker
    // labels are kis ids: label u < T.N is dual vertex u (a cubic face
    // center); label T.N + U is cubic vertex U (T.triangles() /
    // T.dual_graph() order).  Cone i is vertex i of D.
    DelaunayTriangulation D;
    std::vector<coord3d> cone_pos;         // pentagon-incident cubic vertices (20..60)
    std::vector<tri_t>   cone_triangle;    // their dual triangles, T ORIGINAL labels (CCW)
    std::vector<int>     cone_kis_vertex;  // kis id of cone i (= T.N + its cubic vertex U)
};

// AlexandrovIDTCubic realization of T's cubic metric, with every flat
// kis vertex tracked onto the polytope.  Throws PaintError mirroring
// realize_dual's staging (kis/iDT trouble -> IDT_COMPUTE, non-simplicial
// cubic iDT -> NON_SIMPLICIAL, any other non-validating outcome ->
// ALEXANDROV).
CubicPolytope realize_cubic(const TriangulationView& T);

// Realization of every kis vertex on the cubic polytope: cones at their
// solver positions; every tracked point by barycentric combination
// against its flat cell's three cone-corner positions -- exactly ON the
// polytope surface, since the kappa = 0 cells are flat pieces of it.
// Split by kis label: label >= Nv is cubic vertex (label - Nv), label
// < Nv is dual vertex label (a cubic face center).
struct CubicRealization {
    std::vector<coord3d> cubic_coords;   // size 2*Nv - 4, T.triangles() order
    std::vector<coord3d> dual_coords;    // size Nv, dual-vertex labels
};

// The cubic counterpart of evaluate(): evaluate the tracked kis complex
// on the polytope.  Nv = the dual triangulation's vertex count (= T.N).
// @pre  C from realize_cubic(T): every tracked label in range, cones'
//       kis ids in the cubic-vertex range
// @post coverage: every cubic vertex and every dual vertex written
//       exactly once (cones + tracked points partition the kis vertices)
// @throws PaintError(COVERAGE) when some vertex is unwritten; deep
//       invariant violations (double-written or out-of-range labels)
//       throw std::runtime_error
CubicRealization evaluate_tracked_complex(const CubicPolytope& C, int Nv);

// =====================================================================
// Facades.  Never throw; check .status.
// =====================================================================

// Exact dual-metric geometry: every dual vertex placed exactly on the
// Alexandrov dual polytope's surface.  coords in T's ORIGINAL labels.
// @post status.code in { OK, IDT_COMPUTE, NON_SIMPLICIAL, ALEXANDROV,
//       EMBED, COVERAGE, INTERPOLATE, UNEXPECTED }  (the dual facade's
//       closed failure-mode set)
struct DualGeometry {
    std::vector<coord3d> coords;      // size T.N on OK; empty otherwise
    Status               status;
};
DualGeometry dual_geometry(const TriangulationView& T);

// Cubic-metric geometry for the carbon atoms, EXACT on the cubic
// Alexandrov polytope (see banner).  cubic_coords[U] is the position of
// cubic vertex U in the T.dual_graph() labelling (U indexes
// T.triangles()); dual_coords[u] is dual vertex u (= cubic face center)
// located on the SAME polytope, kept for deltahedron-side consumers.
//
// bond_* are the cubic-graph edge-length statistics of cubic_coords: a
// bond is an intrinsic unit segment, so bond_max <= 1 up to solver
// tolerance (crease-crossing bonds are strictly shorter chords); a
// larger value indicates a bug, not geometry.
// @post status.code in { OK, IDT_COMPUTE, NON_SIMPLICIAL, ALEXANDROV,
//       COVERAGE, UNEXPECTED }  (the cubic facade's closed failure-mode
//       set; EMBED/INTERPOLATE are dual-path-only)
struct CubicGeometry {
    std::vector<coord3d> cubic_coords;   // size 2*T.N - 4 on OK
    std::vector<coord3d> dual_coords;    // size T.N on OK, original labels
    Status               status;
    double bond_min = 0, bond_mean = 0, bond_max = 0, bond_cv = 0;
};
CubicGeometry cubic_geometry(const TriangulationView& T);

}  // namespace eisenstein_paint
