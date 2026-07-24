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
//   evaluate(A, ...)   -- integer barycentric combination of the three
//                         corner anchors per cell, gcd-reduced on edge
//                         points so adjacent cells produce bit-identical
//                         output across shared edges (idempotent paint).
//
// Cubic-metric warm start (cubic_geometry): the carbon-atom geometry
// with the global shape of the CUBIC polyhedral metric (flat unit
// pentagons/hexagons; see AlexandrovIDTCubic, delaunay_alexandrov.hh).
// M_cubic = sqrt(3) * M_dual outside 12 constant-size pentagon caps, so
// the dual-metric barycentric weights remain geometrically faithful
// there.  Construction:
//
//   c1. realize_cubic(T): AlexandrovIDTCubic -> exact 3D positions of
//       the 20..60 pentagon-incident cubic vertices (= dual triangles).
//   c2. pentagon_anchors: anchor each of the 12 dual cones at the
//       centroid of its pentagon's 5 surrounding cone vertices.
//   c3. realize_dual(S) + parametrize + evaluate with the c2 anchors:
//       the charts are evaluated on the FLAT cells of the dual
//       polytope, so the only approximation is the (small, cap-local)
//       dual-vs-cubic shape difference -- never a chord across a
//       non-flat raw-iDT cell.  (Painting over the raw iDT is what
//       mangled elongated isomers, e.g. C392-[1,8,10,12,14,120,132,
//       157,181,189,191,198]-fullerene; the API no longer expresses it.)
//   c4. assemble_cubic: cubic vertex U = centroid of its painted dual
//       corners; pentagon-incident cubic vertices are then overwritten
//       with their exact polytope positions from c1.
//
// Because the cubic anchors are a heuristic substitution (unlike the
// exact dual paint), CubicGeometry carries bond-length statistics so
// callers can see the warm-start quality instead of trusting stage OK.
//
// Failure handling: composable functions throw StageError; the facades
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
// Throws StageError: IDT_COMPUTE / NON_SIMPLICIAL / ALEXANDROV per the
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
//   - F.ok == true
//   - anchors.size() == n_cones, finite values
//   - pos3d.size() >= the charted complex's vertex count
void interpolate_cell(const Cell& F,
                      const LatticeMap& lmap,
                      const std::vector<coord3d>& anchors,
                      int n_cones,
                      std::vector<coord3d>& pos3d);

// Evaluate every chart of A against `anchors` (position of cone c at
// anchors[c], c < A.n_cones) and back-permute to ORIGINAL labels via
// `perm` (= SortedDual::perm).  The complex A was charted on must have
// flat cells w.r.t. the surface the anchors live on -- i.e. A must be
// parametrize(P.D, S) for a realized polytope P whose cone positions
// (or a same-surface substitute like pentagon_anchors) are the anchors.
// Throws StageError(INTERPOLATE).
std::vector<coord3d> evaluate(const SurfaceParametrization& A,
                              const std::vector<coord3d>& anchors,
                              const Permutation& perm);

// =====================================================================
// The realized cubic polytope + cubic assembly.
// =====================================================================

struct CubicPolytope {
    std::vector<coord3d> cone_pos;        // pentagon-incident cubic vertices (20..60)
    std::vector<tri_t>   cone_triangle;   // their dual triangles, T ORIGINAL labels (CCW)
};

// AlexandrovIDTCubic realization of T's cubic metric.  Throws
// StageError mirroring realize_dual's staging (kis/iDT trouble ->
// IDT_COMPUTE, non-simplicial cubic iDT -> NON_SIMPLICIAL, any other
// non-validating outcome -> ALEXANDROV).
CubicPolytope realize_cubic(const TriangulationView& T);

// Anchor each of the 12 dual cones at the centroid of its pentagon's 5
// surrounding cubic cone vertices.  Every dual triangle containing a
// deg-5 vertex is pentagon-incident and hence a cone of C, so each
// pentagon accumulates exactly 5 contributions; anything else throws
// (broken cone bookkeeping).  anchors[perm[p]] for pentagon p.
std::vector<coord3d> pentagon_anchors(const TriangulationView& T,
                                      const Permutation& perm,
                                      const CubicPolytope& C);

// Cubic vertex U = dual triangle T.triangles()[U] (the T.dual_graph()
// labelling).  Position: centroid of the painted dual corners, then
// exact polytope positions for the pentagon-incident cubic vertices.
std::vector<coord3d> assemble_cubic(const TriangulationView& T,
                                    const std::vector<coord3d>& dual_coords,
                                    const CubicPolytope& C);

// =====================================================================
// Facades.  Never throw; check .status.
// =====================================================================

// Exact dual-metric geometry: every dual vertex placed exactly on the
// Alexandrov dual polytope's surface.  coords in T's ORIGINAL labels.
struct DualGeometry {
    std::vector<coord3d> coords;      // size T.N on OK; empty otherwise
    Status               status;
};
DualGeometry dual_geometry(const TriangulationView& T);

// Cubic-metric warm start for the carbon atoms.  cubic_coords[U] is the
// position of cubic vertex U in the T.dual_graph() labelling (U indexes
// T.triangles()); dual_coords is the painted dual (T original labels),
// kept for diagnostics and deltahedron-side consumers.
//
// bond_* are the cubic-graph edge-length statistics of cubic_coords
// (nominal bond = 1): the warm start is a heuristic (see banner), so
// its quality is reported rather than silently trusted.
struct CubicGeometry {
    std::vector<coord3d> cubic_coords;   // size 2*T.N - 4 on OK
    std::vector<coord3d> dual_coords;    // size T.N on OK, original labels
    Status               status;
    double bond_min = 0, bond_mean = 0, bond_max = 0, bond_cv = 0;
};
CubicGeometry cubic_geometry(const TriangulationView& T);

}  // namespace eisenstein_paint
