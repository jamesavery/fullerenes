#pragma once

// =====================================================================
// Eisenstein-paint: deterministic warm-start geometry for fullerene
// dual triangulations.
//
// Algorithm:
//   1. Sort the dual triangulation T cones-first (deg != 6 first) to
//      get T_sorted; record the permutation for back-mapping.
//   2. Compute the intrinsic Delaunay triangulation D = iDT(T_sorted).
//   3. Run the Bobenko-Izmestiev Alexandrov solver on D to obtain the
//      12 cone-vertex 3D positions; D is mutated to its post-flip
//      simplicial form.
//   4. For every iDT 2-cell F (= live face of D), embed F into its
//      own Eisenstein-lattice frame: place corners (P0, P1, P2) from
//      iDT lengths + interior angle, then walk T_sorted in F's frame
//      along each CCW arc to fill F's edge strips.
//   5. Enumerate every Z[w] lattice point inside F's polygon and map
//      it to its T_sorted vertex id (Abrash polygon scan + per-row
//      east walker seeded from the edge strips).
//   6. For each lattice point, evaluate its 3D position by integer
//      barycentric interpolation of the 3 cone corner positions
//      (with gcd-reduction on edge points so the two adjacent cells
//      produce bit-identical 3D output across the shared edge).
//   7. Back-permute the per-vertex 3D coords to T's original labels.
//
// Tier 1   `eisenstein_paint::run(T)`            full pipeline -> Result.
// Tier 1'  `eisenstein_paint::run_cubic(T)`      cubic-metric variant -> CubicResult.
// Tier 2   `eisenstein_paint::prepare(T)`        prelude (sort + iDT + Alexandrov)
// Tier 2'  `eisenstein_paint::prepare_iDT(T)`    iDT only, no Alexandrov (~20x faster)
// Tier 2'' `eisenstein_paint::prepare_inputs(T)` slim return for diagnostic tools
// Tier 2c  `eisenstein_paint::prepare_cubic(T)`  cubic prelude (sort + iDT + IDTCubic)
//
// Cubic-metric variant (run_cubic): same integer paint machinery, but
// the three cell-corner 3D anchors come from the AlexandrovIDTCubic
// polytope (the convex realization of the CUBIC polyhedral metric --
// flat unit pentagons/hexagons, kappa = k*pi/15 at the 20..60
// pentagon-incident cubic vertices) instead of the 12-cone dual-metric
// polytope.  M_cubic = sqrt(3) * M_dual outside 12 constant-size
// pentagon caps (see delaunay_alexandrov.hh), so the dual-metric
// barycentric weights remain geometrically faithful there, and the
// weights are scale-free per cell.  Construction:
//
//   c1. AlexandrovIDTCubic::solve(T) -> exact 3D positions of the
//       20..60 pentagon-incident cubic vertices (= dual triangles).
//   c2. Anchor each of the 12 dual cones at the centroid of its
//       pentagon's 5 surrounding cone vertices (the pentagon is
//       intrinsically flat; the centroid is its exact center whenever
//       the pentagon is realized flat on the polytope).
//   c3. Paint the dual vertices with the standard integer pipeline
//       (raw iDT of the dual metric; anchors from c2).
//   c4. Cubic vertex U (= dual triangle U in T.triangles() /
//       T.dual_graph() order) = centroid of its painted dual corners;
//       pentagon-incident cubic vertices are then overwritten with
//       their exact polytope positions from c1.
//
// The result is a warm start for the cubic-graph force-field optimizer
// with the global shape of the cubic metric (edge lengths ~1 bond).
//
// Per-cell primitives (Tier 3, for callers that want to inspect or
// drive individual phases): `embed_cell`, `embed_all_cells`,
// `enumerate_cell_lattice`, `interpolate_cell`.
//
// Failure handling: Result.stage names where the pipeline stopped,
// with Result.why carrying a free-form detail string.  Internal stage
// helpers throw StageError(stage, msg); run() converts to Result.
// Callers that want to catch directly may invoke the prepare* /
// per-cell primitives and let StageError propagate.
// =====================================================================

#include "fullerenes/triangulation.hh"
#include "fullerenes/delaunay.hh"
#include "fullerenes/permutation.hh"
#include "fullerenes/geometry.hh"
#include "fullerenes/eisenstein.hh"
#include "fullerenes/delaunay_strip.hh"      // StripVertex, Strip

#include <iosfwd>
#include <stdexcept>
#include <string>
#include <vector>

namespace eisenstein_paint {

// =====================================================================
// Result of run(T): coords + a stage discriminator if anything failed.
// =====================================================================

struct Result {
    enum class Stage {
        OK,
        IDT_COMPUTE,           // DelaunayTriangulation::compute throw or D.nv != 12
        ALEXANDROV,            // AlexandrovSolver: any non-validating outcome
        NON_SIMPLICIAL,        // iDT carries multi-edges / self-loops / non-triangle face
        EMBED,                 // some cell's embed_cell returned ok=false (or lmap throw)
        COVERAGE,              // some hex vertex unclaimed or claimed >= 3 times
        INTERPOLATE,           // interpolate_cell threw
        UNEXPECTED,            // unrecognised std::exception escaped the pipeline
    };

    std::vector<coord3d> coords;             // size T.N (original labels) on OK; empty otherwise
    Stage                stage = Stage::OK;
    std::string          why;

    bool        ok()         const { return stage == Stage::OK; }
    const char* stage_name() const;
};

// Convenience alias: callers write `eisenstein_paint::Stage::EMBED`
// instead of the full `Result::Stage::EMBED`.
using Stage = Result::Stage;

// Stage discriminator -> lowercase name ("ok", "alexandrov", ...).
const char* stage_name(Stage s);

// Result of run_cubic(T): warm-start coordinates for the CUBIC graph.
// cubic_coords[U] is the 3D position of cubic vertex U in the
// T.dual_graph() labelling, i.e. U indexes the dual triangle
// T.triangles()[U] (compute_faces_oriented order).  dual_coords is the
// painted dual (size T.N, original labels), kept for diagnostics and
// for deltahedron-side consumers.
struct CubicResult {
    std::vector<coord3d> cubic_coords;   // size 2*T.N - 4 on OK; empty otherwise
    std::vector<coord3d> dual_coords;    // size T.N (original labels) on OK
    Stage                stage = Stage::OK;
    std::string          why;

    bool        ok()         const { return stage == Stage::OK; }
    const char* stage_name() const;
};

// =====================================================================
// Prelude state, shared by run() and the prepare* tier-2 entry points.
// =====================================================================

struct Pipeline {
    Triangulation         T_sorted;          // input with cones first (degree != 6 first)
    Permutation           perm;              // perm[u_orig] = u_sorted; .inverse() back-permutes
    DelaunayTriangulation D;                  // iDT on T_sorted (post-Alexandrov, simplicial)
    std::vector<coord3d>  cone_positions;    // 12 entries from Alexandrov, indexed in T_sorted
};

// Slim return for callers needing only T_sorted and the iDT.
// Designed for one-line structured binding:
//
//     auto [T_sorted, D] = eisenstein_paint::prepare_inputs(T);
struct Inputs {
    Triangulation         T_sorted;
    DelaunayTriangulation D;
};

// Single discriminated failure type for the prelude / paint stages.
class StageError : public std::runtime_error {
public:
    Stage stage;
    StageError(Stage s, std::string msg)
        : std::runtime_error(std::move(msg)), stage(s) {}
};

// Tier 1 -- full pipeline; never throws.
Result run(const Triangulation& T);

// Tier 1' -- cubic-metric variant (see banner); never throws.
CubicResult run_cubic(const Triangulation& T);

// Tier 2 -- prelude only; throws StageError on failure.
Pipeline prepare       (const Triangulation& T);   // sort + iDT + Alexandrov
Pipeline prepare_iDT   (const Triangulation& T);   // sort + iDT only (no Alexandrov)
Inputs   prepare_inputs(const Triangulation& T);   // structured-binding shortcut

// Tier 2c -- cubic prelude: sort + raw iDT + AlexandrovIDTCubic +
// pentagon-center anchor derivation.  P.cone_positions holds the 12
// anchors (T_sorted labels); the exact polytope data is kept alongside
// so run_cubic / diagnostics can overwrite the pentagon-incident cubic
// vertices with their true positions.  Throws StageError on failure.
struct CubicPipeline {
    Pipeline             P;                       // D is the RAW dual iDT (no dual Alexandrov)
    std::vector<coord3d> cone_vertex_positions;   // one per AlexandrovIDTCubic cone (20..60)
    std::vector<tri_t>   cone_triangle;           // per cone: dual triangle, T ORIGINAL labels (CCW)
};
CubicPipeline prepare_cubic(const Triangulation& T);

// =====================================================================
// Per-cell types and primitives (Tier 3).
// =====================================================================

// Half of one edge's Strip, scanline-indexed in a Cell's frame.
//
// `by_scanline[b - b_min]` lists every StripVertex whose F-frame
// position has second-coord = b.  Inner vectors may contain multiple
// entries (the strip is up to 1 lattice wide perpendicular to the
// arc).  Empty iff b_max < b_min.
struct EdgeStrip {
    int b_min = 0;
    int b_max = -1;
    std::vector<std::vector<StripVertex>> by_scanline;

    bool empty() const { return b_max < b_min; }
};

// One iDT 2-cell embedded into its own Eisenstein-lattice frame.
struct Cell {
    int        cell_id = -1;
    int        c0 = -1, c1 = -1, c2 = -1;       // CCW corner T_sorted vertex ids
    Eisenstein P0{0, 0};                         // = (0, 0)
    Eisenstein P1{0, 0};
    Eisenstein P2{0, 0};

    // Per-edge strips, in the cell's own frame, scanline-indexed.
    // edge_01 covers the arc c0 -> c1, etc.
    EdgeStrip edge_01;
    EdgeStrip edge_12;
    EdgeStrip edge_20;

    bool ok = false;
};

// A corner placement (P0, P1, P2) of one iDT cell in its own Eisenstein frame.
struct CornerCandidate { Eisenstein P0, P1, P2; };

// Enumerate valid corner placements from the three integer edge length-squares
// (N01, N12, N20) + the interior angle alpha_0 at c0 (L20 = length of edge
// c2->c0).  Returns 1 candidate for single-orbit N01, up to 2 for a split-prime
// N01 (its two mirror geodesics).  Choosing which candidate matches the actual
// surface is the caller's job: embed_cell (with T) or build_intrinsic_atlas
// (T-free, by cross-edge transition consistency).
std::vector<CornerCandidate>
enumerate_corner_candidates(double L20, double alpha_0,
                            long N01, long N12, long N20);

// Embed one live iDT cell into its own Eisenstein frame by walking
// T_sorted.  F's corners are placed from iDT lengths + the interior
// angle at c0 (= D.he_angle of the f_he half-edge); each CCW arc's
// strip is filled by an F-frame walker targeting the next corner.
//
// Preconditions:
//   - 0 <= cell_id < D.nf  and  D.f_he[cell_id] >= 0 (live face).
//
// Throws on contract violations (iDT angle/length not realisable in
// the Eisenstein lattice; walker cannot reach a target from its start).
Cell embed_cell(const DelaunayTriangulation& D,
                const Triangulation& T_sorted,
                int cell_id);

// Embed every live iDT cell in parallel.  Returns a vector of size
// D.nf; dead slots have ok == false.
std::vector<Cell> embed_all_cells(const DelaunayTriangulation& D,
                                   const Triangulation& T_sorted);

// (F-frame lattice pos) -> T_sorted vertex_id, for every Z[w] lattice
// point inside F's polygon.  Computed by polygon scan + per-scanline
// horizontal east-walker seeded from the edge strips.
struct LatticeMap {
    int cell_id = -1;
    // Ordered by scanline then by a (scan order).  Each entry is
    // ((a, b), vertex_id).
    std::vector<std::pair<Eisenstein, int>> entries;
};

LatticeMap enumerate_cell_lattice(const Cell& F,
                                  const Triangulation& T_sorted);

// Interpolate cone positions over F's lattice via integer barycentric,
// writing into pos3d[vertex_id] for every non-cone entry.  Idempotent
// on shared on-edge hex (one barycentric weight is 0 -> reduce by gcd
// of the other two before computing, eliminating the cell-wedge
// denominator from the FP arithmetic so both adjacent cells produce
// bit-identical 3D output across the shared edge).
//
// Preconditions:
//   - F.ok == true
//   - cone_positions.size() == 12 with finite values
//   - pos3d.size() >= T_sorted.N
void interpolate_cell(const Cell& F,
                      const LatticeMap& lmap,
                      const std::vector<coord3d>& cone_positions,
                      std::vector<coord3d>& pos3d);

// =====================================================================
// TikZ visualisations (diagnostic only).
// =====================================================================

// Self-contained TikZ standalone document of one Cell.  Draws F's
// lattice triangle outline (iDT-geodesic arcs P0->P1->P2->P0 in red),
// T_sorted edges between strip vertices from any of F's three edge
// strips (gray), and the strip vertices themselves -- corners as red
// squares, interior strip vertices as circles colour-coded by source
// edge (edge_01 blue, edge_12 green, edge_20 orange).
void dump_cell_tikz(const Cell& F,
                    const Triangulation& T_sorted,
                    std::ostream& os,
                    double scale = 0.7);

// Self-contained TikZ standalone document of F's lattice triangle +
// every lattice point in F's polygon labelled with its T_sorted
// vertex_id (from the LatticeMap).  Hex vertices are coloured by
// whether they're on F's boundary (green = chain pixel on an iDT edge
// of F; orange = on-edge hex with gcd > 1) or strict interior (blue).
// Cones at the corners are drawn as red squares.  If
// `highlight_vertex_id >= 0` and that vertex appears in the lattice
// map, draw it with a thick red outline.
void dump_lattice_map_tikz(const Cell& F,
                           const LatticeMap& lmap,
                           const Triangulation& T_sorted,
                           std::ostream& os,
                           double scale = 0.7,
                           int highlight_vertex_id = -1);

}  // namespace eisenstein_paint
