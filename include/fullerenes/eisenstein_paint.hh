#pragma once

// =====================================================================
// Eisenstein-paint, Layer I: INTRINSIC surface parametrization of a
// fullerene dual triangulation.  Everything in this header lives in the
// flat cone metric -- integer/combinatorial data only.  No function or
// type here mentions a 3D coordinate; realizing the parametrization in
// R^3 is Layer II (eisenstein_paint_geometry.hh).
//
// Pipeline (intrinsic part):
//   1. sorted_dual(T): relabel the dual triangulation cones-first
//      (deg != 6 first); record the permutation for back-mapping.
//   2. dual_idt(S): the intrinsic Delaunay triangulation of the cone
//      metric.  A delta-complex: multi-edges between one cone pair are
//      legitimate and accepted (the walker's cone-interiority + cross-
//      edge discriminators develop them correctly); what is rejected is
//      what the cell machinery genuinely cannot chart -- a self-loop
//      (a face with a repeated corner; see refactor-debt entry
//      2026-07-24-paint-repeated-corner-cells for lifting this).
//   3. Per iDT 2-cell F: embed F's three corners into F's own
//      Eisenstein-lattice frame from iDT lengths + interior angle, then
//      walk T_sorted in F's frame along each CCW arc to fill F's edge
//      strips (embed_cell).
//   4. Enumerate every Z[w] lattice point inside F's triangle and map
//      it to its T_sorted vertex id -- Abrash polygon scan + per-row
//      east walker seeded from the edge strips (enumerate_cell_lattice).
//   5. parametrize(D, S) bundles 3-4 over all live cells, checks
//      coverage (every non-cone vertex claimed by 1 or 2 cells), and
//      records every (cell, position) occurrence per vertex.
//
// The result is a chart-wise integer parametrization of the surface:
// each vertex has exact integer barycentric coordinates w.r.t. its
// cell's corners (via integer_barycentric, barycentric.hh).  Pure
// surface computations (the exact cross-cell atlas eisenstein_atlas.hh,
// unfoldings, surface sampling) build on this with no embedding at all.
//
// The intrinsic layer works on any intrinsically-Delaunay cone complex
// whose cells are Eisenstein-realisable; it is curvature-sign-agnostic
// (nothing assumes cone degree < 6).  Cone vertices are the sorted
// labels < S.n_cones (= 12 for fullerene duals).
//
// Failure handling: the composable functions throw PaintError(code,
// msg) -- one channel, typed at the throw site.  The Layer-II facades
// convert to a Status value; callers of the composable functions let
// PaintError propagate.
//
// GPU note: every read-only input and product of this layer has a
// trivially-copyable span form -- TriangulationView for the mesh,
// DelaunayView for the cone iDT (delaunay_view.hh), and the flat
// tables of eisenstein_paint_tables.hh for the parametrization.
// =====================================================================

#include "fullerenes/triangulation.hh"
#include "fullerenes/delaunay.hh"
#include "fullerenes/permutation.hh"
#include "fullerenes/eisenstein.hh"
#include "fullerenes/delaunay_strip.hh"      // StripVertex, Strip
#include "fullerenes/eisenstein_paint_tables.hh"  // the flat table form

#include <array>
#include <iosfwd>
#include <stdexcept>
#include <string>
#include <vector>

namespace eisenstein_paint {

// =====================================================================
// Failure vocabulary, shared by both layers: a closed set of named
// failure REASONS (style-failures.md's Code/metadata shape).
// =====================================================================

enum class Code {
    OK,
    IDT_COMPUTE,           // DelaunayTriangulation::compute threw
    ALEXANDROV,            // Layer II: any non-validating solver outcome
    NON_SIMPLICIAL,        // non-embeddable iDT: self-loop / non-triangle face, or
                           // Alexandrov FAIL_NOT_SIMPLE (metric not realizable as a
                           // simple convex polytope). Multi-edges are ACCEPTED.
    EMBED,                 // some cell admits no consistent chart (or lmap threw)
    COVERAGE,              // some non-cone vertex claimed 0 or >= 3 times
    INTERPOLATE,           // Layer II: interpolate_cell threw
    UNEXPECTED,            // unrecognised std::exception escaped a facade
};

// Code discriminator -> lowercase name ("ok", "alexandrov", ...).
// (The emitted strings are sweep-output-stable; do not change them.)
const char* code_name(Code c);

// Failure value carried by the Layer-II facade results.
struct Status {
    Code        code = Code::OK;
    std::string why;

    bool        ok()   const { return code == Code::OK; }
    const char* name() const { return code_name(code); }
};

// Single discriminated failure type thrown by the composable functions.
class PaintError : public std::runtime_error {
public:
    Code code;
    PaintError(Code c, std::string msg)
        : std::runtime_error(std::move(msg)), code(c) {}
};

// =====================================================================
// sorted_dual: the cones-first relabelling of the input dual.
// =====================================================================

// Owns the relabelled triangulation (the allocation point of the
// pipeline); everything downstream reads it through TriangulationView.
struct SortedDual {
    Triangulation T;         // input with cones first (degree != 6 first)
    Permutation   perm;      // perm[u_orig] = u_sorted; .inverse() back-permutes
    int           n_cones = 0;   // # vertices with degree != 6 (12 for fullerene duals)
};

SortedDual sorted_dual(const TriangulationView& T);

// The intrinsic Delaunay triangulation of S's cone metric -- a
// delta-complex in general (multi-edges accepted).  Throws PaintError
// (IDT_COMPUTE on compute failure; NON_SIMPLICIAL on a self-loop or a
// malformed face, which the cell machinery cannot chart).
DelaunayTriangulation dual_idt(const SortedDual& S);

// =====================================================================
// Per-cell types and primitives.
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

// One iDT 2-cell embedded into its own Eisenstein-lattice frame -- the
// cell's CHART.  Corner k sits at P[k]; arc k runs corners[k] ->
// corners[(k+1)%3] and carries edge strip edge[k].  (Same k-indexed
// corner/position form as ParamTablesView's corner_ids/frame_points,
// the flat product this construction fills.)
struct Cell {
    int                         cell_id = -1;
    std::array<int, 3>          corners{ -1, -1, -1 };  // CCW T_sorted vertex ids
    std::array<Eisenstein, 3>   P{};                    // corner positions; P[0] == (0,0)
    std::array<EdgeStrip, 3>    edge{};                 // strip of arc k: corners[k] -> corners[k+1]

    bool ok = false;
};

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
                const TriangulationView& T_sorted,
                int cell_id);

// Embed every live iDT cell in parallel.  Returns a vector of size
// D.nf; dead slots have ok == false.
std::vector<Cell> embed_all_cells(const DelaunayTriangulation& D,
                                  const TriangulationView& T_sorted);

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
                                  const TriangulationView& T_sorted);

// =====================================================================
// SurfaceParametrization: every cell charted, every vertex located.
// =====================================================================

// The intrinsic product of the pipeline, CSR-backed: all live cells
// embedded and scanned, coverage-checked, every product a flat table
// (element structs + extents in eisenstein_paint_tables.hh).
// parametrize no longer STORES Cell/LatticeMap -- they remain the
// public construction vocabulary (embed_cell, enumerate_cell_lattice,
// the strip-level tikz diagnostics, which draw data the tables
// intentionally omit); readers go through view()'s accessors.
//
// LIFETIME/STABILITY: D and T alias the inputs' buffers, not their
// owners.  Valid while the owners live AND do not grow -- every
// owner-level growth op (ensure_*/alloc_*/split_face/
// bisect_multi_edges) repoints the owner and dangles these views; a
// copy of the owner does not re-bind them.  Re-run parametrize after
// any growth.  D carries the view's mutation API; holders of a
// SurfaceParametrization must treat it as read-only.
struct SurfaceParametrization {
    DelaunayView      D;                // the charted complex (view)
    TriangulationView T;                // T_sorted
    int               n_cones = 0;      // = SortedDual::n_cones

    // Owning vectors backing ParamTablesView, in its span order.
    std::vector<CellCorners>      cells;        // [nf]
    std::vector<CellFrame>        frames;       // [nf]  the winning chart
    std::vector<ScanBlock>        scans;        // [nf]
    std::vector<ScanRow>          rows;         // CSR row blocks
    std::vector<int32_t>          entry_first;  // [nf+1]
    std::vector<LatticePoint>     entries;      // scanline-major per cell
    std::vector<int32_t>          perm;         // [N_orig] orig -> sorted
    std::vector<int32_t>          occ_first;    // [N_sorted+1]
    std::vector<VertexOccurrence> occ;

    // The span mirror of this parametrization; valid while *this is
    // alive and unmodified.
    ParamTablesView view() const {
        ParamTablesView v;
        v.nf = (int32_t)cells.size();
        v.N_sorted = (int32_t)T.N;
        v.N_orig = (int32_t)perm.size();
        v.n_cones = n_cones;
        v.entries_total = (int64_t)entries.size();
        v.cells = cells; v.frames = frames; v.scans = scans; v.rows = rows;
        v.entry_first = entry_first; v.entries = entries; v.perm = perm;
        v.occ_first = occ_first; v.occ = occ;
        return v;
    }
};

// Chart all live cells of D and check coverage: every non-cone vertex
// (id >= S.n_cones) must appear in exactly 1 or 2 charts (2 = on a
// shared iDT edge).  Throws PaintError(EMBED) when a live cell fails to
// embed or scan, PaintError(COVERAGE) on a coverage violation.
//
// Works on any intrinsically-Delaunay D over S.T -- the raw dual_idt(S)
// or Layer II's realized DualPolytope::D alike.
SurfaceParametrization parametrize(const DelaunayTriangulation& D,
                                   const SortedDual& S);

// A T_sorted vertex is itself a lattice point: locate it via its first
// recorded occurrence (the flat VertexOccurrence: host cell + chart
// position (a, b)).  Throws std::logic_error if v is unclaimed
// (impossible on a parametrize() output: coverage checked).
VertexOccurrence locate_vertex(const SurfaceParametrization& P, int v);

// =====================================================================
// TikZ visualisations (diagnostic only; intrinsic -- lattice pictures).
// =====================================================================

// Self-contained TikZ standalone document of one Cell.  Draws F's
// lattice triangle outline (iDT-geodesic arcs P0->P1->P2->P0 in red),
// T_sorted edges between strip vertices from any of F's three edge
// strips (gray), and the strip vertices themselves -- corners as red
// squares, interior strip vertices as circles colour-coded by source
// edge (edge_01 blue, edge_12 green, edge_20 orange).
void dump_cell_tikz(const Cell& F,
                    const TriangulationView& T_sorted,
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
                           const TriangulationView& T_sorted,
                           std::ostream& os,
                           double scale = 0.7,
                           int highlight_vertex_id = -1);

}  // namespace eisenstein_paint
