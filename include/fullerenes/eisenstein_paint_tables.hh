#pragma once

// =====================================================================
// Eisenstein-paint, flat table form: the trivially-copyable, span-based
// vocabulary of the paint pipeline's products, completing the family
//   eisenstein_paint.hh          -- intrinsic construction (Layer I)
//   eisenstein_paint_geometry.hh -- extrinsic realization (Layer II)
//   eisenstein_paint_tables.hh   -- the flat data form (this header)
//
// Two table shapes share one vocabulary, matching the two pipeline
// stages:
//
//   ParamTablesView  -- the POST-EMBED product: one winning frame per
//     live cell, entries filled, occurrence CSR built.  This is what
//     SurfaceParametrization (eisenstein_paint.hh) stores and what the
//     atlas consumes.  Non-owning, trivially copyable, device-legal.
//
//   CandidateTables  -- the PRE-EMBED prepare form: up to two candidate
//     frames per cell with per-candidate triangle scans, plus the
//     scratch capacities an executor needs to size embed/enumerate
//     workspaces.  Host-built (build_candidate_tables); owning.
//
// Every element struct is trivially copyable; every extent is named.
// Cone positions stay OUT of both (Layer-I discipline: the tables are
// intrinsic; Layer II passes anchors alongside).
//
// The claim lookup (ParamTablesView::claim) is the form
// eisenstein_atlas reads cell charts through: a cell is a CONVEX
// lattice triangle, so each scanline's claimed points are one
// contiguous a-range and (row, a-offset) is a bijection onto the
// cell's entries -- the lookup is pure row arithmetic into the CSR
// (@ref paint-cell-convexity, enforced at flatten time).
//
// "Device-legal" here means: trivially copyable, no allocation and no
// exceptions on the READ path (the element structs and ParamTablesView).
// CandidateTables and its builder are host construction (owning
// vectors, throwing).  The accessors are UNCHECKED projections
// (@pre 0 <= f < nf, 0 <= v < N_sorted); claim is the only lookup
// whose domain includes misses.
// =====================================================================

#include "fullerenes/eisenstein.hh"
#include "fullerenes/delaunay_view.hh"   // DelaunayView (the builder's input)

#include <cstdint>
#include <span>
#include <type_traits>
#include <vector>

namespace eisenstein_paint {

// ---------------------------------------------------------------------
// Element structs (all trivially copyable).
// ---------------------------------------------------------------------

// The three CCW T_sorted corner ids of a cell; c0 < 0 marks a dead slot.
struct CellCorners { int32_t c0 = -1, c1 = -1, c2 = -1; };

// A cell's chart frame: corner positions P1, P2 with P0 = (0,0).
struct CellFrame { int32_t p1a = 0, p1b = 0, p2a = 0, p2b = 0; };

// One triangle scan: rows [rows_first, rows_first + (b_max - b_min + 1))
// in the rows table; n_entries = total lattice points over its rows.
// (Same empty conventions as the raster twins ScanLines/ScanLine.)
struct ScanBlock {
    int32_t b_min = 0, b_max = -1, rows_first = 0, n_entries = 0;
    bool empty() const { return b_max < b_min; }
};

// One scanline: inclusive a-range and its exclusive prefix offset into
// the owning block's entry range.
struct ScanRow {
    int32_t a_left = 1, a_right = 0, entry_off = 0;
    bool empty() const { return a_left > a_right; }
};

// One enumerated lattice point: F-frame position (a, b) and the
// T_sorted vertex id claimed there (the LatticeMap entry, flattened;
// same per-cell order: scanline-major).
struct LatticePoint { int32_t a = 0, b = 0, vid = -1; };

// One appearance of a T_sorted vertex in some cell's chart (the
// Occurrence of eisenstein_paint.hh, flattened).
struct VertexOccurrence { int32_t cell = -1, a = 0, b = 0; };

static_assert(std::is_trivially_copyable_v<CellCorners>);
static_assert(std::is_trivially_copyable_v<CellFrame>);
static_assert(std::is_trivially_copyable_v<ScanBlock>);
static_assert(std::is_trivially_copyable_v<ScanRow>);
static_assert(std::is_trivially_copyable_v<LatticePoint>);
static_assert(std::is_trivially_copyable_v<VertexOccurrence>);

// ---------------------------------------------------------------------
// ParamTablesView: the post-embed product as spans.
// ---------------------------------------------------------------------
struct ParamTablesView {
    int32_t nf = 0;            // cell slots (D face ids; dead slots included)
    int32_t N_sorted = 0;      // T_sorted vertex count
    int32_t N_orig = 0;        // input vertex count (= N_sorted)
    int32_t n_cones = 0;       // sorted labels < n_cones are cones
    int64_t entries_total = 0; // = entry_first[nf]

    std::span<const CellCorners>      cells;        // [nf]
    std::span<const CellFrame>        frames;       // [nf]   the WINNING chart
    std::span<const ScanBlock>        scans;        // [nf]
    std::span<const ScanRow>          rows;         // CSR row blocks
    std::span<const int32_t>          entry_first;  // [nf+1] exclusive prefix
    std::span<const LatticePoint>     entries;      // [entries_total]
    std::span<const int32_t>          perm;         // [N_orig] orig -> sorted
    std::span<const int32_t>          occ_first;    // [N_sorted+1] occurrence CSR
    std::span<const VertexOccurrence> occ;          // 1-2 per non-cone vertex

    bool cell_live(int f) const { return cells[f].c0 >= 0; }

    // The cell's chart corners as lattice points (P0 at the origin) and
    // its CCW corner ids -- the k-indexed forms every "arc k runs corner
    // k -> k+1" loop wants.
    std::array<Eisenstein, 3> frame_points(int f) const {
        return { Eisenstein(0, 0),
                 Eisenstein(frames[f].p1a, frames[f].p1b),
                 Eisenstein(frames[f].p2a, frames[f].p2b) };
    }
    std::array<int, 3> corner_ids(int f) const {
        return { cells[f].c0, cells[f].c1, cells[f].c2 };
    }

    // All lattice points of cell f, scanline-major.
    std::span<const LatticePoint> cell_entries(int f) const {
        return entries.subspan(entry_first[f],
                               entry_first[f + 1] - entry_first[f]);
    }

    // All chart appearances of T_sorted vertex v (1 or 2 for a non-cone
    // vertex on a coverage-checked parametrization).
    std::span<const VertexOccurrence> occurrences(int v) const {
        return occ.subspan(occ_first[v], occ_first[v + 1] - occ_first[v]);
    }

    // THE claim lookup: is lattice point p claimed by cell f?  Pure row
    // arithmetic (see banner); any bounds miss = unclaimed.  Returns
    // nullptr when unclaimed.
    const LatticePoint* claim(int f, Eisenstein p) const {
        if (f < 0 || f >= nf || !cell_live(f)) return nullptr;
        const ScanBlock& sb = scans[f];
        const int b = p.second;
        if (b < sb.b_min || b > sb.b_max) return nullptr;
        const ScanRow& row = rows[sb.rows_first + (b - sb.b_min)];
        const int a = p.first;
        if (a < row.a_left || a > row.a_right) return nullptr;
        const LatticePoint* e =
            &entries[entry_first[f] + row.entry_off + (a - row.a_left)];
        return e;
    }
};

static_assert(std::is_trivially_copyable_v<ParamTablesView>);

// ---------------------------------------------------------------------
// CandidateTables: the pre-embed prepare form (host-built, owning).
// Candidate k of cell f lives at index f*2 + k in frames/scans; a cell
// has n_cand[f] in {0, 1, 2} candidates (<= 2 by the split-prime
// argument), and the embed phase picks the accepted one -- after which
// the accepted slice has exactly ParamTablesView's shape.
// ---------------------------------------------------------------------
struct CandidateTables {
    int32_t nf = 0;
    int32_t N_sorted = 0;

    std::vector<CellCorners> cells;        // [nf]
    std::vector<int32_t>     n_cand;       // [nf]
    std::vector<CellFrame>   frames;       // [nf*2]
    std::vector<ScanBlock>   scans;        // [nf*2]
    std::vector<ScanRow>     rows;         // flattened row blocks

    // CAPACITY prefix, not the post-embed count: cell f's width is the
    // MAX entry count over its candidates, so the accepted slice may be
    // narrower than the slot (ScanBlock::n_entries is the truth there).
    // The name differs from ParamTablesView::entry_first deliberately --
    // cell_entries()-style indexing does not apply to candidate tables.
    std::vector<int32_t> entry_capacity_first;  // [nf+1]
    int64_t entries_capacity = 0;               // entry_capacity_first[nf]

    // Executor sizing: batch-domain widths + embed/enumerate scratch
    // capacities (properties of the DATA -- any executor needs them).
    int32_t rows_cap = 0;          // max row count over (cell, cand)
    int32_t max_cell_entries = 0;  // max per-cell entry width
    int64_t wcap = 0;              // walk-registration scratch (pow2)
    int64_t bcap = 0;              // boundary-map scratch (pow2)
};

// Per-edge candidate-development bound (loud CAPACITY on overflow) --
// a property of the data like the scratch capacities above; executors
// size their per-edge strip slots by it.
inline constexpr int max_edge_candidates = 8;

struct SortedDual;   // eisenstein_paint.hh

// Build the candidate tables of an intrinsically-Delaunay complex D
// over S: corner data, up to two candidate frames per live cell
// (sector-0 representatives of the squared corner distance, filtered
// by the two remaining side lengths and orientation), per-candidate
// triangle scans, and the scratch capacity formulas.  Host, allocating.
// Only S.T.N is read from S today; the SortedDualView facade is the
// eventual parameter.  Throws PaintError(EMBED) when a live cell
// admits no candidate frame, or more than two (the split-prime bound);
// a degenerate frame surfaces as scan_triangle's runtime_error.
CandidateTables build_candidate_tables(const ::DelaunayView& D,
                                       const SortedDual& S);

}  // namespace eisenstein_paint
