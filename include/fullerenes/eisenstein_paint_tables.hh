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
//   ParamTablesView  -- the POST-SELECTION product: one chart per live
//     cell (the accepted development), entries filled, occurrence CSR
//     built.  This is what SurfaceParametrization (eisenstein_paint.hh)
//     stores and what the atlas consumes.  Non-owning, trivially
//     copyable, device-legal.  Per-cell access goes through chart(f)
//     -- the ChartView projection -- so algorithms speak in charts,
//     not in table offsets.
//
//   CellDevelopments -- the PRE-SELECTION form: every admissible
//     lattice development of every live cell, in the canonical
//     enumeration order, with per-development triangle scans and the
//     scratch capacities an executor needs.  Host-built
//     (cell_developments); owning.  Per-cell access goes through
//     n_developments(f) / development(f, k) -- the DevelopmentView
//     projection.
//
// A LATTICE DEVELOPMENT of a cell is an isometric placement of its
// geodesic triangle into Z[w] in the canonical gauge (corner 0 at the
// origin, corner 1 a sector-0 representative of the base norm N01).
// THE COUNT (the one full statement; other sites reference here): the
// apex direction tau = P2/P1 is fixed by the cell metric, so a
// sector-0 representative z realizes a development iff tau*z is a
// lattice point, i.e. iff delta | z with delta = P1/gcd(P1, P2)
// (metric-determined up to units).  Hence
//   #developments = #ideals of Z[w] of norm N(gcd(P1, P2))
//                   (+1 exactly when delta | sqrt(N01)),
// which is at most #ideals(N01) + [N01 a perfect square]: typically
// 1 or 2, but UNBOUNDED over the isomer space -- a GC(7,0)-scaled
// cell (side norm 49) admits 4 -- so storage is CSR over cells,
// never a fixed slot count.  The accepted development -- the first
// whose OWN three boundary strips glue into one consistent
// position -> vertex map; no neighbour data enters the selection --
// becomes the cell's CHART, and it is unique (the cell-development
// lemma at first_consistent_triple, eisenstein_paint.cc).
//
// Every element struct is trivially copyable; every extent is named.
// Cone positions stay OUT of both (Layer-I discipline: the tables are
// intrinsic; Layer II passes anchors alongside).
//
// The claim lookup (ChartView::claim) is the form eisenstein_atlas
// reads cell charts through: a cell is a CONVEX lattice triangle, so
// each scanline's claimed points are one contiguous a-range and
// (row, a-offset) is a bijection onto the chart's entries -- the
// lookup is pure row arithmetic into the CSR
// (@ref paint-cell-convexity, enforced at flatten time).
//
// "Device-legal" here means: trivially copyable, no allocation and no
// exceptions on the READ path (the element structs, the views, the
// projections).  CellDevelopments and its builder are host
// construction (owning vectors, throwing).  The accessors are
// UNCHECKED projections (@pre 0 <= f < nf, 0 <= v < N_sorted); claim
// is the only lookup whose domain includes misses.
// =====================================================================

#include "fullerenes/eisenstein.hh"
#include "fullerenes/delaunay_view.hh"   // DelaunayView (the builder's input)

#include <array>
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

// A development's frame: corner positions P1, P2 with P0 = (0,0).
struct CellFrame {
    int32_t p1a = 0, p1b = 0, p2a = 0, p2b = 0;
    // The gauge frame as corner points (P0 at the origin) -- the one
    // body behind every frame_points() accessor.
    std::array<Eisenstein, 3> points() const {
        return { Eisenstein(0, 0), Eisenstein(p1a, p1b),
                 Eisenstein(p2a, p2b) };
    }
};

// One triangle scan: rows [rows_first, rows_first + (b_max - b_min + 1))
// in the rows table; n_entries = total lattice points over its rows.
// (Same empty conventions as the raster twins ScanLines/ScanLine.)
struct ScanBlock {
    int32_t b_min = 0, b_max = -1, rows_first = 0, n_entries = 0;
    bool    empty() const { return b_max < b_min; }
    int32_t n_rows() const { return empty() ? 0 : b_max - b_min + 1; }
};

// One scanline: inclusive a-range and its exclusive prefix offset into
// the owning block's entry range.
struct ScanRow {
    int32_t a_left = 1, a_right = 0, entry_off = 0;
    bool empty() const { return a_left > a_right; }
};

// One enumerated lattice point: chart position (a, b) and the T_sorted
// vertex id claimed there (scanline-major per cell).
struct LatticePoint {
    int32_t a = 0, b = 0, vid = -1;
    Eisenstein pos() const { return Eisenstein(a, b); }
};

// One appearance of a T_sorted vertex in some cell's chart: the host
// cell and the chart position (a, b).
struct VertexOccurrence {
    int32_t cell = -1, a = 0, b = 0;
    Eisenstein pos() const { return Eisenstein(a, b); }
};

static_assert(std::is_trivially_copyable_v<CellCorners>);
static_assert(std::is_trivially_copyable_v<CellFrame>);
static_assert(std::is_trivially_copyable_v<ScanBlock>);
static_assert(std::is_trivially_copyable_v<ScanRow>);
static_assert(std::is_trivially_copyable_v<LatticePoint>);
static_assert(std::is_trivially_copyable_v<VertexOccurrence>);

// ---------------------------------------------------------------------
// ChartView: the chart of ONE cell -- its CCW corner ids, the
// development phi on the corners, and the claim map lambda presented
// as the claimed lattice points plus the row-arithmetic lookup.
// ---------------------------------------------------------------------
struct ChartView {
    std::array<int32_t, 3>    corners{ -1, -1, -1 };
    CellFrame                 frame;
    ScanBlock                 scan;
    std::span<const ScanRow>      rows;      // the scan's own rows
    std::span<const LatticePoint> entries;   // L_f, scanline-major

    bool live() const { return corners[0] >= 0; }

    // phi on the corners (P0 at the origin).
    std::array<Eisenstein, 3> frame_points() const { return frame.points(); }

    // lambda at p: the claimed lattice point, or nullptr when p is not
    // claimed.  Pure row arithmetic (see the banner); any bounds miss
    // = unclaimed.
    const LatticePoint* claim(Eisenstein p) const {
        if (!live()) return nullptr;
        const int b = p.second;
        if (b < scan.b_min || b > scan.b_max) return nullptr;
        const ScanRow& row = rows[b - scan.b_min];
        const int a = p.first;
        if (a < row.a_left || a > row.a_right) return nullptr;
        return &entries[row.entry_off + (a - row.a_left)];
    }
};

static_assert(std::is_trivially_copyable_v<ChartView>);

// ---------------------------------------------------------------------
// ParamTablesView: the post-selection product as spans.
// ---------------------------------------------------------------------
struct ParamTablesView {
    int32_t nf = 0;            // cell slots (D face ids; dead slots included)
    int32_t N_sorted = 0;      // T_sorted vertex count
    int32_t N_orig = 0;        // input vertex count (= N_sorted)
    int32_t n_cones = 0;       // sorted labels < n_cones are cones
    int64_t entries_total = 0; // = entry_first[nf]

    std::span<const CellCorners>      cells;        // [nf]
    std::span<const CellFrame>        frames;       // [nf]   the accepted chart
    std::span<const ScanBlock>        scans;        // [nf]
    std::span<const ScanRow>          rows;         // CSR row blocks
    std::span<const int32_t>          entry_first;  // [nf+1] exclusive prefix
    std::span<const LatticePoint>     entries;      // [entries_total]
    std::span<const int32_t>          perm;         // [N_orig] orig -> sorted
    std::span<const int32_t>          occ_first;    // [N_sorted+1] occurrence CSR
    std::span<const VertexOccurrence> occ;          // 1-2 per non-cone vertex

    bool cell_live(int f) const { return cells[f].c0 >= 0; }

    // THE per-cell projection: cell f's chart.  Algorithms read charts,
    // never table offsets.
    ChartView chart(int f) const {
        ChartView c;
        c.corners = { cells[f].c0, cells[f].c1, cells[f].c2 };
        c.frame   = frames[f];
        c.scan    = scans[f];
        if (c.scan.n_rows() > 0)
            c.rows = rows.subspan(c.scan.rows_first, c.scan.n_rows());
        c.entries = entries.subspan(entry_first[f],
                                    entry_first[f + 1] - entry_first[f]);
        return c;
    }

    // The k-indexed corner/position forms every "arc k runs corner
    // k -> k+1" loop wants (thin forwards of chart(f)).
    std::array<Eisenstein, 3> frame_points(int f) const {
        return chart(f).frame_points();
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

    // lambda_f at p (thin forward of chart(f).claim; nullptr when f is
    // out of range, dead, or p unclaimed).
    const LatticePoint* claim(int f, Eisenstein p) const {
        if (f < 0 || f >= nf || !cell_live(f)) return nullptr;
        return chart(f).claim(p);
    }
};

static_assert(std::is_trivially_copyable_v<ParamTablesView>);

// ---------------------------------------------------------------------
// DevelopmentView: one lattice development of one cell -- the frame
// phi on the corners and the scanline fibration of phi(F) n Z[w].
// ---------------------------------------------------------------------
struct DevelopmentView {
    CellFrame                frame;
    ScanBlock                scan;
    std::span<const ScanRow> rows;   // the scan's own rows

    std::array<Eisenstein, 3> frame_points() const { return frame.points(); }
    int32_t n_lattice_points() const { return scan.n_entries; }
    ScanRow row(int r) const { return rows[r]; }
};

static_assert(std::is_trivially_copyable_v<DevelopmentView>);

// ---------------------------------------------------------------------
// CellDevelopmentsView: the span mirror of CellDevelopments (below) --
// trivially copyable, device-legal; the capture unit of any executor
// that selects charts.  Same projections as the owner.
// ---------------------------------------------------------------------
struct CellDevelopmentsView {
    int32_t nf = 0;
    int32_t N_sorted = 0;
    int32_t n_cones = 0;       // sorted labels < n_cones are cones

    std::span<const CellCorners> cells;       // [nf]
    std::span<const int32_t>     dev_first;   // [nf+1] CSR prefix
    std::span<const CellFrame>   frames;      // [dev_first[nf]]
    std::span<const ScanBlock>   scans;       // [dev_first[nf]]
    std::span<const ScanRow>     rows;        // flattened row blocks

    int n_developments(int f) const { return dev_first[f + 1] - dev_first[f]; }
    CellCorners corners(int f) const { return cells[f]; }
    DevelopmentView development(int f, int k) const {
        DevelopmentView d;
        const int i = dev_first[f] + k;
        d.frame = frames[i];
        d.scan  = scans[i];
        if (d.scan.n_rows() > 0)
            d.rows = rows.subspan(d.scan.rows_first, d.scan.n_rows());
        return d;
    }
};

static_assert(std::is_trivially_copyable_v<CellDevelopmentsView>);

// ---------------------------------------------------------------------
// CellDevelopments: every admissible lattice development of every live
// cell, in the canonical enumeration order (host-built, owning).
// Storage is CSR over cells (dev_first); per-cell access goes through
// n_developments(f) / development(f, k).  The accepted development's
// slice has exactly ParamTablesView's per-cell shape.
// ---------------------------------------------------------------------
struct CellDevelopments {
    int32_t nf = 0;
    int32_t N_sorted = 0;
    int32_t n_cones = 0;       // sorted labels < n_cones are cones

    std::vector<CellCorners> cells;       // [nf]
    std::vector<int32_t>     dev_first;   // [nf+1] CSR prefix over developments
    std::vector<CellFrame>   frames;      // [dev_first[nf]]
    std::vector<ScanBlock>   scans;       // [dev_first[nf]]
    std::vector<ScanRow>     rows;        // flattened row blocks

    CellDevelopmentsView view() const {
        CellDevelopmentsView v;
        v.nf = nf; v.N_sorted = N_sorted; v.n_cones = n_cones;
        v.cells = cells; v.dev_first = dev_first;
        v.frames = frames; v.scans = scans; v.rows = rows;
        return v;
    }
    int n_developments(int f) const { return view().n_developments(f); }
    CellCorners corners(int f) const { return cells[f]; }
    DevelopmentView development(int f, int k) const { return view().development(f, k); }

    // --- executor sizing (properties of the DATA, not of the ---------
    // --- mathematics; fenced here so the projections above stay ------
    // --- purely geometric) -------------------------------------------

    // CAPACITY prefix, not a post-selection count: cell f's width is
    // the MAX entry count over its developments, so the accepted slice
    // may be narrower than the slot (ScanBlock::n_entries is the truth
    // there).  The name differs from ParamTablesView::entry_first
    // deliberately -- cell_entries()-style indexing does not apply
    // pre-selection.
    std::vector<int32_t> entry_capacity_first;  // [nf+1]
    int64_t entries_capacity = 0;               // entry_capacity_first[nf]

    // The ONE sanctioned use of the capacity prefix as an ADDRESS:
    // an executor's entries arena places cell f's slots at
    // [entry_base(f), entry_base(f + 1)).  Containment holds because
    // f's width is the max over its developments, so the accepted
    // development's n_entries never crosses into cell f+1's block.
    int32_t entry_base(int f) const { return entry_capacity_first[f]; }
    // Device-capturable span form of entry_base; same contract.
    std::span<const int32_t> entry_bases() const { return entry_capacity_first; }

    int32_t rows_cap = 0;          // max row count over (cell, development)
    int32_t max_cell_entries = 0;  // max per-cell entry width
    int64_t wcap = 0;              // walk-registration scratch (pow2)
    int64_t bcap = 0;              // boundary-map scratch (pow2)
};

// Per-edge strip-development bound (loud CAPACITY on overflow) -- a
// property of the data like the scratch capacities above; executors
// size their per-edge strip slots by it.  This bounds the WALKER'S
// admissible developments of one boundary strip, not a cell's frame
// count (which is CSR, above) -- and unlike that count it does NOT
// grow with the norm: the walker only ever tries the raw displacement
// and its sector-0 ROTATION representative (<= 2 endpoints), from
// <= 6 arcs x 6 directions, so distinct strip developments are
// hard-bounded by 72 for EVERY norm.  8 is CALIBRATED headroom over
// the observed maximum (1 for the common strip, 2 for split-prime /
// thin-cell ambiguities); exceeding it is a loud capacity trip,
// never truncation.
inline constexpr int max_strip_developments = 8;

struct SortedDual;   // eisenstein_paint.hh

// Every admissible lattice development of every live cell of an
// intrinsically-Delaunay complex D over S: corner data, frames in the
// canonical enumeration order, per-development triangle scans, and the
// scratch capacity formulas.  Host, allocating.  Only S.T.N and
// S.n_cones are read from S today; the SortedDualView facade is the
// eventual parameter.  Throws PaintError(EMBED) when a live cell
// admits NO lattice development (its metric is not realizable -- the
// loud refusal).
CellDevelopments cell_developments(const ::DelaunayView& D,
                                   const SortedDual& S);

}  // namespace eisenstein_paint
