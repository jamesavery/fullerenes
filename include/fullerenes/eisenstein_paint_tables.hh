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
// never a fixed slot count (a bounded-storage executor slots cells
// at max_cell_developments and refuses LOUDLY on overflow, never
// truncates).  The accepted development -- the first
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
// projections) -- and, for the BUILDER vocabulary at the end of this
// header (cell_metric_total ... cell_developments_into), on the write
// path too: span bodies over caller storage, failure by sentinel.
// The owning CellDevelopments and cell_developments stay host
// construction (owning vectors, throwing) wrapping those bodies --
// ONE body each, the scan_triangle_into / append_scan_rows_into
// pattern.  The accessors are UNCHECKED projections (@pre 0 <= f < nf,
// 0 <= v < N_sorted); claim is the only lookup whose domain includes
// misses.
// =====================================================================

#include "fullerenes/eisenstein.hh"
#include "fullerenes/eisenstein_raster.hh"   // ScanLine, scan_triangle_into
#include "fullerenes/delaunay_view.hh"   // DelaunayView (the builder's input)
#include "fullerenes/delaunay_geometry.hh"   // lsq_integrality_band (the float->exact entry)

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <optional>
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

// append_scan_rows_into — the ONE row-table derivation over CALLER
// storage (device-legal: no allocation, no throw): append one triangle
// scan's rows (a ScanLine span from scan_triangle_into) at rows[nrows..]
// and return its block, advancing nrows.  Empty scanlines keep the
// scan's own (a_left > a_right) values, so every derivation produces
// identical bytes.  rows_first == -1 signals capacity (rows too small).
// The owning append_scan_rows (eisenstein_paint.cc) wraps this.
inline ScanBlock append_scan_rows_into(std::span<ScanRow> rows, long& nrows,
                                       std::span<const ScanLine> lines,
                                       int32_t b_min, int32_t b_max) {
    ScanBlock sb;
    sb.b_min = b_min;
    sb.b_max = b_max;
    sb.rows_first = (int32_t)nrows;
    const long n = (long)b_max - b_min + 1;
    if (n > 0 && nrows + n > (long)rows.size()) { sb.rows_first = -1; return sb; }
    int32_t running = 0;
    for (long i = 0; i < n; ++i) {
        rows[(std::size_t)nrows++] =
            ScanRow{lines[(std::size_t)i].a_left, lines[(std::size_t)i].a_right, running};
        if (!lines[(std::size_t)i].empty())
            running += lines[(std::size_t)i].a_right - lines[(std::size_t)i].a_left + 1;
    }
    sb.n_entries = running;
    return sb;
}

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

struct DevelopmentSizing;   // the builder's sizing block (defined below)

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

    // The owner's executor sizing in the builder's block form — THE one
    // owner -> DevelopmentSizing projection (byte gates and span
    // surfaces compare/carry this form).  Defined after the block.
    DevelopmentSizing sizing() const;
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

// Per-cell development cap for BOUNDED-STORAGE executors.  Unlike the
// strip cap above, the per-cell development COUNT is unbounded over
// the isomer space (the banner's formula; no fixed per-norm cap is
// sound), so the owning cell_developments never caps it -- it grows.
// An executor with fixed slabs slots each cell at this many
// developments and refuses LOUDLY (the builder's capacity sentinel)
// when exceeded, never truncating.  8 is calibrated headroom over the
// observed corpus maximum (C980's cells admit 4); the first norm past
// it is N01 = (7*13)^2 = 8281 with 10 sector-0 reps (~C165,000 in the
// GC family).
inline constexpr int max_cell_developments = 8;

// =====================================================================
// The builder vocabulary (device-legal: no allocation, no throw): the
// span/sentinel bodies the owning cell_developments wraps -- ONE body
// for the lib pipeline and any device executor, completing the
// scan_triangle_into / append_scan_rows_into pattern.
// =====================================================================

// One cell's iDT metric datum: CCW corner ids + the squared side norms
// -- the input of the chart-frame construction.  THE one derivation;
// embed_cell and cell_developments both open with it.
struct CellMetric {
    std::array<int, 3> corners{ -1, -1, -1 };
    long N01 = 0, N12 = 0, N20 = 0;
};

// The float->exact entry of the paint pipeline, total form: the iDT's
// post-flip edge lengths (doubles -- the Alexandrov homotopy's banded
// flips discard the exact Lsq carry) re-enter the integer regime here,
// through the SAME named integrality trust boundary every other such
// conversion uses (delaunay_detail::lsq_integrality_band; cf.
// Diamond::squared and the exact-reduction entry in delaunay.cc).
// nullopt = not an integer within the band (a non-Loeschian metric --
// without the guard such a value could mis-round to a NEIGHBOURING
// valid norm and chart the wrong cell silently), or non-finite /
// >= 2^52 (a broken metric, and past the exact-cast domain below).
// The owning checked_integer_norm (eisenstein_paint.cc) throws the
// per-edge diagnosis.
inline std::optional<long> checked_integer_norm_total(double L) {
    const double sq = L * L;
    if (!(sq >= 0.0 && sq < 4503599627370496.0))   // 2^52
        return std::nullopt;
    // (long)(sq + 0.5) == lround(sq) unconditionally on [0, 2^52): for
    // sq >= 1, 0.5 is a multiple of ulp(sq) so sq + 0.5 is exact except
    // in the exponent-bump region, where the floor still lands on the
    // bumped power; below 1 the two forms differ only ~0.5 away from an
    // integer, which the band rejects either way.  Spelled without
    // lround because device SSCP passes lack a builtin for it.
    const long N = (long)(sq + 0.5);
    if (std::fabs(sq - (double)N) >
        delaunay_detail::lsq_integrality_band * std::max(1.0, sq))
        return std::nullopt;
    return N;
}

// One cell's metric, total form: nullopt when any edge's squared length
// fails the trust boundary above (the owning cell_metric throws with
// the per-edge diagnosis).  @pre D.f_he[f] >= 0.
inline std::optional<CellMetric> cell_metric_total(const ::DelaunayView& D,
                                                   int f) {
    const int h0 = D.f_he[f];
    const int h1 = D.he_next[h0];
    const int h2 = D.he_next[h1];
    const auto N01 = checked_integer_norm_total(D.he_length[h0]);
    const auto N12 = checked_integer_norm_total(D.he_length[h1]);
    const auto N20 = checked_integer_norm_total(D.he_length[h2]);
    if (!N01 || !N12 || !N20) return std::nullopt;
    CellMetric m;
    m.corners = { D.he_origin[h0], D.he_origin[h1], D.he_origin[h2] };
    m.N01 = *N01; m.N12 = *N12; m.N20 = *N20;
    return m;
}

// for_each_development -- THE enumeration body (the banner carries THE
// count statement): visit every admissible development frame (P1, P2)
// of the metric (N01, N12, N20) in the canonical sector-0 order, P0 at
// the origin.  Exact in Z[w]: for each sector-0 base P1 the apex is
// the UNIQUE lattice point at squared distances N20 from P0 and N12
// from P1 on the CCW side (place_third_eis_total; nullopt = this base
// admits no lattice apex).  visit returns false to stop early (a
// bounded executor's capacity trip).  The owning
// enumerate_developments (eisenstein_paint.cc) wraps this.
template <class Visit>
inline void for_each_development(long N01, long N12, long N20,
                                 Visit&& visit) {
    const Eisenstein P0(0, 0);
    for (Eisenstein P1 : Sector0Reps((int)N01)) {
        const auto P2 =
            place_third_eis_total(P0, P1, (int)N01, (int)N20, (int)N12, +1);
        if (!P2) continue;
        if (!visit(P1, *P2)) return;
    }
}

// Scratch capacity formulas (properties of the data; any executor
// sizing embed/enumerate scratch reads these):
//  * walk registration: distinct keys are T_sorted vertex ids, and a
//    walk registers <= 3 vertices per pushed face with <=
//    walk_max_steps faces per primitive sub-walk, so
//    min(N_sorted, 3*(walk_max_steps+2)) bounds the live count; x2
//    keeps the open-addressing load factor comfortable.
//  * boundary map: a qualifying walk to displacement s = l1_norm
//    pushes <= 3s faces, so <= 3*(3s+1)+3 registrations per edge,
//    three edges; x2 headroom.  s_max = max l1_norm over development
//    frame edges.
struct ScratchCapacities { int64_t wcap = 0, bcap = 0; };

inline ScratchCapacities scratch_capacities(int64_t N_sorted, long s_max) {
    // Smallest power of two >= x, floored at 64 (open-addressing
    // scratch wants pow2 extents with load-factor headroom).
    const auto pow2_ceil64 = [](int64_t x) {
        int64_t p = 64;
        while (p < x) p <<= 1;
        return p;
    };
    ScratchCapacities c;
    c.wcap = pow2_ceil64(
        2 * (std::min<int64_t>(N_sorted, 3 * (int64_t)(walk_max_steps + 2)) + 8));
    c.bcap = pow2_ceil64(2 * (3 * (3 * s_max + 4) * 3 + 16));
    return c;
}

// The builder's sizing block -- the span form of the owner's fenced
// executor-sizing fields (defaulted == so both legs of a byte gate
// compare it as one value).
struct DevelopmentSizing {
    int32_t n_frames = 0;          // = dev_first[nf]
    int32_t rows_cap = 0;          // max row count over (cell, development)
    int32_t max_cell_entries = 0;  // max per-cell entry width
    int64_t n_rows = 0;            // rows written (= the owner's rows.size())
    int64_t entries_capacity = 0;  // = entry_capacity_first[nf]
    int64_t wcap = 0, bcap = 0;    // scratch_capacities of the data
    friend bool operator==(const DevelopmentSizing&,
                           const DevelopmentSizing&) = default;
};

static_assert(std::is_trivially_copyable_v<DevelopmentSizing>);

inline DevelopmentSizing CellDevelopments::sizing() const {
    DevelopmentSizing s;
    s.n_frames         = (int32_t)frames.size();
    s.rows_cap         = rows_cap;
    s.max_cell_entries = max_cell_entries;
    s.n_rows           = (int64_t)rows.size();
    s.entries_capacity = entries_capacity;
    s.wcap             = wcap;
    s.bcap             = bcap;
    return s;
}

// CellDevelopmentsOut: CellDevelopments over CALLER storage -- the
// builder's out-form.  The identity triple + the mutable table slabs +
// the sizing block; view() projects the filled prefix as the
// device-legal CellDevelopmentsView.  fail_cell names the refusing
// cell on a degenerate (-1) return (the owner's diagnosis handle).
struct CellDevelopmentsOut {
    int32_t nf = 0;            // written by the builder (= D.nf)
    int32_t N_sorted = 0;      // written by the builder (caller identity)
    int32_t n_cones = 0;       // written by the builder (caller identity)

    std::span<CellCorners> cells;                 // [>= nf]
    std::span<int32_t>     dev_first;             // [>= nf+1]
    std::span<CellFrame>   frames;                // executor-sized (max_cell_developments slots)
    std::span<ScanBlock>   scans;                 // same extent as frames
    std::span<ScanRow>     rows;                  // executor-sized
    std::span<int32_t>     entry_capacity_first;  // [>= nf+1]

    DevelopmentSizing sizing;
    int32_t fail_cell = -1;

    CellDevelopmentsView view() const {
        CellDevelopmentsView v;
        v.nf = nf; v.N_sorted = N_sorted; v.n_cones = n_cones;
        v.cells     = cells.first((std::size_t)nf);
        v.dev_first = dev_first.first((std::size_t)nf + 1);
        v.frames    = frames.first((std::size_t)sizing.n_frames);
        v.scans     = scans.first((std::size_t)sizing.n_frames);
        v.rows      = rows.first((std::size_t)sizing.n_rows);
        return v;
    }
};

// cell_developments_into -- THE builder body over CALLER storage
// (device-legal: no allocation, no throw): every admissible lattice
// development of every live cell of D, CSR by cell, in the canonical
// enumeration order, with per-development triangle scans and the
// scratch capacities.  The owning cell_developments below wraps it --
// ONE body.
//
// @pre  D's live faces are triangles (f_he[f] < 0 marks a dead slot);
//       N_sorted / n_cones are the sorted dual's identity (stored on
//       out, never read by the builder).
// @pre  lines is the per-development scanline scratch (its size bounds
//       the rows of ONE development; scan_triangle_into refuses past
//       it).
// @error -2 (capacity) when cells/dev_first/entry_capacity_first are
//       smaller than nf/nf+1/nf+1, or frames/scans/rows/lines are
//       exhausted -- grow and re-run (the owner grows; a bounded
//       executor refuses loudly).
// @error -1 (degenerate, out.fail_cell = the cell) when a live cell's
//       metric fails the integrality boundary, admits no lattice
//       development, or only a zero-area one (lattice tau == 0, where
//       scan_triangle_into reports the collinear frame) -- the owner's
//       PaintError(EMBED).
// @post on success (return = n_frames >= 0): dev_first is the
//       development CSR with dev_first[nf] == n_frames;
//       entry_capacity_first is the exclusive prefix of per-cell max
//       entry widths; sizing is filled; the filled prefixes equal the
//       owning cell_developments' tables byte-for-byte.  On any
//       negative return the tables are partial: read only fail_cell.
inline long cell_developments_into(const ::DelaunayView& D,
                                   int32_t N_sorted, int32_t n_cones,
                                   CellDevelopmentsOut& out,
                                   std::span<ScanLine> lines) {
    const int nf = D.nf;
    if ((long)out.cells.size() < nf ||
        (long)out.dev_first.size() < nf + 1 ||
        (long)out.entry_capacity_first.size() < nf + 1)
        return -2;
    out.nf = nf;
    out.N_sorted = N_sorted;
    out.n_cones = n_cones;
    out.sizing = DevelopmentSizing{};
    out.fail_cell = -1;
    out.entry_capacity_first[0] = 0;
    long n_rows = 0;     // append_scan_rows_into's running row count
    long s_max = 0;      // max l1_norm over development frame edges
    const Eisenstein P0(0, 0);

    for (int f = 0; f < nf; ++f) {
        out.dev_first[f] = out.sizing.n_frames;
        out.entry_capacity_first[f + 1] = 0;
        out.cells[f] = CellCorners{};
        if (D.f_he[f] < 0) continue;    // dead slot: 0 developments

        const auto m = cell_metric_total(D, f);
        if (!m) { out.fail_cell = f; return -1; }
        out.cells[f] = CellCorners{ m->corners[0], m->corners[1],
                                    m->corners[2] };

        // Per-development frame + triangle scan (scan_triangle_into +
        // append_scan_rows_into: the same derivations the
        // post-selection flatten uses).
        int32_t n_dev = 0, cap_f = 0;
        long rc = 0;
        for_each_development(m->N01, m->N12, m->N20,
                             [&](Eisenstein P1, Eisenstein P2) -> bool {
            const long slot = out.sizing.n_frames;
            if (slot >= (long)out.frames.size() ||
                slot >= (long)out.scans.size()) { rc = -2; return false; }
            out.frames[slot] = CellFrame{ P1.first, P1.second,
                                          P2.first, P2.second };
            for (Eisenstein e : { P1, P2 - P1, P0 - P2 }) {
                const long sd = e.l1_norm();
                if (sd > s_max) s_max = sd;
            }
            int b_min = 0, b_max = -1;
            const long nl = scan_triangle_into(P0, P1, P2, lines,
                                               b_min, b_max);
            if (nl == -1) { rc = -1; return false; }   // zero-area frame
            if (nl == -2) { rc = -2; return false; }   // lines too small
            const ScanBlock sb = append_scan_rows_into(
                out.rows, n_rows, lines.first((std::size_t)nl),
                (int32_t)b_min, (int32_t)b_max);
            if (sb.rows_first < 0) { rc = -2; return false; }
            out.scans[slot] = sb;
            const int32_t nr = sb.n_rows();
            if (nr > out.sizing.rows_cap) out.sizing.rows_cap = nr;
            if (sb.n_entries > cap_f) cap_f = sb.n_entries;
            ++out.sizing.n_frames;
            ++n_dev;
            return true;
        });
        if (rc != 0) { if (rc == -1) out.fail_cell = f; return rc; }
        if (n_dev == 0) { out.fail_cell = f; return -1; }
        if (cap_f > out.sizing.max_cell_entries)
            out.sizing.max_cell_entries = cap_f;
        out.entry_capacity_first[f + 1] = cap_f;
    }
    out.dev_first[nf] = out.sizing.n_frames;
    for (int f = 0; f < nf; ++f)    // exclusive prefix
        out.entry_capacity_first[f + 1] += out.entry_capacity_first[f];
    out.sizing.n_rows = n_rows;
    out.sizing.entries_capacity = out.entry_capacity_first[nf];
    const ScratchCapacities sc = scratch_capacities(N_sorted, s_max);
    out.sizing.wcap = sc.wcap;
    out.sizing.bcap = sc.bcap;
    return out.sizing.n_frames;
}

struct SortedDual;   // eisenstein_paint.hh

// Every admissible lattice development of every live cell of an
// intrinsically-Delaunay complex D over S: corner data, frames in the
// canonical enumeration order, per-development triangle scans, and the
// scratch capacity formulas.  Host, allocating; owning wrapper of
// cell_developments_into -- one body (it grows the slabs on the
// capacity sentinel, so it never refuses on capacity).  Only S.T.N and
// S.n_cones are read from S today; the SortedDualView facade is the
// eventual parameter.  Throws PaintError(EMBED) on the degenerate
// sentinel: a live cell whose metric fails the integrality boundary
// (the per-edge diagnosis) or admits no non-degenerate lattice
// development.
CellDevelopments cell_developments(const ::DelaunayView& D,
                                   const SortedDual& S);

// The T-free form -- THE body, which the SortedDual overload forwards
// to.  N_sorted / n_cones are carried on the product and enter only the
// scratch-capacity formulas, never the developments themselves, so a
// pure flat cone metric (where every D vertex IS a cone, and no dual
// triangulation is known) passes its own vertex count for both.  Same
// failure contract.
CellDevelopments cell_developments(const ::DelaunayView& D,
                                   int32_t N_sorted, int32_t n_cones);

}  // namespace eisenstein_paint
