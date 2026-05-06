#include "fullerenes/eisenstein_paint.hh"

#include "fullerenes/eisenstein.hh"
#include "fullerenes/eisenstein_raster.hh"
#include "fullerenes/eisenstein_tikz.hh"
#include "fullerenes/delaunay_strip.hh"
#include "fullerenes/delaunay_alexandrov12.hh"
#include "fullerenes/barycentric.hh"

#include <algorithm>
#include <array>
#include <cassert>
#include <climits>
#include <cmath>
#include <cstdarg>
#include <cstdio>
#include <cstdlib>
#include <optional>
#include <ostream>
#include <set>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

#ifdef _OPENMP
#include <omp.h>
#endif

namespace eisenstein_paint {

namespace {

// =====================================================================
// Internal: printf-style throw helper.
// =====================================================================

[[noreturn]] void paint_throw(const char* fmt, ...) {
    char buf[512];
    va_list ap; va_start(ap, fmt);
    std::vsnprintf(buf, sizeof buf, fmt, ap);
    va_end(ap);
    throw std::runtime_error(buf);
}

}  // namespace

// =====================================================================
// Result::stage_name
// =====================================================================

const char* Result::stage_name() const {
    switch (stage) {
        case Stage::OK:              return "ok";
        case Stage::IDT_COMPUTE:     return "idt_compute";
        case Stage::ALEXANDROV:      return "alexandrov";
        case Stage::NON_SIMPLICIAL:  return "non_simplicial";
        case Stage::EMBED:           return "embed";
        case Stage::COVERAGE:        return "coverage";
        case Stage::INTERPOLATE:     return "interpolate";
        case Stage::UNEXPECTED:      return "unexpected";
    }
    return "unknown";
}

// =====================================================================
// Prelude: sort, iDT compute, Alexandrov solve.
// =====================================================================

namespace {

void require_simplicial(const DelaunayTriangulation& D, const char* who) {
    if (!D.is_simplicial())
        throw StageError(Stage::NON_SIMPLICIAL,
            std::string(who) + ": iDT is not simplicial "
            "(self-loops or multi-edges present)");
    if (!D.is_well_formed())
        throw StageError(Stage::NON_SIMPLICIAL,
            std::string(who) + ": iDT is not well-formed "
            "(some live half-edge is not in a length-3 face cycle)");
}

void sort_and_compute_idt(const Triangulation& T, Pipeline& P) {
    P.perm     = T.sort_flat_last();
    P.T_sorted = T;
    P.T_sorted.apply_permutation(P.perm);
    try {
        P.D = DelaunayTriangulation::compute(P.T_sorted);
    } catch (const std::exception& e) {
        throw StageError(Stage::IDT_COMPUTE,
            std::string("eisenstein_paint::prepare: iDT compute: ") + e.what());
    }
}

void run_alexandrov(Pipeline& P) {
    AlexandrovSolver A;
    A.D              = std::move(P.D);
    P.cone_positions = A.solve();
    P.D              = std::move(A.D);
    if (A.valid()) return;

    using VS = AlexandrovSolver::ValidationStatus;
    const std::string detail = std::string("AlexandrovSolver: ")
                             + AlexandrovSolver::status_str(A.stats_status);
    const Stage st = (A.stats_status == VS::FAIL_NOT_SIMPLE)
                       ? Stage::NON_SIMPLICIAL
                       : Stage::ALEXANDROV;
    throw StageError(st, "eisenstein_paint::prepare: " + detail);
}

}  // namespace

Pipeline prepare_iDT(const Triangulation& T) {
    Pipeline P;
    sort_and_compute_idt(T, P);
    require_simplicial(P.D, "eisenstein_paint::prepare_iDT");
    return P;
}

Pipeline prepare(const Triangulation& T) {
    Pipeline P;
    sort_and_compute_idt(T, P);
    run_alexandrov(P);
    return P;
}

Inputs prepare_inputs(const Triangulation& T) {
    Pipeline P = prepare(T);
    return { std::move(P.T_sorted), std::move(P.D) };
}

// =====================================================================
// embed_cell: place one iDT cell into its own Eisenstein frame.
// =====================================================================

namespace {

// Lattice direction (0..5) of a unit Eisenstein vector, -1 otherwise.
int direction_of_unit(Eisenstein z) {
    for (int d = 0; d < 6; ++d) if (z == unit_direction(d)) return d;
    return -1;
}

// Derive k0 (frame offset) at vertex w given a neighbour `nbr` of w
// whose position relative to w is along lattice direction `d`.
int derive_k0(const Triangulation& T, int w, int nbr, int d) {
    const auto& nbrs = T[w];
    const int deg = (int)nbrs.size();
    for (int j = 0; j < deg; ++j)
        if (nbrs[j] == nbr) return ((d - j) % 6 + 6) % 6;
    return -1;
}

// Register every vertex in `wr.face_path` into `out` (first registration
// wins so we keep the earliest position/k0 on the walk).  Derives each
// vertex's k0 from one of its face-neighbours via its unit-direction
// displacement.
bool extract_walk_vertices(const Triangulation& T,
                           const WalkResult& wr,
                           std::unordered_map<int, std::pair<Eisenstein, int>>& out)
{
    if (wr.face_path.size() != wr.face_pos.size()) return false;
    for (size_t fi = 0; fi < wr.face_path.size(); ++fi) {
        const auto& F = wr.face_path[fi];
        const auto& P = wr.face_pos[fi];
        for (int k = 0; k < 3; ++k) {
            const int a = F[k];
            const int b = F[(k + 1) % 3];
            const Eisenstein dvec = P[(k + 1) % 3] - P[k];
            const int d = direction_of_unit(dvec);
            if (d < 0) return false;
            const int k0 = derive_k0(T, a, b, d);
            if (k0 < 0) return false;
            if (out.find(a) == out.end()) out[a] = {P[k], k0};
        }
    }
    return true;
}

// Walker output, indexed by T_sorted vertex id, in the START vertex's
// LOCAL frame (start_vertex sits at (0, 0)).
struct WalkVertices {
    std::unordered_map<int, std::pair<Eisenstein, int>> by_id;
};

// A cone vertex (id < 12) at any STRICT interior step of wr.walk means
// the directed line passes through a cone -- and at a cone the
// angular structure is undefined, so a geodesic cannot pass through
// one.  Trivially false when gcd(target_rel) == 1.
bool walks_through_cone(const WalkResult& wr) {
    for (size_t k = 1; k + 1 < wr.walk.size(); ++k)
        if (wr.walk[k].second < 12) return true;
    return false;
}

// Try every (face_arc i, dir_uv) start combination.  Accept the first
// chain whose walker terminates at target_vertex with a clean walk
// (start at (0,0), end at target_rel, no cone on strict interior).
std::optional<WalkVertices>
find_chain(const Triangulation& T_sorted,
           int start_vertex, int target_vertex, Eisenstein target_rel)
{
    const int deg = (int)T_sorted[start_vertex].size();
    for (int i = 0; i < deg; ++i) {
        const int v_a = T_sorted[start_vertex][i];
        const int v_b = T_sorted[start_vertex][(i + 1) % deg];
        for (int dir_uv = 0; dir_uv < 6; ++dir_uv) {
            WalkResult wr = walk_line(T_sorted, start_vertex, v_a, v_b,
                                      dir_uv, target_rel);
            if (wr.final_vertex != target_vertex) continue;

            WalkVertices V;
            if (!extract_walk_vertices(T_sorted, wr, V.by_id)) continue;
            const auto itu = V.by_id.find(start_vertex);
            const auto itv = V.by_id.find(target_vertex);
            if (itu == V.by_id.end() || itu->second.first != Eisenstein(0, 0)) continue;
            if (itv == V.by_id.end() || itv->second.first != target_rel)       continue;
            if (walks_through_cone(wr)) continue;

            return V;
        }
    }
    return std::nullopt;
}

// Translate verts to F's frame (add start_pos), clip to F's lattice
// triangle, bucket by scanline (b coordinate), write into `out`.
// Returns false if no vertex survives the clip.
bool clip_to_meta_by_scanline(const WalkVertices& V,
                              Eisenstein start_pos,
                              Eisenstein FP0, Eisenstein FP1, Eisenstein FP2,
                              EdgeStrip& out)
{
    const auto inside_meta = [&](Eisenstein p) {
        return wedge(FP1 - FP0, p - FP0) >= 0
            && wedge(FP2 - FP1, p - FP1) >= 0
            && wedge(FP0 - FP2, p - FP2) >= 0;
    };

    int b_min = INT_MAX, b_max = INT_MIN;
    std::vector<StripVertex> kept;
    kept.reserve(V.by_id.size());
    for (const auto& [vid, pk] : V.by_id) {
        const Eisenstein pos_F = pk.first + start_pos;
        if (!inside_meta(pos_F)) continue;
        StripVertex sv;
        sv.position     = pos_F;
        sv.vertex_id    = vid;
        sv.frame_offset = pk.second;        // k0 is translation-invariant
        if (pos_F.second < b_min) b_min = pos_F.second;
        if (pos_F.second > b_max) b_max = pos_F.second;
        kept.push_back(sv);
    }
    if (kept.empty()) return false;

    out.b_min = b_min;
    out.b_max = b_max;
    out.by_scanline.assign(b_max - b_min + 1, {});
    for (const auto& sv : kept)
        out.by_scanline[sv.position.second - b_min].push_back(sv);
    return true;
}

// Run the walker from start_vertex at start_pos in F's frame, aiming
// at target_vertex at target_pos.  On success, fills `out` with the
// F-frame strip (clipped to FP0/FP1/FP2 lattice triangle) and returns
// true.  Returns false if no chain qualifies or no vertex survives.
bool try_frame_walker(const Triangulation& T_sorted,
                      int start_vertex, Eisenstein start_pos,
                      int target_vertex, Eisenstein target_pos,
                      Eisenstein FP0, Eisenstein FP1, Eisenstein FP2,
                      EdgeStrip& out)
{
    const auto chain = find_chain(T_sorted, start_vertex, target_vertex,
                                  target_pos - start_pos);
    if (!chain) return false;
    return clip_to_meta_by_scanline(*chain, start_pos, FP0, FP1, FP2, out);
}

// Enumerate valid (P0, P1, P2) corner placements for F in F's frame.
// Returns up to 2 candidates: 1 if N01 is single-orbit, up to 2 for
// split-prime N01.  Each satisfies all three norm constraints and CCW
// orientation; picking among them which one matches the SURFACE
// geodesics is the walker's job (in embed_cell).
struct CornerCandidate { Eisenstein P0, P1, P2; };

std::vector<CornerCandidate>
enumerate_corner_candidates(double L20, double alpha_0,
                            long N01, long N12, long N20)
{
    std::vector<CornerCandidate> out;
    const Eisenstein P0(0, 0);
    std::vector<Eisenstein> P1_candidates = sector0_reps_of_norm((int)N01);

    for (Eisenstein P1 : P1_candidates) {
        auto [P1x, P1y] = P1.coord();
        const double theta_01 = std::atan2(P1y, P1x);
        const double theta_02 = theta_01 + alpha_0;
        const double P2x = L20 * std::cos(theta_02);
        const double P2y = L20 * std::sin(theta_02);
        const Eisenstein P2(std::pair<double,double>{P2x, P2y});

        if ((long)P2.norm2() == N20
            && (long)(P2 - P1).norm2() == N12
            && wedge(P1, P2) > 0)
        {
            out.push_back({P0, P1, P2});
        }
    }
    return out;
}

}  // namespace

Cell embed_cell(const DelaunayTriangulation& D,
                const Triangulation& T_sorted,
                int cell_id)
{
    Cell F;
    F.cell_id = cell_id;
    if (cell_id < 0 || cell_id >= D.nf || D.f_he[cell_id] < 0) return F;

    const int h0 = D.f_he[cell_id];
    const int h1 = D.he_next[h0];
    const int h2 = D.he_next[h1];
    F.c0 = D.he_origin[h0];
    F.c1 = D.he_origin[h1];
    F.c2 = D.he_origin[h2];

    const double L01 = D.he_length[h0];
    const double L20 = D.he_length[h2];
    const long N01 = (long)std::lround(L01 * L01);
    const long N12 = (long)std::lround(D.he_length[h1] * D.he_length[h1]);
    const long N20 = (long)std::lround(L20 * L20);
    const double alpha_0 = D.he_angle[h0];   // F's interior angle at c0

    std::vector<CornerCandidate> candidates =
        enumerate_corner_candidates(L20, alpha_0, N01, N12, N20);
    if (candidates.empty()) return F;

    // Try each corner candidate; accept the one where all three F-frame
    // walkers succeed.  For split-prime N01 there can be 2 candidates;
    // only one matches the surface geodesics' actual unfolding direction.
    for (const CornerCandidate& C : candidates) {
        EdgeStrip e01, e12, e20;
        if (!try_frame_walker(T_sorted, F.c0, C.P0, F.c1, C.P1,
                              C.P0, C.P1, C.P2, e01)) continue;
        if (!try_frame_walker(T_sorted, F.c1, C.P1, F.c2, C.P2,
                              C.P0, C.P1, C.P2, e12)) continue;
        if (!try_frame_walker(T_sorted, F.c2, C.P2, F.c0, C.P0,
                              C.P0, C.P1, C.P2, e20)) continue;
        F.P0 = C.P0;
        F.P1 = C.P1;
        F.P2 = C.P2;
        F.edge_01 = std::move(e01);
        F.edge_12 = std::move(e12);
        F.edge_20 = std::move(e20);
        F.ok = true;
        break;
    }
    return F;
}

std::vector<Cell> embed_all_cells(const DelaunayTriangulation& D,
                                   const Triangulation& T_sorted)
{
    std::vector<Cell> cells(D.nf);
#ifdef _OPENMP
    #pragma omp parallel for schedule(dynamic)
#endif
    for (int f = 0; f < D.nf; ++f)
        cells[f] = embed_cell(D, T_sorted, f);
    return cells;
}

// =====================================================================
// enumerate_cell_lattice: per-scanline east-walker over F's polygon.
// =====================================================================

namespace {

// Find the StripVertex with given F-frame lattice position in any of
// F's edge strips.  Returns nullptr if not found.
const StripVertex* find_seed(const Cell& F, Eisenstein p) {
    auto search = [&](const EdgeStrip& B) -> const StripVertex* {
        if (p.second < B.b_min || p.second > B.b_max) return nullptr;
        const auto& sl = B.by_scanline[p.second - B.b_min];
        for (const auto& sv : sl) if (sv.position == p) return &sv;
        return nullptr;
    };
    const StripVertex* sv = search(F.edge_01);
    if (sv) return sv;
    sv = search(F.edge_12);
    if (sv) return sv;
    sv = search(F.edge_20);
    return sv;
}

// Per-scanline cursor state for the east walk.
struct VertexK0 { int vertex; int k0; };

// One east-direction (dir 0) lattice step along T_sorted.  Throws on
// cone-gap (deg-5 vertex with nbr_idx == 5) or asymmetric T_sorted
// neighbour list.
VertexK0 east_step(const Triangulation& T_sorted,
                   VertexK0 cur, int cell_id, int b, int a)
{
    const int nbr_idx = ((0 - cur.k0) % 6 + 6) % 6;
    const int deg     = (int)T_sorted[cur.vertex].size();
    if (deg == 5 && nbr_idx == 5)
        paint_throw("eisenstein_paint::enumerate_cell_lattice(cell %d, b=%d, a=%d): "
                    "cone-gap step at vertex %d (deg=5, k0=%d, dir=0 -> nbr_idx=5)",
                    cell_id, b, a, cur.vertex, cur.k0);
    if (nbr_idx >= deg)
        paint_throw("eisenstein_paint::enumerate_cell_lattice(cell %d, b=%d, a=%d): "
                    "nbr_idx=%d >= deg=%d at vertex %d (k0=%d)",
                    cell_id, b, a, nbr_idx, deg, cur.vertex, cur.k0);
    const int nbr = T_sorted[cur.vertex][nbr_idx];
    int j = -1;
    const auto& nbr_list = T_sorted[nbr];
    for (int i = 0; i < (int)nbr_list.size(); ++i)
        if (nbr_list[i] == cur.vertex) { j = i; break; }
    if (j < 0)
        paint_throw("eisenstein_paint::enumerate_cell_lattice(cell %d): "
                    "T_sorted asymmetric (vertex %d not in T[%d])",
                    cell_id, cur.vertex, nbr);
    return { nbr, ((0 + 3 - j) % 6 + 6) % 6 };
}

}  // namespace

LatticeMap enumerate_cell_lattice(const Cell& F,
                                   const Triangulation& T_sorted)
{
    LatticeMap out;
    out.cell_id = F.cell_id;
    if (!F.ok)
        paint_throw("eisenstein_paint::enumerate_cell_lattice: F.ok=false (cell %d)",
                    F.cell_id);
    const ScanLines scan = scan_triangle(F.P0, F.P1, F.P2);

    for (int b = scan.b_min; b <= scan.b_max; ++b) {
        const ScanLine& sl = scan.lines[b - scan.b_min];
        if (sl.empty()) continue;

        Eisenstein pos(sl.a_left, b);
        const StripVertex* sv = find_seed(F, pos);
        if (!sv)
            paint_throw("eisenstein_paint::enumerate_cell_lattice(cell %d): no seed for "
                        "left chain pixel (%d, %d) in any edge",
                        F.cell_id, sl.a_left, b);
        VertexK0 cur{ sv->vertex_id, sv->frame_offset };
        out.entries.push_back({ pos, cur.vertex });

        for (int a = sl.a_left + 1; a <= sl.a_right; ++a) {
            pos = pos + unit_direction(0);
            // Edge strip is authoritative on F's boundary -- it carries
            // the canonical k0 derived from the strip walker.  Otherwise
            // step interior via T.
            if (const StripVertex* bv = find_seed(F, pos))
                cur = { bv->vertex_id, bv->frame_offset };
            else
                cur = east_step(T_sorted, cur, F.cell_id, b, a);
            out.entries.push_back({ pos, cur.vertex });
        }
    }
    return out;
}

// =====================================================================
// interpolate_cell: barycentric paint into pos3d.
// =====================================================================

void interpolate_cell(const Cell& F,
                      const LatticeMap& lmap,
                      const std::vector<coord3d>& cone_positions,
                      std::vector<coord3d>& pos3d)
{
    if (!F.ok)
        paint_throw("eisenstein_paint::interpolate_cell: F.ok=false");
    if (cone_positions.size() != 12)
        paint_throw("eisenstein_paint::interpolate_cell: cone_positions size %zu != 12",
                    cone_positions.size());
    const long g = wedge(F.P1 - F.P0, F.P2 - F.P0);
    if (g <= 0)
        paint_throw("eisenstein_paint::interpolate_cell(cell %d): g=%ld <= 0",
                    F.cell_id, g);
    const coord3d& C0 = cone_positions[F.c0];
    const coord3d& C1 = cone_positions[F.c1];
    const coord3d& C2 = cone_positions[F.c2];

    for (const auto& [p, vid] : lmap.entries) {
        if (vid < 12) continue;       // cone -- pre-painted
        const IntBary bw = integer_barycentric(p, F.P0, F.P1, F.P2);
        if (bw.n0 < 0 || bw.n1 < 0 || bw.n2 < 0 || bw.denom != g)
            paint_throw("eisenstein_paint::interpolate_cell(cell %d): bad "
                        "barycentric for vertex %d at pos=(%d,%d): "
                        "(n0=%ld n1=%ld n2=%ld, denom=%ld, g=%ld)",
                        F.cell_id, vid, p.first, p.second,
                        bw.n0, bw.n1, bw.n2, bw.denom, g);
        pos3d[vid] = barycentric_combine(reduce_to_lowest_terms(bw), C0, C1, C2);
    }
}

// =====================================================================
// run: full pipeline orchestrator.
// =====================================================================

namespace {

std::vector<Cell>
embed_all_or_stage(const DelaunayTriangulation& D, const Triangulation& T_sorted)
{
    std::vector<Cell> cells;
    try {
        cells = embed_all_cells(D, T_sorted);
    } catch (const std::exception& e) {
        throw StageError(Stage::EMBED,
            std::string("embed_all_throw: ") + e.what());
    }
    int n_live = 0, n_ok = 0;
    for (int f = 0; f < D.nf; ++f) {
        if (D.f_he[f] < 0) continue;
        ++n_live;
        if (cells[f].ok) ++n_ok;
    }
    if (n_ok != n_live)
        throw StageError(Stage::EMBED,
            "embed_cell partial: " + std::to_string(n_ok) + "/" +
            std::to_string(n_live));
    return cells;
}

std::vector<LatticeMap>
lmaps_or_stage(const std::vector<Cell>& cells, const Triangulation& T_sorted)
{
    std::vector<LatticeMap> lmaps(cells.size());
    try {
        for (size_t fi = 0; fi < cells.size(); ++fi)
            if (cells[fi].ok)
                lmaps[fi] = enumerate_cell_lattice(cells[fi], T_sorted);
    } catch (const std::exception& e) {
        throw StageError(Stage::EMBED,
            std::string("lmap_throw: ") + e.what());
    }
    return lmaps;
}

struct Coverage { int unclaimed = 0; int three_plus = 0; };

Coverage hex_coverage(const std::vector<Cell>& cells,
                      const std::vector<LatticeMap>& lmaps,
                      int Nv_sorted)
{
    std::vector<int> claim(Nv_sorted, 0);
    for (size_t fi = 0; fi < cells.size(); ++fi) {
        if (!cells[fi].ok) continue;
        for (const auto& [p, vid] : lmaps[fi].entries) {
            if (vid < 12) continue;
            ++claim[vid];
        }
    }
    Coverage c;
    for (int v = 12; v < Nv_sorted; ++v) {
        const int n = claim[v];
        if      (n == 0) ++c.unclaimed;
        else if (n >= 3) ++c.three_plus;
    }
    return c;
}

void check_coverage_or_stage(const std::vector<Cell>& cells,
                             const std::vector<LatticeMap>& lmaps,
                             int Nv_sorted)
{
    const Coverage c = hex_coverage(cells, lmaps, Nv_sorted);
    if (c.unclaimed == 0 && c.three_plus == 0) return;
    throw StageError(Stage::COVERAGE,
        "coverage(unclaimed=" + std::to_string(c.unclaimed) +
        ",three_plus="        + std::to_string(c.three_plus) + ")");
}

std::vector<coord3d>
interpolate_pos3d(const std::vector<Cell>& cells,
                  const std::vector<LatticeMap>& lmaps,
                  const std::vector<coord3d>& cone_positions,
                  int Nv_sorted)
{
    // Cones are written directly; interpolate_cell skips them
    // ("pre-painted" contract).  Hex slots are written by interpolation,
    // on-edge hex idempotently from both adjacent cells.
    std::vector<coord3d> pos3d(Nv_sorted);
    for (int c = 0; c < 12; ++c) pos3d[c] = cone_positions[c];
    try {
        for (size_t fi = 0; fi < cells.size(); ++fi)
            if (cells[fi].ok)
                interpolate_cell(cells[fi], lmaps[fi], cone_positions, pos3d);
    } catch (const std::exception& e) {
        throw StageError(Stage::INTERPOLATE,
            std::string("interpolate_throw: ") + e.what());
    }
    return pos3d;
}

// perm[u_orig] = u_sorted, so the position written at sorted vertex
// perm[u_orig] is u_orig's 3D coordinate.
std::vector<coord3d>
back_permute(const std::vector<coord3d>& pos_sorted,
             const Permutation&          perm,
             int                         Nv_orig)
{
    std::vector<coord3d> out(Nv_orig);
    for (int u = 0; u < Nv_orig; ++u)
        out[u] = pos_sorted[perm[u]];
    return out;
}

}  // namespace

Result run(const Triangulation& T) {
    Result R;
    try {
        Pipeline P     = prepare(T);
        auto     cells = embed_all_or_stage(P.D, P.T_sorted);
        auto     lmaps = lmaps_or_stage    (cells, P.T_sorted);
        check_coverage_or_stage(cells, lmaps, P.T_sorted.N);
        auto pos3d_sorted = interpolate_pos3d(cells, lmaps,
                                              P.cone_positions, P.T_sorted.N);
        R.coords = back_permute(pos3d_sorted, P.perm, T.N);
        R.stage  = Stage::OK;
    } catch (const StageError& e) {
        R.stage = e.stage;
        R.why   = e.what();
    } catch (const std::exception& e) {
        R.stage = Stage::UNEXPECTED;
        R.why   = std::string("unexpected: ") + e.what();
    }
    return R;
}

// =====================================================================
// TikZ visualisations.
// =====================================================================

namespace {
using tikz::cart;
using tikz::BBox;
}  // namespace

void dump_cell_tikz(const Cell& F,
                    const Triangulation& T_sorted,
                    std::ostream& os,
                    double scale)
{
    // Collect unique vertices with source-edge index.
    // 0 = edge_01, 1 = edge_12, 2 = edge_20.  First occurrence wins.
    std::unordered_map<int, std::pair<Eisenstein, int>> verts;
    auto add_edge = [&](const EdgeStrip& B, int idx) {
        for (const auto& sl : B.by_scanline)
            for (const auto& sv : sl)
                if (verts.find(sv.vertex_id) == verts.end())
                    verts[sv.vertex_id] = {sv.position, idx};
    };
    add_edge(F.edge_01, 0);
    add_edge(F.edge_12, 1);
    add_edge(F.edge_20, 2);

    os << "\\documentclass[tikz,border=4pt]{standalone}\n"
          "\\usepackage{tikz}\n"
          "\\begin{document}\n"
          "\\begin{tikzpicture}[scale=" << scale
       << ", every node/.style={font=\\tiny}]\n";

    os << "% cell " << F.cell_id
       << "  c0=" << F.c0 << " c1=" << F.c1 << " c2=" << F.c2
       << "  P0=(0,0) P1=(" << F.P1.first << "," << F.P1.second
       << ") P2=(" << F.P2.first << "," << F.P2.second
       << ")  edge sizes = " << verts.size() << " unique\n";

    BBox bb;
    bb.bump(F.P0); bb.bump(F.P1); bb.bump(F.P2);
    for (const auto& kv : verts) bb.bump(kv.second.first);
    const double margin = std::max(5.0, 0.5 * bb.span());
    bb.xmin -= margin; bb.xmax += margin;
    bb.ymin -= margin; bb.ymax += margin;
    os << "% Eisenstein lattice background\n";
    tikz::emit_grid(os, bb);

    // T_sorted edges among strip vertices.
    os << "% T_sorted edges among strip vertices\n";
    std::set<std::pair<int, int>> drawn;
    for (const auto& kv : verts) {
        const int w = kv.first;
        const Eisenstein pw = kv.second.first;
        for (int n : T_sorted[w]) {
            auto it = verts.find(n);
            if (it == verts.end()) continue;
            std::pair<int, int> e = std::minmax(w, n);
            if (!drawn.insert(e).second) continue;
            auto [x1, y1] = cart(pw);
            auto [x2, y2] = cart(it->second.first);
            os << "  \\draw[gray!55,line width=0.3pt] ("
               << x1 << "," << y1 << ") -- (" << x2 << "," << y2 << ");\n";
        }
    }

    // F's lattice triangle in red.
    os << "% iDT geodesics P0 -> P1 -> P2 -> P0\n";
    tikz::emit_lattice_triangle(os, F.P0, F.P1, F.P2);

    // Strip vertices.  Corners as red squares; others as circles
    // colour-coded by source edge.  Vertices that fall OUTSIDE F's
    // lattice triangle (a bug) get a thick red outline.
    auto inside_meta = [&](Eisenstein p) {
        return wedge(F.P1 - F.P0, p - F.P0) >= 0
            && wedge(F.P2 - F.P1, p - F.P1) >= 0
            && wedge(F.P0 - F.P2, p - F.P2) >= 0;
    };
    os << "% strip vertices (red ring = OUTSIDE F's lattice triangle)\n";
    for (const auto& kv : verts) {
        const int v_id = kv.first;
        const Eisenstein p = kv.second.first;
        const int bi = kv.second.second;
        auto [x, y] = cart(p);
        const bool corner  = (v_id == F.c0 || v_id == F.c1 || v_id == F.c2);
        const bool outside = !inside_meta(p);
        const char* fill = corner       ? "red!75!black"
                          : (bi == 0)   ? "blue!70!black"
                          : (bi == 1)   ? "green!55!black"
                                        : "orange!75!black";
        const char* shape = corner ? "rectangle" : "circle";
        const char* draw_spec = outside ? "draw=red,line width=0.9pt"
                                        : "draw=black";
        os << "  \\node[" << shape << "," << draw_spec << ",fill=" << fill
           << ",inner sep=1.6pt,label=below right:{\\tiny " << v_id
           << "}] at (" << x << "," << y << ") {};\n";
    }

    os << "\\end{tikzpicture}\n\\end{document}\n";
}

void dump_lattice_map_tikz(const Cell& F,
                           const LatticeMap& lmap,
                           const Triangulation& T_sorted,
                           std::ostream& os,
                           double scale,
                           int highlight_vertex_id)
{
    (void)T_sorted;     // not currently used; kept for signature symmetry

    os << "\\documentclass[tikz,border=4pt]{standalone}\n"
          "\\usepackage{tikz}\n"
          "\\begin{document}\n"
          "\\begin{tikzpicture}[scale=" << scale
       << ", every node/.style={font=\\tiny}]\n";
    os << "% cell " << F.cell_id
       << "  c0=" << F.c0 << " c1=" << F.c1 << " c2=" << F.c2
       << "  P0=(0,0) P1=(" << F.P1.first << "," << F.P1.second
       << ") P2=(" << F.P2.first << "," << F.P2.second << ")  "
       << lmap.entries.size() << " lattice pts\n";

    BBox bb;
    bb.bump(F.P0); bb.bump(F.P1); bb.bump(F.P2);
    for (const auto& [p, v] : lmap.entries) bb.bump(p);
    const double margin = std::max(3.0, 0.3 * bb.span());
    bb.xmin -= margin; bb.xmax += margin;
    bb.ymin -= margin; bb.ymax += margin;
    tikz::emit_grid(os, bb);
    tikz::emit_lattice_triangle(os, F.P0, F.P1, F.P2);

    // Classify each lattice point: corner cone, on-edge hex (in any
    // edge strip), or strict interior.
    auto is_on_bdry = [&](Eisenstein p) -> bool {
        auto check = [&](const EdgeStrip& B) {
            if (p.second < B.b_min || p.second > B.b_max) return false;
            const auto& sl = B.by_scanline[p.second - B.b_min];
            for (const auto& sv : sl) if (sv.position == p) return true;
            return false;
        };
        return check(F.edge_01) || check(F.edge_12) || check(F.edge_20);
    };

    for (const auto& [p, v] : lmap.entries) {
        auto [x, y] = cart(p);
        const bool is_corner = (v == F.c0 || v == F.c1 || v == F.c2);
        const bool on_bdry   = is_on_bdry(p);
        const bool highlight = (v == highlight_vertex_id);
        const char* shape = is_corner ? "rectangle" : "circle";
        const char* fill  = is_corner ? "red!75!black"
                          : on_bdry   ? "green!55!black"
                                      : "blue!70!black";
        const char* draw  = highlight ? "draw=red,line width=1.2pt"
                                      : "draw=black";
        const double sep = highlight ? 2.6 : 1.6;
        os << "  \\node[" << shape << "," << draw << ",fill=" << fill
           << ",inner sep=" << sep << "pt,label=below right:{\\tiny " << v
           << "}] at (" << x << "," << y << ") {};\n";
    }
    os << "\\end{tikzpicture}\n\\end{document}\n";
}

}  // namespace eisenstein_paint
