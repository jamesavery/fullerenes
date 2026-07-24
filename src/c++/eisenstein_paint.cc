// Eisenstein-paint, Layer I: intrinsic surface parametrization (see the
// header banner).  Everything here is integer/combinatorial -- no coord3d
// appears in this translation unit.  Layer II (the Alexandrov realization
// and chart evaluation) lives in eisenstein_paint_geometry.cc.

#include "fullerenes/eisenstein_paint.hh"

#include "fullerenes/eisenstein.hh"
#include "fullerenes/eisenstein_raster.hh"
#include "fullerenes/eisenstein_tikz.hh"
#include "fullerenes/delaunay_strip.hh"

#include <algorithm>
#include <array>
#include <cassert>
#include <climits>
#include <cmath>
#include <cstdarg>
#include <cstdio>
#include <cstdlib>
#include <map>
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
// stage_name
// =====================================================================

const char* stage_name(Stage s) {
    switch (s) {
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
// sorted_dual / dual_idt.
// =====================================================================

SortedDual sorted_dual(const TriangulationView& T) {
    SortedDual S;
    S.perm = T.sort_flat_last();
    S.T    = Triangulation(T);
    S.T.apply_permutation(S.perm);
    S.n_cones = 0;
    for (node_t u = 0; u < S.T.N; ++u)
        if (S.T.degree(u) != 6) ++S.n_cones;
    return S;
}

namespace {

// The raw iDT of S's cone metric; shared by dual_idt (guarded) and
// Layer II's realize_dual (unguarded -- the Alexandrov flips supersede).
DelaunayTriangulation compute_dual_idt(const SortedDual& S, const char* who) {
    try {
        return DelaunayTriangulation::compute(S.T);
    } catch (const std::exception& e) {
        throw StageError(Stage::IDT_COMPUTE,
            std::string(who) + ": iDT compute: " + e.what());
    }
}

// Guard the iDT before the chart machinery runs.  Multi-edges (two
// same-length geodesics between one cone pair) are legitimate delta-complex
// features that develop correctly (find_chains' discriminators), so they are
// ACCEPTED.  What stays a hard reject is what embed_cell genuinely cannot
// chart: a self-loop (a cone edge (v, v)) makes a face with a repeated
// corner, and is_well_formed catches bigons / non-length-3 faces.  (Lifting
// the repeated-corner limitation is queued: refactor-debt entry
// 2026-07-24-paint-repeated-corner-cells.)
void require_chartable_idt(const DelaunayTriangulation& D, const char* who) {
    for (int h = 0; h < D.nh; ++h)
        if (D.alive(h) && D.he_origin[h] == D.dest(h))
            throw StageError(Stage::NON_SIMPLICIAL,
                std::string(who) + ": iDT has a self-loop (a cone edge to "
                "itself); a face with a repeated corner cannot be charted");
    if (!D.is_well_formed())
        throw StageError(Stage::NON_SIMPLICIAL,
            std::string(who) + ": iDT is not well-formed "
            "(some live half-edge is not in a length-3 face cycle)");
}

}  // namespace

// Exposed to Layer II (eisenstein_paint_geometry.cc) via this internal
// declaration; not part of the public header.
DelaunayTriangulation detail_compute_dual_idt(const SortedDual& S, const char* who) {
    return compute_dual_idt(S, who);
}

DelaunayTriangulation dual_idt(const SortedDual& S) {
    DelaunayTriangulation D = compute_dual_idt(S, "eisenstein_paint::dual_idt");
    require_chartable_idt(D, "eisenstein_paint::dual_idt");
    return D;
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
int derive_k0(const TriangulationView& T, int w, int nbr, int d) {
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
bool extract_walk_vertices(const TriangulationView& T,
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

// A cone vertex (id < n_cones) at any STRICT interior step of wr.walk means
// the directed line passes through a cone -- and at a cone the
// angular structure is undefined, so a geodesic cannot pass through
// one.  Trivially false when gcd(target_rel) == 1.
bool walks_through_cone(const WalkResult& wr, int n_cones) {
    for (size_t k = 1; k + 1 < wr.walk.size(); ++k)
        if (wr.walk[k].second < n_cones) return true;
    return false;
}

// A developed vertex `p` (cell frame) is STRICTLY inside the meta triangle
// (FP0, FP1, FP2) iff all three edge wedges are > 0.  Used to detect a
// FOLDED development: an embedded iDT face carries cones only at its three
// corners, so any OTHER cone landing strictly interior means the walk
// developed the wrong same-length geodesic (a split-prime multi-edge whose
// wrong branch wraps a neighbouring cone into the strip) or a reflected
// corridor -- both non-embeddings that the charts later reject.
bool strictly_inside_meta(Eisenstein p, Eisenstein FP0, Eisenstein FP1, Eisenstein FP2) {
    return wedge(FP1 - FP0, p - FP0) > 0
        && wedge(FP2 - FP1, p - FP1) > 0
        && wedge(FP0 - FP2, p - FP2) > 0;
}

// Enumerate ALL developments (deduplicated) whose walker terminates at
// target_vertex with a clean, orientation-preserving, NON-FOLDING development:
// start at (0,0), end at target_rel, no cone on the strict interior of the
// walk, and no non-corner cone strictly inside the meta triangle (FP0,FP1,FP2).
//
// Why ALL, not the first.  A single cone pair can be joined by MULTIPLE
// same-length lattice geodesics that all terminate at target_vertex and all
// pass the per-edge filters above, yet develop DIFFERENT mesh face corridors.
// Two distinct phenomena produce this:
//   (i)  Split-prime edge lengths (norm N with two sector-0 reps in mirror
//        orbits, e.g. N=19,28,37): two geodesics whose wrong branch wraps a
//        NEIGHBOURING cone -- caught by the meta-triangle cone-interiority test.
//   (ii) Obtuse/thin cells (one interior angle near 0): two long edge geodesics
//        run nearly parallel through the sliver, and walk_line, trying start
//        fans in list order, can reach target_vertex via a HEX-only corridor on
//        the wrong side of the true edge (e.g. C102 cell 13, corners (8,5,3),
//        L=(27,1,28): edge c0->c1 develops the false corridor {8,41,42,33,34,5}
//        instead of the true {8,38,27,26,25,5}).  Both corridors are cone-free,
//        so the cone-interiority test cannot separate them.
// The discriminator for (ii) is CROSS-EDGE consistency, resolved one level up in
// embed_cell: each F-frame lattice position is one surface point = one mesh
// vertex, so the correct triple of edge developments must agree wherever their
// clipped strips share a position.  find_chains therefore returns every
// qualifying development (deduped by its (vertex_id -> position, k0) map) and
// lets embed_cell pick the mutually-consistent triple.  For the common,
// non-thin cell each edge yields exactly one development and the search is a
// single trivially-consistent triple.
std::vector<WalkVertices>
find_chains(const TriangulationView& T_sorted, int n_cones,
           int start_vertex, int target_vertex, Eisenstein target_rel,
           Eisenstein start_pos, Eisenstein FP0, Eisenstein FP1, Eisenstein FP2)
{
    std::vector<WalkVertices> results;
    std::set<std::vector<std::array<long,4>>> seen;   // dedup by development map
    // Candidate endpoints, all ORIENTATION-PRESERVING (`back` is a rotation,
    // never a reflection): the raw cell-frame displacement, and its sector0
    // rotation-representative (a >= 0, b >= 0), where walk_line's corridor is
    // developed.  Rotation keeps the SAME chirality orbit -- the same geodesic
    // -- so the neighbour cycles stay CCW, which the CCW-only east-walker in
    // enumerate_cell_lattice requires.  (A reflection would swap to the mirror
    // geodesic and reverse every cycle, which folds; we never do it.)  The six
    // unit rotations of any nonzero vector are 60 degrees apart, so exactly one
    // lands in the 60-degree sector0.
    Eisenstein rotrep{0, 0};
    bool have_rot = false;
    for (int k = 0; k < 6; ++k) {
        const Eisenstein E = target_rel * unit_direction(k);   // rotate by w^k
        if (E.first >= 0 && E.second >= 0) { rotrep = E; have_rot = true; break; }
    }
    std::vector<Eisenstein> endpoints;
    endpoints.push_back(target_rel);
    if (have_rot && rotrep != target_rel) endpoints.push_back(rotrep);

    const int deg = (int)T_sorted[start_vertex].size();
    for (const Eisenstein endpoint : endpoints) {
        const LatticeIsometry back = align(endpoint, target_rel);   // walk frame -> cell frame (a rotation)
        for (int i = 0; i < deg; ++i) {
            const int v_a = T_sorted[start_vertex][i];
            const int v_b = T_sorted[start_vertex][(i + 1) % deg];
            for (int dir_uv = 0; dir_uv < 6; ++dir_uv) {
                WalkResult wr = walk_line(T_sorted, start_vertex, v_a, v_b, dir_uv, endpoint);
                if (wr.final_vertex != target_vertex) continue;

                WalkVertices V;
                if (!extract_walk_vertices(T_sorted, wr, V.by_id)) continue;
                const auto itu = V.by_id.find(start_vertex);
                const auto itv = V.by_id.find(target_vertex);
                if (itu == V.by_id.end() || itu->second.first != Eisenstein(0, 0)) continue;
                if (itv == V.by_id.end() || itv->second.first != endpoint)         continue;
                if (walks_through_cone(wr, n_cones)) continue;

                for (auto& [vid, pk] : V.by_id) {
                    pk.first = back.apply(pk.first);
                    // k0 is a lattice DIRECTION (T[v][0]'s developed direction), so
                    // it rotates with the frame: map it through back's linear part.
                    const Eisenstein d  = unit_direction(pk.second);
                    const Eisenstein td = back.u * (back.reflect ? d.complex_conj() : d);
                    pk.second = direction_of_unit(td);
                }

                // Reject a FOLDED development: any non-corner cone strictly
                // inside the meta triangle means this walk took the wrong
                // same-length geodesic branch.  Corners (start/target and the
                // third corner) sit ON the boundary, so they are never strictly
                // interior and need no special-casing.
                bool folds_cone = false;
                for (const auto& [vid, pk] : V.by_id) {
                    if (vid >= n_cones) continue;                                // hex
                    if (vid == start_vertex || vid == target_vertex) continue;   // this edge's cones
                    if (strictly_inside_meta(pk.first + start_pos, FP0, FP1, FP2)) {
                        folds_cone = true; break;
                    }
                }
                if (folds_cone) continue;

                // Deduplicate: many (arc, dir_uv) start combinations develop
                // the identical corridor.  Signature is the sorted set of
                // (vertex_id, pos.a, pos.b, k0) tuples in the cell frame.
                std::vector<std::array<long,4>> sig;
                sig.reserve(V.by_id.size());
                for (const auto& [vid, pk] : V.by_id)
                    sig.push_back({vid, pk.first.first, pk.first.second, pk.second});
                std::sort(sig.begin(), sig.end());
                if (seen.insert(sig).second) results.push_back(std::move(V));
            }
        }
    }
    return results;
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

// Run the walker from start_vertex at start_pos in F's frame, aiming at
// target_vertex at target_pos.  Returns EVERY qualifying development as an
// F-frame strip (clipped to the FP0/FP1/FP2 lattice triangle).  Empty if no
// chain qualifies or none survives the clip.  embed_cell picks the one strip
// per edge that is mutually consistent across the three edges.
std::vector<EdgeStrip>
frame_walker_candidates(const TriangulationView& T_sorted, int n_cones,
                        int start_vertex, Eisenstein start_pos,
                        int target_vertex, Eisenstein target_pos,
                        Eisenstein FP0, Eisenstein FP1, Eisenstein FP2)
{
    std::vector<EdgeStrip> out;
    const auto chains = find_chains(T_sorted, n_cones, start_vertex, target_vertex,
                                    target_pos - start_pos,
                                    start_pos, FP0, FP1, FP2);
    for (const WalkVertices& chain : chains) {
        EdgeStrip strip;
        if (clip_to_meta_by_scanline(chain, start_pos, FP0, FP1, FP2, strip))
            out.push_back(std::move(strip));
    }
    return out;
}

// Add every (position -> vertex_id) of `B` into `claim`.  Returns false the
// moment a position is already claimed by a DIFFERENT vertex -- the cross-edge
// fold that means the strips developed inconsistent corridors.
bool merge_strip_claims(const EdgeStrip& B,
                        std::map<std::pair<int,int>, int>& claim)
{
    for (const auto& row : B.by_scanline)
        for (const StripVertex& sv : row) {
            const std::pair<int,int> key{sv.position.first, sv.position.second};
            const auto [it, ins] = claim.emplace(key, sv.vertex_id);
            if (!ins && it->second != sv.vertex_id) return false;
        }
    return true;
}

// True iff the three clipped edge strips form ONE consistent position->vertex
// map: every F-frame lattice position claimed by more than one strip carries
// the same mesh vertex.  This is the cross-edge discriminator that rejects a
// thin/obtuse cell's wrong-branch edge development (see find_chains).
bool strips_consistent(const EdgeStrip& e01, const EdgeStrip& e12, const EdgeStrip& e20)
{
    std::map<std::pair<int,int>, int> claim;
    return merge_strip_claims(e01, claim)
        && merge_strip_claims(e12, claim)
        && merge_strip_claims(e20, claim);
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

// The cone count of the charted complex: cones are the sorted labels
// whose T_sorted degree != 6.  They are sorted FIRST, so the count is
// the first index with degree 6 (or N if none).  O(n_cones) and purely
// combinatorial -- keeps the per-cell primitives' signatures free of a
// redundant n_cones parameter.
int count_cones(const TriangulationView& T_sorted) {
    int n = 0;
    while (n < (int)T_sorted.N && T_sorted.degree(n) != 6) ++n;
    return n;
}

}  // namespace

Cell embed_cell(const DelaunayTriangulation& D,
                const TriangulationView& T_sorted,
                int cell_id)
{
    Cell F;
    F.cell_id = cell_id;
    if (cell_id < 0 || cell_id >= D.nf || D.f_he[cell_id] < 0) return F;
    const int n_cones = count_cones(T_sorted);

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

    // Try each corner candidate; accept the one where the three F-frame edge
    // walkers admit a MUTUALLY CONSISTENT development.  Each edge can yield
    // several qualifying developments (split-prime multi-edges, or the
    // wrong-side corridor of a thin/obtuse cell -- see find_chains).  The
    // correct cell development is unique: a shared F-frame lattice position is
    // one surface point = one mesh vertex, so the true triple agrees wherever
    // two strips overlap.  We therefore search the (small, deduped) product of
    // the three edges' candidates for the first pairwise-consistent triple.
    // For the common non-thin cell each edge has exactly one candidate and this
    // is a single trivially-consistent check.  For split-prime N01 there can be
    // 2 corner candidates; only one matches the surface geodesics.
    for (const CornerCandidate& C : candidates) {
        const std::vector<EdgeStrip> cand01 =
            frame_walker_candidates(T_sorted, n_cones, F.c0, C.P0, F.c1, C.P1, C.P0, C.P1, C.P2);
        const std::vector<EdgeStrip> cand12 =
            frame_walker_candidates(T_sorted, n_cones, F.c1, C.P1, F.c2, C.P2, C.P0, C.P1, C.P2);
        const std::vector<EdgeStrip> cand20 =
            frame_walker_candidates(T_sorted, n_cones, F.c2, C.P2, F.c0, C.P0, C.P0, C.P1, C.P2);
        if (cand01.empty() || cand12.empty() || cand20.empty()) continue;

        bool placed = false;
        for (const EdgeStrip& e01 : cand01) {
            for (const EdgeStrip& e12 : cand12) {
                for (const EdgeStrip& e20 : cand20) {
                    if (!strips_consistent(e01, e12, e20)) continue;
                    F.P0 = C.P0;
                    F.P1 = C.P1;
                    F.P2 = C.P2;
                    F.edge_01 = e01;
                    F.edge_12 = e12;
                    F.edge_20 = e20;
                    F.ok = true;
                    placed = true;
                    break;
                }
                if (placed) break;
            }
            if (placed) break;
        }
        if (placed) break;
    }
    return F;
}

std::vector<Cell> embed_all_cells(const DelaunayTriangulation& D,
                                  const TriangulationView& T_sorted)
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
VertexK0 east_step(const TriangulationView& T_sorted,
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
                                  const TriangulationView& T_sorted)
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
// parametrize: charts + coverage + occurrences.
// =====================================================================

namespace {

std::vector<Cell>
embed_all_or_stage(const DelaunayTriangulation& D, const TriangulationView& T_sorted)
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
lmaps_or_stage(const std::vector<Cell>& cells, const TriangulationView& T_sorted)
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

}  // namespace

SurfaceParametrization parametrize(const DelaunayTriangulation& D,
                                   const SortedDual& S)
{
    SurfaceParametrization P;
    P.D       = &D;
    P.T       = S.T;
    P.n_cones = S.n_cones;
    P.cells   = embed_all_or_stage(D, S.T);
    P.lmaps   = lmaps_or_stage(P.cells, S.T);

    // Occurrences: one entry per (cell, lattice claim).  Coverage on the
    // non-cone vertices: claimed by 1 cell (interior) or 2 (on a shared
    // iDT edge, idempotent); 0 or >= 3 is a chart fold.
    P.occurrences.assign(S.T.N, {});
    for (size_t fi = 0; fi < P.cells.size(); ++fi) {
        if (!P.cells[fi].ok) continue;
        for (const auto& [pos, vid] : P.lmaps[fi].entries)
            P.occurrences[vid].push_back({ (int)fi, pos });
    }
    int unclaimed = 0, three_plus = 0;
    for (int v = S.n_cones; v < (int)S.T.N; ++v) {
        const int n = (int)P.occurrences[v].size();
        if      (n == 0) ++unclaimed;
        else if (n >= 3) ++three_plus;
    }
    if (unclaimed != 0 || three_plus != 0)
        throw StageError(Stage::COVERAGE,
            "coverage(unclaimed=" + std::to_string(unclaimed) +
            ",three_plus="        + std::to_string(three_plus) + ")");
    return P;
}

Occurrence locate_vertex(const SurfaceParametrization& P, int v) {
    if (v < 0 || v >= (int)P.occurrences.size() || P.occurrences[v].empty())
        throw std::logic_error("eisenstein_paint::locate_vertex: vertex "
                               + std::to_string(v) + " unclaimed by any cell");
    return P.occurrences[v].front();
}

// =====================================================================
// TikZ visualisations.
// =====================================================================

namespace {
using tikz::cart;
using tikz::BBox;
}  // namespace

void dump_cell_tikz(const Cell& F,
                    const TriangulationView& T_sorted,
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
                           const TriangulationView& T_sorted,
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
