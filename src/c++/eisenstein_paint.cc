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
// Internal helpers shared across the intrinsic layer.
// =====================================================================

// Layer I has ONE failure channel: every modeled failure throws
// PaintError with its stage typed AT THE THROW SITE (the facades
// translate to Status; nothing re-wraps or relabels).
[[noreturn]] void embed_fail(std::string msg) {
    throw PaintError(Code::EMBED, std::move(msg));
}

// The cone count of a cones-first-sorted triangulation: cones (degree
// != 6) are sorted FIRST, so the count is the first index of degree 6.
// The single definition of "how many cones" -- sorted_dual and the
// per-cell primitives both call it.
// @pre cones-first order (sort_flat_last applied)
int count_cones(const TriangulationView& T_sorted) {
    int n = 0;
    while (n < (int)T_sorted.N && T_sorted.degree(n) != 6) ++n;
    return n;
}

}  // namespace

// =====================================================================
// code_name
// =====================================================================

const char* code_name(Code s) {
    switch (s) {
        case Code::OK:              return "ok";
        case Code::IDT_COMPUTE:     return "idt_compute";
        case Code::ALEXANDROV:      return "alexandrov";
        case Code::NON_SIMPLICIAL:  return "non_simplicial";
        case Code::EMBED:           return "embed";
        case Code::COVERAGE:        return "coverage";
        case Code::INTERPOLATE:     return "interpolate";
        case Code::UNEXPECTED:      return "unexpected";
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
    S.n_cones = count_cones(S.T);
    return S;
}

namespace {

// The raw iDT of S's cone metric; shared by dual_idt (guarded) and
// Layer II's realize_dual (unguarded -- the Alexandrov flips supersede).
DelaunayTriangulation compute_dual_idt(const SortedDual& S, const char* who) {
    try {
        return DelaunayTriangulation::compute(S.T);
    } catch (const std::exception& e) {
        throw PaintError(Code::IDT_COMPUTE,
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
            throw PaintError(Code::NON_SIMPLICIAL,
                std::string(who) + ": iDT has a self-loop (a cone edge to "
                "itself); a face with a repeated corner cannot be charted");
    if (!D.is_well_formed())
        throw PaintError(Code::NON_SIMPLICIAL,
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

// The cell's corner triangle in some Eisenstein frame, with the two
// exact containment predicates of the chart machinery (integer wedges).
struct MetaTriangle {
    Eisenstein P0, P1, P2;

    bool contains_strict(Eisenstein p) const {
        return wedge(P1 - P0, p - P0) > 0
            && wedge(P2 - P1, p - P1) > 0
            && wedge(P0 - P2, p - P2) > 0;
    }
    bool contains_closed(Eisenstein p) const {
        return wedge(P1 - P0, p - P0) >= 0
            && wedge(P2 - P1, p - P1) >= 0
            && wedge(P0 - P2, p - P2) >= 0;
    }
    // The same triangle expressed relative to a new origin (wedges are
    // translation-invariant, so containment answers are unchanged).
    MetaTriangle relative_to(Eisenstein origin) const {
        return { P0 - origin, P1 - origin, P2 - origin };
    }
};

// The isometric development of one iDT arc's strip into the START
// vertex's local frame (start at (0, 0)): every developed T_sorted
// vertex mapped to its position and frame offset k0.
struct EdgeDevelopment {
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

// ALL flat, orientation-preserving developments of one arc's strip
// (deduplicated): every development whose walker terminates at
// target_vertex with start at (0,0), end at target_rel, no cone on the
// strict interior of the walk, and no non-corner cone strictly inside the
// cell triangle M_rel (given in the start-relative frame).
//
// Why ALL, not the first.  A single cone pair can be joined by MULTIPLE
// same-length lattice geodesics that all terminate at target_vertex and all
// pass the per-edge filters above, yet develop DIFFERENT mesh face corridors.
// Two distinct phenomena produce this:
//   (i)  Split-prime edge lengths (norm N with two sector-0 reps in mirror
//        orbits, e.g. N=19,28,37): two geodesics whose wrong branch wraps a
//        NEIGHBOURING cone -- caught by the cell-triangle cone-interiority test.
//   (ii) Obtuse/thin cells (one interior angle near 0): two long edge geodesics
//        run nearly parallel through the sliver, and walk_line, trying start
//        fans in list order, can reach target_vertex via a HEX-only corridor on
//        the wrong side of the true edge (e.g. C102 cell 13, corners (8,5,3),
//        L=(27,1,28): edge c0->c1 develops the false corridor {8,41,42,33,34,5}
//        instead of the true {8,38,27,26,25,5}).  Both corridors are cone-free,
//        so the cone-interiority test cannot separate them.
// The discriminator for (ii) is CROSS-EDGE consistency, resolved one level up
// in embed_cell: each F-frame lattice position is one surface point = one mesh
// vertex, so the correct triple of edge developments must agree wherever their
// clipped strips share a position.  edge_developments therefore returns every
// qualifying development (deduped by its (vertex_id -> position, k0) map) and
// lets embed_cell pick the mutually-consistent triple.  For the common,
// non-thin cell each edge yields exactly one development and the search is a
// single trivially-consistent triple.
//
// DETERMINISM: the enumeration order below -- endpoints [target_rel,
// rotation-representative], start fans in neighbour-list order, dir_uv
// 0..5, first-wins dedup -- is part of the algorithm (downstream selection
// takes the FIRST consistent triple); do not reorder.
std::vector<EdgeDevelopment>
edge_developments(const TriangulationView& T_sorted, int n_cones,
                  int start_vertex, int target_vertex, Eisenstein target_rel,
                  const MetaTriangle& M_rel)
{
    std::vector<EdgeDevelopment> results;
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

                EdgeDevelopment V;
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
                // inside the cell triangle means this walk took the wrong
                // same-length geodesic branch.  Corners (start/target and the
                // third corner) sit ON the boundary, so they are never strictly
                // interior and need no special-casing.
                bool folds_cone = false;
                for (const auto& [vid, pk] : V.by_id) {
                    if (vid >= n_cones) continue;                                // hex
                    if (vid == start_vertex || vid == target_vertex) continue;   // this edge's cones
                    if (M_rel.contains_strict(pk.first)) {
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

// Restrict a development to the cell triangle: translate to F's frame
// (add start_pos), keep the vertices inside the closed cell, bucket by
// scanline (b coordinate).  Returns false if no vertex survives.
bool clip_to_cell_by_scanline(const EdgeDevelopment& V,
                              Eisenstein start_pos,
                              const MetaTriangle& M_cell,
                              EdgeStrip& out)
{
    int b_min = INT_MAX, b_max = INT_MIN;
    std::vector<StripVertex> kept;
    kept.reserve(V.by_id.size());
    for (const auto& [vid, pk] : V.by_id) {
        const Eisenstein pos_F = pk.first + start_pos;
        if (!M_cell.contains_closed(pos_F)) continue;
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

// All strip developments of one arc: every qualifying development of
// the arc start->target, clipped to the cell triangle and scanline-
// bucketed.  Empty if no development qualifies or none survives the
// clip; embed_cell picks the one strip per arc that is mutually
// consistent across the three arcs.
std::vector<EdgeStrip>
arc_strip_developments(const TriangulationView& T_sorted, int n_cones,
                       int start_vertex, Eisenstein start_pos,
                       int target_vertex, Eisenstein target_pos,
                       const MetaTriangle& M_cell)
{
    std::vector<EdgeStrip> out;
    for (const EdgeDevelopment& dev :
         edge_developments(T_sorted, n_cones, start_vertex, target_vertex,
                           target_pos - start_pos,
                           M_cell.relative_to(start_pos))) {
        EdgeStrip strip;
        if (clip_to_cell_by_scanline(dev, start_pos, M_cell, strip))
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

// A lattice development of the cell: a valid (P0, P1, P2) corner
// placement in F's frame -- one per sector-0 representative of N01
// whose induced apex P2 (its direction is fixed by the metric) is a
// lattice point meeting the remaining norm constraints and CCW
// orientation.  The count is a(N(gcd(P1, P2))) + [delta | sqrt(N01)]
// -- see the eisenstein_paint_tables.hh banner for THE statement --
// typically 1 or 2 but UNBOUNDED over the isomer space (C980's cells
// admit 4), so callers must never cap it.  Picking which development
// matches the SURFACE geodesics is the walker's job (in embed_cell).
struct Development { Eisenstein P0, P1, P2; };

// One cell's iDT metric datum: CCW corner ids, the squared side norms
// and the interior angle at corner 0 -- the input of the chart-frame
// construction.  THE one derivation; embed_cell and
// cell_developments both open with it.  @pre D.f_he[f] >= 0.
struct CellMetric {
    std::array<int, 3> corners;
    double L20 = 0, alpha_0 = 0;
    long   N01 = 0, N12 = 0, N20 = 0;
};

CellMetric cell_metric(const DelaunayView& D, int f)
{
    const int h0 = D.f_he[f];
    const int h1 = D.he_next[h0];
    const int h2 = D.he_next[h1];
    CellMetric m;
    m.corners = { D.he_origin[h0], D.he_origin[h1], D.he_origin[h2] };
    const double L01 = D.he_length[h0];
    m.L20     = D.he_length[h2];
    m.N01     = (long)std::lround(L01 * L01);
    m.N12     = (long)std::lround(D.he_length[h1] * D.he_length[h1]);
    m.N20     = (long)std::lround(m.L20 * m.L20);
    m.alpha_0 = D.he_angle[h0];   // the interior angle at corner 0
    return m;
}

std::vector<Development>
enumerate_developments(double L20, double alpha_0,
                       long N01, long N12, long N20)
{
    std::vector<Development> out;
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

// The cross-edge consistency selection: the FIRST triple of arc strips
// (in the given candidate orders -- the order is part of the algorithm's
// determinism) that glue to ONE position -> vertex map wherever they
// overlap (strips_consistent).  First-match IS the search for the
// unique solution (the cell-development lemma): the true geodesic
// c0 -> c1 has one developed displacement in the surface's own frame,
// exactly one of its six rotations lies in sector 0, and every OTHER
// sector-0 representative of N01 sits in a different unit orbit -- no
// walk along the surface can terminate on it -- so at most one
// development admits a consistent strip triple.
bool first_consistent_triple(const std::array<std::vector<EdgeStrip>, 3>& cand,
                             std::array<const EdgeStrip*, 3>& out)
{
    for (const EdgeStrip& e01 : cand[0])
        for (const EdgeStrip& e12 : cand[1])
            for (const EdgeStrip& e20 : cand[2])
                if (strips_consistent(e01, e12, e20)) {
                    out = { &e01, &e12, &e20 };
                    return true;
                }
    return false;
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

    const CellMetric m = cell_metric(D, cell_id);
    F.corners = m.corners;

    // Each development is a lattice-realisable corner placement (see
    // Development for the count -- unbounded, never capped); each arc
    // can admit several strip developments (see edge_developments).
    // The cell's chart is the FIRST placement + strip triple that is
    // cross-edge consistent -- unique by the cell-development lemma
    // (first_consistent_triple) -- searched in deterministic order.
    for (const Development& C :
         enumerate_developments(m.L20, m.alpha_0, m.N01, m.N12, m.N20)) {
        const MetaTriangle M = { C.P0, C.P1, C.P2 };
        const std::array<Eisenstein, 3> P = { C.P0, C.P1, C.P2 };
        std::array<std::vector<EdgeStrip>, 3> cand;
        for (int k = 0; k < 3; ++k)
            cand[k] = arc_strip_developments(T_sorted, n_cones,
                                             F.corners[k],           P[k],
                                             F.corners[(k + 1) % 3], P[(k + 1) % 3],
                                             M);
        if (cand[0].empty() || cand[1].empty() || cand[2].empty()) continue;

        std::array<const EdgeStrip*, 3> strips;
        if (first_consistent_triple(cand, strips)) {
            F.P = P;
            for (int k = 0; k < 3; ++k) F.edge[k] = *strips[k];
            F.ok = true;
            break;
        }
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
    for (int k = 0; k < 3; ++k) {
        const EdgeStrip& B = F.edge[k];
        if (p.second < B.b_min || p.second > B.b_max) continue;
        for (const auto& sv : B.by_scanline[p.second - B.b_min])
            if (sv.position == p) return &sv;
    }
    return nullptr;
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
        embed_fail("eisenstein_paint::enumerate_cell_lattice(cell " +
                   std::to_string(cell_id) + ", b=" + std::to_string(b) +
                   ", a=" + std::to_string(a) + "): cone-gap step at vertex " +
                   std::to_string(cur.vertex) + " (deg=5, k0=" +
                   std::to_string(cur.k0) + ", dir=0 -> nbr_idx=5)");
    if (nbr_idx >= deg)
        embed_fail("eisenstein_paint::enumerate_cell_lattice(cell " +
                   std::to_string(cell_id) + ", b=" + std::to_string(b) +
                   ", a=" + std::to_string(a) + "): nbr_idx=" +
                   std::to_string(nbr_idx) + " >= deg=" + std::to_string(deg) +
                   " at vertex " + std::to_string(cur.vertex) +
                   " (k0=" + std::to_string(cur.k0) + ")");
    const int nbr = T_sorted[cur.vertex][nbr_idx];
    int j = -1;
    const auto& nbr_list = T_sorted[nbr];
    for (int i = 0; i < (int)nbr_list.size(); ++i)
        if (nbr_list[i] == cur.vertex) { j = i; break; }
    if (j < 0)
        embed_fail("eisenstein_paint::enumerate_cell_lattice(cell " +
                   std::to_string(cell_id) + "): T_sorted asymmetric (vertex " +
                   std::to_string(cur.vertex) + " not in T[" +
                   std::to_string(nbr) + "])");
    return { nbr, ((0 + 3 - j) % 6 + 6) % 6 };
}

}  // namespace

LatticeMap enumerate_cell_lattice(const Cell& F,
                                  const TriangulationView& T_sorted)
{
    LatticeMap out;
    out.cell_id = F.cell_id;
    if (!F.ok)
        embed_fail("eisenstein_paint::enumerate_cell_lattice: F.ok=false (cell " +
                   std::to_string(F.cell_id) + ")");
    const ScanLines scan = scan_triangle(F.P[0], F.P[1], F.P[2]);

    for (int b = scan.b_min; b <= scan.b_max; ++b) {
        const ScanLine& sl = scan.lines[b - scan.b_min];
        if (sl.empty()) continue;

        Eisenstein pos(sl.a_left, b);
        const StripVertex* sv = find_seed(F, pos);
        if (!sv)
            embed_fail("eisenstein_paint::enumerate_cell_lattice(cell " +
                       std::to_string(F.cell_id) + "): no seed for left chain "
                       "pixel (" + std::to_string(sl.a_left) + ", " +
                       std::to_string(b) + ") in any edge");
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

// The cell-development lemma as a guard: every live cell of an
// intrinsically-Delaunay complex over the cone metric admits a flat
// isometric Eisenstein chart, so a live cell with no consistent chart is
// an EMBED failure, reported by cell id.
void require_all_charted(const DelaunayTriangulation& D,
                         const std::vector<Cell>& cells)
{
    for (int f = 0; f < D.nf; ++f)
        if (D.f_he[f] >= 0 && !cells[f].ok)
            embed_fail("eisenstein_paint::parametrize: live cell " +
                       std::to_string(f) + " admits no consistent chart");
}

// append_scan_rows -- the ONE row-table derivation, shared by the
// post-embed flatten and the candidate builder: append one triangle
// scan's rows and return its block.  Empty scanlines keep the scan's
// own (a_left > a_right) values, so both derivations produce identical
// bytes.
ScanBlock append_scan_rows(std::vector<ScanRow>& rows, const ScanLines& scan)
{
    ScanBlock sb;
    sb.b_min = scan.b_min;
    sb.b_max = scan.b_max;
    sb.rows_first = (int32_t)rows.size();
    int32_t running = 0;
    for (const ScanLine& sl : scan.lines) {
        rows.push_back(ScanRow{sl.a_left, sl.a_right, running});
        if (!sl.empty()) running += sl.a_right - sl.a_left + 1;
    }
    sb.n_entries = running;
    return sb;
}

// Flatten one cell's lattice map into the CSR: the row table comes from
// the winning frame's OWN triangle scan (append_scan_rows -- the same
// derivation the candidate tables use), and the entries are verified to
// land exactly where the rows say as they are appended.
//
// The convexity lemma (@anchor paint-cell-convexity, PROMOTION-DESIGN
// sec 2.3): a cell is a convex lattice triangle, so each scanline's
// claimed points are ONE contiguous a-range -- claim()'s (row, a-offset)
// bijection rests on it.  enumerate_cell_lattice produces exactly the
// scan's runs by construction; this check is that postcondition made
// loud (a divergence is a chart fold or a scan drift, never silently
// mis-indexed later).
void flatten_cell(SurfaceParametrization& P, int f, const LatticeMap& lmap)
{
    const Eisenstein P0(0, 0);
    const Eisenstein P1(P.frames[f].p1a, P.frames[f].p1b);
    const Eisenstein P2(P.frames[f].p2a, P.frames[f].p2b);
    P.scans[f] = append_scan_rows(P.rows, scan_triangle(P0, P1, P2));
    const ScanBlock& sb = P.scans[f];

    auto fold = [&](int b, const char* what) {
        throw PaintError(Code::EMBED,
            "flatten_cell: cell " + std::to_string(f) + " scanline b=" +
            std::to_string(b) + ": " + what + " (chart fold)");
    };
    size_t i = 0;
    for (int b = sb.b_min; b <= sb.b_max; ++b) {
        const ScanRow& row = P.rows[(size_t)sb.rows_first + (b - sb.b_min)];
        for (int32_t a = row.a_left; a <= row.a_right; ++a) {
            if (i >= lmap.entries.size()) fold(b, "entries end before the scan");
            const auto& [pos, vid] = lmap.entries[i];
            if (pos.first != a || pos.second != b)
                fold(b, "entry does not match the scan run");
            P.entries.push_back(LatticePoint{pos.first, pos.second, vid});
            ++i;
        }
    }
    if (i != lmap.entries.size())
        fold(sb.b_max, "entries outlast the scan");
}

// Per-vertex chart appearances as a CSR: count, prefix, fill -- the
// fill iterates cells in id order, so each vertex's occurrences keep
// the old push_back order (cell-major).
void build_occurrence_csr(SurfaceParametrization& P)
{
    const int Nv = (int)P.T.N;
    P.occ_first.assign(Nv + 1, 0);
    for (const LatticePoint& e : P.entries) P.occ_first[e.vid + 1]++;
    for (int v = 0; v < Nv; ++v) P.occ_first[v + 1] += P.occ_first[v];
    P.occ.resize(P.entries.size());
    std::vector<int32_t> cursor(P.occ_first.begin(), P.occ_first.end() - 1);
    const ParamTablesView V = P.view();
    for (int f = 0; f < (int)P.cells.size(); ++f)
        for (const LatticePoint& e : V.cell_entries(f))
            P.occ[cursor[e.vid]++] = VertexOccurrence{f, e.a, e.b};
}

// The coverage lemma as a guard: every non-cone vertex appears in exactly
// 1 chart (interior) or 2 (on a shared iDT edge, where the paint is
// idempotent); 0 or >= 3 appearances is a chart fold.
// @anchor paint-coverage
void check_coverage(const SurfaceParametrization& P)
{
    const ParamTablesView V = P.view();
    int unclaimed = 0, three_plus = 0;
    for (int v = P.n_cones; v < (int)P.T.N; ++v) {
        const int n = (int)V.occurrences(v).size();
        if      (n == 0) ++unclaimed;
        else if (n >= 3) ++three_plus;
    }
    if (unclaimed != 0 || three_plus != 0)
        throw PaintError(Code::COVERAGE,
            "coverage(unclaimed=" + std::to_string(unclaimed) +
            ",three_plus="        + std::to_string(three_plus) + ")");
}

}  // namespace

SurfaceParametrization parametrize(const DelaunayTriangulation& D,
                                   const SortedDual& S)
{
    SurfaceParametrization P;
    P.D       = D;                       // owner IS-A view: alias, caller keeps D alive
    P.T       = S.T;
    P.n_cones = S.n_cones;
    P.perm.resize(S.perm.size());
    for (size_t u = 0; u < S.perm.size(); ++u) P.perm[u] = S.perm[u];

    const std::vector<Cell> construction = embed_all_cells(D, S.T);
    require_all_charted(D, construction);

    const int nf = (int)construction.size();
    P.cells.assign(nf, CellCorners{});
    P.frames.assign(nf, CellFrame{});
    P.scans.assign(nf, ScanBlock{});
    P.entry_first.assign(nf + 1, 0);

    for (int f = 0; f < nf; ++f) {
        P.entry_first[f] = (int32_t)P.entries.size();
        if (!construction[f].ok) {
            // dead slot: empty block, rows_first kept monotone
            P.scans[f].rows_first = (int32_t)P.rows.size();
            continue;
        }
        const Cell& F = construction[f];
        P.cells[f]  = CellCorners{F.corners[0], F.corners[1], F.corners[2]};
        P.frames[f] = CellFrame{F.P[1].first, F.P[1].second,
                                F.P[2].first, F.P[2].second};
        flatten_cell(P, f, enumerate_cell_lattice(F, S.T));
    }
    P.entry_first[nf] = (int32_t)P.entries.size();

    build_occurrence_csr(P);
    check_coverage(P);
    return P;
}

VertexOccurrence locate_vertex(const SurfaceParametrization& P, int v) {
    if (v < 0 || v >= (int)P.T.N)
        throw std::logic_error("eisenstein_paint::locate_vertex: vertex "
                               + std::to_string(v) + " out of range");
    const auto os = P.view().occurrences(v);
    if (os.empty())
        throw std::logic_error("eisenstein_paint::locate_vertex: vertex "
                               + std::to_string(v) + " unclaimed by any cell");
    return os[0];
}

// =====================================================================
// cell_developments: every admissible lattice development, CSR by cell.
// =====================================================================

namespace {

// Smallest power of two >= x, floored at 64 (open-addressing scratch
// wants pow2 extents with load-factor headroom).
int64_t pow2_ceil64(int64_t x) { int64_t p = 64; while (p < x) p <<= 1; return p; }

}  // namespace

CellDevelopments cell_developments(const ::DelaunayView& D,
                                   const SortedDual& S)
{
    CellDevelopments cd;
    cd.nf       = D.nf;
    cd.N_sorted = (int32_t)S.T.N;
    cd.n_cones  = S.n_cones;
    cd.cells.assign(cd.nf, CellCorners{});
    cd.dev_first.assign(cd.nf + 1, 0);
    cd.entry_capacity_first.assign(cd.nf + 1, 0);

    long s_max = 0;    // max |a|+|b| over development frame edges

    for (int f = 0; f < cd.nf; ++f) {
        cd.dev_first[f] = (int32_t)cd.frames.size();
        if (D.f_he[f] < 0) continue;    // dead slot: 0 developments

        const CellMetric m = cell_metric(D, f);
        cd.cells[f] = CellCorners{ m.corners[0], m.corners[1], m.corners[2] };

        // The development count is unbounded over the isomer space
        // (see the header banner for the formula; C980's cells admit
        // 4), hence the CSR.  A live cell with NO development has an
        // unrealizable metric: refuse loudly.
        const auto devs =
            enumerate_developments(m.L20, m.alpha_0, m.N01, m.N12, m.N20);
        if (devs.empty())
            throw PaintError(Code::EMBED,
                "cell_developments: live cell " + std::to_string(f) +
                " admits no lattice development");

        // Per-development frame + triangle scan (append_scan_rows: the
        // same row derivation the post-selection flatten uses).
        int32_t cap_f = 0;
        for (const Development& C : devs) {
            const Eisenstein P1 = C.P1, P2 = C.P2;
            cd.frames.push_back(CellFrame{ P1.first, P1.second, P2.first, P2.second });
            for (Eisenstein e : { P1, P2 - P1, C.P0 - P2 }) {
                const long sd = std::labs(e.first) + std::labs(e.second);
                if (sd > s_max) s_max = sd;
            }
            const ScanBlock sb =
                append_scan_rows(cd.rows, scan_triangle(C.P0, P1, P2));
            cd.scans.push_back(sb);
            const int32_t n_rows = sb.b_max - sb.b_min + 1;
            if (n_rows > cd.rows_cap) cd.rows_cap = n_rows;
            if (sb.n_entries > cap_f) cap_f = sb.n_entries;
        }
        if (cap_f > cd.max_cell_entries) cd.max_cell_entries = cap_f;
        cd.entry_capacity_first[f + 1] = cap_f;
    }
    cd.dev_first[cd.nf] = (int32_t)cd.frames.size();
    for (int f = 0; f < cd.nf; ++f)    // exclusive prefix
        cd.entry_capacity_first[f + 1] += cd.entry_capacity_first[f];
    cd.entries_capacity = cd.entry_capacity_first[cd.nf];

    // Scratch capacity formulas (properties of the data; any executor
    // sizing embed/enumerate scratch reads these):
    //  * walk registration: distinct keys are T_sorted vertex ids, and a
    //    walk registers <= 3 vertices per pushed face with <=
    //    walk_max_steps faces per primitive sub-walk, so
    //    min(N_sorted, 3*(walk_max_steps+2)) bounds the live count; x2
    //    keeps the open-addressing load factor comfortable.
    //  * boundary map: a qualifying walk to displacement s = |a|+|b|
    //    pushes <= 3s faces, so <= 3*(3s+1)+3 registrations per edge,
    //    three edges; x2 headroom.
    cd.wcap = pow2_ceil64(
        2 * (std::min<int64_t>(cd.N_sorted, 3 * (int64_t)(walk_max_steps + 2)) + 8));
    cd.bcap = pow2_ceil64(2 * (3 * (3 * s_max + 4) * 3 + 16));

    return cd;
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
    // Collect unique vertices with source-arc index k (strip of arc
    // corners[k] -> corners[k+1]).  First occurrence wins.
    std::unordered_map<int, std::pair<Eisenstein, int>> verts;
    for (int k = 0; k < 3; ++k)
        for (const auto& sl : F.edge[k].by_scanline)
            for (const auto& sv : sl)
                if (verts.find(sv.vertex_id) == verts.end())
                    verts[sv.vertex_id] = {sv.position, k};

    os << "\\documentclass[tikz,border=4pt]{standalone}\n"
          "\\usepackage{tikz}\n"
          "\\begin{document}\n"
          "\\begin{tikzpicture}[scale=" << scale
       << ", every node/.style={font=\\tiny}]\n";

    os << "% cell " << F.cell_id
       << "  c0=" << F.corners[0] << " c1=" << F.corners[1] << " c2=" << F.corners[2]
       << "  P0=(0,0) P1=(" << F.P[1].first << "," << F.P[1].second
       << ") P2=(" << F.P[2].first << "," << F.P[2].second
       << ")  edge sizes = " << verts.size() << " unique\n";

    BBox bb;
    bb.bump(F.P[0]); bb.bump(F.P[1]); bb.bump(F.P[2]);
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
    tikz::emit_lattice_triangle(os, F.P[0], F.P[1], F.P[2]);

    // Strip vertices.  Corners as red squares; others as circles
    // colour-coded by source arc.  Vertices that fall OUTSIDE F's
    // lattice triangle (a bug) get a thick red outline.
    const MetaTriangle M = { F.P[0], F.P[1], F.P[2] };
    os << "% strip vertices (red ring = OUTSIDE F's lattice triangle)\n";
    for (const auto& kv : verts) {
        const int v_id = kv.first;
        const Eisenstein p = kv.second.first;
        const int bi = kv.second.second;
        auto [x, y] = cart(p);
        const bool corner  = (v_id == F.corners[0] || v_id == F.corners[1] ||
                              v_id == F.corners[2]);
        const bool outside = !M.contains_closed(p);
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
       << "  c0=" << F.corners[0] << " c1=" << F.corners[1] << " c2=" << F.corners[2]
       << "  P0=(0,0) P1=(" << F.P[1].first << "," << F.P[1].second
       << ") P2=(" << F.P[2].first << "," << F.P[2].second << ")  "
       << lmap.entries.size() << " lattice pts\n";

    BBox bb;
    bb.bump(F.P[0]); bb.bump(F.P[1]); bb.bump(F.P[2]);
    for (const auto& [p, v] : lmap.entries) bb.bump(p);
    const double margin = std::max(3.0, 0.3 * bb.span());
    bb.xmin -= margin; bb.xmax += margin;
    bb.ymin -= margin; bb.ymax += margin;
    tikz::emit_grid(os, bb);
    tikz::emit_lattice_triangle(os, F.P[0], F.P[1], F.P[2]);

    // Classify each lattice point: corner cone, on-edge hex (in any
    // edge strip -- find_seed), or strict interior.
    for (const auto& [p, v] : lmap.entries) {
        auto [x, y] = cart(p);
        const bool is_corner = (v == F.corners[0] || v == F.corners[1] ||
                                v == F.corners[2]);
        const bool on_bdry   = (find_seed(F, p) != nullptr);
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
