#include "fullerenes/buckinverse.hh"
#include <algorithm>
#include <set>

namespace buckinverse {

// =====================================================================
// Seed graph construction from RSPI (ring spiral pentagon indices)
// RSPIs looked up from the isomer database via IsomerDB::readPDB.
// The DB stores 1-indexed values; these are converted to 0-indexed.
// =====================================================================

// C20 (Ih): 1 isomer. DB RSPI = {1..12}, 0-indexed = {0..11}.
Triangulation makeSeedC20() {
    return FullereneDual(20, {0,1,2,3,4,5,6,7,8,9,10,11});
}

// C28 (Td): isomer #2 of 2. DB RSPI = {1,2,3,5,7,9,10,11,12,13,14,15}.
Triangulation makeSeedC28() {
    return FullereneDual(28, {0,1,2,4,6,8,9,10,11,12,13,14});
}

// C30 (D5h): isomer #1 of 3. DB RSPI = {1,2,3,4,5,6,12,13,14,15,16,17}.
Triangulation makeSeedC30() {
    return FullereneDual(30, {0,1,2,3,4,5,11,12,13,14,15,16});
}

SeedType identifySeed(const Triangulation& t) {
    // Seeds are uniquely identified by dual vertex count among fullerenes.
    // C20: 12 dual vertices (all degree-5)
    // C28: 16 dual vertices (12 deg-5, 4 deg-6)
    // C30: 17 dual vertices (12 deg-5, 5 deg-6)
    switch (t.N) {
        case 12: return SeedType::C20;
        case 16: return SeedType::C28;
        case 17: return SeedType::C30;
        default: return SeedType::NotASeed;
    }
}

// =====================================================================
// Path computation
// =====================================================================

PathInfo computeStraightPath(const Graph& g, node_t u, node_t v,
                             Dir dir, int numEntries) {
    PathInfo pi;
    pi.path.resize(numEntries);
    pi.path[0] = u;
    pi.path[1] = v;
    for (int k = 2; k < numEntries; ++k)
        pi.path[k] = straightAhead(g, dir, pi.path[k-1], pi.path[k-2]);

    // Parallel path: numEntries entries (matching Haskell convention).
    // For k < numEntries-1: sideNbr at edge path[k] -> path[k+1].
    // For k = numEntries-1: sideNbr at path[k] -> continuation direction.
    // The last entry is the "true parallel" at the far endpoint (not a strip vertex).
    pi.parallel.resize(numEntries);
    for (int k = 0; k < numEntries - 1; ++k)
        pi.parallel[k] = sideNbr(g, dir, pi.path[k], pi.path[k+1]);
    {
        int k = numEntries - 1;
        node_t cont = straightAhead(g, dir, pi.path[k], pi.path[k-1]);
        pi.parallel[k] = sideNbr(g, dir, pi.path[k], cont);
    }

    // Check for self-intersection
    std::set<node_t> seen;
    for (node_t x : pi.path)     if (!seen.insert(x).second) return pi;
    for (node_t x : pi.parallel) if (!seen.insert(x).second) return pi;
    pi.valid = true;
    return pi;
}

PathInfo computeBentZeroPath(const Graph& g, node_t u, node_t v, Dir dir) {
    PathInfo pi;
    pi.path.resize(5);
    pi.path[0] = u;
    pi.path[1] = v;
    pi.path[2] = straightAhead(g, dir, v, u);
    pi.path[3] = turnAhead(g, dir, pi.path[2], pi.path[1]);
    pi.path[4] = straightAhead(g, dir, pi.path[3], pi.path[2]);

    // Parallel path: skip edge INTO the turn vertex (edge 1->2).
    // par[0] = sideNbr at edge 0->1 (before turn)
    // par[1] = sideNbr at edge 2->3 (after turn)
    // par[2] = sideNbr at edge 3->4 (after turn)
    pi.parallel.resize(3);
    pi.parallel[0] = sideNbr(g, dir, pi.path[0], pi.path[1]);
    pi.parallel[1] = sideNbr(g, dir, pi.path[2], pi.path[3]);
    pi.parallel[2] = sideNbr(g, dir, pi.path[3], pi.path[4]);

    std::set<node_t> seen;
    for (node_t x : pi.path)     if (!seen.insert(x).second) return pi;
    for (node_t x : pi.parallel) if (!seen.insert(x).second) return pi;
    pi.valid = true;
    return pi;
}

PathInfo computeBentPath(const Graph& g, node_t u, node_t v,
                         Dir dir, int bi, int bj) {
    PathInfo pi;
    int pathLen = bi + bj + 5;
    pi.path.resize(pathLen);
    pi.path[0] = u;
    pi.path[1] = v;

    // Pre-bend phase: bi+2 straight steps
    for (int k = 2; k <= bi + 2; ++k)
        pi.path[k] = straightAhead(g, dir, pi.path[k-1], pi.path[k-2]);

    // Turn at bend point
    pi.path[bi + 3] = turnAhead(g, dir, pi.path[bi+2], pi.path[bi+1]);

    // Post-bend phase: bj+1 straight steps
    for (int k = bi + 4; k < pathLen; ++k)
        pi.path[k] = straightAhead(g, dir, pi.path[k-1], pi.path[k-2]);

    // Parallel path: skip edge INTO the turn vertex (edge bi+1 -> bi+2).
    // Pre-bend: sideNbr at edges k -> k+1 for k = 0..bi (bi+1 entries)
    // Post-bend: sideNbr at edges k -> k+1 for k = bi+2..pathLen-2 (bj+2 entries)
    // Total: bi+1 + bj+2 = bi+bj+3
    int parLen = bi + bj + 3;
    pi.parallel.resize(parLen);
    int pi_idx = 0;
    for (int k = 0; k <= bi; ++k)
        pi.parallel[pi_idx++] = sideNbr(g, dir, pi.path[k], pi.path[k+1]);
    for (int k = bi + 2; k < pathLen - 1; ++k)
        pi.parallel[pi_idx++] = sideNbr(g, dir, pi.path[k], pi.path[k+1]);

    std::set<node_t> seen;
    for (node_t x : pi.path)     if (!seen.insert(x).second) return pi;
    for (node_t x : pi.parallel) if (!seen.insert(x).second) return pi;
    pi.valid = true;
    return pi;
}

// =====================================================================
// Reduction surgery (Phase 3)
// =====================================================================

// Adjacency list modification helpers (operate on a mutable copy)
using adj_t = std::vector<std::vector<node_t>>;

// Returns false if old_nbr is not found in adj[u].
static bool replaceNbrAdj(adj_t& adj, node_t u, node_t old_nbr, node_t new_nbr) {
    for (auto& n : adj[u]) {
        if (n == old_nbr) { n = new_nbr; return true; }
    }
    return false;
}

static void removeNbrAdj(adj_t& adj, node_t u, node_t v) {
    auto& ns = adj[u];
    ns.erase(std::remove(ns.begin(), ns.end(), v), ns.end());
}

// Find non-path, non-strip neighbours of a vertex
static std::vector<node_t> nonPathStripNbrs(const Graph& g, node_t u,
        const std::set<node_t>& pathSet, const std::set<node_t>& stripSet) {
    std::vector<node_t> result;
    for (int i = 0; i < g.degree(u); ++i) {
        node_t n = g.nbrs(u)[i];
        if (!pathSet.count(n) && !stripSet.count(n))
            result.push_back(n);
    }
    return result;
}

// Find true parallel vertices for a straight reduction (L_i).
// strip has pathlength entries, true_par has pathlength+1 entries.
// The last true_par entry comes from the extended parallel (far_par).
// Returns empty vector if true parallel can't be found (invalid topology).
static std::vector<node_t> findTrueParStraight(const Graph& g,
        const std::vector<node_t>& path, const std::vector<node_t>& strip,
        node_t far_par) {
    int pathlength = (int)strip.size();
    std::set<node_t> pathSet(path.begin(), path.end());
    std::set<node_t> stripSet(strip.begin(), strip.end());

    std::vector<node_t> tp(pathlength + 1);

    // tp[0]: non-path, non-strip nbr of strip[0] adjacent to path[0]
    auto nps0 = nonPathStripNbrs(g, strip[0], pathSet, stripSet);
    tp[0] = -1;
    for (node_t n : nps0) {
        if (g.arc_ix(path[0], n) >= 0) { tp[0] = n; break; }
    }
    if (tp[0] < 0) return {};

    // tp[k] for k=1..pathlength-1: common non-path, non-strip nbr of strip[k-1] and strip[k]
    for (int k = 1; k < pathlength; ++k) {
        auto nps_prev = nonPathStripNbrs(g, strip[k-1], pathSet, stripSet);
        auto nps_curr = nonPathStripNbrs(g, strip[k], pathSet, stripSet);
        std::set<node_t> prev_set(nps_prev.begin(), nps_prev.end());
        tp[k] = -1;
        for (node_t n : nps_curr) {
            if (prev_set.count(n)) { tp[k] = n; break; }
        }
        if (tp[k] < 0) return {};
    }

    // tp[pathlength]: the far parallel vertex (from extended parallel)
    tp[pathlength] = far_par;

    return tp;
}

// Find true parallel vertices for a bent reduction (B_{i,j}).
// strip has numNew entries, true_par has numNew entries.
// Returns empty vector if true parallel can't be found (invalid topology).
static std::vector<node_t> findTrueParBent(const Graph& g,
        const std::vector<node_t>& path, const std::vector<node_t>& strip,
        int bentPos) {
    int numNew = (int)strip.size();
    std::set<node_t> pathSet(path.begin(), path.end());
    std::set<node_t> stripSet(strip.begin(), strip.end());

    std::vector<node_t> tp(numNew);

    // tp[0]: non-path, non-strip nbr of strip[0] adjacent to path[0]
    auto nps0 = nonPathStripNbrs(g, strip[0], pathSet, stripSet);
    tp[0] = -1;
    for (node_t n : nps0) {
        if (g.arc_ix(path[0], n) >= 0) { tp[0] = n; break; }
    }
    if (tp[0] < 0) return {};

    // Pre-bend: tp[k] for k=1..bentPos: common of strip[k-1] and strip[k]
    for (int k = 1; k <= bentPos; ++k) {
        auto nps_prev = nonPathStripNbrs(g, strip[k-1], pathSet, stripSet);
        auto nps_curr = nonPathStripNbrs(g, strip[k], pathSet, stripSet);
        std::set<node_t> prev_set(nps_prev.begin(), nps_prev.end());
        tp[k] = -1;
        for (node_t n : nps_curr) {
            if (prev_set.count(n)) { tp[k] = n; break; }
        }
        if (tp[k] < 0) return {};
    }

    // Bend vertex (strip[bentPos+1]): sole non-path, non-strip nbr
    int bendIdx = bentPos + 1;
    auto npsBend = nonPathStripNbrs(g, strip[bendIdx], pathSet, stripSet);
    if (npsBend.size() != 1) return {};
    tp[bendIdx] = npsBend[0];

    // Post-bend: use "other than previous" approach.
    // strip[k] has two non-path, non-strip nbrs: tp[k-1] and tp[k].
    // We know tp[k-1] (from bend or previous iteration), so tp[k] is the other.
    for (int k = bendIdx + 1; k < numNew; ++k) {
        auto nps = nonPathStripNbrs(g, strip[k], pathSet, stripSet);
        tp[k] = -1;
        for (node_t n : nps) {
            if (n != tp[k-1]) { tp[k] = n; break; }
        }
        if (tp[k] < 0) return {};
    }

    return tp;
}

// Validate that all strip vertex neighbors are in {path, strip, tp}.
// Returns false if any strip vertex has an "extra" neighbor, which means
// the reduction site doesn't have clean expansion topology and can't
// be surgically inverted.
static bool stripTopologyClean(const Graph& g,
        const std::vector<node_t>& path, const std::vector<node_t>& strip,
        const std::vector<node_t>& tp) {
    std::set<node_t> validSet;
    for (node_t x : path) validSet.insert(x);
    for (node_t x : strip) validSet.insert(x);
    for (node_t x : tp) validSet.insert(x);

    for (node_t s : strip) {
        for (int i = 0; i < g.degree(s); ++i) {
            node_t n = g.nbrs(s)[i];
            if (!validSet.count(n)) return false;
        }
    }
    return true;
}

// Forward declarations for reconnection functions
static bool reconnectStraight(adj_t& adj,
        const std::vector<node_t>& path, const std::vector<node_t>& strip,
        const std::vector<node_t>& tp);
static bool reconnectBent(adj_t& adj,
        const std::vector<node_t>& path, const std::vector<node_t>& strip,
        const std::vector<node_t>& tp, int bentPos, int bentLen);

// =====================================================================
// Inversion: reduce arbitrary graphs (Phase 3)
// Strip-based approach: identify strip vertices (degree-5 pair for L0)
// and reconstruct path/tp from CW ordering around them.
// =====================================================================

std::vector<InvSite> findL0InvSites(const Graph& g) {
    std::vector<InvSite> result;
    auto d5 = deg5vertices(g);

    for (node_t s0 : d5) {
        for (int ni = 0; ni < g.degree(s0); ++ni) {
            node_t s1 = g.nbrs(s0)[ni];
            if (g.degree(s1) != 5) continue;
            if (s1 <= s0) continue;  // each pair once

            // Try both directions
            for (Dir dir : {Dir::DRight, Dir::DLeft}) {
                // Identify vertices from CW ordering around s0.
                // DRight: s0's CW = [path[0], tp[0], tp[1], s1, path[1]]
                // DLeft:  s0's CW = [path[0], path[1], s1, tp[1], tp[0]]
                node_t path_1, tp_1, path_0, tp_0;
                if (dir == Dir::DRight) {
                    path_1 = g.next(s0, s1);
                    tp_1   = g.prev(s0, s1);
                    path_0 = g.next(s0, path_1);
                    tp_0   = g.prev(s0, tp_1);
                } else {
                    path_1 = g.prev(s0, s1);
                    tp_1   = g.next(s0, s1);
                    path_0 = g.prev(s0, path_1);
                    tp_0   = g.next(s0, tp_1);
                }

                // Identify from s1's CW ordering.
                // DRight: s1's CW = [s0, tp[1], tp[2], path[2], path[1]]
                // DLeft:  s1's CW = [s0, path[1], path[2], tp[2], tp[1]]
                node_t path_2, tp_2;
                if (dir == Dir::DRight) {
                    path_2 = g.prev(s1, path_1);
                    tp_2   = g.prev(s1, path_2);
                } else {
                    path_2 = g.next(s1, path_1);
                    tp_2   = g.next(s1, path_2);
                }

                // Validate
                // path[0] should be degree 6 (promoted by expansion)
                if (g.degree(path_0) != 6) continue;
                // tp[2] should be degree 6 (promoted by expansion)
                if (g.degree(tp_2) != 6) continue;
                // path[0] must be adj to path[1]
                if (g.arc_ix(path_0, path_1) < 0) continue;
                // path[2] must be adj to path[1]
                if (g.arc_ix(path_2, path_1) < 0) continue;

                // All 8 vertices must be distinct
                std::set<node_t> all = {s0, s1, path_0, path_1, path_2,
                                         tp_0, tp_1, tp_2};
                if ((int)all.size() != 8) continue;

                InvSite site;
                site.kind = Lk(0);
                site.dir = dir;
                site.strip = {s0, s1};
                site.path = {path_0, path_1, path_2};
                site.tp = {tp_0, tp_1, tp_2};
                result.push_back(site);
            }
        }
    }
    return result;
}

// Find all valid L_i (i >= 1) inversion sites.
// Strip chain: s0 (deg-5), s1..s_i (deg-6), s_{i+1} (deg-5).
// Chain walk through degree-6 interior is direction-independent (advanceCW 3).
static std::vector<InvSite> findLiInvSites(const Graph& g, int maxRedLen) {
    std::vector<InvSite> result;
    auto d5 = deg5vertices(g);

    for (node_t s0 : d5) {
        for (int ni = 0; ni < g.degree(s0); ++ni) {
            node_t s1 = g.nbrs(s0)[ni];
            if (g.degree(s1) != 6) continue;

            // Walk strip chain through degree-6 vertices
            std::vector<node_t> sc = {s0, s1};  // strip chain
            node_t prv = s0, cur = s1;
            bool chainOK = true;
            while (g.degree(cur) == 6) {
                node_t nxt = advanceCW(g, cur, prv, 3);
                for (node_t x : sc)
                    if (x == nxt) { chainOK = false; break; }
                if (!chainOK) break;
                sc.push_back(nxt);
                prv = cur; cur = nxt;
                if ((int)sc.size() > maxRedLen + 1) { chainOK = false; break; }
            }
            if (!chainOK || g.degree(cur) != 5) continue;

            int i_val = (int)sc.size() - 2;
            if (i_val < 1 || i_val + 1 > maxRedLen) continue;

            for (Dir dir : {Dir::DRight, Dir::DLeft}) {
                int pLen = i_val + 3;
                std::vector<node_t> path(pLen), tp(pLen);

                // From s0 (degree 5)
                if (dir == Dir::DRight) {
                    path[1] = g.next(s0, s1);
                    path[0] = g.next(s0, path[1]);
                    tp[0]   = g.next(s0, path[0]);
                    tp[1]   = g.next(s0, tp[0]);
                } else {
                    path[1] = g.prev(s0, s1);
                    path[0] = g.prev(s0, path[1]);
                    tp[1]   = g.next(s0, s1);
                    tp[0]   = g.next(s0, tp[1]);
                }

                // Interior s_j (degree 6, j=1..i_val)
                bool ok = true;
                for (int j = 1; j <= i_val && ok; ++j) {
                    node_t sj = sc[j], sp = sc[j-1], sn = sc[j+1];
                    if (dir == Dir::DRight) {
                        if (g.next(sj, sp) != tp[j]) { ok = false; break; }
                        tp[j+1] = g.next(sj, tp[j]);
                        if (g.next(sj, tp[j+1]) != sn) { ok = false; break; }
                        path[j+1] = g.next(sj, sn);
                        if (g.next(sj, path[j+1]) != path[j]) { ok = false; break; }
                    } else {
                        if (g.next(sj, sp) != path[j]) { ok = false; break; }
                        path[j+1] = g.next(sj, path[j]);
                        if (g.next(sj, path[j+1]) != sn) { ok = false; break; }
                        tp[j+1] = g.next(sj, sn);
                        if (g.next(sj, tp[j+1]) != tp[j]) { ok = false; break; }
                    }
                }
                if (!ok) continue;

                // Far endpoint s_{i+1} (degree 5)
                node_t sf = sc[i_val+1], sfp = sc[i_val];
                if (dir == Dir::DRight) {
                    if (g.next(sf, sfp) != tp[i_val+1]) continue;
                    tp[i_val+2] = g.next(sf, tp[i_val+1]);
                    path[i_val+2] = g.next(sf, tp[i_val+2]);
                    if (g.next(sf, path[i_val+2]) != path[i_val+1]) continue;
                } else {
                    if (g.next(sf, sfp) != path[i_val+1]) continue;
                    path[i_val+2] = g.next(sf, path[i_val+1]);
                    tp[i_val+2] = g.next(sf, path[i_val+2]);
                    if (g.next(sf, tp[i_val+2]) != tp[i_val+1]) continue;
                }

                // Validate endpoint degrees
                if (g.degree(path[0]) != 6) continue;
                if (g.degree(tp[i_val+2]) != 6) continue;

                // All vertices distinct
                std::set<node_t> all;
                for (node_t x : sc) all.insert(x);
                for (node_t x : path) all.insert(x);
                for (node_t x : tp) all.insert(x);
                if ((int)all.size() != (int)(sc.size() + path.size() + tp.size()))
                    continue;

                if (!stripTopologyClean(g, path, sc, tp)) continue;

                InvSite site;
                site.kind = Lk(i_val);
                site.dir = dir;
                site.strip = sc;
                site.path = path;
                site.tp = tp;
                result.push_back(site);
            }
        }
    }
    return result;
}

// Find B(0,0) inversion sites.
// Strip = [a(deg5), b(deg6), c(deg5)] where a adj b, b adj c.
// DRight CW: a=[p0,q0,q1,b,p1], b=[a,q1,c,p3,p2,p1], c=[b,q1,q2,p4,p3]
// DLeft CW:  a=[p0,p1,b,q1,q0], b=[a,p1,p2,p3,c,q1], c=[b,p3,p4,q2,q1]
static std::vector<InvSite> findB00InvSites(const Graph& g) {
    std::vector<InvSite> result;
    auto d5 = deg5vertices(g);

    for (node_t a : d5) {
        for (int ni = 0; ni < g.degree(a); ++ni) {
            node_t b = g.nbrs(a)[ni];
            if (g.degree(b) != 6) continue;

            // Find degree-5 neighbors of b (candidates for c)
            for (int bi = 0; bi < g.degree(b); ++bi) {
                node_t c = g.nbrs(b)[bi];
                if (c == a) continue;
                if (g.degree(c) != 5) continue;
                // a must NOT be adjacent to c (otherwise it's L-type)
                if (g.arc_ix(a, c) >= 0) continue;

                for (Dir dir : {Dir::DRight, Dir::DLeft}) {
                    node_t p0, p1, q0, q1_a;
                    // Extract from a's CW (degree 5)
                    if (dir == Dir::DRight) {
                        // a CW: [p0, q0, q1, b, p1]
                        p1   = g.next(a, b);
                        p0   = g.next(a, p1);
                        q0   = g.next(a, p0);
                        q1_a = g.next(a, q0);
                        if (g.next(a, q1_a) != b) continue;
                    } else {
                        // a CW: [p0, p1, b, q1, q0]
                        q1_a = g.next(a, b);
                        q0   = g.next(a, q1_a);
                        p0   = g.next(a, q0);
                        p1   = g.next(a, p0);
                        if (g.next(a, p1) != b) continue;
                    }

                    // Extract from b's CW (degree 6)
                    node_t p2, p3, q1_b;
                    if (dir == Dir::DRight) {
                        // b CW: [a, q1, c, p3, p2, p1]
                        q1_b = g.next(b, a);
                        if (g.next(b, q1_b) != c) continue;
                        p3 = g.next(b, c);
                        p2 = g.next(b, p3);
                        if (g.next(b, p2) != p1) continue;
                    } else {
                        // b CW: [a, p1, p2, p3, c, q1]
                        if (g.next(b, a) != p1) continue;
                        p2 = g.next(b, p1);
                        p3 = g.next(b, p2);
                        if (g.next(b, p3) != c) continue;
                        q1_b = g.next(b, c);
                    }
                    if (q1_a != q1_b) continue;  // q1 must be consistent

                    // Extract from c's CW (degree 5)
                    node_t p4, q2;
                    if (dir == Dir::DRight) {
                        // c CW: [b, q1, q2, p4, p3]
                        if (g.next(c, b) != q1_a) continue;
                        q2 = g.next(c, q1_a);
                        p4 = g.next(c, q2);
                        if (g.next(c, p4) != p3) continue;
                    } else {
                        // c CW: [b, p3, p4, q2, q1]
                        if (g.next(c, b) != p3) continue;
                        p4 = g.next(c, p3);
                        q2 = g.next(c, p4);
                        if (g.next(c, q2) != q1_a) continue;
                    }

                    // Validate: p0 and p4 must be degree 6 (promoted by expansion)
                    if (g.degree(p0) != 6) continue;
                    if (g.degree(p4) != 6) continue;

                    // All 11 vertices must be distinct
                    std::set<node_t> all = {a, b, c, p0, p1, p2, p3, p4,
                                             q0, q1_a, q2};
                    if ((int)all.size() != 11) continue;

                    // Strip topology check
                    std::vector<node_t> strip = {a, b, c};
                    std::vector<node_t> path = {p0, p1, p2, p3, p4};
                    std::vector<node_t> tp = {q0, q1_a, q2};
                    if (!stripTopologyClean(g, path, strip, tp)) continue;

                    InvSite site;
                    site.kind = Bk(0, 0);
                    site.dir = dir;
                    site.strip = strip;
                    site.path = path;
                    site.tp = tp;
                    result.push_back(site);
                }
            }
        }
    }
    return result;
}

// Find F (nanotube ring) inversion sites.
// Detects two parallel 5-cycles of degree-6 vertices (ring and strip).
// strip[i] sits between ring[i] and ring[(i+1)%5].
// outer[i] is the common non-ring non-strip neighbor of strip[i-1] and strip[i].
std::vector<InvSite> findFRingInvSites(const Graph& g) {
    std::vector<InvSite> result;

    // Collect degree-6 vertices
    std::vector<node_t> deg6;
    for (int v = 0; v < g.N; ++v)
        if (g.degree(v) == 6) deg6.push_back(v);
    if ((int)deg6.size() < 10) return result;  // need at least 2 rings of 5

    std::set<std::set<node_t>> seen_rings;

    for (node_t u : deg6) {
        for (int ni = 0; ni < g.degree(u); ++ni) {
            node_t v = g.nbrs(u)[ni];
            if (g.degree(v) != 6) continue;

            // Trace 5-cycle via straight-ahead (advanceCW 3 at degree 6)
            std::vector<node_t> ring = {u, v};
            node_t prv = u, cur = v;
            bool valid = true;
            for (int step = 2; step < 5; ++step) {
                node_t nxt = advanceCW(g, cur, prv, 3);
                if (g.degree(nxt) != 6) { valid = false; break; }
                for (node_t x : ring) if (x == nxt) { valid = false; break; }
                if (!valid) break;
                ring.push_back(nxt);
                prv = cur; cur = nxt;
            }
            if (!valid) continue;

            // Check closure
            if (advanceCW(g, cur, prv, 3) != u) continue;

            // Deduplicate (same vertex set, different starting edge)
            std::set<node_t> ringSet(ring.begin(), ring.end());
            if (!seen_rings.insert(ringSet).second) continue;

            // Find strip: for each ring edge (r_i, r_{i+1}), find common
            // non-ring degree-6 neighbor
            std::vector<node_t> strip(5, -1);
            bool stripOK = true;
            for (int i = 0; i < 5 && stripOK; ++i) {
                int j = (i + 1) % 5;
                for (int k = 0; k < g.degree(ring[i]); ++k) {
                    node_t w = g.nbrs(ring[i])[k];
                    if (ringSet.count(w)) continue;
                    if (g.degree(w) == 6 && g.arc_ix(ring[j], w) >= 0) {
                        strip[i] = w; break;
                    }
                }
                if (strip[i] < 0) stripOK = false;
            }
            if (!stripOK) continue;

            // Verify strip is a 5-cycle of 5 distinct vertices
            std::set<node_t> stripSet(strip.begin(), strip.end());
            if ((int)stripSet.size() != 5) continue;
            for (int i = 0; i < 5; ++i)
                if (g.arc_ix(strip[i], strip[(i+1)%5]) < 0) { stripOK = false; break; }
            if (!stripOK) continue;

            // Find outer vertices: for each strip[i], its 2 non-ring non-strip
            // neighbors are {outer[i], outer[(i+1)%5]}.
            // outer[i] = common element of strip[i-1]'s extras and strip[i]'s extras.
            std::vector<std::vector<node_t>> extras(5);
            bool extrasOK = true;
            for (int i = 0; i < 5 && extrasOK; ++i) {
                for (int k = 0; k < g.degree(strip[i]); ++k) {
                    node_t w = g.nbrs(strip[i])[k];
                    if (!ringSet.count(w) && !stripSet.count(w))
                        extras[i].push_back(w);
                }
                if ((int)extras[i].size() != 2) extrasOK = false;
            }
            if (!extrasOK) continue;

            std::vector<node_t> outer(5, -1);
            bool outerOK = true;
            for (int i = 0; i < 5 && outerOK; ++i) {
                int im1 = (i - 1 + 5) % 5;
                // outer[i] is in both extras[im1] and extras[i]
                outer[i] = -1;
                for (node_t x : extras[im1])
                    for (node_t y : extras[i])
                        if (x == y) { outer[i] = x; goto found_outer; }
                found_outer:
                if (outer[i] < 0) outerOK = false;
            }
            if (!outerOK) continue;

            // All outer must be distinct and disjoint from ring/strip
            std::set<node_t> outerSet(outer.begin(), outer.end());
            if ((int)outerSet.size() != 5) continue;
            for (node_t o : outer)
                if (ringSet.count(o) || stripSet.count(o)) { outerOK = false; break; }
            if (!outerOK) continue;

            // Strip topology: all strip neighbors in {ring, strip, outer}
            if (!stripTopologyClean(g, ring, strip, outer)) continue;

            // Note: we use path=ring, tp=outer for F-type InvSite
            InvSite site;
            site.kind = Fk();
            site.dir = Dir::DRight;  // unused for F-type
            site.strip = strip;
            site.path = ring;
            site.tp = outer;
            result.push_back(site);
        }
    }
    return result;
}

std::vector<InvSite> allInvSites(const Graph& g, int maxRedLen) {
    std::vector<InvSite> result;
    if (maxRedLen >= 1) {
        auto l0 = findL0InvSites(g);
        result.insert(result.end(), l0.begin(), l0.end());
    }
    if (maxRedLen >= 2) {
        auto li = findLiInvSites(g, maxRedLen);
        result.insert(result.end(), li.begin(), li.end());
    }
    if (maxRedLen >= 2) {
        auto b00 = findB00InvSites(g);
        result.insert(result.end(), b00.begin(), b00.end());
    }
    {
        auto fr = findFRingInvSites(g);
        result.insert(result.end(), fr.begin(), fr.end());
    }
    return result;
}

Graph invertReduction(const Graph& g, const InvSite& site) {
    const auto& strip = site.strip;
    const auto& path = site.path;
    const auto& tp = site.tp;
    int numStrip = (int)strip.size();

    // Copy adjacency lists
    adj_t adj(g.N);
    for (int i = 0; i < g.N; ++i) adj[i] = g.neighbours[i];

    if (site.kind.type == ExpKind::L_type) {
        if (!reconnectStraight(adj, path, strip, tp))
            return Graph();
    } else if (site.kind.type == ExpKind::B_type) {
        int bentPos = site.kind.i;
        int bentLen = site.kind.i + site.kind.j;
        if (!reconnectBent(adj, path, strip, tp, bentPos, bentLen))
            return Graph();
    } else if (site.kind.type == ExpKind::F_type) {
        // F-ring reconnection: path=ring, tp=outer
        // ring[i]: replace strip[(i-1+5)%5] → outer[(i-1+5)%5],
        //          replace strip[i] → outer[i]
        // outer[i]: replace strip[(i-1+5)%5] → ring[i],
        //           replace strip[i] → ring[(i+1)%5]
        for (int i = 0; i < 5; ++i) {
            int im1 = (i - 1 + 5) % 5;
            int ip1 = (i + 1) % 5;
            if (!replaceNbrAdj(adj, path[i], strip[im1], tp[im1])) return Graph();
            if (!replaceNbrAdj(adj, path[i], strip[i], tp[i])) return Graph();
            if (!replaceNbrAdj(adj, tp[i], strip[im1], path[i])) return Graph();
            if (!replaceNbrAdj(adj, tp[i], strip[i], path[ip1])) return Graph();
        }
    } else {
        return Graph();
    }

    // Validate post-surgery degrees
    std::set<node_t> stripSet(strip.begin(), strip.end());
    for (int i = 0; i < g.N; ++i) {
        if (stripSet.count(i)) continue;
        int d = (int)adj[i].size();
        if (d != 5 && d != 6) return Graph();
        for (node_t nbr : adj[i])
            if (stripSet.count(nbr)) return Graph();
    }

    // Renumber
    std::vector<int> old_to_new(g.N, -1);
    int new_id = 0;
    for (int i = 0; i < g.N; ++i)
        if (!stripSet.count(i)) old_to_new[i] = new_id++;

    int new_N = g.N - numStrip;
    neighbours_t new_adj(new_N, GRAPH_DMAX);
    for (int i = 0; i < g.N; ++i) {
        if (old_to_new[i] < 0) continue;
        for (node_t nbr : adj[i])
            new_adj[old_to_new[i]].push_back(old_to_new[nbr]);
    }
    return Graph(new_adj);
}

// Apply reconnection for a straight reduction (L_0 and L_i).
// pathlength = strip.size() = i+2
// Returns false if any replace operation fails (topology mismatch).
static bool reconnectStraight(adj_t& adj,
        const std::vector<node_t>& path, const std::vector<node_t>& strip,
        const std::vector<node_t>& tp) {
    int pathlength = (int)strip.size();

    // path[0]: remove strip[0] (degree 6 -> 5)
    removeNbrAdj(adj, path[0], strip[0]);

    // path[k] for k=1..pathlength: replace strip neighbours with true_par
    for (int k = 1; k <= pathlength; ++k) {
        if (!replaceNbrAdj(adj, path[k], strip[k-1], tp[k-1])) return false;
        if (k < pathlength)
            if (!replaceNbrAdj(adj, path[k], strip[k], tp[k])) return false;
    }

    // true_par[k] for k=0..pathlength-1: replace strip with path
    for (int k = 0; k < pathlength; ++k) {
        if (k > 0)
            if (!replaceNbrAdj(adj, tp[k], strip[k-1], path[k])) return false;
        if (!replaceNbrAdj(adj, tp[k], strip[k], path[k+1])) return false;
    }

    // true_par[pathlength]: remove strip[pathlength-1] (degree 6 -> 5)
    removeNbrAdj(adj, tp[pathlength], strip[pathlength-1]);
    return true;
}

// Apply reconnection for a bent reduction (B_{0,0} and B_{i,j}).
// bentPos = i, bentLen = i+j, numNew = strip.size() = i+j+3
// Returns false if any replace operation fails (topology mismatch).
static bool reconnectBent(adj_t& adj,
        const std::vector<node_t>& path, const std::vector<node_t>& strip,
        const std::vector<node_t>& tp, int bentPos, int bentLen) {
    int numNew = (int)strip.size();

    // Endpoint removals (degree 6 -> 5)
    removeNbrAdj(adj, path[0], strip[0]);
    removeNbrAdj(adj, path[bentLen + 4], strip[numNew - 1]);

    // Before-bend path restoration: path[k] for k=1..bentPos+1
    for (int k = 1; k <= bentPos + 1; ++k) {
        if (!replaceNbrAdj(adj, path[k], strip[k-1], tp[k-1])) return false;
        if (!replaceNbrAdj(adj, path[k], strip[k], tp[k])) return false;
    }

    // Before-bend parallel restoration: tp[k] for k=0..bentPos
    for (int k = 0; k <= bentPos; ++k) {
        if (k > 0)
            if (!replaceNbrAdj(adj, tp[k], strip[k-1], path[k])) return false;
        if (!replaceNbrAdj(adj, tp[k], strip[k], path[k+1])) return false;
    }

    // Bend area: bendI = bentPos + 2
    int bendI = bentPos + 2;
    if (!replaceNbrAdj(adj, path[bendI], strip[bendI-1], tp[bendI-1])) return false;
    if (!replaceNbrAdj(adj, tp[bendI-1], strip[bendI-2], path[bendI-1])) return false;
    if (!replaceNbrAdj(adj, tp[bendI-1], strip[bendI-1], path[bendI])) return false;
    if (!replaceNbrAdj(adj, tp[bendI-1], strip[bendI], path[bendI+1])) return false;

    // After-bend path restoration: path[k] for k=bendI+1..bentLen+3
    for (int k = bendI + 1; k <= bentLen + 3; ++k) {
        if (!replaceNbrAdj(adj, path[k], strip[k-2], tp[k-2])) return false;
        if (!replaceNbrAdj(adj, path[k], strip[k-1], tp[k-1])) return false;
    }

    // After-bend parallel restoration: tp[k] for k=bentPos+2..bentLen+2
    for (int k = bentPos + 2; k <= bentLen + 2; ++k) {
        if (!replaceNbrAdj(adj, tp[k], strip[k], path[k+1])) return false;
        if (k < bentLen + 2)
            if (!replaceNbrAdj(adj, tp[k], strip[k+1], path[k+2])) return false;
    }
    return true;
}

Graph applyReduction(const Graph& g, const Reduction& red) {
    const auto& kind = red.kind;
    node_t u = red.u, v = red.v;
    Dir dir = red.dir;

    std::vector<node_t> path, strip, tp;
    int numStrip;

    if (kind.type == ExpKind::L_type) {
        // Straight reduction L_i: pathlength = i+2, numEntries = i+3
        int numEntries = kind.i + 3;
        int pathlength = kind.i + 2;
        auto pi = computeStraightPath(g, u, v, dir, numEntries);
        if (!pi.valid) return Graph();

        path = std::move(pi.path);
        // strip = par[0..pathlength-1], far_par = par[pathlength]
        strip.assign(pi.parallel.begin(), pi.parallel.begin() + pathlength);
        node_t far_par = pi.parallel[pathlength];
        numStrip = pathlength;

        tp = findTrueParStraight(g, path, strip, far_par);
        if (tp.empty()) return Graph();

        // Validate strip topology: all strip vertex neighbors must be in
        // {path, strip, tp}. If not, this site doesn't have clean expansion
        // topology and can't be surgically inverted.
        if (!stripTopologyClean(g, path, strip, tp))
            return Graph();

        // Copy adjacency and reconnect
        adj_t adj(g.N);
        for (int i = 0; i < g.N; ++i) adj[i] = g.neighbours[i];
        if (!reconnectStraight(adj, path, strip, tp))
            return Graph();

        // Validate post-surgery degrees and build result
        std::set<node_t> stripSet(strip.begin(), strip.end());
        for (int i = 0; i < g.N; ++i) {
            if (stripSet.count(i)) continue;
            int d = (int)adj[i].size();
            if (d != 5 && d != 6) return Graph();
            // Check no remaining strip references
            for (node_t nbr : adj[i])
                if (stripSet.count(nbr)) return Graph();
        }

        // Renumber
        std::vector<int> old_to_new(g.N, -1);
        int new_id = 0;
        for (int i = 0; i < g.N; ++i)
            if (!stripSet.count(i)) old_to_new[i] = new_id++;

        int new_N = g.N - numStrip;
        neighbours_t new_adj(new_N, GRAPH_DMAX);
        for (int i = 0; i < g.N; ++i) {
            if (old_to_new[i] < 0) continue;
            for (node_t nbr : adj[i])
                new_adj[old_to_new[i]].push_back(old_to_new[nbr]);
        }
        return Graph(new_adj);

    } else if (kind.type == ExpKind::B_type) {
        // Bent reduction B_{i,j}
        int bi = kind.i, bj = kind.j;
        int bentPos = bi;
        int bentLen = bi + bj;
        int numNew = bentLen + 3;

        PathInfo pi;
        if (bi == 0 && bj == 0)
            pi = computeBentZeroPath(g, u, v, dir);
        else
            pi = computeBentPath(g, u, v, dir, bi, bj);
        if (!pi.valid) return Graph();

        path = std::move(pi.path);
        strip = std::move(pi.parallel); // all parallel entries are strip vertices
        numStrip = numNew;
        assert((int)strip.size() == numNew);

        tp = findTrueParBent(g, path, strip, bentPos);
        if (tp.empty()) return Graph();

        // Validate strip topology
        if (!stripTopologyClean(g, path, strip, tp))
            return Graph();

        // Copy adjacency and reconnect
        adj_t adj(g.N);
        for (int i = 0; i < g.N; ++i) adj[i] = g.neighbours[i];
        if (!reconnectBent(adj, path, strip, tp, bentPos, bentLen))
            return Graph();

        // Validate post-surgery degrees and build result
        std::set<node_t> stripSet(strip.begin(), strip.end());
        for (int i = 0; i < g.N; ++i) {
            if (stripSet.count(i)) continue;
            int d = (int)adj[i].size();
            if (d != 5 && d != 6) return Graph();
            for (node_t nbr : adj[i])
                if (stripSet.count(nbr)) return Graph();
        }

        // Renumber
        std::vector<int> old_to_new(g.N, -1);
        int new_id = 0;
        for (int i = 0; i < g.N; ++i)
            if (!stripSet.count(i)) old_to_new[i] = new_id++;

        int new_N = g.N - numStrip;
        neighbours_t new_adj(new_N, GRAPH_DMAX);
        for (int i = 0; i < g.N; ++i) {
            if (old_to_new[i] < 0) continue;
            for (node_t nbr : adj[i])
                new_adj[old_to_new[i]].push_back(old_to_new[nbr]);
        }
        return Graph(new_adj);

    } else {
        assert(false && "F-type reductions not supported");
        return Graph();
    }
}

// =====================================================================
// Reduction enumeration (Phase 2)
// =====================================================================

// Follow straight-ahead from u through v until a degree-5 vertex is found.
// Direction is irrelevant at degree-6 interior vertices (straightAhead gives
// the same result for both directions), so we use DLeft arbitrarily.
std::optional<StraightEndpoint> followStraightToFive(
        const Graph& g, node_t u, node_t v, int maxDist) {
    std::set<node_t> visited;
    visited.insert(u);
    node_t from = u, to = v;
    for (int dist = 1; dist <= maxDist; ++dist) {
        if (visited.count(to)) return std::nullopt;  // cycle
        if (g.degree(to) == 5) return StraightEndpoint{to, from, dist};
        visited.insert(to);
        node_t next = straightAhead(g, Dir::DLeft, to, from);
        from = to;
        to = next;
    }
    return std::nullopt;
}

// L0 validity: flanking vertices 2 positions from the edge at each endpoint
// must be degree-6.
static bool isValidL0Direction(const Graph& g, node_t u, node_t v, Dir dir) {
    if (dir == Dir::DRight)
        return g.degree(advanceCW(g, u, v, 2)) == 6
            && g.degree(advanceCW(g, v, u, 2)) == 6;
    else
        return g.degree(advanceCW(g, u, v, g.degree(u) - 2)) == 6
            && g.degree(advanceCW(g, v, u, g.degree(v) - 2)) == 6;
}

// Straight reduction validity: same flanking check as L0, at both endpoints.
static bool isValidStraightReduction(const Graph& g,
        node_t u, node_t v, node_t w, node_t prevW, Dir dir) {
    if (dir == Dir::DRight)
        return g.degree(advanceCW(g, u, v, 2)) == 6
            && g.degree(advanceCW(g, w, prevW, 2)) == 6;
    else
        return g.degree(advanceCW(g, u, v, g.degree(u) - 2)) == 6
            && g.degree(advanceCW(g, w, prevW, g.degree(w) - 2)) == 6;
}

// Bent endpoint flanking: OPPOSITE side from straight reductions.
static node_t bentEndpointFlank(const Graph& g, node_t w, node_t prevW, Dir dir) {
    if (dir == Dir::DLeft)
        return advanceCW(g, w, prevW, 2);
    else
        return advanceCW(g, w, prevW, g.degree(w) - 2);
}

// Safe bent path computation for validity checking.
// Returns nullopt on self-intersection. Computes uniform parallel path
// (no skipping at the bend), unlike computeBentPath which is for surgery.
static std::optional<PathInfo> computeBentPathSafe(
        const Graph& g, node_t u, node_t v, Dir dir, int bentPos, int bentLen) {
    int pathLen = bentLen + 5;
    std::vector<node_t> path(pathLen);
    path[0] = u;
    path[1] = v;

    // Pre-bend: straight steps to path[bentPos+2]
    for (int k = 2; k <= bentPos + 2; ++k)
        path[k] = straightAhead(g, dir, path[k-1], path[k-2]);

    // Turn at the bend vertex path[bentPos+2]
    path[bentPos + 3] = turnAhead(g, dir, path[bentPos+2], path[bentPos+1]);

    // Post-bend: straight steps
    for (int k = bentPos + 4; k < pathLen; ++k)
        path[k] = straightAhead(g, dir, path[k-1], path[k-2]);

    // Check for self-intersection in path
    std::set<node_t> seen(path.begin(), path.end());
    if ((int)seen.size() != pathLen) return std::nullopt;

    // Uniform parallel: sideNbr for ALL consecutive pairs k=0..bentLen+2
    int parLen = bentLen + 3;
    std::vector<node_t> par(parLen);
    for (int k = 0; k < parLen; ++k)
        par[k] = sideNbr(g, dir, path[k], path[k+1]);

    PathInfo pi;
    pi.path = std::move(path);
    pi.parallel = std::move(par);
    pi.valid = true;
    return pi;
}

// B00 validity: path non-self-intersecting, far endpoint degree-5,
// first != last, path and parallel disjoint.
static bool isValidB00Reduction(const Graph& g, node_t u, node_t v, Dir dir) {
    PathInfo pi = computeBentZeroPath(g, u, v, dir);
    if (pi.path.empty()) return false;

    // Check path has no duplicates
    std::set<node_t> pathSet(pi.path.begin(), pi.path.end());
    if ((int)pathSet.size() != (int)pi.path.size()) return false;

    // Far endpoint must be degree-5
    if (g.degree(pi.path.back()) != 5) return false;

    // Must not be a loop
    if (pi.path.front() == pi.path.back()) return false;

    // Path and parallel must be disjoint
    for (node_t x : pi.parallel)
        if (pathSet.count(x)) return false;

    return true;
}

// Bent reduction validity (i+j > 0): path non-self-intersecting, all
// intermediate vertices degree-6, far endpoint degree-5, flanking degree-6.
static bool isValidBentReduction(const Graph& g, node_t u, node_t v,
                                  Dir dir, int i, int j) {
    auto opt = computeBentPathSafe(g, u, v, dir, i, i + j);
    if (!opt) return false;

    const auto& path = opt->path;
    int n = (int)path.size();

    // Far endpoint must be degree-5
    if (g.degree(path[n-1]) != 5) return false;

    // Must not be a loop
    if (path[0] == path[n-1]) return false;

    // All intermediate vertices must be degree-6
    for (int k = 1; k < n - 1; ++k)
        if (g.degree(path[k]) != 6) return false;

    // Flanking vertex at far endpoint must be degree-6 and not on path
    node_t flank = bentEndpointFlank(g, path[n-1], path[n-2], dir);
    if (g.degree(flank) != 6) return false;
    for (int k = 0; k < n; ++k)
        if (path[k] == flank) return false;

    return true;
}

std::vector<Reduction> allReductions(const Graph& g, int maxRedLen) {
    std::vector<Reduction> result;
    auto d5 = deg5vertices(g);

    // L0 reductions
    for (node_t u : d5) {
        for (int ni = 0; ni < g.degree(u); ++ni) {
            node_t v = g.nbrs(u)[ni];
            if (g.degree(v) != 5) continue;
            for (Dir d : {Dir::DLeft, Dir::DRight}) {
                if (!isValidL0Direction(g, u, v, d)) continue;
                // Check path/parallel disjointness (matching Haskell isValidStraightSite)
                auto pi = computeStraightPath(g, u, v, d, 3);
                if (!pi.valid) continue;
                result.push_back({Lk(0), u, v, d});
            }
        }
    }

    // Straight reductions L_i (i >= 1)
    if (maxRedLen >= 2) {
        for (node_t u : d5) {
            for (int ni = 0; ni < g.degree(u); ++ni) {
                node_t v = g.nbrs(u)[ni];
                if (g.degree(v) == 5) continue;  // L0 handled above
                auto ep = followStraightToFive(g, u, v, maxRedLen);
                if (!ep) continue;
                for (Dir d : {Dir::DLeft, Dir::DRight}) {
                    if (!isValidStraightReduction(g, u, v, ep->endpoint, ep->prev, d))
                        continue;
                    int i = ep->distance - 1;
                    // Check path/parallel disjointness
                    auto pi = computeStraightPath(g, u, v, d, i + 3);
                    if (!pi.valid) continue;
                    result.push_back({Lk(i), u, v, d});
                }
            }
        }
    }

    // B(0,0) reductions
    if (maxRedLen >= 2) {
        for (node_t u : d5) {
            for (int ni = 0; ni < g.degree(u); ++ni) {
                node_t v = g.nbrs(u)[ni];
                for (Dir d : {Dir::DLeft, Dir::DRight}) {
                    if (isValidB00Reduction(g, u, v, d))
                        result.push_back({Bk(0, 0), u, v, d});
                }
            }
        }
    }

    // Bent reductions B(i,j) with i+j > 0
    if (maxRedLen >= 3) {
        for (node_t u : d5) {
            for (int ni = 0; ni < g.degree(u); ++ni) {
                node_t v = g.nbrs(u)[ni];
                for (Dir d : {Dir::DLeft, Dir::DRight}) {
                    for (int bentLen = 1; bentLen <= maxRedLen - 2; ++bentLen) {
                        for (int bentPos = 0; bentPos <= bentLen; ++bentPos) {
                            if (isValidBentReduction(g, u, v, d, bentPos, bentLen - bentPos))
                                result.push_back({Bk(bentPos, bentLen - bentPos), u, v, d});
                        }
                    }
                }
            }
        }
    }

    return result;
}

// =====================================================================
// ReducibleDual implementation
// =====================================================================

ReducibleDual::ReducibleDual(const Graph& g)
    : V(g.N), n_alive(g.N)
{
    for (int u = 0; u < g.N; u++) {
        int d = g.degree(u);
        assert(d <= D_MAX);
        V[u].active = (1u << d) - 1;
        for (int i = 0; i < d; i++)
            V[u].nbr[i] = g.nbrs(u)[i];
        if (d == 5) deg5.insert(u);
    }
}

ReducibleDual::ReducibleDual(int capacity)
    : V(capacity), n_alive(0)
{
    for (auto& v : V) v.active = 0;
}

int ReducibleDual::arc_ix(node_t u, node_t v) const {
    uint8_t m = V[u].active;
    for (; m; m &= m - 1) {
        int i = __builtin_ctz(m);
        if (V[u].nbr[i] == v) return i;
    }
    return -1;
}

node_t ReducibleDual::next(node_t u, node_t v) const {
    int pos = arc_ix(u, v);
    if (pos < 0) return -1;
    for (int step = 1; step <= D_MAX; step++) {
        int p = (pos + step) % D_MAX;
        if (V[u].active & (1u << p)) return V[u].nbr[p];
    }
    return -1;
}

node_t ReducibleDual::prev(node_t u, node_t v) const {
    int pos = arc_ix(u, v);
    if (pos < 0) return -1;
    for (int step = 1; step <= D_MAX; step++) {
        int p = (pos - step + D_MAX) % D_MAX;
        if (V[u].active & (1u << p)) return V[u].nbr[p];
    }
    return -1;
}

node_t ReducibleDual::advanceCW(node_t u, node_t v, int k) const {
    node_t cur = v;
    for (int i = 0; i < k; i++) {
        cur = next(u, cur);
        if (cur < 0) return -1;
    }
    return cur;
}

node_t ReducibleDual::straightAhead(Dir dir, node_t u, node_t from) const {
    int d = degree(u);
    if (d == 6) return advanceCW(u, from, 3);
    if (dir == Dir::DRight) return advanceCW(u, from, 3);
    else return advanceCW(u, from, 2);
}

void ReducibleDual::splice(node_t u, node_t old_v, node_t new_v) {
    uint8_t m = V[u].active;
    for (; m; m &= m - 1) {
        int i = __builtin_ctz(m);
        if (V[u].nbr[i] == old_v) { V[u].nbr[i] = new_v; return; }
    }
    assert(false && "splice: old_v not found in active neighbours of u");
}

void ReducibleDual::snip(node_t u, node_t v) {
    uint8_t m = V[u].active;
    for (; m; m &= m - 1) {
        int i = __builtin_ctz(m);
        if (V[u].nbr[i] == v) { V[u].active &= ~(1u << i); return; }
    }
    assert(false && "snip: v not found in active neighbours of u");
}

void ReducibleDual::kill(node_t v) {
    V[v].active = 0;
    n_alive--;
    deg5.erase(v);
}

void ReducibleDual::unsnip(node_t u, node_t v) {
    for (int i = 0; i < D_MAX; i++) {
        if (!(V[u].active & (1u << i)) && V[u].nbr[i] == v) {
            V[u].active |= (1u << i);
            return;
        }
    }
    assert(false && "unsnip: v not found in inactive positions of u");
}

void ReducibleDual::insertAfter(node_t u, node_t v, node_t after) {
    // Insert v into u's CW ordering right after 'after'.
    // u must be degree 5 (one inactive position). Promotes u to degree 6.
    assert(degree(u) == 5);

    // Extract current 5 active neighbors in CW order (by physical position)
    node_t cw[5]; int ci = 0;
    for (int p = 0; p < D_MAX; p++)
        if (V[u].active & (1u << p))
            cw[ci++] = V[u].nbr[p];
    assert(ci == 5);

    // Find 'after' in CW sequence
    int ins = -1;
    for (int i = 0; i < 5; i++)
        if (cw[i] == after) { ins = i; break; }
    assert(ins >= 0 && "insertAfter: 'after' not found in CW of u");

    // Build new 6-element CW: insert v right after 'after'
    node_t new_cw[6];
    for (int i = 0; i <= ins; i++) new_cw[i] = cw[i];
    new_cw[ins + 1] = v;
    for (int i = ins + 1; i < 5; i++) new_cw[i + 1] = cw[i];

    for (int i = 0; i < D_MAX; i++) V[u].nbr[i] = new_cw[i];
    V[u].active = 0x3f;  // all 6 bits active
}

void ReducibleDual::unkill(node_t v, const Vertex& saved) {
    V[v] = saved;
    n_alive++;
    if (degree(v) == 5) deg5.insert(v);
}

// --- Reconnection ---

void ReducibleDual::reconnectStraight(
        const std::vector<node_t>& path,
        const std::vector<node_t>& strip,
        const std::vector<node_t>& tp) {
    int n = (int)strip.size();

    snip(path[0], strip[0]);

    for (int k = 1; k <= n; k++) {
        splice(path[k], strip[k-1], tp[k-1]);
        if (k < n) splice(path[k], strip[k], tp[k]);
    }

    for (int k = 0; k < n; k++) {
        if (k > 0) splice(tp[k], strip[k-1], path[k]);
        splice(tp[k], strip[k], path[k+1]);
    }

    snip(tp[n], strip[n-1]);
}

void ReducibleDual::reconnectBent(
        const std::vector<node_t>& path,
        const std::vector<node_t>& strip,
        const std::vector<node_t>& tp,
        int bentPos, int bentLen) {
    int numNew = (int)strip.size();

    snip(path[0], strip[0]);
    snip(path[bentLen + 4], strip[numNew - 1]);

    // Before bend
    for (int k = 1; k <= bentPos + 1; k++) {
        splice(path[k], strip[k-1], tp[k-1]);
        splice(path[k], strip[k], tp[k]);
    }
    for (int k = 0; k <= bentPos; k++) {
        if (k > 0) splice(tp[k], strip[k-1], path[k]);
        splice(tp[k], strip[k], path[k+1]);
    }

    // Bend
    int bi = bentPos + 2;
    splice(path[bi], strip[bi-1], tp[bi-1]);
    splice(tp[bi-1], strip[bi-2], path[bi-1]);
    splice(tp[bi-1], strip[bi-1], path[bi]);
    splice(tp[bi-1], strip[bi], path[bi+1]);

    // After bend
    for (int k = bi + 1; k <= bentLen + 3; k++) {
        splice(path[k], strip[k-2], tp[k-2]);
        splice(path[k], strip[k-1], tp[k-1]);
    }
    for (int k = bentPos + 2; k <= bentLen + 2; k++) {
        splice(tp[k], strip[k], path[k+1]);
        if (k < bentLen + 2) splice(tp[k], strip[k+1], path[k+2]);
    }
}

void ReducibleDual::reconnectRing(
        const std::vector<node_t>& ring,
        const std::vector<node_t>& strip,
        const std::vector<node_t>& outer) {
    for (int i = 0; i < 5; i++) {
        int im1 = (i + 4) % 5, ip1 = (i + 1) % 5;
        splice(ring[i], strip[im1], outer[im1]);
        splice(ring[i], strip[i], outer[i]);
        splice(outer[i], strip[im1], ring[i]);
        splice(outer[i], strip[i], ring[ip1]);
    }
}

// --- Expansion reconnection (reverse of reduction) ---

void ReducibleDual::expandStraight(
        const std::vector<node_t>& path,
        const std::vector<node_t>& strip,
        const std::vector<node_t>& tp) {
    int n = (int)strip.size();

    // Reverse of reconnectStraight, undoing each operation in reverse order.
    // reconnectStraight does:
    //   snip(path[0], strip[0])
    //   for k=1..n: splice(path[k], strip[k-1], tp[k-1])
    //              splice(path[k], strip[k], tp[k])     [k<n]
    //   for k=0..n-1: splice(tp[k], strip[k-1], path[k])  [k>0]
    //                 splice(tp[k], strip[k], path[k+1])
    //   snip(tp[n], strip[n-1])

    unsnip(tp[n], strip[n-1]);

    for (int k = n - 1; k >= 0; k--) {
        splice(tp[k], path[k+1], strip[k]);
        if (k > 0) splice(tp[k], path[k], strip[k-1]);
    }

    for (int k = n; k >= 1; k--) {
        if (k < n) splice(path[k], tp[k], strip[k]);
        splice(path[k], tp[k-1], strip[k-1]);
    }

    unsnip(path[0], strip[0]);
}

void ReducibleDual::expandBent(
        const std::vector<node_t>& path,
        const std::vector<node_t>& strip,
        const std::vector<node_t>& tp,
        int bentPos, int bentLen) {
    int numNew = (int)strip.size();

    // Reverse of reconnectBent, undoing each operation in reverse order.
    // reconnectBent does:
    //   snip(path[0], strip[0])
    //   snip(path[bentLen+4], strip[numNew-1])
    //   before-bend path: k=1..bentPos+1
    //   before-bend tp: k=0..bentPos
    //   bend: bi=bentPos+2
    //   after-bend path: k=bi+1..bentLen+3
    //   after-bend tp: k=bentPos+2..bentLen+2

    int bi = bentPos + 2;

    // Undo after-bend tp (last operations first)
    for (int k = bentLen + 2; k >= bentPos + 2; k--) {
        if (k < bentLen + 2) splice(tp[k], path[k+2], strip[k+1]);
        splice(tp[k], path[k+1], strip[k]);
    }

    // Undo after-bend path
    for (int k = bentLen + 3; k >= bi + 1; k--) {
        splice(path[k], tp[k-1], strip[k-1]);
        splice(path[k], tp[k-2], strip[k-2]);
    }

    // Undo bend
    splice(tp[bi-1], path[bi+1], strip[bi]);
    splice(tp[bi-1], path[bi], strip[bi-1]);
    splice(tp[bi-1], path[bi-1], strip[bi-2]);
    splice(path[bi], tp[bi-1], strip[bi-1]);

    // Undo before-bend tp
    for (int k = bentPos; k >= 0; k--) {
        splice(tp[k], path[k+1], strip[k]);
        if (k > 0) splice(tp[k], path[k], strip[k-1]);
    }

    // Undo before-bend path
    for (int k = bentPos + 1; k >= 1; k--) {
        splice(path[k], tp[k], strip[k]);
        splice(path[k], tp[k-1], strip[k-1]);
    }

    // Undo endpoint snips
    unsnip(path[bentLen + 4], strip[numNew - 1]);
    unsnip(path[0], strip[0]);
}

void ReducibleDual::expandRing(
        const std::vector<node_t>& ring,
        const std::vector<node_t>& strip,
        const std::vector<node_t>& outer) {
    // Reverse of reconnectRing. reconnectRing does for i=0..4:
    //   splice(ring[i], strip[im1], outer[im1])
    //   splice(ring[i], strip[i], outer[i])
    //   splice(outer[i], strip[im1], ring[i])
    //   splice(outer[i], strip[i], ring[ip1])
    // Undo in reverse order (i=4..0, operations reversed within each i):
    for (int i = 4; i >= 0; i--) {
        int im1 = (i + 4) % 5, ip1 = (i + 1) % 5;
        splice(outer[i], ring[ip1], strip[i]);
        splice(outer[i], ring[i], strip[im1]);
        splice(ring[i], outer[i], strip[i]);
        splice(ring[i], outer[im1], strip[im1]);
    }
}

// --- Expand (inverse of reduce) ---

void ReducibleDual::expand(const ReductionStep& step) {
    const auto& site = step.site;

    // 1. Undo deg5 updates from reduce()
    if (site.kind.type != ExpKind::F_type) {
        deg5.erase(site.path[0]);
        deg5.erase(site.tp.back());
    }

    // 2. Unkill strip vertices (restore saved state)
    for (auto& [id, saved] : step.saved)
        unkill(id, saved);

    // 3. Undo reconnection
    if (site.kind.type == ExpKind::L_type)
        expandStraight(site.path, site.strip, site.tp);
    else if (site.kind.type == ExpKind::B_type)
        expandBent(site.path, site.strip, site.tp,
                   site.kind.i, site.kind.i + site.kind.j);
    else
        expandRing(site.path, site.strip, site.tp);
}

// --- Standalone expand (from ExtensionStep, no saved state) ---
// Uses insertAfter instead of unsnip for snipped endpoints, because
// strip vertices created from CW formulas don't have inactive positions.

void ReducibleDual::expand(const ExtensionStep& step) {
    const auto& kind = step.kind;
    const auto& dir = step.dir;
    const auto& strip = step.strip;
    const auto& path = step.path;
    const auto& tp = step.tp;

    // 1. Undo deg5 updates from the original reduce()
    if (kind.type != ExpKind::F_type) {
        deg5.erase(path[0]);
        deg5.erase(tp.back());
    }

    // 2. Create strip vertices with computed CW adjacency
    int n = (int)strip.size();

    if (kind.type == ExpKind::L_type) {
        for (int si = 0; si < n; si++) {
            node_t v = strip[si];
            if (si == 0) {
                if (dir == Dir::DRight) {
                    V[v].nbr[0]=path[0]; V[v].nbr[1]=tp[0]; V[v].nbr[2]=tp[1];
                    V[v].nbr[3]=strip[1]; V[v].nbr[4]=path[1];
                } else {
                    V[v].nbr[0]=path[0]; V[v].nbr[1]=path[1]; V[v].nbr[2]=strip[1];
                    V[v].nbr[3]=tp[1]; V[v].nbr[4]=tp[0];
                }
                V[v].active = 0x1f; n_alive++; deg5.insert(v);
            } else if (si == n - 1) {
                if (dir == Dir::DRight) {
                    V[v].nbr[0]=strip[si-1]; V[v].nbr[1]=tp[si]; V[v].nbr[2]=tp[si+1];
                    V[v].nbr[3]=path[si+1]; V[v].nbr[4]=path[si];
                } else {
                    V[v].nbr[0]=strip[si-1]; V[v].nbr[1]=path[si]; V[v].nbr[2]=path[si+1];
                    V[v].nbr[3]=tp[si+1]; V[v].nbr[4]=tp[si];
                }
                V[v].active = 0x1f; n_alive++; deg5.insert(v);
            } else {
                if (dir == Dir::DRight) {
                    V[v].nbr[0]=strip[si-1]; V[v].nbr[1]=tp[si]; V[v].nbr[2]=tp[si+1];
                    V[v].nbr[3]=strip[si+1]; V[v].nbr[4]=path[si+1]; V[v].nbr[5]=path[si];
                } else {
                    V[v].nbr[0]=strip[si-1]; V[v].nbr[1]=path[si]; V[v].nbr[2]=path[si+1];
                    V[v].nbr[3]=strip[si+1]; V[v].nbr[4]=tp[si+1]; V[v].nbr[5]=tp[si];
                }
                V[v].active = 0x3f; n_alive++;
            }
        }

        // 3. Fix path/tp vertices: insertAfter for endpoints, splice for interior.
        // Insert strip[n-1] into tp[n] (5→6)
        if (dir == Dir::DRight)
            insertAfter(tp[n], strip[n-1], path[n]);
        else
            insertAfter(tp[n], strip[n-1], tp[n-1]);

        // Reverse tp splices (same as expandStraight but without unsnip)
        for (int k = n - 1; k >= 0; k--) {
            splice(tp[k], path[k+1], strip[k]);
            if (k > 0) splice(tp[k], path[k], strip[k-1]);
        }

        // Reverse path splices
        for (int k = n; k >= 1; k--) {
            if (k < n) splice(path[k], tp[k], strip[k]);
            splice(path[k], tp[k-1], strip[k-1]);
        }

        // Insert strip[0] into path[0] (5→6)
        if (dir == Dir::DRight)
            insertAfter(path[0], strip[0], tp[0]);
        else
            insertAfter(path[0], strip[0], path[1]);

    } else if (kind.type == ExpKind::B_type) {
        assert(kind.i == 0 && kind.j == 0 && n == 3);
        node_t a = strip[0], b = strip[1], c = strip[2];
        int bentPos = 0, bentLen = 0, numNew = 3;

        // Create a, b, c
        if (dir == Dir::DRight) {
            V[a].nbr[0]=path[0]; V[a].nbr[1]=tp[0]; V[a].nbr[2]=tp[1];
            V[a].nbr[3]=b; V[a].nbr[4]=path[1];
        } else {
            V[a].nbr[0]=path[0]; V[a].nbr[1]=path[1]; V[a].nbr[2]=b;
            V[a].nbr[3]=tp[1]; V[a].nbr[4]=tp[0];
        }
        V[a].active = 0x1f; n_alive++; deg5.insert(a);

        if (dir == Dir::DRight) {
            V[b].nbr[0]=a; V[b].nbr[1]=tp[1]; V[b].nbr[2]=c;
            V[b].nbr[3]=path[3]; V[b].nbr[4]=path[2]; V[b].nbr[5]=path[1];
        } else {
            V[b].nbr[0]=a; V[b].nbr[1]=path[1]; V[b].nbr[2]=path[2];
            V[b].nbr[3]=path[3]; V[b].nbr[4]=c; V[b].nbr[5]=tp[1];
        }
        V[b].active = 0x3f; n_alive++;

        if (dir == Dir::DRight) {
            V[c].nbr[0]=b; V[c].nbr[1]=tp[1]; V[c].nbr[2]=tp[2];
            V[c].nbr[3]=path[4]; V[c].nbr[4]=path[3];
        } else {
            V[c].nbr[0]=b; V[c].nbr[1]=path[3]; V[c].nbr[2]=path[4];
            V[c].nbr[3]=tp[2]; V[c].nbr[4]=tp[1];
        }
        V[c].active = 0x1f; n_alive++; deg5.insert(c);

        // 3. Fix path/tp: insertAfter for endpoints, splice for interior.
        // Insert strip[numNew-1] = c into path[bentLen+4] = path[4] (5→6)
        if (dir == Dir::DRight)
            insertAfter(path[bentLen + 4], strip[numNew - 1], path[bentLen + 3]);
        else
            insertAfter(path[bentLen + 4], strip[numNew - 1], tp[numNew - 1]);

        // Reverse after-bend tp splices
        int bi = bentPos + 2;
        for (int k = bentLen + 2; k >= bentPos + 2; k--) {
            if (k < bentLen + 2) splice(tp[k], path[k+2], strip[k+1]);
            splice(tp[k], path[k+1], strip[k]);
        }

        // Reverse after-bend path splices
        for (int k = bentLen + 3; k >= bi + 1; k--) {
            splice(path[k], tp[k-1], strip[k-1]);
            splice(path[k], tp[k-2], strip[k-2]);
        }

        // Reverse bend
        splice(tp[bi-1], path[bi+1], strip[bi]);
        splice(tp[bi-1], path[bi], strip[bi-1]);
        splice(tp[bi-1], path[bi-1], strip[bi-2]);
        splice(path[bi], tp[bi-1], strip[bi-1]);

        // Reverse before-bend tp splices
        for (int k = bentPos; k >= 0; k--) {
            splice(tp[k], path[k+1], strip[k]);
            if (k > 0) splice(tp[k], path[k], strip[k-1]);
        }

        // Reverse before-bend path splices
        for (int k = bentPos + 1; k >= 1; k--) {
            splice(path[k], tp[k], strip[k]);
            splice(path[k], tp[k-1], strip[k-1]);
        }

        // Insert strip[0] = a into path[0] (5→6)
        if (dir == Dir::DRight)
            insertAfter(path[0], strip[0], tp[0]);
        else
            insertAfter(path[0], strip[0], path[1]);

    } else {
        // F-ring: detect chirality, create strip, use expandRing (no unsnip needed)
        node_t nxt = next(path[0], tp[0]);
        bool orderingA = (nxt == tp[4]);

        for (int i = 0; i < 5; i++) {
            int im1 = (i + 4) % 5, ip1 = (i + 1) % 5;
            node_t v = strip[i];
            if (orderingA) {
                V[v].nbr[0]=path[i]; V[v].nbr[1]=path[ip1]; V[v].nbr[2]=strip[ip1];
                V[v].nbr[3]=tp[ip1]; V[v].nbr[4]=tp[i]; V[v].nbr[5]=strip[im1];
            } else {
                V[v].nbr[0]=path[i]; V[v].nbr[1]=strip[im1]; V[v].nbr[2]=tp[i];
                V[v].nbr[3]=tp[ip1]; V[v].nbr[4]=strip[ip1]; V[v].nbr[5]=path[ip1];
            }
            V[v].active = 0x3f; n_alive++;
        }

        expandRing(path, strip, tp);
    }
}

// --- Reduce ---

void ReducibleDual::reduce(const InvSite& site) {
    if (site.kind.type == ExpKind::L_type)
        reconnectStraight(site.path, site.strip, site.tp);
    else if (site.kind.type == ExpKind::B_type)
        reconnectBent(site.path, site.strip, site.tp,
                      site.kind.i, site.kind.i + site.kind.j);
    else
        reconnectRing(site.path, site.strip, site.tp);

    for (node_t s : site.strip) kill(s);

    // Update deg5: snip reduces degree by 1 (6→5 at path endpoints)
    if (site.kind.type != ExpKind::F_type) {
        deg5.insert(site.path[0]);
        deg5.insert(site.tp.back());
    }
}

// --- Site Finding ---

static bool rdStripTopologyClean(const ReducibleDual& g,
        const std::vector<node_t>& path,
        const std::vector<node_t>& strip,
        const std::vector<node_t>& tp) {
    std::set<node_t> valid;
    for (node_t x : path) valid.insert(x);
    for (node_t x : strip) valid.insert(x);
    for (node_t x : tp) valid.insert(x);

    for (node_t s : strip) {
        uint8_t m = g.V[s].active;
        for (; m; m &= m - 1) {
            node_t n = g.V[s].nbr[__builtin_ctz(m)];
            if (!valid.count(n)) return false;
        }
    }
    return true;
}

std::optional<InvSite> ReducibleDual::findL0Site() const {
    for (node_t s0 : deg5) {
        uint8_t m = V[s0].active;
        for (; m; m &= m - 1) {
            node_t s1 = V[s0].nbr[__builtin_ctz(m)];
            if (degree(s1) != 5 || s1 <= s0) continue;

            for (Dir dir : {Dir::DRight, Dir::DLeft}) {
                node_t path_1, tp_1, path_0, tp_0;
                if (dir == Dir::DRight) {
                    path_1 = next(s0, s1);
                    tp_1   = prev(s0, s1);
                    path_0 = next(s0, path_1);
                    tp_0   = prev(s0, tp_1);
                } else {
                    path_1 = prev(s0, s1);
                    tp_1   = next(s0, s1);
                    path_0 = prev(s0, path_1);
                    tp_0   = next(s0, tp_1);
                }

                node_t path_2, tp_2;
                if (dir == Dir::DRight) {
                    path_2 = prev(s1, path_1);
                    tp_2   = prev(s1, path_2);
                } else {
                    path_2 = next(s1, path_1);
                    tp_2   = next(s1, path_2);
                }

                if (degree(path_0) != 6 || degree(tp_2) != 6) continue;
                if (arc_ix(path_0, path_1) < 0 || arc_ix(path_2, path_1) < 0) continue;

                std::set<node_t> all = {s0, s1, path_0, path_1, path_2, tp_0, tp_1, tp_2};
                if ((int)all.size() != 8) continue;

                InvSite site;
                site.kind = Lk(0); site.dir = dir;
                site.strip = {s0, s1};
                site.path = {path_0, path_1, path_2};
                site.tp = {tp_0, tp_1, tp_2};
                if (!rdStripTopologyClean(*this, site.path, site.strip, site.tp)) continue;
                return site;
            }
        }
    }
    return std::nullopt;
}

std::optional<InvSite> ReducibleDual::findLiSite(int maxRedLen) const {
    for (node_t s0 : deg5) {
        uint8_t m0 = V[s0].active;
        for (; m0; m0 &= m0 - 1) {
            node_t s1 = V[s0].nbr[__builtin_ctz(m0)];
            if (degree(s1) != 6) continue;

            // Walk strip chain through degree-6 vertices
            std::vector<node_t> sc = {s0, s1};
            node_t prv = s0, cur = s1;
            bool chainOK = true;
            while (degree(cur) == 6) {
                node_t nxt = advanceCW(cur, prv, 3);
                for (node_t x : sc) if (x == nxt) { chainOK = false; break; }
                if (!chainOK) break;
                sc.push_back(nxt);
                prv = cur; cur = nxt;
                if ((int)sc.size() > maxRedLen + 1) { chainOK = false; break; }
            }
            if (!chainOK || degree(cur) != 5) continue;

            int i_val = (int)sc.size() - 2;
            if (i_val < 1 || i_val + 1 > maxRedLen) continue;

            for (Dir dir : {Dir::DRight, Dir::DLeft}) {
                int pLen = i_val + 3;
                std::vector<node_t> path(pLen), tp(pLen);

                if (dir == Dir::DRight) {
                    path[1] = next(s0, s1);
                    path[0] = next(s0, path[1]);
                    tp[0]   = next(s0, path[0]);
                    tp[1]   = next(s0, tp[0]);
                } else {
                    path[1] = prev(s0, s1);
                    path[0] = prev(s0, path[1]);
                    tp[1]   = next(s0, s1);
                    tp[0]   = next(s0, tp[1]);
                }

                bool ok = true;
                for (int j = 1; j <= i_val && ok; j++) {
                    node_t sj = sc[j], sp = sc[j-1], sn = sc[j+1];
                    if (dir == Dir::DRight) {
                        if (next(sj, sp) != tp[j]) { ok = false; break; }
                        tp[j+1] = next(sj, tp[j]);
                        if (next(sj, tp[j+1]) != sn) { ok = false; break; }
                        path[j+1] = next(sj, sn);
                        if (next(sj, path[j+1]) != path[j]) { ok = false; break; }
                    } else {
                        if (next(sj, sp) != path[j]) { ok = false; break; }
                        path[j+1] = next(sj, path[j]);
                        if (next(sj, path[j+1]) != sn) { ok = false; break; }
                        tp[j+1] = next(sj, sn);
                        if (next(sj, tp[j+1]) != tp[j]) { ok = false; break; }
                    }
                }
                if (!ok) continue;

                node_t sf = sc[i_val+1], sfp = sc[i_val];
                if (dir == Dir::DRight) {
                    if (next(sf, sfp) != tp[i_val+1]) continue;
                    tp[i_val+2] = next(sf, tp[i_val+1]);
                    path[i_val+2] = next(sf, tp[i_val+2]);
                    if (next(sf, path[i_val+2]) != path[i_val+1]) continue;
                } else {
                    if (next(sf, sfp) != path[i_val+1]) continue;
                    path[i_val+2] = next(sf, path[i_val+1]);
                    tp[i_val+2] = next(sf, path[i_val+2]);
                    if (next(sf, tp[i_val+2]) != tp[i_val+1]) continue;
                }

                if (degree(path[0]) != 6 || degree(tp[i_val+2]) != 6) continue;

                std::set<node_t> all;
                for (node_t x : sc) all.insert(x);
                for (node_t x : path) all.insert(x);
                for (node_t x : tp) all.insert(x);
                if ((int)all.size() != (int)(sc.size() + path.size() + tp.size()))
                    continue;

                if (!rdStripTopologyClean(*this, path, sc, tp)) continue;

                InvSite site;
                site.kind = Lk(i_val); site.dir = dir;
                site.strip = sc; site.path = path; site.tp = tp;
                return site;
            }
        }
    }
    return std::nullopt;
}

std::optional<InvSite> ReducibleDual::findB00Site() const {
    for (node_t a : deg5) {
        uint8_t ma = V[a].active;
        for (; ma; ma &= ma - 1) {
            node_t b = V[a].nbr[__builtin_ctz(ma)];
            if (degree(b) != 6) continue;

            uint8_t mb = V[b].active;
            for (; mb; mb &= mb - 1) {
                node_t c = V[b].nbr[__builtin_ctz(mb)];
                if (c == a || degree(c) != 5) continue;
                if (arc_ix(a, c) >= 0) continue;  // a NOT adj c

                for (Dir dir : {Dir::DRight, Dir::DLeft}) {
                    node_t p0, q0, q1, p1;
                    if (dir == Dir::DRight) {
                        p1 = next(a, b); p0 = next(a, p1);
                        q0 = next(a, p0); q1 = next(a, q0);
                    } else {
                        p1 = prev(a, b); p0 = prev(a, p1);
                        q1 = next(a, b); q0 = next(a, q1);
                    }

                    // Validate q1 from b's perspective
                    node_t bq1;
                    if (dir == Dir::DRight) {
                        bq1 = next(b, a);
                    } else {
                        bq1 = next(b, c);
                    }
                    if (bq1 != q1) continue;

                    node_t p2, p3;
                    if (dir == Dir::DRight) {
                        p2 = prev(b, p1); p3 = prev(b, p2);
                    } else {
                        p2 = next(b, p1); p3 = next(b, p2);
                    }

                    node_t q2, p4;
                    if (dir == Dir::DRight) {
                        q2 = next(c, q1); p4 = next(c, q2);
                    } else {
                        q2 = prev(c, q1); p4 = prev(c, q2);
                    }

                    if (degree(p0) != 6 || degree(p4) != 6) continue;

                    std::set<node_t> all = {a, b, c, p0, p1, p2, p3, p4, q0, q1, q2};
                    if ((int)all.size() != 11) continue;

                    std::vector<node_t> strip = {a, b, c};
                    std::vector<node_t> path = {p0, p1, p2, p3, p4};
                    std::vector<node_t> tp = {q0, q1, q2};

                    if (!rdStripTopologyClean(*this, path, strip, tp)) continue;

                    InvSite site;
                    site.kind = Bk(0, 0); site.dir = dir;
                    site.strip = strip; site.path = path; site.tp = tp;
                    return site;
                }
            }
        }
    }
    return std::nullopt;
}

std::optional<InvSite> ReducibleDual::findFRingSite() const {
    // Detect nanotube pole: a degree-5 vertex all of whose neighbours are degree-5
    node_t pole = -1;
    for (node_t v : deg5) {
        bool all_deg5 = true;
        uint8_t m = V[v].active;
        for (uint8_t t = m; t; t &= t - 1) {
            if (degree(V[v].nbr[__builtin_ctz(t)]) != 5) { all_deg5 = false; break; }
        }
        if (all_deg5) { pole = v; break; }
    }
    if (pole < 0) return std::nullopt;

    // Collect pole-neighbours (all degree-5)
    std::set<node_t> poleNbrs;
    uint8_t mp = V[pole].active;
    for (; mp; mp &= mp - 1)
        poleNbrs.insert(V[pole].nbr[__builtin_ctz(mp)]);

    // From a pole-neighbour q, find ring_start (degree-6 neighbour of q)
    node_t q = *poleNbrs.begin();
    node_t ring_start = -1;
    uint8_t mq = V[q].active;
    for (; mq; mq &= mq - 1) {
        node_t w = V[q].nbr[__builtin_ctz(mq)];
        if (degree(w) == 6) { ring_start = w; break; }
    }
    if (ring_start < 0) return std::nullopt;

    // Find ring_second: a degree-6 neighbour of ring_start that is ALSO adjacent
    // to a pole-neighbour (distinguishes ring vertices from strip vertices)
    node_t ring_second = -1;
    uint8_t mr = V[ring_start].active;
    for (; mr; mr &= mr - 1) {
        node_t w = V[ring_start].nbr[__builtin_ctz(mr)];
        if (degree(w) != 6) continue;
        // Check if w is adjacent to any pole-neighbour
        uint8_t mw = V[w].active;
        bool adj_pole_nbr = false;
        for (uint8_t tw = mw; tw; tw &= tw - 1) {
            if (poleNbrs.count(V[w].nbr[__builtin_ctz(tw)])) { adj_pole_nbr = true; break; }
        }
        if (adj_pole_nbr) { ring_second = w; break; }
    }
    if (ring_second < 0) return std::nullopt;

    // Trace 5-cycle via straight-ahead (advanceCW 3 at degree 6)
    std::vector<node_t> ring = {ring_start, ring_second};
    node_t prv = ring_start, cur = ring_second;
    for (int step = 2; step < 5; step++) {
        node_t nxt = advanceCW(cur, prv, 3);
        if (degree(nxt) != 6) return std::nullopt;
        ring.push_back(nxt);
        prv = cur; cur = nxt;
    }
    if (advanceCW(cur, prv, 3) != ring_start) return std::nullopt;

    std::set<node_t> ringSet(ring.begin(), ring.end());
    if ((int)ringSet.size() != 5) return std::nullopt;

    // Find strip: common non-ring degree-6 neighbour of consecutive ring vertices
    std::vector<node_t> strip(5, -1);
    for (int i = 0; i < 5; i++) {
        int j = (i + 1) % 5;
        uint8_t mi = V[ring[i]].active;
        for (; mi; mi &= mi - 1) {
            node_t w = V[ring[i]].nbr[__builtin_ctz(mi)];
            if (ringSet.count(w)) continue;
            if (degree(w) == 6 && arc_ix(ring[j], w) >= 0) { strip[i] = w; break; }
        }
        if (strip[i] < 0) return std::nullopt;
    }

    std::set<node_t> stripSet(strip.begin(), strip.end());
    if ((int)stripSet.size() != 5) return std::nullopt;

    // Verify strip 5-cycle
    for (int i = 0; i < 5; i++)
        if (arc_ix(strip[i], strip[(i+1)%5]) < 0) return std::nullopt;

    // Find outer vertices
    std::vector<std::vector<node_t>> extras(5);
    for (int i = 0; i < 5; i++) {
        uint8_t mi = V[strip[i]].active;
        for (; mi; mi &= mi - 1) {
            node_t w = V[strip[i]].nbr[__builtin_ctz(mi)];
            if (!ringSet.count(w) && !stripSet.count(w)) extras[i].push_back(w);
        }
        if ((int)extras[i].size() != 2) return std::nullopt;
    }

    std::vector<node_t> outer(5, -1);
    for (int i = 0; i < 5; i++) {
        int im1 = (i + 4) % 5;
        for (node_t x : extras[im1])
            for (node_t y : extras[i])
                if (x == y) { outer[i] = x; goto rd_found_outer; }
        rd_found_outer:
        if (outer[i] < 0) return std::nullopt;
    }

    std::set<node_t> outerSet(outer.begin(), outer.end());
    if ((int)outerSet.size() != 5) return std::nullopt;

    if (!rdStripTopologyClean(*this, ring, strip, outer)) return std::nullopt;

    InvSite site;
    site.kind = Fk(); site.dir = Dir::DRight;
    site.strip = strip; site.path = ring; site.tp = outer;
    return site;
}

std::optional<InvSite> ReducibleDual::findSite(int maxRedLen) const {
    // Only accept reductions that produce valid fullerene dual vertex counts:
    // 12 (C20) or >= 14 (C24+). C22 (13v) doesn't exist; below 12 is impossible.
    auto ok = [&](const std::optional<InvSite>& s) -> bool {
        if (!s) return false;
        int result = n_alive - (int)s->strip.size();
        return result == 12 || result >= 14;
    };

    if (auto s = findL0Site(); ok(s)) return s;
    if (maxRedLen >= 2) {
        if (auto s = findLiSite(maxRedLen); ok(s)) return s;
        if (auto s = findB00Site(); ok(s)) return s;
    }
    if (auto s = findFRingSite(); ok(s)) return s;
    return std::nullopt;
}

SeedType ReducibleDual::reduceToSeed(int maxRedLen) {
    int max_steps = (int)V.size();  // At most N reductions possible
    while (auto site = findSite(maxRedLen)) {
        if (--max_steps < 0) return SeedType::NotASeed;  // Safety limit
        reduce(*site);
    }

    switch (n_alive) {
        case 12: return SeedType::C20;
        case 16: return SeedType::C28;
        case 17: return SeedType::C30;
        default: return SeedType::NotASeed;
    }
}

SeedType ReducibleDual::reduceToSeed(std::vector<ReductionStep>& path, int maxRedLen) {
    path.clear();
    int max_steps = (int)V.size();
    while (auto site = findSite(maxRedLen)) {
        if (--max_steps < 0) return SeedType::NotASeed;

        // Save strip vertex states before reduction destroys them
        ReductionStep step;
        step.site = *site;
        for (node_t s : site->strip)
            step.saved.push_back({s, V[s]});

        reduce(*site);
        path.push_back(std::move(step));
    }

    switch (n_alive) {
        case 12: return SeedType::C20;
        case 16: return SeedType::C28;
        case 17: return SeedType::C30;
        default: return SeedType::NotASeed;
    }
}

ExtensionPath ReducibleDual::reduceToExtensionPath(int maxRedLen) {
    std::vector<ReductionStep> steps;
    SeedType seed = reduceToSeed(steps, maxRedLen);

    ExtensionPath ep;
    ep.seed = seed;
    ep.full_N = (int)V.size();

    // Record seed state: full vertex data (including inactive positions for unsnip)
    for (int u = 0; u < (int)V.size(); u++) {
        if (!alive(u)) continue;
        ExtensionPath::SeedVertex sv;
        sv.id = u;
        for (int i = 0; i < D_MAX; i++) sv.nbr[i] = V[u].nbr[i];
        sv.active = V[u].active;
        ep.seed_state.push_back(sv);
    }

    // Convert reduction steps to expansion steps in reverse order
    for (int i = (int)steps.size() - 1; i >= 0; i--) {
        const auto& site = steps[i].site;
        ep.steps.push_back({site.kind, site.dir, site.strip, site.path, site.tp});
    }
    return ep;
}

Graph ReducibleDual::toGraph() const {
    std::vector<int> remap(V.size(), -1);
    int id = 0;
    for (int u = 0; u < (int)V.size(); u++)
        if (alive(u)) remap[u] = id++;

    neighbours_t adj(id, GRAPH_DMAX);
    for (int u = 0; u < (int)V.size(); u++) {
        if (!alive(u)) continue;
        uint8_t m = V[u].active;
        for (; m; m &= m - 1)
            adj[remap[u]].push_back(remap[V[u].nbr[__builtin_ctz(m)]]);
    }
    return Graph(adj);
}

Graph graphFromExtensionPath(const ExtensionPath& ep) {
    // Create empty ReducibleDual with full capacity
    ReducibleDual rd(ep.full_N);

    // Place seed vertices (full state including inactive positions)
    for (auto& sv : ep.seed_state) {
        for (int i = 0; i < ReducibleDual::D_MAX; i++)
            rd.V[sv.id].nbr[i] = sv.nbr[i];
        rd.V[sv.id].active = sv.active;
        rd.n_alive++;
        if (rd.degree(sv.id) == 5) rd.deg5.insert(sv.id);
    }

    // Apply expansion steps
    for (auto& step : ep.steps)
        rd.expand(step);

    return rd.toGraph();
}

} // namespace buckinverse
