#include "buckinverse.hh"
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
        node_t n = g.neighbours[u][i];
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
            node_t n = g.neighbours[s][i];
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
            node_t s1 = g.neighbours[s0][ni];
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
            node_t s1 = g.neighbours[s0][ni];
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
            node_t b = g.neighbours[a][ni];
            if (g.degree(b) != 6) continue;

            // Find degree-5 neighbors of b (candidates for c)
            for (int bi = 0; bi < g.degree(b); ++bi) {
                node_t c = g.neighbours[b][bi];
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
            node_t v = g.neighbours[u][ni];
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
                    node_t w = g.neighbours[ring[i]][k];
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
                    node_t w = g.neighbours[strip[i]][k];
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
    neighbours_t new_adj(new_N);
    for (int i = 0; i < g.N; ++i) {
        if (old_to_new[i] < 0) continue;
        for (node_t nbr : adj[i])
            new_adj[old_to_new[i]].push_back(old_to_new[nbr]);
    }
    return Graph(new_adj, true);
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
        neighbours_t new_adj(new_N);
        for (int i = 0; i < g.N; ++i) {
            if (old_to_new[i] < 0) continue;
            for (node_t nbr : adj[i])
                new_adj[old_to_new[i]].push_back(old_to_new[nbr]);
        }
        return Graph(new_adj, true);

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
        neighbours_t new_adj(new_N);
        for (int i = 0; i < g.N; ++i) {
            if (old_to_new[i] < 0) continue;
            for (node_t nbr : adj[i])
                new_adj[old_to_new[i]].push_back(old_to_new[nbr]);
        }
        return Graph(new_adj, true);

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
            node_t v = g.neighbours[u][ni];
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
                node_t v = g.neighbours[u][ni];
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
                node_t v = g.neighbours[u][ni];
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
                node_t v = g.neighbours[u][ni];
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

} // namespace buckinverse
