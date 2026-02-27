#include "buckinverse.hh"
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

    pi.parallel.resize(numEntries - 1);
    for (int k = 0; k < numEntries - 1; ++k)
        pi.parallel[k] = sideNbr(g, dir, pi.path[k], pi.path[k+1]);

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
                if (isValidL0Direction(g, u, v, d))
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
                    if (isValidStraightReduction(g, u, v, ep->endpoint, ep->prev, d))
                        result.push_back({Lk(ep->distance - 1), u, v, d});
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
