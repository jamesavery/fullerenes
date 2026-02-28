// Phase 2 tests: reduction enumeration.

#include <gtest/gtest.h>
#include "fullerenes/buckinverse.hh"
#include "fullerenes/buckygen-wrapper.hh"
#include <map>
#include <set>

using namespace buckinverse;
using namespace std;

// =====================================================================
// Test: allReductions on seed graphs
// =====================================================================

TEST(BuckinversePhase2, ReductionsC20) {
    auto c20 = makeSeedC20();
    auto reds = allReductions(c20);

    map<string, int> counts;
    for (const auto& r : reds)
        counts[r.kind.toString()]++;

    EXPECT_GT(reds.size(), 0u) << "C20 has reductions (all B(0,0))";
    EXPECT_EQ(counts.count("L0"), 0u) << "C20 has no L0 reductions (no degree-6 flanking)";
    EXPECT_GT(counts.count("B(0,0)"), 0u) << "C20 has B(0,0) reductions";

    for (const auto& r : reds) {
        EXPECT_EQ(c20.degree(r.u), 5) << "start vertex is degree-5";
        EXPECT_GE(c20.arc_ix(r.u, r.v), 0) << "edge exists";
    }
}

TEST(BuckinversePhase2, ReductionsC28) {
    auto c28 = makeSeedC28();
    auto reds = allReductions(c28);

    EXPECT_GT(reds.size(), 0u) << "C28 has reductions";

    for (const auto& r : reds) {
        EXPECT_EQ(c28.degree(r.u), 5) << r.toString() << " start vertex is degree-5";
        EXPECT_GE(c28.arc_ix(r.u, r.v), 0) << r.toString() << " edge exists";
    }
}

TEST(BuckinversePhase2, ReductionsC30) {
    auto c30 = makeSeedC30();
    auto reds = allReductions(c30);

    EXPECT_GT(reds.size(), 0u) << "C30 has reductions";

    for (const auto& r : reds) {
        EXPECT_EQ(c30.degree(r.u), 5) << r.toString() << " start vertex is degree-5";
        EXPECT_GE(c30.arc_ix(r.u, r.v), 0) << r.toString() << " edge exists";
    }
}

// =====================================================================
// Test: followStraightToFive
// =====================================================================

TEST(BuckinversePhase2, FollowStraight) {
    auto c28 = makeSeedC28();

    int found = 0;
    for (node_t u : deg5vertices(c28)) {
        for (int ni = 0; ni < c28.degree(u); ++ni) {
            node_t v = c28.neighbours[u][ni];
            if (c28.degree(v) == 5) continue;

            auto ep = followStraightToFive(c28, u, v, 5);
            if (ep) {
                EXPECT_EQ(c28.degree(ep->endpoint), 5)
                    << "followStraight endpoint is degree-5";
                EXPECT_GE(ep->distance, 2)
                    << "followStraight distance >= 2";
                EXPECT_LE(ep->distance, 5)
                    << "followStraight distance <= maxDist";
                EXPECT_GE(c28.arc_ix(ep->prev, ep->endpoint), 0)
                    << "followStraight prev is adjacent to endpoint";
                found++;
            }
        }
    }
    EXPECT_GT(found, 0) << "Found at least one straight walk endpoint";
}

// =====================================================================
// Test: B00 symmetry
// =====================================================================

TEST(BuckinversePhase2, B00Symmetry) {
    auto c20 = makeSeedC20();
    auto reds = allReductions(c20, 2);

    for (const auto& r : reds) {
        EXPECT_EQ(c20.degree(r.u), 5) << "B00 start is degree-5";
    }
    EXPECT_GT(reds.size(), 0u) << "B00 reductions on C20";
}

// =====================================================================
// Test: Reduction path validity
// =====================================================================

TEST(BuckinversePhase2, ReductionPaths) {
    auto c28 = makeSeedC28();
    auto reds = allReductions(c28);

    for (const auto& r : reds) {
        if (r.kind.type == ExpKind::L_type) {
            int numEntries = r.kind.i + 3;
            auto pi = computeStraightPath(c28, r.u, r.v, r.dir, numEntries);
            EXPECT_TRUE(pi.valid) << r.toString() << " straight path is valid";
            if (pi.valid) {
                EXPECT_EQ(c28.degree(pi.path.front()), 5)
                    << r.toString() << " path starts at degree-5";
                EXPECT_EQ(c28.degree(pi.path.back()), 5)
                    << r.toString() << " path ends at degree-5";
            }
        } else if (r.kind == Bk(0, 0)) {
            auto pi = computeBentZeroPath(c28, r.u, r.v, r.dir);
            set<node_t> pathSet(pi.path.begin(), pi.path.end());
            EXPECT_EQ((int)pathSet.size(), (int)pi.path.size())
                << r.toString() << " B00 path has no duplicates";
            EXPECT_EQ(c28.degree(pi.path.back()), 5)
                << r.toString() << " B00 path ends at degree-5";
        } else if (r.kind.type == ExpKind::B_type) {
            auto pi = computeBentPath(c28, r.u, r.v, r.dir, r.kind.i, r.kind.j);
            EXPECT_EQ((int)pi.path.size(), r.kind.i + r.kind.j + 5)
                << r.toString() << " bent path has correct length";
        }
    }
}

// =====================================================================
// Test: No duplicate reductions
// =====================================================================

TEST(BuckinversePhase2, NoDuplicates) {
    auto c28 = makeSeedC28();
    auto reds = allReductions(c28);

    set<tuple<int, int, int, int, node_t, node_t, int>> seen;
    int dups = 0;
    for (const auto& r : reds) {
        auto key = make_tuple(r.kind.type, r.kind.i, r.kind.j,
                              (int)r.dir, r.u, r.v, 0);
        if (!seen.insert(key).second) dups++;
    }
    EXPECT_EQ(dups, 0) << "No duplicate reductions in C28";
}

// =====================================================================
// Test: maxRedLen parameter
// =====================================================================

TEST(BuckinversePhase2, MaxRedLen) {
    auto c28 = makeSeedC28();
    auto reds1 = allReductions(c28, 1);
    auto reds2 = allReductions(c28, 2);
    auto reds5 = allReductions(c28, 5);

    EXPECT_LE(reds1.size(), reds2.size()) << "maxRedLen=1 <= maxRedLen=2";
    EXPECT_LE(reds2.size(), reds5.size()) << "maxRedLen=2 <= maxRedLen=5";

    for (const auto& r : reds1) {
        EXPECT_EQ(r.kind, Lk(0)) << "maxRedLen=1 produces only L0";
    }
}

// =====================================================================
// Test: Reductions on generated fullerenes (non-seeds)
// =====================================================================

TEST(BuckinversePhase2, GeneratedFullerenes) {
    for (int N = 32; N <= 40; N += 2) {
        auto Q = BuckyGen::start(N, false);
        Graph G;
        int isomer_count = 0;

        while (BuckyGen::next_fullerene(Q, G)) {
            isomer_count++;
            auto reds = allReductions(G, 5);

            EXPECT_FALSE(reds.empty())
                << "C" << N << " isomer " << isomer_count
                << " has at least one reduction";

            for (const auto& r : reds) {
                EXPECT_EQ(G.degree(r.u), 5) << "start vertex is degree-5";
                EXPECT_GE(G.arc_ix(r.u, r.v), 0) << "edge exists";
            }
        }
        BuckyGen::stop(Q);

        EXPECT_GT(isomer_count, 0) << "C" << N << " has isomers";
    }
}
