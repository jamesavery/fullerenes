// Phase 1 tests: navigation primitives, seed construction, path computation.

#include <gtest/gtest.h>
#include "fullerenes/buckinverse.hh"

using namespace buckinverse;
using namespace std;

// =====================================================================
// Test helpers
// =====================================================================

static bool isValidTriangulation(const Triangulation& t) {
    // For a triangulation of the sphere: F = 2V - 4, E = 3V - 6
    int V = t.N;
    int E = 0;
    for (int v = 0; v < V; ++v) E += t.degree(v);
    E /= 2;
    int F = E - V + 2;
    return (F == 2 * V - 4) && (E == 3 * V - 6);
}

static bool allDegreesValid(const Triangulation& t) {
    for (int v = 0; v < t.N; ++v) {
        int d = t.degree(v);
        if (d != 5 && d != 6) return false;
    }
    return true;
}

static int countDeg5(const Triangulation& t) {
    int count = 0;
    for (int v = 0; v < t.N; ++v)
        if (t.degree(v) == 5) count++;
    return count;
}

// =====================================================================
// Test: Seed construction
// =====================================================================

TEST(BuckinversePhase1, Seeds) {
    // C20
    {
        auto c20 = makeSeedC20();
        EXPECT_EQ(c20.N, 12) << "C20 has 12 vertices";
        EXPECT_TRUE(allDegreesValid(c20)) << "C20 all degrees 5-6";
        EXPECT_EQ(countDeg5(c20), 12) << "C20 has 12 degree-5 vertices (all)";
        EXPECT_TRUE(isValidTriangulation(c20)) << "C20 is a valid triangulation";
        EXPECT_TRUE(c20.is_consistently_oriented()) << "C20 is oriented";
        EXPECT_EQ(identifySeed(c20), SeedType::C20) << "C20 identified as C20";
    }

    // C28
    {
        auto c28 = makeSeedC28();
        EXPECT_EQ(c28.N, 16) << "C28 has 16 vertices";
        EXPECT_TRUE(allDegreesValid(c28)) << "C28 all degrees 5-6";
        EXPECT_EQ(countDeg5(c28), 12) << "C28 has 12 degree-5 vertices";
        EXPECT_TRUE(isValidTriangulation(c28)) << "C28 is a valid triangulation";
        EXPECT_TRUE(c28.is_consistently_oriented()) << "C28 is oriented";
        EXPECT_EQ(identifySeed(c28), SeedType::C28) << "C28 identified as C28";
    }

    // C30
    {
        auto c30 = makeSeedC30();
        EXPECT_EQ(c30.N, 17) << "C30 has 17 vertices";
        EXPECT_TRUE(allDegreesValid(c30)) << "C30 all degrees 5-6";
        EXPECT_EQ(countDeg5(c30), 12) << "C30 has 12 degree-5 vertices";
        EXPECT_TRUE(isValidTriangulation(c30)) << "C30 is a valid triangulation";
        EXPECT_TRUE(c30.is_consistently_oriented()) << "C30 is oriented";
        EXPECT_EQ(identifySeed(c30), SeedType::C30) << "C30 identified as C30";
    }
}

// =====================================================================
// Test: Navigation primitives
// =====================================================================

TEST(BuckinversePhase1, Navigation) {
    auto c20 = makeSeedC20();

    // In C20, all vertices are degree-5.
    for (int v = 0; v < c20.N; ++v) {
        EXPECT_EQ(c20.degree(v), 5) << "C20 vertex " << v << " degree 5";

        node_t nb0 = c20.nbrs(v)[0];
        EXPECT_EQ(advanceCW(c20, v, nb0, 0), nb0)
            << "advanceCW(v," << nb0 << ",0) == " << nb0;

        EXPECT_EQ(advanceCW(c20, v, nb0, 5), nb0)
            << "advanceCW wraps at degree";

        EXPECT_EQ(c20.next(v, nb0), advanceCW(c20, v, nb0, 1))
            << "next == advanceCW(1)";

        EXPECT_EQ(c20.prev(v, nb0), advanceCW(c20, v, nb0, 4))
            << "prev == advanceCW(4)";
    }

    // Test straightAhead on C20 (all degree-5).
    {
        node_t u = 0;
        node_t v = c20.nbrs(u)[0];
        EXPECT_EQ(straightAhead(c20, Dir::DRight, v, u), advanceCW(c20, v, u, 3))
            << "straightAhead DRight at deg-5 == advance 3";
        EXPECT_EQ(straightAhead(c20, Dir::DLeft, v, u), advanceCW(c20, v, u, 2))
            << "straightAhead DLeft at deg-5 == advance 2";
    }

    // Test sideNbr
    {
        node_t u = 0;
        node_t v = c20.nbrs(u)[0];
        EXPECT_EQ(sideNbr(c20, Dir::DRight, u, v), c20.prev(u, v))
            << "sideNbr DRight == prev";
        EXPECT_EQ(sideNbr(c20, Dir::DLeft, u, v), c20.next(u, v))
            << "sideNbr DLeft == next";
    }

    // Test on a graph with degree-6 vertices (C28)
    auto c28 = makeSeedC28();
    {
        node_t v6 = -1;
        for (int v = 0; v < c28.N; ++v)
            if (c28.degree(v) == 6) { v6 = v; break; }
        ASSERT_GE(v6, 0) << "C28 has a degree-6 vertex";

        node_t nb = c28.nbrs(v6)[0];
        EXPECT_EQ(straightAhead(c28, Dir::DRight, v6, nb), advanceCW(c28, v6, nb, 3))
            << "straightAhead DRight at deg-6 == advance 3";
        EXPECT_EQ(straightAhead(c28, Dir::DLeft, v6, nb), advanceCW(c28, v6, nb, 3))
            << "straightAhead DLeft at deg-6 == advance 3";

        EXPECT_EQ(turnAhead(c28, Dir::DRight, v6, nb), advanceCW(c28, v6, nb, 2))
            << "turnAhead DRight at deg-6 == advance 2";
        EXPECT_EQ(turnAhead(c28, Dir::DLeft, v6, nb), advanceCW(c28, v6, nb, 4))
            << "turnAhead DLeft at deg-6 == advance 4";
    }
}

// =====================================================================
// Test: Path computation
// =====================================================================

TEST(BuckinversePhase1, Paths) {
    auto c28 = makeSeedC28();
    auto d5 = deg5vertices(c28);
    EXPECT_EQ((int)d5.size(), 12) << "C28 deg5vertices returns 12";

    // Test straight path computation from each degree-5 vertex
    int valid_L0_count = 0;
    for (node_t u : d5) {
        for (int ni = 0; ni < c28.degree(u); ++ni) {
            node_t v = c28.nbrs(u)[ni];
            for (Dir d : {Dir::DRight, Dir::DLeft}) {
                auto pi = computeStraightPath(c28, u, v, d, 3);
                if (pi.valid) {
                    EXPECT_EQ((int)pi.path.size(), 3) << "L0 path has 3 entries";
                    EXPECT_EQ((int)pi.parallel.size(), 3) << "L0 parallel has 3 entries";
                    EXPECT_GE(c28.arc_ix(pi.path[0], pi.path[1]), 0)
                        << "path[0]-path[1] are adjacent";
                    EXPECT_GE(c28.arc_ix(pi.path[1], pi.path[2]), 0)
                        << "path[1]-path[2] are adjacent";
                    valid_L0_count++;
                }
            }
        }
    }
    EXPECT_GT(valid_L0_count, 0) << "Found valid L0 paths in C28";

    // Test longer straight paths (L1: 4 entries)
    int valid_L1_count = 0;
    for (node_t u : d5) {
        for (int ni = 0; ni < c28.degree(u); ++ni) {
            node_t v = c28.nbrs(u)[ni];
            for (Dir d : {Dir::DRight, Dir::DLeft}) {
                auto pi = computeStraightPath(c28, u, v, d, 4);
                if (pi.valid) {
                    EXPECT_EQ((int)pi.path.size(), 4) << "L1 path has 4 entries";
                    EXPECT_EQ((int)pi.parallel.size(), 4) << "L1 parallel has 4 entries";
                    valid_L1_count++;
                }
            }
        }
    }
    EXPECT_GT(valid_L1_count, 0) << "Found valid L1 paths in C28";

    // Test bent zero path
    int valid_B00_count = 0;
    for (node_t u : d5) {
        for (int ni = 0; ni < c28.degree(u); ++ni) {
            node_t v = c28.nbrs(u)[ni];
            for (Dir d : {Dir::DRight, Dir::DLeft}) {
                auto pi = computeBentZeroPath(c28, u, v, d);
                if (pi.valid) {
                    EXPECT_EQ((int)pi.path.size(), 5) << "B00 path has 5 entries";
                    EXPECT_EQ((int)pi.parallel.size(), 3) << "B00 parallel has 3 entries";
                    valid_B00_count++;
                }
            }
        }
    }
    EXPECT_GT(valid_B00_count, 0) << "Found valid B(0,0) paths in C28";
}

// =====================================================================
// Test: Degree-5 enumeration helper
// =====================================================================

TEST(BuckinversePhase1, Deg5) {
    EXPECT_EQ((int)deg5vertices(makeSeedC20()).size(), 12) << "C20 has 12 degree-5 vertices";
    EXPECT_EQ((int)deg5vertices(makeSeedC28()).size(), 12) << "C28 has 12 degree-5 vertices";
    EXPECT_EQ((int)deg5vertices(makeSeedC30()).size(), 12) << "C30 has 12 degree-5 vertices";
}

// =====================================================================
// Test: Navigation consistency
// =====================================================================

TEST(BuckinversePhase1, NavigationConsistency) {
    auto c28 = makeSeedC28();

    for (int u = 0; u < c28.N; ++u) {
        for (node_t v : c28.nbrs(u)) {
            for (Dir d : {Dir::DRight, Dir::DLeft}) {
                node_t sa = straightAhead(c28, d, u, v);
                EXPECT_GE(sa, 0);
                EXPECT_LT(sa, c28.N);
                EXPECT_GE(c28.arc_ix(u, sa), 0)
                    << "straightAhead result is a neighbor of u";

                node_t ta = turnAhead(c28, d, u, v);
                EXPECT_GE(ta, 0);
                EXPECT_LT(ta, c28.N);
                EXPECT_GE(c28.arc_ix(u, ta), 0)
                    << "turnAhead result is a neighbor of u";

                node_t sn = sideNbr(c28, d, u, v);
                EXPECT_GE(sn, 0);
                EXPECT_LT(sn, c28.N);
                EXPECT_GE(c28.arc_ix(u, sn), 0)
                    << "sideNbr result is a neighbor of u";
            }
        }
    }
}
