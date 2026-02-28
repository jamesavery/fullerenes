// Round-trip and standalone expansion tests.
// Tests both round-trip (ReductionStep-based) and standalone (ExtensionPath-based) expansion.

#include <gtest/gtest.h>
#include "fullerenes/buckinverse.hh"
#include "fullerenes/buckygen-wrapper.hh"
#include <map>

using namespace buckinverse;
using namespace std;

// Exact element-by-element comparison (same physical layout)
static bool graphsExact(const Graph& a, const Graph& b) {
    if (a.N != b.N) return false;
    for (int u = 0; u < a.N; u++) {
        if (a.degree(u) != b.degree(u)) return false;
        for (int j = 0; j < a.degree(u); j++)
            if (a.neighbours[u][j] != b.neighbours[u][j]) return false;
    }
    return true;
}

// CW-rotation-aware comparison (same circular neighbor ordering)
static bool graphsCW(const Graph& a, const Graph& b) {
    if (a.N != b.N) return false;
    for (int u = 0; u < a.N; u++) {
        int d = a.degree(u);
        if (b.degree(u) != d) return false;
        int off = -1;
        for (int k = 0; k < d; k++)
            if (b.neighbours[u][k] == a.neighbours[u][0]) { off = k; break; }
        if (off < 0) return false;
        for (int j = 0; j < d; j++)
            if (a.neighbours[u][j] != b.neighbours[u][(off + j) % d]) return false;
    }
    return true;
}

// =====================================================================
// Test: Seed graphs (trivial round-trip: 0 steps)
// =====================================================================

TEST(BuckinverseRoundtrip, SeedGraphs) {
    for (int seedN : {20, 28, 30}) {
        Graph G;
        if (seedN == 20) G = makeSeedC20();
        else if (seedN == 28) G = makeSeedC28();
        else G = makeSeedC30();

        // Round-trip test
        ReducibleDual rd(G);
        vector<ReducibleDual::ReductionStep> steps;
        SeedType seed = rd.reduceToSeed(steps);

        ASSERT_NE(seed, SeedType::NotASeed) << "C" << seedN << " seed recognized";
        EXPECT_TRUE(steps.empty()) << "C" << seedN << " seed: 0 reduction steps";

        Graph back = rd.toGraph();
        EXPECT_TRUE(graphsExact(G, back))
            << "C" << seedN << " seed: identity round-trip mismatch";

        // Standalone test
        ReducibleDual rd2(G);
        ExtensionPath ep = rd2.reduceToExtensionPath();
        Graph sa = graphFromExtensionPath(ep);
        EXPECT_TRUE(graphsCW(G, sa))
            << "C" << seedN << " seed: standalone mismatch";
    }
}

// =====================================================================
// Test: All isomers C32-C40 (both round-trip and standalone)
// =====================================================================

TEST(BuckinverseRoundtrip, AllIsomersC40) {
    int Nmax = 40;
    int total = 0, rt_failures = 0, sa_failures = 0;

    for (int N = 32; N <= Nmax; N += 2) {
        auto Q = BuckyGen::start(N, false);
        Graph G;
        int isomer_count = 0;

        while (BuckyGen::next_fullerene(Q, G)) {
            isomer_count++;
            total++;

            // --- Round-trip test ---
            ReducibleDual rd(G);
            vector<ReducibleDual::ReductionStep> steps;
            SeedType seed = rd.reduceToSeed(steps);

            if (seed == SeedType::NotASeed) {
                ADD_FAILURE() << "C" << N << " #" << isomer_count
                     << ": stuck at " << rd.N() << "v";
                rt_failures++;
                sa_failures++;
                continue;
            }

            for (int i = (int)steps.size() - 1; i >= 0; i--)
                rd.expand(steps[i]);

            Graph back = rd.toGraph();
            if (!graphsExact(G, back)) {
                ADD_FAILURE() << "C" << N << " #" << isomer_count
                     << ": round-trip mismatch";
                rt_failures++;
            }

            // --- Standalone expansion test ---
            ReducibleDual rd2(G);
            ExtensionPath ep = rd2.reduceToExtensionPath();
            Graph sa = graphFromExtensionPath(ep);

            if (!graphsCW(G, sa)) {
                ADD_FAILURE() << "C" << N << " #" << isomer_count
                     << ": standalone mismatch";
                sa_failures++;
            }
        }
        BuckyGen::stop(Q);
    }

    EXPECT_EQ(rt_failures, 0) << "Round-trip: " << rt_failures << " failures";
    EXPECT_EQ(sa_failures, 0) << "Standalone: " << sa_failures << " failures";
    EXPECT_EQ(total, 84) << "Tested all 84 isomers C32-C40";
}
