// Phase 3 tests: reduction surgery via strip-based inversion.

#include <gtest/gtest.h>
#include "fullerenes/buckinverse.hh"
#include "fullerenes/buckygen-wrapper.hh"
#include <map>
#include <set>

using namespace buckinverse;
using namespace std;

// =====================================================================
// Test helpers
// =====================================================================

static bool isValidTriangulation(const Graph& g, string& err) {
    int deg5_count = 0;
    for (int u = 0; u < g.N; ++u) {
        int d = g.degree(u);
        if (d != 5 && d != 6) {
            err = "vertex " + to_string(u) + " has degree " + to_string(d);
            return false;
        }
        if (d == 5) deg5_count++;

        set<node_t> seen;
        for (int i = 0; i < d; ++i) {
            node_t v = g.neighbours[u][i];
            if (v == u) { err = "self-loop at " + to_string(u); return false; }
            if (v < 0 || v >= g.N) {
                err = "invalid neighbour " + to_string(v) + " at " + to_string(u);
                return false;
            }
            if (!seen.insert(v).second) {
                err = "duplicate neighbour " + to_string(v) + " at " + to_string(u);
                return false;
            }
        }

        for (int i = 0; i < d; ++i) {
            node_t v = g.neighbours[u][i];
            if (g.arc_ix(v, u) < 0) {
                err = "asymmetric: " + to_string(u) + " adj " + to_string(v)
                    + " but not vice versa";
                return false;
            }
        }
    }
    if (deg5_count != 12) {
        err = "has " + to_string(deg5_count) + " degree-5 vertices (expected 12)";
        return false;
    }
    return true;
}

// =====================================================================
// Test: L0 inversion sites are found
// =====================================================================

TEST(BuckinversePhase3, L0Sites) {
    for (int N = 32; N <= 40; N += 2) {
        auto Q = BuckyGen::start(N, false);
        Graph G;
        int isomer_idx = 0;
        int with_sites = 0;

        while (BuckyGen::next_fullerene(Q, G)) {
            isomer_idx++;
            auto sites = findL0InvSites(G);
            if (!sites.empty()) with_sites++;
        }
        BuckyGen::stop(Q);

        EXPECT_GT(with_sites, 0) << "C" << N << " has isomers with L0 sites";
    }
}

// =====================================================================
// Test: single-step inversion produces valid triangulation
// =====================================================================

TEST(BuckinversePhase3, SingleInversion) {
    for (int N = 32; N <= 40; N += 2) {
        auto Q = BuckyGen::start(N, false);
        Graph G;
        int isomer_idx = 0;
        int ok = 0, fail = 0;

        while (BuckyGen::next_fullerene(Q, G)) {
            isomer_idx++;
            auto sites = allInvSites(G, 5);
            if (sites.empty()) continue;

            for (const auto& site : sites) {
                Graph reduced = invertReduction(G, site);
                if (reduced.N == 0) continue;

                int expected_N = G.N - site.kind.newVertices();
                EXPECT_EQ(reduced.N, expected_N)
                    << "C" << N << " #" << isomer_idx << " " << site.kind.toString()
                    << " reduced from " << G.N << " to " << reduced.N;

                string err;
                bool valid = isValidTriangulation(reduced, err);
                EXPECT_TRUE(valid)
                    << "C" << N << " #" << isomer_idx << " " << site.kind.toString()
                    << " result: " << err;

                if (valid && reduced.N == expected_N) ok++;
                else fail++;
                break;
            }
        }
        BuckyGen::stop(Q);

        EXPECT_EQ(fail, 0) << "C" << N << ": " << fail << " inversion failures";
    }
}

// =====================================================================
// Test: iterative reduction to seed
// =====================================================================

TEST(BuckinversePhase3, ReduceToSeed) {
    map<string, int> seed_tally;
    int total = 0, failures = 0;

    for (int N = 32; N <= 40; N += 2) {
        auto Q = BuckyGen::start(N, false);
        Graph G;
        int isomer_idx = 0;

        while (BuckyGen::next_fullerene(Q, G)) {
            isomer_idx++;
            total++;

            Graph current = G;
            int steps = 0;
            bool success = true;

            while (current.N > 17) {
                auto sites = allInvSites(current, 5);
                if (sites.empty()) {
                    success = false;
                    failures++;
                    ADD_FAILURE() << "C" << N << " #" << isomer_idx
                         << " stuck at " << current.N << " vertices after "
                         << steps << " steps";
                    break;
                }

                Graph next;
                bool found = false;
                for (const auto& site : sites) {
                    next = invertReduction(current, site);
                    if (next.N > 0) {
                        string err;
                        if (isValidTriangulation(next, err)) {
                            found = true;
                            break;
                        }
                    }
                }

                if (!found) {
                    success = false;
                    failures++;
                    ADD_FAILURE() << "C" << N << " #" << isomer_idx
                         << " stuck at " << current.N << " vertices after "
                         << steps << " steps (all " << sites.size()
                         << " sites failed)";
                    break;
                }
                current = next;
                steps++;
            }

            if (success) {
                Triangulation t(current.neighbours);
                SeedType seed = identifySeed(t);
                if (seed == SeedType::NotASeed) {
                    auto sites = allInvSites(current, 5);
                    if (!sites.empty()) {
                        for (const auto& site : sites) {
                            Graph next = invertReduction(current, site);
                            if (next.N > 0) {
                                string err;
                                if (isValidTriangulation(next, err)) {
                                    current = next;
                                    break;
                                }
                            }
                        }
                        t = Triangulation(current.neighbours);
                        seed = identifySeed(t);
                    }
                }
                EXPECT_NE(seed, SeedType::NotASeed)
                    << "C" << N << " #" << isomer_idx << " reduced to seed in "
                    << steps << " steps";
                seed_tally[seedName(seed)]++;
            }
        }
        BuckyGen::stop(Q);
    }

    EXPECT_EQ(failures, 0) << "All isomers reduced to seeds";
    EXPECT_EQ(total, 84) << "Tested all 84 isomers C32-C40";
}
