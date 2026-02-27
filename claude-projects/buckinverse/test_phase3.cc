// Phase 3 tests: reduction surgery via strip-based inversion.

#include "buckinverse.hh"
#include "fullerenes/buckygen-wrapper.hh"
#include <iostream>
#include <map>
#include <set>
#include <cassert>

using namespace buckinverse;
using namespace std;

// =====================================================================
// Test helpers
// =====================================================================

static int tests_passed = 0;
static int tests_failed = 0;

#define CHECK(cond, msg) do { \
    if (cond) { tests_passed++; } \
    else { tests_failed++; cerr << "FAIL: " << msg << " (line " << __LINE__ << ")\n"; } \
} while(0)

// Validate that a graph is a valid fullerene dual triangulation:
// - All vertices have degree 5 or 6
// - Exactly 12 degree-5 vertices
// - All neighbour lists are consistent (u adj v iff v adj u)
// - No self-loops or duplicate neighbours
static bool isValidTriangulation(const Graph& g, string& err) {
    int deg5_count = 0;
    for (int u = 0; u < g.N; ++u) {
        int d = g.degree(u);
        if (d != 5 && d != 6) {
            err = "vertex " + to_string(u) + " has degree " + to_string(d);
            return false;
        }
        if (d == 5) deg5_count++;

        // Check for self-loops and duplicates
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

        // Check symmetry: u adj v implies v adj u
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

void test_l0_sites() {
    cout << "--- L0 inversion sites ---\n";

    for (int N = 32; N <= 40; N += 2) {
        cout << flush; cerr << flush;
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

        cout << "  C" << N << ": " << isomer_idx << " isomers, "
             << with_sites << " have L0 inversion sites\n";
        CHECK(with_sites > 0, "C" + to_string(N) + " has isomers with L0 sites");
    }
}

// =====================================================================
// Test: single-step L0 inversion produces valid triangulation
// =====================================================================

void test_single_inversion() {
    cout << "--- Single-step L0 inversion ---\n";

    for (int N = 32; N <= 40; N += 2) {
        cout << flush; cerr << flush;
        auto Q = BuckyGen::start(N, false);
        Graph G;
        int isomer_idx = 0;
        int ok = 0, fail = 0, no_site = 0;

        while (BuckyGen::next_fullerene(Q, G)) {
            isomer_idx++;
            auto sites = allInvSites(G, 5);
            if (sites.empty()) {
                no_site++;
                continue;
            }

            // Try first site
            bool found_valid = false;
            for (const auto& site : sites) {
                Graph reduced = invertReduction(G, site);
                if (reduced.N == 0) continue;

                int expected_N = G.N - site.kind.newVertices();
                CHECK(reduced.N == expected_N,
                      "C" + to_string(N) + " #" + to_string(isomer_idx)
                      + " " + site.kind.toString()
                      + " reduced from " + to_string(G.N) + " to "
                      + to_string(reduced.N) + " (expected "
                      + to_string(expected_N) + ")");

                string err;
                bool valid = isValidTriangulation(reduced, err);
                CHECK(valid,
                      "C" + to_string(N) + " #" + to_string(isomer_idx)
                      + " " + site.kind.toString() + " result: " + err);

                if (valid && reduced.N == expected_N) {
                    ok++;
                    found_valid = true;
                } else {
                    fail++;
                }
                break;
            }
            if (!found_valid && !sites.empty()) {
                fail++;
            }
        }
        BuckyGen::stop(Q);

        cout << "  C" << N << ": " << isomer_idx << " isomers, "
             << ok << " ok, " << fail << " failed, "
             << no_site << " no site\n";
    }
}

// =====================================================================
// Test: iterative reduction to seed
// =====================================================================

void test_reduce_to_seed() {
    cout << "--- Iterative reduction to seed ---\n";

    map<string, int> seed_tally;
    int total = 0, failures = 0;

    for (int N = 32; N <= 40; N += 2) {
        cout << flush; cerr << flush;
        auto Q = BuckyGen::start(N, false);
        Graph G;
        int isomer_idx = 0;

        while (BuckyGen::next_fullerene(Q, G)) {
            isomer_idx++;
            total++;

            Graph current = G;
            int steps = 0;
            bool success = true;

            while (current.N > 17) { // largest seed has 17 vertices
                auto sites = allInvSites(current, 5);
                if (sites.empty()) {
                    cerr << "FAIL: C" << N << " #" << isomer_idx
                         << " stuck at " << current.N << " vertices after "
                         << steps << " steps (no inversion sites)\n";
                    success = false;
                    failures++;
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
                    cerr << "FAIL: C" << N << " #" << isomer_idx
                         << " stuck at " << current.N << " vertices after "
                         << steps << " steps (all " << sites.size()
                         << " sites failed)\n";
                    success = false;
                    failures++;
                    break;
                }
                current = next;
                steps++;
            }

            if (success) {
                Triangulation t(current.neighbours, true);
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
                        t = Triangulation(current.neighbours, true);
                        seed = identifySeed(t);
                    }
                }
                string sname = seedName(seed);
                seed_tally[sname]++;
                CHECK(seed != SeedType::NotASeed,
                      "C" + to_string(N) + " #" + to_string(isomer_idx)
                      + " reduced to " + sname + " in " + to_string(steps)
                      + " steps");
            }
        }
        BuckyGen::stop(Q);
    }

    cout << "  Total: " << total << " isomers, " << failures << " failures\n";
    cout << "  Seed tally:";
    for (auto& [s, c] : seed_tally)
        cout << " " << s << "=" << c;
    cout << "\n";
}

// =====================================================================
// Main
// =====================================================================

int main() {
    test_l0_sites();
    test_single_inversion();
    test_reduce_to_seed();

    cout << "\n=== Results: " << tests_passed << " passed, "
         << tests_failed << " failed ===\n";
    return tests_failed > 0 ? 1 : 0;
}
