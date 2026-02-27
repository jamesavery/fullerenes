// Scalable reduction-to-seed test for larger fullerenes.
// Usage: ./test_scale [Nmax]   (default: 60)

#include "buckinverse.hh"
#include "fullerenes/buckygen-wrapper.hh"
#include <iostream>
#include <map>
#include <set>
#include <cassert>
#include <chrono>

using namespace buckinverse;
using namespace std;

// Validate that a graph is a valid fullerene dual triangulation.
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

// Diagnose a stuck graph: what site types exist?
static void diagnoseStuck(const Graph& g, int N, int isomer_idx, int steps) {
    auto inv = allInvSites(g, 5);
    auto red = allReductions(g, 5);

    // Count by type
    map<string, int> inv_types, red_types;
    for (auto& s : inv) inv_types[s.kind.toString()]++;
    for (auto& r : red) red_types[r.kind.toString()]++;

    cerr << "  STUCK C" << N << " #" << isomer_idx
         << " at " << g.N << " dual vertices after " << steps << " steps\n";
    cerr << "    allInvSites found " << inv.size() << " sites:";
    for (auto& [t, c] : inv_types) cerr << " " << t << "=" << c;
    if (inv.empty()) cerr << " (none)";
    cerr << "\n";
    cerr << "    allReductions found " << red.size() << " sites:";
    for (auto& [t, c] : red_types) cerr << " " << t << "=" << c;
    if (red.empty()) cerr << " (none)";
    cerr << "\n";

    // Show types in allReductions but not in allInvSites
    for (auto& [t, c] : red_types) {
        if (inv_types.find(t) == inv_types.end())
            cerr << "    MISSING inversion type: " << t
                 << " (" << c << " reduction sites)\n";
    }
}

int main(int argc, char* argv[]) {
    int Nmax = 60;
    if (argc > 1) Nmax = atoi(argv[1]);
    if (Nmax < 32 || Nmax % 2 != 0) {
        cerr << "Usage: " << argv[0] << " [Nmax]  (even, >= 32, default 60)\n";
        return 1;
    }

    cout << "Reducing all fullerene isomers C32-C" << Nmax << " to seeds...\n";

    auto t0 = chrono::steady_clock::now();

    map<string, int> seed_tally;
    int total = 0, failures = 0;

    for (int N = 32; N <= Nmax; N += 2) {
        auto tN = chrono::steady_clock::now();
        auto Q = BuckyGen::start(N, false);
        Graph G;
        int isomer_count = 0;
        int ok = 0, fail = 0;
        int max_steps = 0;
        int total_steps = 0;

        while (BuckyGen::next_fullerene(Q, G)) {
            isomer_count++;
            total++;

            Graph current = G;
            int steps = 0;
            bool success = true;

            while (current.N > 17) {
                auto sites = allInvSites(current, 5);
                if (sites.empty()) {
                    diagnoseStuck(current, N, isomer_count, steps);
                    success = false;
                    fail++;
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
                    diagnoseStuck(current, N, isomer_count, steps);
                    success = false;
                    fail++;
                    failures++;
                    break;
                }
                current = next;
                steps++;
            }

            if (success) {
                // May need one more reduction if stopped at N=17 but not a seed
                Triangulation t(current.neighbours, true);
                SeedType seed = identifySeed(t);
                if (seed == SeedType::NotASeed) {
                    auto sites = allInvSites(current, 5);
                    for (const auto& site : sites) {
                        Graph next = invertReduction(current, site);
                        if (next.N > 0) {
                            string err;
                            if (isValidTriangulation(next, err)) {
                                current = next;
                                steps++;
                                break;
                            }
                        }
                    }
                    t = Triangulation(current.neighbours, true);
                    seed = identifySeed(t);
                }

                if (seed == SeedType::NotASeed) {
                    cerr << "  C" << N << " #" << isomer_count
                         << " reduced to " << current.N
                         << " vertices but not a seed\n";
                    fail++;
                    failures++;
                } else {
                    seed_tally[seedName(seed)]++;
                    ok++;
                }
                total_steps += steps;
                if (steps > max_steps) max_steps = steps;
            }
        }
        BuckyGen::stop(Q);

        auto tN_end = chrono::steady_clock::now();
        double dt = chrono::duration<double>(tN_end - tN).count();

        cout << "  C" << N << ": " << isomer_count << " isomers, "
             << ok << " ok, " << fail << " fail"
             << "  (max " << max_steps << " steps, avg "
             << (isomer_count > 0 ? (double)total_steps / ok : 0.0)
             << ", " << dt << "s)\n";
        cout << flush;
    }

    auto t1 = chrono::steady_clock::now();
    double elapsed = chrono::duration<double>(t1 - t0).count();

    cout << "\nTotal: " << total << " isomers, " << failures << " failures\n";
    cout << "Seeds:";
    for (auto& [s, c] : seed_tally) cout << " " << s << "=" << c;
    cout << "\n";
    cout << "Time: " << elapsed << "s\n";
    return failures > 0 ? 1 : 0;
}
