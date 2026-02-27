// Round-trip test: reduce fullerene dual to seed, expand back, verify exact match.
// Usage: ./test_roundtrip [Nmax]   (default: 100)
// For every isomer C20-Nmax: reduce → seed → expand → compare adjacency.

#include "buckinverse.hh"
#include "fullerenes/buckygen-wrapper.hh"
#include <iostream>
#include <map>
#include <chrono>

using namespace buckinverse;
using namespace std;
using clk = chrono::steady_clock;

static bool graphsEqual(const Graph& a, const Graph& b) {
    if (a.N != b.N) return false;
    for (int u = 0; u < a.N; u++) {
        if (a.degree(u) != b.degree(u)) return false;
        for (int j = 0; j < a.degree(u); j++)
            if (a.neighbours[u][j] != b.neighbours[u][j]) return false;
    }
    return true;
}

int main(int argc, char* argv[]) {
    int Nmax = 100;
    if (argc > 1) Nmax = atoi(argv[1]);
    if (Nmax < 20 || Nmax % 2 != 0) {
        cerr << "Usage: " << argv[0] << " [Nmax]  (even, >= 20, default 100)\n";
        return 1;
    }

    cout << "Round-trip test: C20-C" << Nmax << "\n\n";

    map<string, int> seed_tally;
    int total = 0, failures = 0;
    double total_red = 0, total_exp = 0;

    // Test seed graphs themselves (trivial: 0 steps)
    for (int seedN : {20, 28, 30}) {
        if (seedN > Nmax) continue;
        Graph G;
        if (seedN == 20) G = makeSeedC20();
        else if (seedN == 28) G = makeSeedC28();
        else G = makeSeedC30();

        ReducibleDual rd(G);
        vector<ReducibleDual::ReductionStep> steps;
        SeedType seed = rd.reduceToSeed(steps);

        if (seed == SeedType::NotASeed) {
            cerr << "  FAIL C" << seedN << " seed: not recognized\n";
            failures++;
        } else if (!steps.empty()) {
            cerr << "  FAIL C" << seedN << " seed: " << steps.size() << " reduction steps (expected 0)\n";
            failures++;
        } else {
            Graph back = rd.toGraph();
            if (!graphsEqual(G, back)) {
                cerr << "  FAIL C" << seedN << " seed: identity round-trip mismatch\n";
                failures++;
            } else {
                seed_tally[seedName(seed)]++;
            }
        }
        total++;
    }

    for (int N = 32; N <= Nmax; N += 2) {
        if (N == 22) continue;  // C22 doesn't exist

        auto Q = BuckyGen::start(N, false);
        Graph G;
        int isomer_count = 0, ok = 0, fail = 0;
        double dt_red = 0, dt_exp = 0;

        while (BuckyGen::next_fullerene(Q, G)) {
            isomer_count++;
            total++;

            // Reduce, recording steps
            auto t0 = clk::now();
            ReducibleDual rd(G);
            vector<ReducibleDual::ReductionStep> steps;
            SeedType seed = rd.reduceToSeed(steps);
            auto t1 = clk::now();
            dt_red += chrono::duration<double>(t1 - t0).count();

            if (seed == SeedType::NotASeed) {
                cerr << "  FAIL C" << N << " #" << isomer_count
                     << ": stuck at " << rd.N() << "v\n";
                fail++; failures++;
                continue;
            }

            // Expand back
            auto t2 = clk::now();
            for (int i = (int)steps.size() - 1; i >= 0; i--)
                rd.expand(steps[i]);
            auto t3 = clk::now();
            dt_exp += chrono::duration<double>(t3 - t2).count();

            // Verify exact match
            Graph back = rd.toGraph();
            if (!graphsEqual(G, back)) {
                cerr << "  FAIL C" << N << " #" << isomer_count
                     << ": round-trip mismatch (N=" << back.N << " vs " << G.N << ")\n";
                fail++; failures++;
                continue;
            }

            seed_tally[seedName(seed)]++;
            ok++;
        }
        BuckyGen::stop(Q);

        total_red += dt_red;
        total_exp += dt_exp;

        int Nv = N / 2 + 2;
        double us_red = isomer_count > 0 ? 1e6 * dt_red / isomer_count : 0;
        double us_exp = isomer_count > 0 ? 1e6 * dt_exp / isomer_count : 0;
        cout << "  C" << N << ": " << isomer_count << " isomers, "
             << ok << " ok, " << fail << " fail  "
             << "red=" << us_red << "  exp=" << us_exp << " us/isomer"
             << "  (" << (us_red + us_exp) / Nv << " us/vertex)\n" << flush;
    }

    cout << "\nTotal: " << total << " isomers, " << failures << " failures\n";
    cout << "Seeds:";
    for (auto& [s, c] : seed_tally) cout << " " << s << "=" << c;
    cout << "\n";
    cout << "Time: red=" << total_red << "s  exp=" << total_exp << "s\n";
    return failures > 0 ? 1 : 0;
}
