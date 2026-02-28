// Round-trip test: reduce fullerene dual to seed, expand back, verify exact match.
// Tests both round-trip (ReductionStep-based) and standalone (ExtensionPath-based) expansion.
// Usage: ./test_roundtrip [Nmax]   (default: 100)

#include "buckinverse.hh"
#include "fullerenes/buckygen-wrapper.hh"
#include <iostream>
#include <map>
#include <chrono>

using namespace buckinverse;
using namespace std;
using clk = chrono::steady_clock;

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
        // Find rotation offset
        int off = -1;
        for (int k = 0; k < d; k++)
            if (b.neighbours[u][k] == a.neighbours[u][0]) { off = k; break; }
        if (off < 0) return false;
        for (int j = 0; j < d; j++)
            if (a.neighbours[u][j] != b.neighbours[u][(off + j) % d]) return false;
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
    int total = 0, rt_failures = 0, sa_failures = 0;
    double total_red = 0, total_exp = 0, total_sa = 0;

    // Test seed graphs (trivial: 0 steps)
    for (int seedN : {20, 28, 30}) {
        if (seedN > Nmax) continue;
        Graph G;
        if (seedN == 20) G = makeSeedC20();
        else if (seedN == 28) G = makeSeedC28();
        else G = makeSeedC30();

        // Round-trip test
        ReducibleDual rd(G);
        vector<ReducibleDual::ReductionStep> steps;
        SeedType seed = rd.reduceToSeed(steps);

        if (seed == SeedType::NotASeed) {
            cerr << "  FAIL C" << seedN << " seed: not recognized\n";
            rt_failures++;
        } else if (!steps.empty()) {
            cerr << "  FAIL C" << seedN << " seed: " << steps.size() << " reduction steps\n";
            rt_failures++;
        } else {
            Graph back = rd.toGraph();
            if (!graphsExact(G, back)) {
                cerr << "  FAIL C" << seedN << " seed: identity round-trip mismatch\n";
                rt_failures++;
            }
        }

        // Standalone test
        ReducibleDual rd2(G);
        ExtensionPath ep = rd2.reduceToExtensionPath();
        Graph sa = graphFromExtensionPath(ep);
        if (!graphsCW(G, sa)) {
            cerr << "  FAIL C" << seedN << " seed: standalone mismatch\n";
            sa_failures++;
        }

        seed_tally[seedName(seed)]++;
        total++;
    }

    for (int N = 32; N <= Nmax; N += 2) {
        if (N == 22) continue;

        auto Q = BuckyGen::start(N, false);
        Graph G;
        int isomer_count = 0, rt_ok = 0, rt_fail = 0, sa_ok = 0, sa_fail = 0;
        double dt_red = 0, dt_exp = 0, dt_sa = 0;

        while (BuckyGen::next_fullerene(Q, G)) {
            isomer_count++;
            total++;

            // --- Round-trip test ---
            auto t0 = clk::now();
            ReducibleDual rd(G);
            vector<ReducibleDual::ReductionStep> steps;
            SeedType seed = rd.reduceToSeed(steps);
            auto t1 = clk::now();
            dt_red += chrono::duration<double>(t1 - t0).count();

            if (seed == SeedType::NotASeed) {
                cerr << "  FAIL C" << N << " #" << isomer_count
                     << ": stuck at " << rd.N() << "v\n";
                rt_fail++; rt_failures++;
                sa_fail++; sa_failures++;
                continue;
            }

            auto t2 = clk::now();
            for (int i = (int)steps.size() - 1; i >= 0; i--)
                rd.expand(steps[i]);
            auto t3 = clk::now();
            dt_exp += chrono::duration<double>(t3 - t2).count();

            Graph back = rd.toGraph();
            if (!graphsExact(G, back)) {
                cerr << "  FAIL C" << N << " #" << isomer_count
                     << ": round-trip mismatch\n";
                rt_fail++; rt_failures++;
            } else {
                rt_ok++;
            }

            // --- Standalone expansion test ---
            auto t4 = clk::now();
            ReducibleDual rd2(G);
            ExtensionPath ep = rd2.reduceToExtensionPath();
            Graph sa = graphFromExtensionPath(ep);
            auto t5 = clk::now();
            dt_sa += chrono::duration<double>(t5 - t4).count();

            if (!graphsCW(G, sa)) {
                cerr << "  FAIL C" << N << " #" << isomer_count
                     << ": standalone mismatch\n";
                sa_fail++; sa_failures++;
            } else {
                sa_ok++;
            }

            seed_tally[seedName(seed)]++;
        }
        BuckyGen::stop(Q);

        total_red += dt_red;
        total_exp += dt_exp;
        total_sa += dt_sa;

        cout << "  C" << N << ": " << isomer_count << " isomers  "
             << "rt=" << rt_ok << "/" << rt_fail << "  "
             << "sa=" << sa_ok << "/" << sa_fail << "\n" << flush;
    }

    cout << "\nTotal: " << total << " isomers\n";
    cout << "Round-trip: " << rt_failures << " failures\n";
    cout << "Standalone: " << sa_failures << " failures\n";
    cout << "Seeds:";
    for (auto& [s, c] : seed_tally) cout << " " << s << "=" << c;
    cout << "\n";
    cout << "Time: red=" << total_red << "s  rt_exp=" << total_exp
         << "s  standalone=" << total_sa << "s\n";
    return (rt_failures + sa_failures) > 0 ? 1 : 0;
}
