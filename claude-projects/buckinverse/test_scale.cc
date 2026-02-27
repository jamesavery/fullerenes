// Scalable reduction-to-seed test for larger fullerenes.
// Usage: ./test_scale [Nmax]   (default: 60)
// Compares O(N^2) reference implementation with O(N) ReducibleDual.

#include "buckinverse.hh"
#include "fullerenes/buckygen-wrapper.hh"
#include <iostream>
#include <map>
#include <set>
#include <cassert>
#include <chrono>

using namespace buckinverse;
using namespace std;

// --- O(N^2) reference: iterative allInvSites + invertReduction ---

static SeedType reduceToSeedRef(const Graph& g) {
    Graph current = g;
    // Keep reducing as long as valid sites exist.
    // Seeds have no valid inversion sites, so this terminates at a seed.
    while (true) {
        auto sites = allInvSites(current, 5);
        bool found = false;
        for (const auto& site : sites) {
            Graph next = invertReduction(current, site);
            if (next.N > 0) { current = next; found = true; break; }
        }
        if (!found) break;
    }
    switch (current.N) {
        case 12: return SeedType::C20;
        case 16: return SeedType::C28;
        case 17: return SeedType::C30;
        default: return SeedType::NotASeed;
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
    cout << "  [ref = O(N^2) allInvSites+invertReduction,  "
         << "fast = O(N) ReducibleDual]\n\n";

    map<string, int> seed_tally;
    int total = 0, failures = 0;
    double total_ref = 0, total_fast = 0;

    for (int N = 32; N <= Nmax; N += 2) {
        auto Q = BuckyGen::start(N, false);
        Graph G;
        int isomer_count = 0, ok = 0, fail = 0;
        double dt_ref = 0, dt_fast = 0;

        while (BuckyGen::next_fullerene(Q, G)) {
            isomer_count++;
            total++;

            // O(N) ReducibleDual
            auto t1 = chrono::steady_clock::now();
            ReducibleDual rd(G);
            SeedType seed_fast = rd.reduceToSeed();
            auto t2 = chrono::steady_clock::now();
            dt_fast += chrono::duration<double>(t2 - t1).count();

            // O(N^2) reference
            auto t3 = chrono::steady_clock::now();
            SeedType seed_ref = reduceToSeedRef(G);
            auto t4 = chrono::steady_clock::now();
            dt_ref += chrono::duration<double>(t4 - t3).count();

            if (seed_fast == SeedType::NotASeed) {
                // Diagnose: extract stuck graph and check what old code finds
                Graph stuck = rd.toGraph();
                auto sites = allInvSites(stuck, 5);
                cerr << "  FAIL C" << N << " #" << isomer_count
                     << ": stuck at " << rd.N() << "v, old finds "
                     << sites.size() << " sites";
                for (auto& s : sites) cerr << " " << s.kind.toString();
                cerr << "\n";
                fail++;
                failures++;
            } else if (seed_ref == SeedType::NotASeed) {
                cerr << "  FAIL C" << N << " #" << isomer_count
                     << ": reference got NotASeed\n";
                fail++;
                failures++;
            } else {
                // Verify C28/C30 are genuine seeds (no further reductions possible)
                if (seed_fast == SeedType::C28 || seed_fast == SeedType::C30) {
                    Graph g = rd.toGraph();
                    auto sites = allInvSites(g, 5);
                    if (!sites.empty()) {
                        cerr << "  FALSE SEED C" << N << " #" << isomer_count
                             << ": " << seedName(seed_fast) << " has "
                             << sites.size() << " sites:";
                        for (auto& s : sites) cerr << " " << s.kind.toString();
                        cerr << "\n";
                        fail++;
                        failures++;
                        continue;
                    }
                }
                seed_tally[seedName(seed_fast)]++;
                ok++;
            }
        }
        BuckyGen::stop(Q);

        total_ref += dt_ref;
        total_fast += dt_fast;

        double speedup = dt_ref > 0 ? dt_ref / dt_fast : 0;
        cout << "  C" << N << ": " << isomer_count << " isomers, "
             << ok << " ok, " << fail << " fail"
             << "  ref=" << dt_ref << "s  fast=" << dt_fast << "s"
             << "  (" << speedup << "x)\n";
        cout << flush;
    }

    cout << "\nTotal: " << total << " isomers, " << failures << " failures\n";
    cout << "Seeds:";
    for (auto& [s, c] : seed_tally) cout << " " << s << "=" << c;
    cout << "\n";
    cout << "Time: ref=" << total_ref << "s  fast=" << total_fast << "s"
         << "  (" << (total_ref > 0 ? total_ref / total_fast : 0) << "x)\n";
    return failures > 0 ? 1 : 0;
}
