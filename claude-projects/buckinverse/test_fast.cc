// Fast-path-only reduction test for large fullerenes.
// Usage: ./test_fast [Nmax]   (default: 100)
// Only runs O(N) ReducibleDual — no reference comparison.

#include "buckinverse.hh"
#include "fullerenes/buckygen-wrapper.hh"
#include <iostream>
#include <map>
#include <chrono>

using namespace buckinverse;
using namespace std;

int main(int argc, char* argv[]) {
    int Nmax = 100;
    if (argc > 1) Nmax = atoi(argv[1]);
    if (Nmax < 32 || Nmax % 2 != 0) {
        cerr << "Usage: " << argv[0] << " [Nmax]  (even, >= 32, default 100)\n";
        return 1;
    }

    cout << "ReducibleDual reduction: C32-C" << Nmax << "\n\n";

    map<string, int> seed_tally;
    int total = 0, failures = 0;
    double total_time = 0;

    for (int N = 32; N <= Nmax; N += 2) {
        auto Q = BuckyGen::start(N, false);
        Graph G;
        int isomer_count = 0, ok = 0, fail = 0;

        auto t0 = chrono::steady_clock::now();
        while (BuckyGen::next_fullerene(Q, G)) {
            isomer_count++;
            total++;

            ReducibleDual rd(G);
            SeedType seed = rd.reduceToSeed();

            if (seed == SeedType::NotASeed) {
                cerr << "  FAIL C" << N << " #" << isomer_count
                     << ": stuck at " << rd.N() << "v\n";
                fail++;
                failures++;
            } else {
                seed_tally[seedName(seed)]++;
                ok++;
            }
        }
        auto t1 = chrono::steady_clock::now();
        BuckyGen::stop(Q);

        double dt = chrono::duration<double>(t1 - t0).count();
        total_time += dt;

        double us_per = isomer_count > 0 ? 1e6 * dt / isomer_count : 0;
        cout << "  C" << N << ": " << isomer_count << " isomers, "
             << ok << " ok, " << fail << " fail  "
             << dt << "s  (" << us_per << " us/isomer)\n" << flush;
    }

    cout << "\nTotal: " << total << " isomers, " << failures << " failures\n";
    cout << "Seeds:";
    for (auto& [s, c] : seed_tally) cout << " " << s << "=" << c;
    cout << "\n";
    cout << "Time: " << total_time << "s\n";
    return failures > 0 ? 1 : 0;
}
