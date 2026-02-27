// Fast-path-only reduction test for large fullerenes.
// Usage: ./test_fast [Nmax]   (default: 100)
// Only runs O(N) ReducibleDual — no reference comparison.
// Reports buckygen and reduction times separately.

#include "buckinverse.hh"
#include "fullerenes/buckygen-wrapper.hh"
#include <iostream>
#include <map>
#include <chrono>

using namespace buckinverse;
using namespace std;
using clk = chrono::steady_clock;

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
    double total_gen = 0, total_red = 0;

    for (int N = 32; N <= Nmax; N += 2) {
        auto Q = BuckyGen::start(N, false);
        Graph G;
        int isomer_count = 0, ok = 0, fail = 0;
        double dt_gen = 0, dt_red = 0;

        while (true) {
            auto tg0 = clk::now();
            bool has_next = BuckyGen::next_fullerene(Q, G);
            auto tg1 = clk::now();
            dt_gen += chrono::duration<double>(tg1 - tg0).count();
            if (!has_next) break;

            isomer_count++;
            total++;

            auto tr0 = clk::now();
            ReducibleDual rd(G);
            SeedType seed = rd.reduceToSeed();
            auto tr1 = clk::now();
            dt_red += chrono::duration<double>(tr1 - tr0).count();

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
        BuckyGen::stop(Q);

        total_gen += dt_gen;
        total_red += dt_red;

        int Nv = N / 2 + 2;  // dual vertex count
        double us_gen = isomer_count > 0 ? 1e6 * dt_gen / isomer_count : 0;
        double us_red = isomer_count > 0 ? 1e6 * dt_red / isomer_count : 0;
        double us_vert = isomer_count > 0 ? 1e6 * dt_red / isomer_count / Nv : 0;
        cout << "  C" << N << ": " << isomer_count << " isomers, "
             << ok << " ok, " << fail << " fail  "
             << "gen=" << us_gen << "  red=" << us_red << " us/isomer"
             << "  (" << us_vert << " us/vertex)\n" << flush;
    }

    cout << "\nTotal: " << total << " isomers, " << failures << " failures\n";
    cout << "Seeds:";
    for (auto& [s, c] : seed_tally) cout << " " << s << "=" << c;
    cout << "\n";
    cout << "Time: gen=" << total_gen << "s  red=" << total_red << "s\n";
    return failures > 0 ? 1 : 0;
}
