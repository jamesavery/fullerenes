#include <chrono>
#include <cstdio>
#include <cstdlib>
#include <map>
#include <vector>
#include "fullerenes/triangulation.hh"
#include "fullerenes/symmetry.hh"
#include "fullerenes/buckygen-wrapper.hh"

using namespace std;
using namespace std::chrono;

int main(int argc, char* argv[])
{
  int max_N = 80;
  if(argc > 1) max_N = atoi(argv[1]);
  if(max_N < 20 || max_N % 2 != 0){ fprintf(stderr, "Usage: %s [max_N]  (even, >= 20)\n", argv[0]); return 1; }

  // Phase 1: Pre-compute all isomers into N-indexed storage.
  map<int, vector<Triangulation>> isomers;
  for(int N = 20; N <= max_N; N += 2){
    if(N == 22) continue; // No fullerenes with 22 vertices
    BuckyGen::buckygen_queue Q = BuckyGen::start(N, false, false);
    Graph g;
    while(BuckyGen::next_fullerene(Q, g))
      isomers[N].emplace_back(g);
    fprintf(stderr, "Generated N=%d: %zu isomers\n", N, isomers[N].size());
  }

  // Phase 2: Time symmetry computation per N.
  // For small N, repeat until we accumulate >= 100ms of wall time.
  const double min_time_ms = 100.0;

  struct Result { int N, n_isomers, repetitions; double total_us, per_isomer_us, per_vertex_us; };
  vector<Result> results;

  for(auto& [N, tris] : isomers){
    int n = tris.size();

    // Warm-up pass (not timed).
    for(auto& T : tris) Symmetry S(T);

    // Timed passes: repeat until total time >= min_time_ms.
    int reps = 0;
    auto t0 = steady_clock::now();
    double elapsed_ms = 0;
    do {
      for(auto& T : tris) Symmetry S(T);
      reps++;
      elapsed_ms = duration<double, milli>(steady_clock::now() - t0).count();
    } while(elapsed_ms < min_time_ms);

    double total_us = duration<double, micro>(steady_clock::now() - t0).count();
    double per_isomer_us = total_us / (reps * n);

    double per_vertex_us = per_isomer_us / N;
    results.push_back({N, n, reps, total_us / reps, per_isomer_us, per_vertex_us});
    fprintf(stderr, "N=%3d: %6d isomers, %4d reps, %.1f us/isomer, %.2f us/vertex\n",
            N, n, reps, per_isomer_us, per_vertex_us);
  }

  // Output JSON.
  printf("{\n  \"max_N\": %d,\n  \"entries\": [\n", max_N);
  for(size_t i = 0; i < results.size(); i++){
    auto& r = results[i];
    printf("    {\"N\":%d, \"n_isomers\":%d, \"repetitions\":%d, "
           "\"total_us\":%.1f, \"per_isomer_us\":%.1f, \"per_vertex_us\":%.2f}%s\n",
           r.N, r.n_isomers, r.repetitions,
           r.total_us, r.per_isomer_us, r.per_vertex_us,
           i + 1 < results.size() ? "," : "");
  }
  printf("  ]\n}\n");

  return 0;
}
