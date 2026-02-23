#include <chrono>
#include <cstdio>
#include <vector>
#include <algorithm>
#include "fullerenes/triangulation.hh"
#include "fullerenes/spiral.hh"

using namespace std;
using namespace std::chrono;

struct Entry { int k, l, T, N_carbon; };

int main()
{
  vector<int> C28_spiral = {5,5,5,6,5,6,5,6,5,5,5,5,5,5,5,6};
  Triangulation C28dual(C28_spiral);
  int V0 = C28dual.N;  // 16
  int max_N = 100800;

  // Collect all (k,l) pairs with N_carbon <= max_N
  vector<Entry> all;
  for(int k = 1; k <= 60; k++){
    int T = k*k;
    int N_carbon = 2*(V0-2)*T;
    if(N_carbon > max_N) break;
    all.push_back({k, 0, T, N_carbon});
  }
  for(int k = 2; k <= 65; k++){
    for(int l = 1; l < k; l++){
      int T = k*k + k*l + l*l;
      int N_carbon = 2*(V0-2)*T;
      if(N_carbon > max_N) continue;
      all.push_back({k, l, T, N_carbon});
    }
  }

  sort(all.begin(), all.end(), [](const Entry& a, const Entry& b){
    return a.N_carbon < b.N_carbon;
  });

  // Pick 100 evenly spaced entries
  int n_total = all.size();
  int n_pick = 100;
  vector<int> indices;
  for(int i = 0; i < n_pick; i++){
    int idx = (int)((long long)i * (n_total - 1) / (n_pick - 1));
    if(indices.empty() || indices.back() != idx)
      indices.push_back(idx);
  }

  fprintf(stderr, "Total (k,l) pairs: %d, picking %d\n", n_total, (int)indices.size());

  // JSON header
  printf("{\n  \"base_spiral\": [5,5,5,6,5,6,5,6,5,5,5,5,5,5,5,6],\n  \"entries\": [\n");

  for(size_t ii = 0; ii < indices.size(); ii++){
    const Entry& e = all[indices[ii]];

    fprintf(stderr, "[%zu/%zu] GC(%d,%d) N=%d ...\n",
            ii+1, indices.size(), e.k, e.l, e.N_carbon);

    // Generate graph via GC transform
    Triangulation graph = C28dual.GCtransform(e.k, e.l);

    // Time spiral computation: graph -> canonical spiral
    auto t0 = steady_clock::now();
    spiral_nomenclature sn(graph, spiral_nomenclature::FULLERENE,
                           spiral_nomenclature::TRIANGULATION, true);
    string spiral_str = sn.to_string();
    auto t1 = steady_clock::now();
    double time_spiral_us = duration<double, micro>(t1 - t0).count();

    // Time spiral windup: spiral -> graph
    auto t2 = steady_clock::now();
    Triangulation roundtrip(sn);
    auto t3 = steady_clock::now();
    double time_windup_us = duration<double, micro>(t3 - t2).count();

    fprintf(stderr, "  spiral: %.1fms  windup: %.1fms  %s\n",
            time_spiral_us/1000.0, time_windup_us/1000.0, spiral_str.c_str());

    printf("    {\"k\":%d, \"l\":%d, \"T\":%d, \"N_carbon\":%d, \"N_dual\":%d, "
           "\"spiral\":\"%s\", \"time_spiral_us\":%.1f, \"time_windup_us\":%.1f}%s\n",
           e.k, e.l, e.T, e.N_carbon, graph.N,
           spiral_str.c_str(), time_spiral_us, time_windup_us,
           ii == indices.size()-1 ? "" : ",");
  }

  printf("  ]\n}\n");
  return 0;
}
