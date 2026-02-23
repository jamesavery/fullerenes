#include <chrono>
#include <cstdio>
#include <vector>
#include "fullerenes/triangulation.hh"

using namespace std;
using namespace std::chrono;

struct Result { int k, l, T, N_dual, N_carbon; double time_us; };

int main()
{
  // C28 dual: 16 nodes, 12 pentagons + 4 hexagons
  vector<int> C28_spiral = {5,5,5,6,5,6,5,6,5,5,5,5,5,5,5,6};
  Triangulation C28dual(C28_spiral);
  int V0 = C28dual.N;  // 16
  int max_N = 100800; // GC(60,0) on C28 = 100,800 carbon atoms

  vector<Result> halma, chiral;

  // --- Halma series: GC(k,0), k=1..60 ---
  for(int k = 1; k <= 60; k++){
    int T = k*k;
    int N_carbon = 2*(V0-2)*T;
    if(N_carbon > max_N) break;

    auto t0 = steady_clock::now();
    Triangulation result = C28dual.GCtransform(k, 0);
    auto t1 = steady_clock::now();

    halma.push_back({k, 0, T, result.N, N_carbon,
                     duration<double, micro>(t1 - t0).count()});
    fprintf(stderr, "halma k=%d N=%d %.1fms\n", k, N_carbon,
            duration<double, milli>(t1 - t0).count());
  }

  // --- Chiral series: GC(k,l), 1 <= l < k ---
  for(int k = 2; k <= 65; k++){
    for(int l = 1; l < k; l++){
      int T = k*k + k*l + l*l;
      int N_carbon = 2*(V0-2)*T;
      if(N_carbon > max_N) continue;

      auto t0 = steady_clock::now();
      Triangulation result = C28dual.GCtransform(k, l);
      auto t1 = steady_clock::now();

      chiral.push_back({k, l, T, result.N, N_carbon,
                        duration<double, micro>(t1 - t0).count()});
      fprintf(stderr, "chiral k=%d l=%d N=%d %.1fms\n", k, l, N_carbon,
              duration<double, milli>(t1 - t0).count());
    }
  }

  // --- Output JSON ---
  auto print_entry = [](const Result& r, bool last){
    printf("    {\"k\":%d, \"l\":%d, \"T\":%d, \"N_dual\":%d, \"N_carbon\":%d, \"time_us\":%.1f}%s\n",
           r.k, r.l, r.T, r.N_dual, r.N_carbon, r.time_us, last ? "" : ",");
  };

  printf("{\n  \"halma\": [\n");
  for(size_t i = 0; i < halma.size(); i++)
    print_entry(halma[i], i == halma.size()-1);
  printf("  ],\n  \"chiral\": [\n");
  for(size_t i = 0; i < chiral.size(); i++)
    print_entry(chiral[i], i == chiral.size()-1);
  printf("  ]\n}\n");

  return 0;
}
