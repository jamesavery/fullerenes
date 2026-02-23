#include <cstdio>
#include <vector>
#include <algorithm>
#include "fullerenes/triangulation.hh"
#include "fullerenes/spiral.hh"

using namespace std;

struct Entry { int k, l, T, N_carbon; };

int main()
{
  vector<int> C28_spiral = {5,5,5,6,5,6,5,6,5,5,5,5,5,5,5,6};
  Triangulation C28dual(C28_spiral);
  int V0 = C28dual.N;  // 16
  int max_N = 100800;

  // Collect all (k,l) pairs with N_carbon <= max_N
  vector<Entry> all;

  // Halma: GC(k,0)
  for(int k = 1; k <= 60; k++){
    int T = k*k;
    int N_carbon = 2*(V0-2)*T;
    if(N_carbon > max_N) break;
    all.push_back({k, 0, T, N_carbon});
  }

  // Chiral: GC(k,l), 1 <= l < k
  for(int k = 2; k <= 65; k++){
    for(int l = 1; l < k; l++){
      int T = k*k + k*l + l*l;
      int N_carbon = 2*(V0-2)*T;
      if(N_carbon > max_N) continue;
      all.push_back({k, l, T, N_carbon});
    }
  }

  // Sort by N_carbon
  sort(all.begin(), all.end(), [](const Entry& a, const Entry& b){
    return a.N_carbon < b.N_carbon;
  });

  // Pick ~100 evenly spaced entries
  int n_total = all.size();
  int n_pick = 100;
  vector<int> indices;
  for(int i = 0; i < n_pick; i++){
    int idx = (int)((long long)i * (n_total - 1) / (n_pick - 1));
    if(indices.empty() || indices.back() != idx)
      indices.push_back(idx);
  }

  fprintf(stderr, "Total (k,l) pairs: %d, picking %d\n", n_total, (int)indices.size());

  printf("[\n");
  for(size_t ii = 0; ii < indices.size(); ii++){
    const Entry& e = all[indices[ii]];

    fprintf(stderr, "[%zu/%zu] GC(%d,%d) N=%d ... ",
            ii+1, indices.size(), e.k, e.l, e.N_carbon);

    Triangulation result = C28dual.GCtransform(e.k, e.l);

    spiral_nomenclature sn(result, spiral_nomenclature::FULLERENE,
                           spiral_nomenclature::TRIANGULATION, true);
    string spiral_str = sn.to_string();

    fprintf(stderr, "%s\n", spiral_str.c_str());

    printf("  {\"k\":%d, \"l\":%d, \"T\":%d, \"N_carbon\":%d, \"spiral\":\"%s\"}%s\n",
           e.k, e.l, e.T, e.N_carbon, spiral_str.c_str(),
           ii == indices.size()-1 ? "" : ",");
  }
  printf("]\n");

  return 0;
}
