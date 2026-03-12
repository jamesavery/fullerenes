#include "fullerenes/planargraph.hh"
#include <cstdio>
using namespace std;

int main(int ac, char** av) {
  int target = ac > 1 ? atoi(av[1]) : 3586;
  string path = string(FULLERENE_ROOT) + "/data/triangulations_15.pl";
  FILE* f = fopen(path.c_str(), "rb");
  fseek(f, 15, SEEK_SET);
  
  for (int idx = 0; idx <= target; idx++) {
    PlanarGraph g = PlanarGraph::read_hog_planarcode(f);
    if (g.N == 0) break;
    if (idx == target) {
      printf("Graph #%d, N=%d\n", idx, g.N);
      printf("Degrees:");
      for (int u = 0; u < g.N; u++) printf(" %d", (int)g[u].size());
      printf("\n");
      for (int u = 0; u < g.N; u++) {
        printf("  v%d (deg %d):", u, (int)g[u].size());
        for (int v : g[u]) printf(" %d", v);
        printf("\n");
      }
      int hist[20] = {};
      for (int u = 0; u < g.N; u++) hist[g[u].size()]++;
      printf("Degree histogram:");
      for (int d = 0; d < 20; d++) if (hist[d]) printf(" deg%d=%d", d, hist[d]);
      printf("\n");
      int gb = 0;
      for (int u = 0; u < g.N; u++) gb += 6 - (int)g[u].size();
      printf("sum(6-d_v) = %d\n", gb);
    }
  }
  fclose(f);
  return 0;
}
