// Dump the input iDT edge-length distribution (before any PALC flips).
// For every alive non-bigon edge in the iDT, count how many distinct
// length values appear.  Surface metric is unit-edge piecewise-flat,
// so iDT length² is integer-valued in the Eisenstein lattice; we round
// to integer ℓ² to report the count of distinct length classes.
//
// Usage:
//   idt_lengths <spiral_name>
//
// Output (one line):
//   <name> n_alive=<k> n_distinct=<m> n_self=<s> n_multi_pair=<p> ℓ²=<a>,<b>,...
//
// where:
//   k = number of alive (oriented) non-bigon undirected edges
//   m = number of distinct ℓ² values
//   s = number of self-loop edges
//   p = number of distinct multi-edge vertex-pairs (count > 1)

#include "fullerenes/triangulation.hh"
#include "fullerenes/delaunay.hh"
#include "fullerenes/spiral.hh"

#include <cstdio>
#include <cmath>
#include <string>
#include <set>
#include <map>
#include <vector>

using namespace std;

int main(int argc, char** argv) {
  if (argc < 2) {
    fprintf(stderr, "Usage: %s <spiral_name>\n", argv[0]);
    return 1;
  }
  string name = argv[1];
  spiral_nomenclature sn(name);
  Triangulation T(sn);
  auto D = DelaunayTriangulation::compute(T);

  set<long long> lsq_set;
  vector<long long> lsq_list;
  int n_alive_nonbigon = 0, n_self = 0;
  map<pair<int,int>, int> by_pair;
  for (int h = 0; h < D.nh; h += 2) {
    if (!D.alive(h)) continue;
    int u = D.he_origin[h], v = D.dest(h);
    if (u == v) { n_self++; continue; }
    if (D.he_face[h] == D.he_face[h ^ 1]) continue;  // bigon
    n_alive_nonbigon++;
    double L = D.he_length[h];
    long long lsq = (long long)std::llround(L * L);
    lsq_set.insert(lsq);
    lsq_list.push_back(lsq);
    by_pair[{min(u,v), max(u,v)}]++;
  }
  int n_multi_pairs = 0;
  for (auto& [k, c] : by_pair) if (c > 1) n_multi_pairs++;

  printf("%s n_alive=%d n_distinct=%zu n_self=%d n_multi_pair=%d  L²=",
         name.c_str(), n_alive_nonbigon, lsq_set.size(),
         n_self, n_multi_pairs);
  bool first = true;
  for (auto v : lsq_set) { printf("%s%lld", first?"":",", v); first = false; }
  printf("\n");
  return 0;
}
