#include <stdlib.h>
#include <iostream>
#include "fullerenes/fullerenegraph.hh"
#include "fullerenes/triangulation.hh"
#include "fullerenes/unfold.hh"

using namespace std;

int main(int ac, char **av)
{
  if(ac < 14){
    cerr << "Usage: " << av[0] << " N rspi[1..12] K L\n"
         << "  N:         number of vertices in fullerene\n"
         << "  rspi[1..12]: 12 ring spiral pentagon indices (1-based)\n"
         << "  K L:       Goldberg-Coxeter parameters\n";
    return 1;
  }

  int N = strtol(av[1], 0, 0);
  vector<int> rspi(12);
  for(int i = 0; i < 12; i++) rspi[i] = strtol(av[2+i], 0, 0) - 1;

  int K = 1, L = 0;
  if(ac >= 16){
    K = strtol(av[14], 0, 0);
    L = strtol(av[15], 0, 0);
  }

  jumplist_t jumps;
  FullereneGraph g(N, rspi, jumps);
  Triangulation dual(g.dual_graph(6));

  cout << "Input: C" << N << " with RSPI = " << (rspi+1) << endl;
  cout << "GC(" << K << "," << L << ") transform: T = " << (K*K + K*L + L*L) << endl;

  Triangulation result = dual.GCtransform(K, L);
  int N_new = (result.N - 2) * 2;

  cout << "Result: C" << N_new << " with " << result.N << " dual nodes" << endl;

  // Get spiral of result
  vector<int> spiral;
  jumplist_t result_jumps;
  FullereneDual Tresult(result);
  if(Tresult.get_rspi(spiral, result_jumps)){
    cout << "RSPI = " << (spiral+1) << endl;
    if(!result_jumps.empty())
      cout << "Jumps = " << result_jumps << endl;
  } else {
    cerr << "Warning: Could not compute spiral for result." << endl;
  }

  return 0;
}
