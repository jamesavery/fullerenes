#include "libgraph/planargraph.hh"
#include "contrib/buckygen-wrapper.hh"

int main(int ac, char **av)
{
  int N = ac>=2? strtol(av[1],0,0) : 20;
  BuckyGen::buckygen_queue Q = BuckyGen::start(N,0);
  
  Graph G;
  while(BuckyGen::next_fullerene(Q,G)){
    //  G.update_from_neighbours();
    cout << G << endl;
  }

  return 0;
}
