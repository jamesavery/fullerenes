#include "libgraph/planargraph.hh"
#include "contrib/buckygen-wrapper.hh"

int main(int ac, char **av)
{
  int N = ac>=2? strtol(av[1],0,0) : 20;
  BuckyGen::buckygen_queue Q = BuckyGen::start(N,0);
  
  Graph G;
  int i=0;
  while(BuckyGen::next_fullerene(Q,G)){
    i++;
    G.update_from_neighbours();
    printf("Graph %d is %sconsistently oriented\n",i,G.is_consistently_oriented()?"":"not ");
  }

  return 0;
}
