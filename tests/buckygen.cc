#include "libgraph/triangulation.hh"
#include "contrib/buckygen-wrapper.hh"

int main(int ac, char **av)
{
  int N = ac>=2? strtol(av[1],0,0) : 20;
  BuckyGen::buckygen_queue Q = BuckyGen::start(N,0);
  
  Triangulation G;
  PlanarGraph dual;
  int i=0;
  while(BuckyGen::next_fullerene(Q,G)){
    i++;
    G = Triangulation(G.neighbours,true);
    dual = G.dual_graph();
    //    G.update_from_neighbours();
    //    if(!G.is_consistently_oriented()) abort();
    //    PlanarGraph dual(G.dual_graph());

    // printf("Graph %d with %d nodes is %sconsistently oriented\n",i,G.N,G.is_consistently_oriented()?"":"not ");
    //	printf("Graph %d is %sconsistently oriented and has %ld perfect matchings\n",i,G.is_consistently_oriented()?"":"not ",G.dual_graph().count_perfect_matchings());
  }
  printf("Generated %d graphs (%d,%d)\n",i,int(G.neighbours.size()),int(dual.neighbours.size()));
  
  return 0;
}
