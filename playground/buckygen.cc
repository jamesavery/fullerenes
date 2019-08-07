#include <limits.h>
#include "libgraph/triangulation.hh"
#include "libgraph/fullerenegraph.hh"
#include "contrib/buckygen-wrapper.hh"

int main(int ac, char **av)
{
  int N = ac>=2? strtol(av[1],0,0) : 20;
  bool  IPR = ac>=3? strtol(av[2],0,0) : 0;
  bool  only_nontrivial = ac>=4? strtol(av[3],0,0) : 0;
  size_t chunk_number   = ac>=5? strtol(av[4],0,0) : 1;
  size_t chunk_index    = ac>=6? strtol(av[5],0,0) : 0;

  char graph_filename[100];
  char dual_filename[100];
  
  sprintf(graph_filename,"graph%03d.py",N);
  sprintf(dual_filename, "dual%03d.py",N);  
  
  // ofstream
  //   graph_output(graph_filename),
  //   dual_output(dual_filename);
  
  size_t i=0;
  Triangulation G;
  PlanarGraph dual;

  BuckyGen::buckygen_queue Q = BuckyGen::start(N,IPR,only_nontrivial,
					       chunk_index,chunk_number);    
  while(BuckyGen::next_fullerene(Q,G)){
    vector<int> rspi;
    jumplist_t jumps;

    i++;
    if(i%100000 == 0) fprintf(stderr,"Reached isomer %ld\n",i);
    //      G = Triangulation(G.neighbours,true);
    //      FullereneDual(G).get_rspi(rspi,jumps,true,true);

    FullereneDual FG(G);
    //    cout << FG.neighbours << "\n";
    FullereneGraph F = FG.dual();
    //    dual = G.dual_graph();
    //G.update_from_neighbours();
    //if(!G.is_consistently_oriented()) abort();
    //    PlanarGraph dual(G.dual_graph());

    // printf("Graph %d with %d nodes is %sconsistently oriented\n",i,G.N,G.is_consistently_oriented()?"":"not ");
    //	printf("Graph %d is %sconsistently oriented and has %ld perfect matchings\n",i,G.is_consistently_oriented()?"":"not ",G.dual_graph().count_perfect_matchings());
  }
  BuckyGen::stop(Q);

  fprintf(stderr,"Generated %ld graphs (%d,%d)\n",i,int(G.neighbours.size()),int(dual.neighbours.size()));
  
  return 0;
}
