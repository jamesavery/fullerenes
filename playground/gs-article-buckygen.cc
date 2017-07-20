#include <limits.h>
#include "libgraph/triangulation.hh"
#include "contrib/buckygen-wrapper.hh"

#include <iostream>
#include <chrono>

vector<int> rspi_from_fullerene_spiral(const vector<int>& spiral)
{
  vector<int> rspi(12);
  for(int i=0,j=0;i<spiral.size();i++){
    if(spiral[i] == 5) rspi[j++] = i+1;
    assert(j<=12);
  }
  return rspi;
}

int main(int ac, char **av)
{
  int N = ac>=2? strtol(av[1],0,0) : 20;
  bool  IPR = ac>=3? strtol(av[2],0,0) : 0;
  bool  only_nontrivial = ac>=4? strtol(av[3],0,0) : 0;
  size_t chunk_number   = ac>=5? strtol(av[4],0,0) : 1;
  size_t chunk_index    = ac>=6? strtol(av[5],0,0) : 0;

  size_t i=0;
  Triangulation G;
  PlanarGraph dual;

  BuckyGen::buckygen_queue Q = BuckyGen::start(N,IPR,only_nontrivial,
					       chunk_index,chunk_number);
    

  while(BuckyGen::next_fullerene(Q,G)){
    vector<int> rspi;
    vector<int> spiral;
    Triangulation::jumplist_t jumps;

    i++;
    //    if(i%100000 == 0)
    if(i==86325){
      G = Triangulation(G.neighbours,true);
      FullereneDual(G).get_spiral(spiral,jumps,true,true,true,true);
      //    if(jumps.size()!=0) fprintf(stderr,"Isomer %ld has jump length %ld\n",i,jumps.size());
      if(jumps.size()!=0) cout << "jumps="<<jumps << endl;
    }
    //   cout << "RSPI: " << rspi_from_fullerene_spiral(spiral) << "; JUMPS: " << jumps << endl;
    //    cout << "SPIRAL: " << spiral << "; JUMPS: " << jumps << endl;
    
    // dual = G.dual_graph();
    //    G.update_from_neighbours();
    //    if(!G.is_consistently_oriented()) abort();
    //    PlanarGraph dual(G.dual_graph());

    // printf("Graph %d with %d nodes is %sconsistently oriented\n",i,G.N,G.is_consistently_oriented()?"":"not ");
    //	printf("Graph %d is %sconsistently oriented and has %ld perfect matchings\n",i,G.is_consistently_oriented()?"":"not ",G.dual_graph().count_perfect_matchings());
  }
  BuckyGen::stop(Q);

  printf("Generated %ld graphs (%d,%d)\n",i,int(G.neighbours.size()),int(dual.neighbours.size()));
  
  return 0;
}
