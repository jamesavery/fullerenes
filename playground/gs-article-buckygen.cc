#include <limits.h>
#include "libgraph/triangulation.hh"
#include "contrib/buckygen-wrapper.hh"

#include <iostream>
#include <chrono>
#include <ctime>

using namespace std;
typedef Triangulation::jumplist_t jumplist_t;

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
  
  cout << "Nv = " << N << ";\n";
  vector<int> isomers_with_jumps;
  vector< vector<int> > rspis_with_jumps;
  vector< jumplist_t > all_jumps;

  bool only_do_buckygen = false;  
  auto cpu_start = clock();
  while(BuckyGen::next_fullerene(Q,G)){
    vector<int> rspi;
    vector<int> spiral;
    Triangulation::jumplist_t jumps;

    i++;
    if(i%100000 == 0) cerr << "Reached isomer " << i << ".\n";
    if(!only_do_buckygen){
      FullereneDual F(G.neighbours);
      F.get_spiral(spiral,jumps,true,true);
      //    if(jumps.size()!=0) fprintf(stderr,"Isomer %ld has jump length %ld\n",i,jumps.size());
      if(jumps.size()!=0){
	vector<int> rspi = rspi_from_fullerene_spiral(spiral);
	for(auto &j: jumps) jumps.first++;
	
	cout << "(* BuckyGen isomer number " << i << " has a jump in its canonical general spiral: *)\n";
	cout << "rspi"<<i<<" = " << rspi << "; jumps"<<i<<" = " << jumps << ";\n";

	isomers_with_jumps.push_back(i);
	rspis_with_jumps.push_back(rspi);
	all_jumps.push_back(jumps);
      }
    }
    //    cout << "SPIRAL: " << spiral << "; JUMPS: " << jumps << endl;
    
    // dual = G.dual_graph();
    //    G.update_from_neighbours();
    //    if(!G.is_consistently_oriented()) abort();
    //    PlanarGraph dual(G.dual_graph());

    // printf("Graph %d with %d nodes is %sconsistently oriented\n",i,G.N,G.is_consistently_oriented()?"":"not ");
    //	printf("Graph %d is %sconsistently oriented and has %ld perfect matchings\n",i,G.is_consistently_oriented()?"":"not ",G.dual_graph().count_perfect_matchings());
  }
  BuckyGen::stop(Q);
  auto cpu_end = clock();
  size_t Nisomers = i;
  double cputime = (cpu_end-cpu_start) * 1.0 / CLOCKS_PER_SEC;
  cout << "Nisomers"<<chunk_index<<" = " << i << ";\n"
       << "cputime"<<chunk_index<<"  = " << cputime <<";\n"
       << "rate"<<chunk_index<<"     = " << cputime/(Nisomers*N) <<";\n"
       << "isomersWithJumps"<<chunk_index<<" = " << isomers_with_jumps << ";\n"
       << "rspisWithJumps"<<chunk_index<<"   = " << rspis_with_jumps << ";\n"
       << "allJumps = " << all_jumps << ";\n";
  
  //  printf("Generated %ld graphs (%d,%d)\n",i,int(G.neighbours.size()),int(dual.neighbours.size()));
  
  return 0;
}
