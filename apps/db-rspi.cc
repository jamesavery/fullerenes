#include <inttypes.h>

#include "libgraph/triangulation.hh"
#include "libgraph/fullerenegraph.hh"
#include "contrib/buckygen-wrapper.hh"

void get_rspi(const vector<int>& spiral, uint8_t rspi[12])
{
  for(int i=0, j=0;i<spiral.size();i++) if(spiral[i] == 5) rspi[j++] = i;
}


int main(int ac, char **av)
{
  int N            = ac>=2? strtol(av[1],0,0) : 20;
  int chunk_number = ac>=3? strtol(av[2],0,0) : 1;
  int chunk_index  = ac>=4? strtol(av[3],0,0) : 0;
  int ipr              = ac>=5? strtol(av[4],0,0) : 0;
  int only_nontrivial  = ac>=6? strtol(av[5],0,0) : 0;

  BuckyGen::buckygen_queue Q = BuckyGen::start(N,ipr,only_nontrivial,chunk_index,chunk_number);
  
  Triangulation G;
  PlanarGraph dual;
  int i=0;

  string outputfile = ("output/c"+to_string(N)+(ipr?"-IPR":"")+(only_nontrivial?"-nontrivial":"")+".rspi"
		       +(chunk_number==1?"":("."+to_string(chunk_number)+"-"+to_string(chunk_index))));
  FILE *output = fopen(outputfile.c_str(),"wb");

  if(!output){
	 perror(outputfile.c_str());
	 abort();
  }

  cerr << "Writing to " << outputfile << "\n";

  while(BuckyGen::next_fullerene(Q,G)){
    uint8_t rspi[12];
    vector<int> spiral;
    Triangulation::jumplist_t jumps;
    
    i++;
    G = Triangulation(G.neighbours,true);
    bool spiral_OK = G.get_spiral(spiral,jumps,true,false,false);
    if(!spiral_OK){
      vector<int> rspi2(12);
      FullereneDual FD(G);

      spiral_OK = FD.get_rspi(rspi2,jumps,true,true);

      cout << jumps << endl;
      cout << rspi2 << endl;
      assert(spiral_OK);
    }

    get_rspi(spiral,rspi);

    if(i % 1000 == 0) cout << "Isomer " << i << ": " << vector<int>(rspi,rspi+12) << endl;

    fwrite(rspi,12,1,output);

  }
  printf("Generated %d graphs (%d,%d)\n",i,int(G.neighbours.size()),int(dual.neighbours.size()));
  fclose(output);
  
  return 0;
}
