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

  BuckyGen::buckygen_queue Q = BuckyGen::start(N,0,chunk_index,chunk_number);
  
  Triangulation G;
  PlanarGraph dual;
  int i=0;

  FILE *output = fopen(("output/c"+to_string(N)+".rspi"
			+(chunk_number==1?"":("."+to_string(chunk_number)+"-"+to_string(chunk_index)))).c_str(),"wb");

  if(!output) abort();

  while(BuckyGen::next_fullerene(Q,G)){
    uint8_t rspi[12];
    vector<int> spiral;
    Triangulation::jumplist_t jumps;
    
    i++;
    G = Triangulation(G.neighbours,true);
    bool spiral_OK = G.get_spiral(spiral,jumps,false);
    if(!spiral_OK){
      vector<int> rspi2(12);
      FullereneGraph F(G.dual_graph());

      F.get_canonical_general_spiral_from_fg(rspi2,jumps);

      cout << jumps << endl;
      cout << rspi2 << endl;
      assert(spiral_OK);
    }

    get_rspi(spiral,rspi);

    if(i % 100 == 0) cout << "Isomer " << i << ": " << vector<int>(rspi,rspi+12) << endl;

    fwrite(rspi,12,1,output);

  }
  printf("Generated %d graphs (%d,%d)\n",i,int(G.neighbours.size()),int(dual.neighbours.size()));
  fclose(output);
  
  return 0;
}
