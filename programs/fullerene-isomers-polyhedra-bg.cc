#include <limits.h>
#include "fullerenes/buckygen-wrapper.hh"
#include "fullerenes/triangulation.hh"
#include "fullerenes/polyhedron.hh"

#include <iostream>

int main(int ac, char **argv)
{
  
  if(ac<3){
    fprintf(stderr,"Syntax: %s <N:int> [output_dir] [IPR:0|1] [only_nontrivial:0|1]\n",argv[0]);
    return -1;
  }
  int N                = strtol(argv[1],0,0);     // Argument 1: Number of vertices N

  string output_dir   = ac>=3? argv[2] : "output";    // Argument 2: directory to output files to
  int IPR             = ac>=4? strtol(argv[3],0,0):0; // Argument 3: Only generate IPR fullerenes?
  int only_nontrivial = ac>=5? strtol(argv[4],0,0):0; // Argument 4: Only generate fullerenes with nontrivial symmetry group?
  
  ofstream failures((output_dir+"/failures.txt").c_str()); // output/failures.txt contains list of any fullerenes that failed optimization

  int i=0;  
  Triangulation dualG;
  BuckyGen::buckygen_queue Q = BuckyGen::start(N,IPR,only_nontrivial);  

  while(BuckyGen::next_fullerene(Q,dualG)){ // Generate all appropriate C_N isomer duals 
    i++;
    dualG.update();		            // Update metadata
    FullereneGraph G = dualG.dual_graph();  // Construct fullerene graph

    Polyhedron P = Polyhedron::fullerene_polyhedron(G); // Generate optimized fullerene geometry

    if(P.is_invalid()){
      vector<int> RSPI;
      jumplist_t  jumps;
      dualG.get_spiral(RSPI,jumps);
      failures << "Isomer " << i << " failed. RSPI="<<(RSPI+1) << ", jumps="<<jumps<<endl;
      continue;
    }

    // Output molecular geometry files
    string filename = output_dir+"/C"+to_string(N)+"-"+to_string(i);        
    Polyhedron::to_file(P,filename+".mol2");
    Polyhedron::to_file(P,filename+".xyz"); 
  }
  failures.close();
  
  return 0;
}
