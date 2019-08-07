#include <limits.h>
#include "fullerenes/buckygen-wrapper.hh"
#include "fullerenes/triangulation.hh"
#include "fullerenes/polyhedron.hh"

#include <iostream>
#include <chrono>
#include <ctime>

int main(int ac, char **av)
{
  int N                = strtol(av[1],0,0);

  if(ac<3 || N<20 || N&1){
    fprintf(stderr,"Syntax: %s <N:int> <RSPIfile:string> <dual:0|1>\n",av[0]);
  }

  const char *RSPIfilename = av[2];
  string RSPIline;  
  FILE *RSPIfile = fopen(RSPIfilename,"r");

  while(getline(RSPIfile,RSPIline)){
    vector<int> RSPI = split<int>(RSPIline,string(" "),string(" \t\r\n"))+(-1);
    cout << RSPI << endl;

    vector<int> spiral(N/2+2,6);
    jumplist_t  jumps;
    
    for(int i=0;i<12;i++) spiral[RSPI[i]] = 5;
    Triangulation T(spiral,jumps);

    string filename;
    stringstream s(filename);
    s << "output/C"<<N << "-" << (RSPI+1);
    cerr << s.str() << endl;
    FullereneGraph g = T.dual_graph();
    g.layout2d = g.tutte_layout();
    Polyhedron P0(g,g.zero_order_geometry(),6);    
    Polyhedron::to_file(P0,s.str()+"-P0.mol2");
      
    Polyhedron P = P0;
    P.points = g.optimized_geometry(P.points);
    Polyhedron::to_file(P,s.str()+"-P.mol2"); 
    
    // P.move_to_origin();
    // // P.align_with_axes();   
    
    // Polyhedron D(P.dual());    
    // D.layout2d = D.tutte_layout();
    // D.faces    = D.compute_faces(3,true);
    // D.face_max = 3;
    // D.optimize();
  }

  
  return 0;
}
