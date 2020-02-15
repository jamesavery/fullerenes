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

  if(ac<3 || N<20 || N==22 || N&1){
    fprintf(stderr,"Syntax: %s <N:int> <RSPIfile:string> [output_dir] [dual:0|1]\n",av[0]);
  }

  const char *RSPIfilename = av[2];
  string output_dir   = ac>=4? av[3] : "output";
  bool  do_dual = ac>=5? strtol(av[4],0,0) : 0;
  
  string RSPIline;  
  FILE *RSPIfile = fopen(RSPIfilename,"r");

  int i=0;
  ofstream failures((output_dir+"/failures.txt").c_str());
  while(getline(RSPIfile,RSPIline)){
    i++;
    vector<int> RSPI = split<int>(RSPIline,string(" "),string(" \t\r\n"))+(-1);
    cout << RSPI << endl;

    vector<int> spiral(N/2+2,6);
    jumplist_t  jumps;
    
    for(int i=0;i<12;i++) spiral[RSPI[i]] = 5;
    Triangulation T(spiral,jumps);

    string filename;
    stringstream s(filename);
    s << output_dir + "/C"<<N << "-" << i;
    cerr << s.str() << endl;
    FullereneGraph g = T.dual_graph();
    g.layout2d = g.tutte_layout();
    Polyhedron P0(g,g.zero_order_geometry(),6);    
    Polyhedron::to_file(P0,s.str()+"-P0.mol2");
      
    Polyhedron P = P0;
    P.points = g.optimized_geometry(P.points);
    bool has_nans = false;
    for(auto p: P.points){
      if(std::isnan(p[0])||std::isnan(p[1])||std::isnan(p[1])) has_nans = true;
    }
    if(has_nans){
      failures << i << ":" << (RSPI+1) << endl;
      continue;
    }
    P.move_to_origin();
    P.align_with_axes();       
    Polyhedron::to_file(P,s.str()+"-P.mol2");
    Polyhedron::to_file(P,s.str()+"-P.xyz"); 
    
    
    
    
    // Polyhedron D(P.dual());    
    // D.layout2d = D.tutte_layout();
    // D.faces    = D.compute_faces(3,true);
    // D.face_max = 3;
    // D.optimize();

    
  }
  failures.close();
  
  return 0;
}
