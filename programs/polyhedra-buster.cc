#include <limits.h>
#include "fullerenes/buckygen-wrapper.hh"
#include "fullerenes/triangulation.hh"
#include "fullerenes/polyhedron.hh"

#include <iostream>
#include <chrono>
#include <ctime>

int face_size(const Graph &g, node_t u, node_t v)
{
  int d = 1;
  node_t u0 = u;
  while(v != u0){
    node_t w = v;
    v = g.next_on_face(u,v);
    u = w;
    d++;
  }
  return d;
}


int main(int ac, char **av)
{
  size_t N                = strtol(av[1],0,0);

  if(ac<3 || N<20 || N&1){
    fprintf(stderr,"Syntax: %s <N:int> <RSPIfile:string> <dual:0|1>\n",av[0]);
  }

  const char *RSPIfilename = av[2];
  string RSPIline;  
  FILE *RSPIfile = fopen(RSPIfilename,"r");

  vector<Polyhedron> P0s, Ps;

  vector<int> all_rspi;
  
  while(getline(RSPIfile,RSPIline)){
    // NB: Only one space between numbers
    vector<int> RSPI = split<int>(RSPIline,string(" "),string(" \t\r\n"))+(-1);

    vector<int> spiral(N/2+2,6);
    jumplist_t  jumps;
    
    for(int i=0;i<12;i++){ spiral[RSPI[i]] = 5; all_rspi.push_back(RSPI[i]+1); }
    Triangulation T(spiral,jumps);

    string filename;
    stringstream s(filename);
   
    FullereneGraph g = T.dual_graph();
    g.layout2d = g.tutte_layout();
    Polyhedron P0(g,g.zero_order_geometry(),6);
    P0s.push_back(P0);    
      
    Polyhedron P = P0;
    P.points = g.optimized_geometry(P.points);

    Ps.push_back(P);
  }

  size_t M = Ps.size();
  size_t Nfaces = N/2+2;
  
  vector<size_t> neighbours(M*N*3), next_on_face(M*N*3),prev_on_face(M*N*3),
                 face_right(M*N*3), neighbours_shape{{M,N,3}}, points_shape{{M,N,3}};
  vector<double> points_start(M*N*3), points_opt(M*N*3);
  vector<vector<face_t>>
    pentagons(M,vector<face_t>(12)),
    hexagons (M,vector<face_t>(Nfaces-12));
  
  for(int i=0;i<M;i++){
    // Everything with neighbours-shape
    for(node_t u=0;u<N;u++)
      for(int j=0;j<3;j++){
	size_t index = i*N*3 + u*3 + j;
	node_t v = Ps[i].neighbours[u][j];
	neighbours[index]   = v;
	next_on_face[index] = Ps[i].next_on_face(u,v);
	prev_on_face[index] = Ps[i].prev_on_face(u,v);	
	points_start[index] = P0s[i].points[u][j];
	points_opt[index]   = Ps[i].points[u][j];
	face_right[index]   = face_size(Ps[i],u,v);
      }

    // Faces
    vector<face_t> faces = Ps[i].compute_faces_oriented();

    for(int j=0,npent=0,nhex=0;j<Nfaces;j++)
      if      (faces[j].size() == 5) pentagons[i][npent++] = faces[j];
      else if (faces[j].size() == 6) hexagons [i][nhex++]  = faces[j];
  }
  cerr << "from numpy import array\n\n";
  cerr << "rspi_shape = " << vector<size_t>{{M,12}} << ";\n";
  cerr << "neighbours_shape = " << neighbours_shape << ";\n";
  
  cerr << "rspi       = array(" << all_rspi << ").reshape(rspi_shape);\n\n";  
  cerr << "neighbours   = array(" << neighbours << ").reshape(neighbours_shape);\n\n";
  cerr << "pentagons    = array(" << pentagons  << ");\n\n";
  cerr << "hexagons     = array(" << hexagons   << ");\n\n";  
  cerr << "next_on_face = array(" << next_on_face << ").reshape(neighbours_shape);\n\n";
  cerr << "prev_on_face = array(" << prev_on_face << ").reshape(neighbours_shape);\n\n";
  cerr << "face_right   = array(" << face_right << ").reshape(neighbours_shape);\n\n"; 
  cerr << "points_start = array(" << points_start << ").reshape(points_shape);\n\n";
  cerr << "points_opt   = array(" << points_opt << ").reshape(points_shape);\n\n";   
  
  return 0;
}
