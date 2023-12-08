#include <limits.h>
#include "fullerenes/buckygen-wrapper.hh"
#include "fullerenes/triangulation.hh"
#include "fullerenes/polyhedron.hh"
#include "fullerenes/unfold.hh"

#include <iostream>
#include <chrono>
#include <ctime>


pair<int,int> middle_hexagon(const Triangulation &dual)
{
  vector<int> pentagons(12);
  for(node_t u=0,p=0;u<dual.N;u++) if(dual.neighbours[u].size()==5) pentagons[p++] = u;
  
  matrix<int> D(12, dual.N);
  for(node_t p=0;p<12;p++)
    dual.single_source_shortest_paths(pentagons[p],&(D[p*dual.N]));
  
  vector<int> Dmin{&D[0], &D[dual.N]};
  for(node_t u=0;u<dual.N;u++) 
    for(node_t p=0;p<12;p++)
      Dmin[u] = min(Dmin[u],D(p,u));
  // Find argmax of Dmin
  int mx=0, u_max = 0;
  for(node_t u=0;u<dual.N;u++) if(Dmin[u]>mx){ mx = Dmin[u]; u_max = u; }
  return {u_max, mx};
}


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
  if(ac<3){
    fprintf(stderr,"Syntax: %s <N:int> <RSPIfile:string> <dual:0|1>\n",av[0]);
    return -1;
  }
  size_t N                = strtol(av[1],0,0);
  assert(N != 22 && N>=20 && !(N&1));

  const char *RSPIfilename = av[2];
  string RSPIline;  
  FILE *RSPIfile = fopen(RSPIfilename,"r");

  vector<Triangulation>   Ts, Tlf;
  vector<Polyhedron> P0s, Ps;
  

  vector<int> all_rspi;
  
  while(getline(RSPIfile,RSPIline)){
    // NB: Only one space between numbers
    vector<int> RSPI = split<int>(RSPIline,string(" "),string(" \t\r\n"))+(-1);

    vector<int> spiral(N/2+2,6);
    jumplist_t  jumps;
    
    for(int i=0;i<12;i++){ spiral[RSPI[i]] = 5; all_rspi.push_back(RSPI[i]+1); }
    Triangulation T(spiral,jumps);
    //    T = T.sort_nodes();		// Move pentagons to dual node id 0,...,11.
    
    string filename;
    stringstream s(filename);
    
    FullereneGraph g = T.dual_graph();
    g.layout2d = g.tutte_layout();
    Polyhedron P0(g,g.zero_order_geometry(),6);
    P0s.push_back(P0);    
      
    Polyhedron P = P0;
    P.points = g.optimized_geometry(P.points);

    Ps.push_back(P);
    Ts.push_back(T);

    Tlf.push_back(g.leapfrog_dual());
  }

  size_t M = Ps.size();
  size_t Nfaces = N/2+2;
  
  vector<int> neighbours(M*N*3), dual_neighbours(M*Nfaces*6,-1), next_on_face(M*N*3),prev_on_face(M*N*3),
    face_right(M*N*3), next_on_tri(M*Nfaces*6,-1);
  vector<size_t> neighbours_shape{{M,N,3}}, dual_neighbours_shape{{M,Nfaces,6}},
					     points_shape{{M,N,3}}, triangles_shape{{M,N,3}};
						    
  vector<double> points_start(M*N*3), points_opt(M*N*3);
  vector<vector<face_t>>
    pentagons(M,vector<face_t>(12)),
    hexagons (M,vector<face_t>(Nfaces-12));
  vector<vector<tri_t>> triangles(M);
  
  for(int i=0;i<M;i++){
    // Everything with neighbours-shape
    for(node_t u=0;u<N;u++){
    
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
    }
    // Faces
    //    vector<face_t> faces = Ps[i].compute_faces_oriented();
    vector<face_t> faces = Ts[i].cubic_faces();

    triangles[i] = Ts[i].compute_faces();

    for(node_t f=0,npent=0,nhex=0;f<Nfaces;f++){
      auto nf = Ts[i].neighbours[f];

     
      for(int j=0;j<nf.size();j++){
	dual_neighbours[i*Nfaces*6+f*6+j] = nf[j];
	next_on_tri[i*Nfaces*6+f*6+j]     = Ts[i].next_on_face(f,nf[j]);
      }
	
      if      (faces[f].size() == 5) pentagons[i][npent++] = faces[f];
      else if (faces[f].size() == 6) hexagons [i][nhex++]  = faces[f];
    }
  }
  cerr << "from numpy import array, nan\n\n";
  cerr << "rspi_shape = " << vector<size_t>{{M,12}} << ";\n";
  cerr << "neighbours_shape = " << neighbours_shape << ";\n";
  cerr << "points_shape     = " << points_shape << ";\n";  
  cerr << "dual_neighbours_shape = " << dual_neighbours_shape << ";\n";
  cerr << "triangles_shape = " << triangles_shape << ";\n";    

  cerr << "rspi         = array(" << all_rspi << ").reshape(rspi_shape);\n\n";
  cerr << "# Cubic graph, its faces, and 3D embedding\n";
  cerr << "cubic_neighbours  = array(" << neighbours << ").reshape(neighbours_shape);\n\n";
  cerr << "pentagons    = array(" << pentagons  << ");\n\n";
  cerr << "hexagons     = array(" << hexagons   << ");\n\n";  
  cerr << "next_on_face = array(" << next_on_face << ").reshape(neighbours_shape);\n\n";
  cerr << "prev_on_face = array(" << prev_on_face << ").reshape(neighbours_shape);\n\n";
  cerr << "face_right   = array(" << face_right << ").reshape(neighbours_shape);\n\n"; 
  cerr << "points_start = array(" << points_start << ").reshape(points_shape);\n\n";
  cerr << "points_opt   = array(" << points_opt << ").reshape(points_shape);\n\n";   

  cerr << "# Dual graph and its faces\n";
  cerr << "dual_neighbours   = array(" << dual_neighbours << ").reshape(dual_neighbours_shape);\n\n";
  cerr << "next_on_tri       = array(" << next_on_tri << ").reshape(dual_neighbours_shape);\n\n";  
  cerr << "triangles         = array(" << triangles << ").reshape(triangles_shape);\n\n";
  cerr << "# prev_on_tri is the same as next_on_tri\n";

  bool do_unfolding = true;
  if(do_unfolding){
    vector<vector<arc_t>> Arcs(Ts.size());
    vector<vector<Unfolding::arccoord_t>> Arcpos(Ts.size());
    vector<vector<vector<Eisenstein>>> Tripos(Ts.size());    
   
    for(int j=0;j<Ts.size();j++){
      int u_max, mx;
      tie(u_max,mx) = middle_hexagon(Tlf[j]);
      arc_t first_arc = {u_max, Tlf[j].neighbours[u_max][0]};    
      
      Unfolding uf(Tlf[j],first_arc);
      Arcs[j]   = get_keys(uf.arc_coords);
      Arcpos[j] = get_values(uf.arc_coords);
      Tripos[j] = uf.tri_coords();
    }
    cerr << "lf_unfolding_arcs   = " << Arcs << ";\n";
    cerr << "lf_unfolding_arcpos = " << Arcpos << ";\n";
    cerr << "lf_unfolding_tripos = " << Tripos << ";\n";
  }
  
  return 0;
}
