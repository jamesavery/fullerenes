
#include "libgraph/fullerenegraph.hh"
#include "libgraph/triangulation.hh"
#include "libgraph/polyhedron.hh"

#include <iostream>
#include <chrono>		// For wall-time
#include <ctime> 		// For CPU-time

typedef Triangulation::jumplist_t jumplist_t;
typedef std::chrono::high_resolution_clock Clock;

int main(int ac, char **av) {
  int N;
  Triangulation::jumplist_t input_jumps;
  vector<int> RSPI(12);
  if(ac>=14) {
    N = strtol(av[1],0,0);
    for(int i=0;i<12;i++) RSPI[i] = strtol(av[i+2],0,0)-1;
  } else{assert(false);} 
  
  if(ac>14){			// General RSPI: RSPI followed by jumps.
    for(int i=14;i<ac;i+=2)
      input_jumps.push_back(make_pair(strtol(av[i],0,0)-1,strtol(av[i+1],0,0)));
  }

  cout << "Nv       = " << N  << ";\n"
       << "jumplist = " << input_jumps << ";\n"
       << "rspi     = " << RSPI << ";\n";


  FullereneDual  g(N, RSPI, input_jumps);
  FullereneGraph F(g.dual_graph());
  F.is_oriented = true;

  vector<face_t> faces = g.dual_faces();
  F.layout2d = F.tutte_layout();
  

  vector<coord3d>  points = F.optimized_geometry(F.zero_order_geometry());
  Polyhedron       P(F,points,6,faces);
  P.move_to_origin();
  P.align_with_axes();
  

  
  vector<coord3d> dpoints(g.N);
  for(int f=0;f<faces.size();f++){
    coord3d sum{0,0,0};
    for(const auto &v: faces[f]) sum += P.points[v] /double(faces[f].size());
    dpoints[f] = sum;
  }

  vector< vector<coord3d> > face_points(faces.size());
  vector< vector<pair<int,int>> >  shared_faces(F.N);
  
  for(node_t f=0;f<faces.size();f++){
    // Calculate points as difference to centroid
    face_points[f].resize(faces[f].size());
    for(int i=0;i<faces[f].size();i++) face_points[f][i] = points[faces[f][i]] - dpoints[f];

    // Identify shared edge with next neighbour
    if(f+1<faces.size()){
      vector<pair<int,int>> shared_nodes;
      
      face_t f0(faces[f]), f1(faces[f+1]);

      for(int i=0;i<f0.size();i++)
	for(int j=0;j<f1.size();j++)
	  if(f0[i] == f1[j]){
	    shared_faces[f0[i]].push_back(make_pair(f,i));
	    shared_faces[f0[i]].push_back(make_pair(f+1,j));
	  }
    }
  }

  // Do the slightly unwound spiral coordinates
  pair<coord3d,coord3d> bbox = P.bounding_box();
  double xmin = bbox.first[0], xmax = bbox.second[0];
  
  vector<coord3d>          spiral_centers(g.N);
  vector<vector<coord3d> > spiral_vertices(g.N);
  // Perturb the face centers along the spiral (t\in [0,1[ )
  for(node_t f=0;f<faces.size();f++){
    double t     = f/double(faces.size());
    double bulge = (t-.5)*(t-.5);
    
    spiral_centers[f] = dpoints[f]- coord3d{0,0,t*t*t*t*t*t*(xmax-xmin)/2.0};
  }

  // Now do the vertices of the faces.
  // First, displace according to spiral
  for(node_t f=0;f<faces.size();f++)
    for(int i=0;i<faces[f].size();i++)
      spiral_vertices[f].push_back( spiral_centers[f] + face_points[f][i] );

  // Now make sure that vertices on shared edges along the spiral have the same position
  for(node_t v=0;v<F.N;v++){
    coord3d x{0,0,0};

    for(const auto &fi: shared_faces[v]){
      int f = fi.first, i = fi.second;
      x += spiral_vertices[f][i];
    }
    x *= 1.0/shared_faces[v].size();

    for(const auto &fi: shared_faces[v]){
      int f = fi.first, i = fi.second;
      spiral_vertices[f][i] = x;
    }
  }
  
  
  Polyhedron     dP(g,dpoints,3,vector<face_t>(g.triangles.begin(),g.triangles.end()));

  cout << "P              = " << P  << ";\n"
       << "dP             = " << dP << ";\n"
       << "dpoints        = " << dpoints << ";\n"
       << "faces          = " << faces << "+1;\n"
       << "facePoints     = " << face_points << ";\n"
       << "spiralCenters  = " << spiral_centers << ";\n"
       << "spiralVertices = " << spiral_vertices << ";\n";

 

  return 0;
}
