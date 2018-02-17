
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

  cerr << "General RSPI: {" << input_jumps << ", " << (RSPI+1) << "};\n";
 

  FullereneDual  g(N, RSPI, input_jumps);
  g.is_oriented = true;		// Should be there already
  FullereneGraph F(g.dual_graph());
  F.is_oriented = true;         // Should be there already

  vector<face_t> faces = g.dual_faces();
  F.layout2d = F.tutte_layout();
  

  vector<coord3d>  points = F.optimized_geometry(F.zero_order_geometry());
  points *= 2;
  
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
    
    spiral_centers[f] = dpoints[f]*(1+t/50)+coord3d{0,0,t*t*
					   t*(xmax-xmin)/5.0};
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

  // cout << "P              = " << P  << ";\n"
  //      << "dP             = " << dP << ";\n"
  //      << "dpoints        = " << dpoints << ";\n"
  //      << "faces          = " << faces << "+1;\n"
  //      << "facePoints     = " << face_points << ";\n"
  //      << "spiralCenters  = " << spiral_centers << ";\n"
  //      << "spiralVertices = " << spiral_vertices << ";\n";

  size_t Npeel = 0;
  for(auto f: faces) Npeel += f.size();

  Graph Gpeel(Npeel);
  Gpeel.is_oriented = true;
  vector<face_t> peel_faces(faces.size());
  vector<coord3d> peel_points(Npeel);

  for(int i=0, v=0;i<faces.size();i++){
    const face_t &f(faces[i]);

    //    fprintf(stderr,"%d %ld\n",v,f.size());
    for(int j=0;j<f.size();j++){
      dedge_t e{v+j,v+((j+1)%f.size())};
      Gpeel.insert_edge(e);
      //      cerr << "insert " << e << endl;
      peel_faces[i].push_back(v+j);
      peel_points[v+j] = spiral_vertices[i][j];
    }
    v += f.size();
  }
  Polyhedron Ppeel(Gpeel, peel_points, 6, peel_faces);

  cout << "Gpeel = " << Gpeel << ";\n";

  string basename("polyhedron-"+to_string(N));
  {
    Graph p(Npeel);
    p.is_oriented = true;
    ofstream mol2(("output/"+basename+"-points.mol2").c_str());
    mol2 << Polyhedron(p,peel_points,6, peel_faces).to_mol2();
    mol2.close();
  }

  {
    ofstream mol2(("output/"+basename+"-peel.mol2").c_str());
    ofstream pov(("output/"+basename+"-peel.pov").c_str());
    mol2 << Ppeel.to_mol2();
    mol2.close();
    pov << Ppeel.to_povray() << endl;
    pov.close();
  }

  {
    ofstream mol2(("output/"+basename+"-dP.mol2").c_str());
    mol2 << dP.to_mol2();
    mol2.close();
  }

  {
    ofstream mol2(("output/"+basename+"-P.mol2").c_str());
    ofstream pov(("output/"+basename+"-P.pov").c_str());
    mol2 << P.to_mol2();
    mol2.close();

    pov << P.to_povray() << endl;
    pov.close();
    
  }
  
  return 0;
}
