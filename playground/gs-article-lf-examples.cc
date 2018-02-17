#include "libgraph/fullerenegraph.hh"
#include "libgraph/polyhedron.hh"
#include "libgraph/isomerdb.hh"
#include "libgraph/triangulation.hh"
#include "libgraph/fullerenegraph.hh"
#include "libgraph/symmetry.hh"

// all lists are CCW

// tetrehedron, testing
Graph example1()
{
  int N = 4;
  neighbours_t neighbours(N,vector<node_t>(3));

  neighbours[0][0] = 1;
  neighbours[0][1] = 2;
  neighbours[0][2] = 3;

  neighbours[1][0] = 2;
  neighbours[1][1] = 0;
  neighbours[1][2] = 3;

  neighbours[2][0] = 3;
  neighbours[2][1] = 0;
  neighbours[2][2] = 1;

  neighbours[3][0] = 1;
  neighbours[3][1] = 0;
  neighbours[3][2] = 2;

  return Graph(neighbours,true);
}

// truncated tetrahedron
Graph example2(){
  const int M=3, N=12;
  neighbours_t neighbours(N,vector<node_t>(3));

  for(int i=0; i!=M; ++i){
    neighbours[i][0] = M+i;
    neighbours[i][1] = (i+1)%M;
    neighbours[i][2] = (i-1+M)%M;
    
    neighbours[1*M+i][0] = i;
    neighbours[1*M+i][1] = 3*M+i;
    neighbours[1*M+i][2] = 2*M+i;
    
    neighbours[2*M+i][0] = M+i;
    neighbours[2*M+i][1] = 3*M+i;
    neighbours[2*M+i][2] = 3*M+(i+1)%M;
    
    neighbours[3*M+i][0] = M+i;
    neighbours[3*M+i][1] = 2*M+(i-1+M)%M;
    neighbours[3*M+i][2] = 2*M+i;
  }
  return Graph(neighbours,true);
}


// triakis tetrahedron
Graph example3(){
  const int M=3, N=8;
  neighbours_t neighbours(N,vector<node_t>());

  for(int i=0; i!=M; ++i){
    neighbours[i].resize(3);
    neighbours[i][0] = 6;
    neighbours[i][1] = M+(i-1+M)%M;
    neighbours[i][2] = M+i;

    neighbours[M+i].resize(6);
    neighbours[M+i][0] = 6;
    neighbours[M+i][1] = i;
    neighbours[M+i][2] = M+(i-1+M)%M;
    neighbours[M+i][3] = 7;
    neighbours[M+i][4] = M+(i+1)%M;
    neighbours[M+i][5] = (i+1)%M;
  }
  neighbours[6].resize(6);
  neighbours[6][0] = 0;
  neighbours[6][1] = 3;
  neighbours[6][2] = 1;
  neighbours[6][3] = 4;
  neighbours[6][4] = 2;
  neighbours[6][5] = 5;

  neighbours[7].resize(3);
  neighbours[7][0] = 3;
  neighbours[7][1] = 5;
  neighbours[7][2] = 4;
  return Graph(neighbours,true);
}


// cuboctahedron
Graph example4(){
  const int M=4, N=12;
  neighbours_t neighbours(N,vector<node_t>(4));

  for(int i=0; i!=M; ++i){
    neighbours[i][0] = (i+1)%M;
    neighbours[i][1] = (i-1+M)%M;
    neighbours[i][2] = M+(i-1+M)%M;
    neighbours[i][3] = M+i;

    neighbours[M+i][0] = (i+1)%M;
    neighbours[M+i][1] = i;
    neighbours[M+i][2] = 2*M+i;
    neighbours[M+i][3] = 2*M+(i+1)%M;

    neighbours[2*M+i][0] = M+i;
    neighbours[2*M+i][1] = M+(i-1+M)%M;
    neighbours[2*M+i][2] = 2*M+(i-1+M)%M;
    neighbours[2*M+i][3] = 2*M+(i+1)%M;
  }
  return Graph(neighbours,true);
}


// rhombicuboctahedron
Graph example5(){
  const int M=4, N=24;
  neighbours_t neighbours(N,vector<node_t>(4));

  for(int i=0; i!=M; ++i){
    neighbours[i][0] = (i+1)%M;
    neighbours[i][1] = (i-1+M)%M;
    neighbours[i][2] = M+i;
    neighbours[i][3] = 2*M+i;

    neighbours[M+i][0] = i;
    neighbours[M+i][1] = 2*M+(i-1+M)%M;
    neighbours[M+i][2] = 3*M+i;
    neighbours[M+i][3] = 2*M+i;

    neighbours[2*M+i][0] = i;
    neighbours[2*M+i][1] = M+i;
    neighbours[2*M+i][2] = 4*M+i;
    neighbours[2*M+i][3] = M+(i+1)%M;

    neighbours[3*M+i][0] = M+i;
    neighbours[3*M+i][1] = 4*M+(i-1+M)%M;
    neighbours[3*M+i][2] = 5*M+i;
    neighbours[3*M+i][3] = 4*M+i;

    neighbours[4*M+i][0] = 2*M+i;
    neighbours[4*M+i][1] = 3*M+i;
    neighbours[4*M+i][2] = 5*M+i;
    neighbours[4*M+i][3] = 3*M+(i+1)%M;

    neighbours[5*M+i][0] = 4*M+i;
    neighbours[5*M+i][1] = 3*M+i;
    neighbours[5*M+i][2] = 5*M+(i-1+M)%M;
    neighbours[5*M+i][3] = 5*M+(i+1)%M;
  }
  return Graph(neighbours,true);
}


// the macromol thing
Graph example6(){
  const int M=4, N=30;
  neighbours_t neighbours(N,vector<node_t>(4));

  for(int i=0; i!=M; ++i){
    neighbours[i][0] = 28;
    neighbours[i][1] = M+(i-1+M)%M;
    neighbours[i][2] = 2*M+i;
    neighbours[i][3] = M+i;

    neighbours[M+i][0] = i;
    neighbours[M+i][1] = 3*M+i;
    neighbours[M+i][2] = 2*M+(i+1)%M;
    neighbours[M+i][3] = (i+1)%M;

    neighbours[2*M+i][0] = i;
    neighbours[2*M+i][1] = M+(i-1+M)%M;
    neighbours[2*M+i][2] = 4*M+i;
    neighbours[2*M+i][3] = 3*M+i;

    neighbours[3*M+i][0] = M+i;
    neighbours[3*M+i][1] = 2*M+i;
    neighbours[3*M+i][2] = 5*M+i;
    neighbours[3*M+i][3] = 4*M+(i+1)%M;

    neighbours[4*M+i][0] = 2*M+i;
    neighbours[4*M+i][1] = 3*M+(i-1+M)%M;
    neighbours[4*M+i][2] = 6*M+i;
    neighbours[4*M+i][3] = 5*M+i;

    neighbours[5*M+i][0] = 3*M+i;
    neighbours[5*M+i][1] = 4*M+i;
    neighbours[5*M+i][2] = 6*M+i;
    neighbours[5*M+i][3] = 6*M+(i+1)%M;

    neighbours[6*M+i][0] = 4*M+i;
    neighbours[6*M+i][1] = 5*M+(i-1+M)%M;
    neighbours[6*M+i][2] = 29;
    neighbours[6*M+i][3] = 5*M+i;
  }
  neighbours[28][0] = 0;
  neighbours[28][1] = 1;
  neighbours[28][2] = 2;
  neighbours[28][3] = 3;

  neighbours[29][0] = 24;
  neighbours[29][1] = 27;
  neighbours[29][2] = 26;
  neighbours[29][3] = 25;
  return Graph(neighbours,true);
}

// elongated square bipyramid (should we want it)
Graph example7(){
  const int M=4, N=10;
  neighbours_t neighbours(N,vector<node_t>(4));

  for(int i=0; i!=M; ++i){
    neighbours[i][0] = 8;
    neighbours[i][1] = (i-1+M)%M;
    neighbours[i][2] = M+i;
    neighbours[i][3] = (i+1)%M;

    neighbours[M+i][0] = i;
    neighbours[M+i][1] = M+(i-1+M)%M;
    neighbours[M+i][2] = 9;
    neighbours[M+i][3] = M+(i+1)%M;
  }
  neighbours[8][0] = 0;
  neighbours[8][1] = 1;
  neighbours[8][2] = 2;
  neighbours[8][3] = 3;

  neighbours[9][0] = 7;
  neighbours[9][1] = 6;
  neighbours[9][2] = 5;
  neighbours[9][3] = 4;
  return Graph(neighbours,true);
}


// tetrahedron
Polyhedron Example1()
{
  // Normalized to unit edge length 
  vector<coord3d> points{{1.05196, 1.75377, 0.513112}, {0.891458, 0.827152, 0.174592}, {0.434181, 1.14097, 1.00679}, {1.41947, 0.969106, 1.01171}};

  return Polyhedron(example1(),points,3);
}

// truncated tetrahedron
Polyhedron Example2()
{
  vector<coord3d> points{{2.41493, 1.37468, 1.46768}, {2.09275, 2.04649, 0.800694}, {1.87136, 2.15801, 1.76916}, {1.92872, 0.246424, 1.48687}, {1.21092, 1.74364, 0.}, {0.717997, 1.99257, 2.15879}, {1.21139, 0., 0.834834}, {0.272162, 1.6085, 0.317471}, {0.322097, 1.07391, 2.17425}, {0.989806, 0.111258, 1.80416}, {0.815354, 0.82527, 0.0154645}, {0., 1.74615, 1.50818}};

  return Polyhedron(example2(),points,7);
}

// triakis tetrahedron
Polyhedron Example3()
{
  vector<coord3d> points{{0., 1.62589, 1.13702}, {0.787645, 0.155512, 0.}, {1.99993, 1.34964, 1.08679}, {0.3575, 0.50278, 1.04533}, {1.34762, 0.364738, 1.02069}, {0.957819, 1.093, 1.58342}, {0.943175, 1.17002, 0.58747}, {0.81648, 0., 2.01168}};
  
  return Polyhedron(example3(),points,3);
}

Polyhedron Example4()
{
  // Extracted from Mathematica's GraphPlot3D
  vector<coord3d> points{{1.34874, 1.14808, 0.}, {1.63636, 0.296786, 0.439735}, {0.682595, 0., 0.489853}, {0.395479, 0.851524, 0.0504444}, {1.90795, 1.17722, 0.828152}, {1.2425, 0.029786, 1.3178}, {0., 0.583842, 0.929148}, {0.666174, 1.73162, 0.440428}, {1.22508, 1.76097, 1.26813}, {1.51408, 0.909684, 1.70714}, {0.559471, 0.613215, 1.75766}, {0.270617, 1.46379, 1.31917}};

  return Polyhedron(example4(),points,4);
}

Polyhedron Example5()
{
  // Extracted from Mathematica's GraphPlot3D
  vector<coord3d> points{{1.9049, 2.78742, 1.87212}, {2.72988, 2.27101, 1.28218}, {2.07476, 2.30732, 0.384266}, {1.24791, 2.83037, 0.976054}, {1.10496, 2.48671, 2.40496}, {2.88235, 1.6184, 2.01215}, {2.29403, 1.35658, 0.198893}, {0.546742, 2.22817, 0.579918}, {2.05013, 2.13251, 2.59767}, {2.97646, 1.32505, 1.08533}, {1.37302, 1.70401, 0.}, {0.451023, 2.51782, 1.50508}, {0.679585, 1.47177, 2.65847}, {2.42336, 0.609225, 2.28318}, {1.85902, 0.349844, 0.453745}, {0.0928809, 1.21083, 0.845492}, {1.59791, 1.12955, 2.85736}, {2.52278, 0.315801, 1.36331}, {0.919415, 0.703607, 0.254353}, {0., 1.50504, 1.77463}, {0.90576, 0.517145, 2.46462}, {1.72458, 0., 1.88074}, {1.0682, 0.0474456, 0.975372}, {0.240998, 0.553828, 1.56541}};

  return Polyhedron(example5(),points,4);
}

Polyhedron Example6()
{
    // Extracted from Mathematica's GraphPlot3D
  vector<coord3d> points{{1.72556, 2.74414, 2.72532}, {2.96765, 2.54854, 1.84286}, {3.13329, 1.15502, 2.36031}, {1.8973, 1.33239, 3.22814}, {2.14506, 3.10691, 1.79109}, {3.36918, 1.57784, 1.50582}, {2.32856, 0.562358, 2.79959}, {1.09231, 2.06214, 3.06661}, {0.785743, 2.70424, 2.40051}, {2.57314, 2.76575, 0.976252}, {3.11153, 0.643044, 1.50912}, {1.33747, 0.574436, 2.93227}, {1.09141, 3.00454, 1.38614}, {2.85836, 1.77772, 0.538115}, {2.22049, 0.0370753, 1.82513}, {0.458455, 1.22466, 2.65056}, {0.062109, 2.00459, 1.96437}, {1.6117, 2.72298, 0.444181}, {2.64039, 0.692725, 0.52523}, {1.11407, 0., 2.03524}, {0.260453, 2.30931, 1.05905}, {1.86688, 1.88681, 0.}, {1.89507, 0.0992983, 0.763338}, {0.295422, 0.500085, 1.81479}, {0., 1.36813, 1.21442}, {0.902812, 2.07085, 0.208289}, {1.75219, 0.809482, 0.0984979}, {0.857399, 0.131403, 1.10422}, {2.61284, 2.05525, 2.76884}, {0.694215, 0.997084, 0.426291}};

  return Polyhedron(example6(),points,4);
}

// elongated square bipyramid
Polyhedron Example7()
{
  vector<coord3d> points{{1.9697, 0.997756, 1.06898}, {1.96952, 0.000364472, 0.996802}, {1.96881, 0.0712316, 0.000124596}, {1.96853, 1.06939, 0.0708048}, {0.682214, 0.998026, 1.06852}, {0.681719, 0., 0.99687}, {0.681605, 0.0722204, 0.}, {0.68192, 1.06942, 0.0714736}, {2.65146, 0.535316, 0.533368}, {0., 0.535071, 0.534176}};

  return Polyhedron(example7(),points,4);
}

// First fullerene without a pentagon-starting spiral is Td-C100
// with classical canonical spiral [CS: 2,8,9,23,24,28,29,37,41,45,46,52]
// and  canonical                  [GS: 43,2; 1,4,5,26,27,31,32,40,43,47,48,52]
Graph example8_TdC100()
{
  vector<int> RSPI{{2,8,9,23,24,28,29,37,41,45,46,52}};
  FullereneDual dF(100,RSPI+(-1));
  return dF.dual_graph();
}

Polyhedron Example8_TdC100()
{
  FullereneGraph F(example8_TdC100());

  vector<face_t>   faces  = F.compute_faces(6);
  cerr << "facesTdC100 = " << faces << ";\n";
  F.layout2d = F.tutte_layout(face_t{{0,71,70,9,14,63}});
  vector<coord3d>  points = F.optimized_geometry(F.zero_order_geometry());

  return Polyhedron(F,points,6,faces);
}

// Random C1-C100 fullerene 
// with classical canonical pentagon indices [GS:  1,2,3,4,5,12,43,46,49,50,51,52]
Graph example9_C1C100()
{
  vector<int> RSPI{{1,2,3,4,5,12,43,46,49,50,51,52}};
  FullereneDual dF(100,RSPI+(-1));
  return dF.dual_graph();
}

Polyhedron Example9_C1C100()
{
  FullereneGraph F(example9_C1C100());

  vector<face_t>   faces  = F.compute_faces(6);
  cerr << "facesC1C100 = " << faces << ";\n";
  F.layout2d = F.tutte_layout(face_t{{54,98,64,88,90}});
  vector<coord3d>  points = F.optimized_geometry(F.zero_order_geometry());

  return Polyhedron(F,points,6,faces);
}


Graph example10_Tutte()
{ 
  Triangulation tutte_dual({5,10,5,5,5,9,5,4,5,4,4,5,4,10,5,5,5,5,5,10,4,5,5,4,5},jumplist_t{{10,1},{16,1}});
  //cerr << "tuttedual = " << tutte_dual << ";\n";
  return tutte_dual.dual_graph();
}

Polyhedron Example10_Tutte()
{
  PlanarGraph g(example10_Tutte());
  g.layout2d = g.tutte_layout();
  vector<coord3d> points = g.zero_order_geometry();
  Polyhedron P(g,points,10,g.compute_faces());

  return P;
}

//[GS: 4, 6^3,(6,4)^5]-24-cage
Graph example11_OmniOct()
{ 
  Triangulation omnioct_dual({{4,6,6,6,6,4,6,4,6,4,6,4,6,4}},jumplist_t{});
  return omnioct_dual.dual_graph();
}

//[GS: 9,3; (8, 3, 6, 3)^2, 3, 8, 3]-18-cage
Graph example12_OmniPrism()
{ 
  Triangulation omniprism_dual({{8,3,6,3,8,3,6,3,3,8,3}},jumplist_t{{8,3}});
  return omniprism_dual.dual_graph();
}

//Smallest non-spiralable polyhedral graph with faces up to hexagons
//[GS: 18,2; 3, 6^11, 3, 6, (6, 3)^2, 6^2]-36-cage
Graph example13_SmallestNS()
{
  Triangulation dual({{3,6,6,6,6,6,6,6,6,6,6,6,3,6,6,3,6,3,6,6}},jumplist_t{{17,2}});
  return dual.dual_graph();
}



Polyhedron Example11_OmniOct()
{
  PlanarGraph g(example11_OmniOct());
  g.layout2d = g.tutte_layout();
  vector<coord3d> points = g.zero_order_geometry();
  Polyhedron P(g,points,10,g.compute_faces());

  return P;
}

Polyhedron Example12_OmniPrism()
{
  PlanarGraph g(example12_OmniPrism());
  g.layout2d = g.tutte_layout();
  vector<coord3d> points = g.zero_order_geometry();
  Polyhedron P(g,points,10,g.compute_faces());

  return P;
}


Polyhedron Example13_SmallestNS()
{
  PlanarGraph g(example13_SmallestNS());
  g.layout2d = g.tutte_layout();
  vector<coord3d> points = g.zero_order_geometry();
  Polyhedron P(g,points,10,g.compute_faces());

  return P;
}


PlanarGraph ExampleGraph(int Nex)
{
  switch(Nex){
  case 1: return example1();
  case 2: return example2();
  case 3: return example3();
  case 4: return example4();
  case 5: return example5();
  case 6: return example6();
  case 7: return example7();
  case 8: return example8_TdC100();
  case 9: return example9_C1C100();
  case 10: return example10_Tutte();
  case 11: return example11_OmniOct();
  case 12: return example12_OmniPrism();
  case 13: return example13_SmallestNS();
  default:
    break;
  }
  cerr << "invalid graph chosen, aborting ..." << endl;
  abort();    
}

Polyhedron ExamplePolyhedron(int Nex)
{
  switch(Nex){
  case 1: return Example1();
  case 2: return Example2();
  case 3: return Example3();
  case 4: return Example4();
  case 5: return Example5();
  case 6: return Example6();
  case 7: return Example7();
  case 8: return Example8_TdC100();
  case 9: return Example9_C1C100();
  case 10: return Example10_Tutte();
  case 11: return Example11_OmniOct();
  case 12: return Example12_OmniPrism();
  case 13: return Example13_SmallestNS();
  default:
    break;
  }
  cerr << "invalid graph chosen, aborting ..." << endl;
  abort();    
}

Polyhedron LFPolyhedron(const Polyhedron& P)
{
  Triangulation LF = P.leapfrog_dual();
  vector<coord3d> points = P.points;
  points.resize(LF.N);

  for(int i=0;i<P.faces.size();i++){
    points[P.N+i] = P.faces[i].centroid(P.points);
  }

  return Polyhedron(LF,points);
}

string jumps_to_string(const jumplist_t &jumps)
{
  string s="";
  for(const auto &j: jumps)
    s += to_string(j.first) +"," + to_string(j.second) + ",";
  s.pop_back();

  return s;
}

string spiral_to_string(const vector<int>& spiral)
{
  string s="";
  for(int i=0;i<spiral.size();i++)
    s += to_string(spiral[i]) + (i+1<spiral.size()? ",":"");
  return s;
}

vector<int> spiral_to_rspi(vector<int> spiral)
{
  vector<int> rspi(12);
  for(int i=0,j=0;i<spiral.size();i++) if(spiral[i] == 5) rspi[j++] = i;
  return rspi;
}


string spiral_to_rspi_string(const vector<int>& spiral)
{
  return spiral_to_string(spiral_to_rspi(spiral)+1);
}

struct name_info {
  enum { CUBIC, TRIANGULATION, GENERAL } graph_type;
  bool is_a_fullerene;
  string atom;
  //  string point_group;

  
  PlanarGraph   graph;
  Triangulation triangulation;
  bool cs;

  general_spiral GS;
  Permutation permutation;

  name_info(const PlanarGraph &g, const string& atom="",bool compatibility=false) : atom(atom), graph(g),cs(compatibility) {
    if(g.is_triangulation()){
      graph_type     = TRIANGULATION;
      is_a_fullerene = g.dual_graph().is_a_fullerene();
      triangulation  = g;
    } else if(g.is_cubic()){
      graph_type     = CUBIC;
      triangulation  = g.dual_graph();
      is_a_fullerene = g.is_a_fullerene();
    } else {
      graph_type     = GENERAL;
      triangulation  = g.leapfrog_dual();
      is_a_fullerene = false;
    }

    permutation.resize(triangulation.N);
    bool spiral_success = triangulation.get_spiral(GS.spiral,GS.jumps,!cs);
    assert(spiral_success);

    // Find spiral order with respect to triangulation vertices
    permutation.resize(triangulation.N);
    spiral_success = triangulation.get_spiral(GS.spiral,GS.jumps,permutation,!cs);

    cout << "original: " << triangulation << endl;

    // See if we can reconstruct the triangulation:
    Triangulation t(GS.spiral,GS.jumps);
    cout << "reproduction: " << t << endl;
    // Nope! Spiral wind-up breaks in the presence of separating triangles; needs to be made oriented.
  }

  friend ostream& operator<<(ostream &s, const name_info &n){
    string graph_type_string[3] = {"","D,","LF,"};
    s << "["<<graph_type_string[n.graph_type] <<(n.cs?"CS":"GS") << ": "
      << (n.GS.jumps.empty()? "": (jumps_to_string(n.GS.jumps)+"; "))
      << (n.is_a_fullerene? spiral_to_rspi_string(n.GS.spiral) : spiral_to_string(n.GS.spiral))
      << "]-" <<n.atom<< n.graph.N <<"-" << (n.is_a_fullerene? "fullerene" : "cage");
    return s;
  }
};

#include "gs-article-fix-examples.cc"

int main(int ac, char **av)
{
  int Nex = ac>=2? strtol(av[1],0,0) : 1;
  
  string basename("gs-ex-"+to_string(Nex));

  PlanarGraph g  = ExampleGraph(Nex);
  PlanarGraph dg = g.dual_graph();

  ofstream goutput(("output/"+basename+"-g.m").c_str());
  goutput << "g = " << g << ";\n"
	  << "dg = " << dg << ";\n";
  goutput.close();
  
  Polyhedron  P = ExamplePolyhedron(Nex);

  P.points *= 3;		// Triangle bond length is 3
  P.move_to_origin();
  P.align_with_axes();
  
  cerr << "Graph has "<<(g.has_separating_triangles()?"":"no ")<<"separating triangles.\n";
  
  Polyhedron LFP = LFPolyhedron(P);

  ofstream output(("output/"+basename+".m").c_str());
  if(Nex == 1){
    g.layout2d = g.tutte_layout();
  } else {
    vector<int> spiral;
    jumplist_t  jumps;
    Triangulation LFT = LFP;

    assert(LFT.is_consistently_oriented());
    LFT.get_spiral(spiral,jumps);

    output << "LFspiral = " << spiral << ";\n"
	   << "LFjumps  = " << jumps  << ";\n";

    name_info name(g,"C");

    name_info compat_name(g,"C",true);

    // End at outer face, or as close to it as possible
    vector<int> outer_faces[13] = {{0,1,2},
				   {6, 9, 8, 11, 7, 10},
				   {3,4,5},{0,1,2,3},{20,21,22,23},{23, 27, 29, 24},{0,3,7,4},{93,95,97,99,94},{54,98,64,88,90},{},
				   {0, 3, 7, 14, 15, 11},{3, 12, 13, 14, 15, 10, 9, 4},{27, 31, 32, 34, 35, 28}};

    face_t outer_face = outer_faces[Nex-1];
    if(outer_face.empty())
      g.layout2d = g.tutte_layout();
    else 
      g.layout2d = g.tutte_layout(outer_faces[Nex-1]);
    
    // Hack:
    if(Nex==8) fix_up_TdC100_layout(g);
    //    if(Nex==9) fix_up_C1C100_layout(g);
    
    vector<face_t> faces = g.compute_faces(12);
    vector<coord2d> dglayout(faces.size()), LFlayout2d(g.layout2d);

    for(int i=0;i<faces.size();i++){
      coord2d x = faces[i].centroid(g.layout2d);
      LFlayout2d.push_back(x);
      dglayout[i] = x;
    }
    LFP.layout2d = LFlayout2d;
    dg.layout2d  = dglayout;

    output << "triangulation = "  << name.triangulation << ";\n"
	   << "permutation   = " << name.permutation << "+1;\n"
	   << "jumps         = " << name.GS.jumps << ";\n"
	   << "spiral        = " << name.GS.spiral << ";\n"
      	   << "cspermutation   = " << compat_name.permutation << "+1;\n"
      	   << "csjumps         = " << compat_name.GS.jumps << ";\n"
    	   << "csspiral        = " << compat_name.GS.spiral << ";\n";

    if(name.is_a_fullerene) output << "rspi = " << spiral_to_rspi(name.GS.spiral) << ";\n";
    output << "gsname = \"" << name << "\";\n";
    output << "csname = \"" << compat_name << "\";\n";

    Symmetry S(name.triangulation);
    cerr << "groupsize = " << S.G.size() << ";\n";
    
  }
  
  output << "g   = " << g << ";\n";
  output << "dg  = " << dg << ";\n";
  output << "LFg = " << PlanarGraph(LFP) << ";\n";
  output << "faces   = " << g.compute_faces() << "+1;\n"
	 << "LFfaces = " << LFP.faces << "+1;\n";
  
  {
    ofstream mol2(("output/"+basename+"-P.mol2").c_str());
    mol2 << P.to_mol2();
    mol2.close();

    ofstream pov(("output/"+basename+"-P.pov").c_str());
    pov << P.to_povray();
    pov.close();
  }

  {
    ofstream mol2(("output/"+basename+"-LFP0.mol2").c_str());
    mol2 << LFP.to_mol2();
    mol2.close();
  }

  // {
  //   ofstream mol2(("output/"+basename+"-LFP.mol2").c_str());
  //   LFP.optimize();		// This doesn't really work well for some reason..?
  //   mol2 << LFP.to_mol2();
  //   mol2.close();
  // }  
  
  {
    ofstream tex(("output/"+basename+"-layout.tex").c_str());
    tex << g.to_latex(10, 10, false, false, true, 0, 0, 0, 0.5, 0.5, 2);
    tex.close();
  }  

  output.close();

  
  return 0;
}
