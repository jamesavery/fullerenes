#include "libgraph/fullerenegraph.hh"
#include "libgraph/polyhedron.hh"
#include "libgraph/isomerdb.hh"
#include "libgraph/triangulation.hh"

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
    neighbours[i][0] = M+1;
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

// tetrahedron
Polyhedron Example1()
{
  // Works with our optimizer
  PlanarGraph g = example1();
  g.layout2d = g.tutte_layout();
  Polyhedron P(g,g.zero_order_geometry(),6);
  P.optimize();

  return P;
}

// truncated tetrahedron
Polyhedron Example2()
{
  // TODO. Lukas?
  PlanarGraph g = example2();
  g.layout2d = g.tutte_layout();
  Polyhedron P(g,g.zero_order_geometry(),6);
  //  P.optimize(); // Doesn't work

  return P;
}

// triakis tetrahedron
Polyhedron Example3()
{
  PlanarGraph g = example2();
  g.layout2d = g.tutte_layout();
  Polyhedron P(g,g.zero_order_geometry(),6);
  //  P.optimize(); // Doesn't work
  
  return P;
}

Polyhedron Example4()
{
  // Extracted from Mathematica's GraphPlot3D
  vector<coord3d> points{{1.15848, 1.73211, 1.42869}, {1.35437, 1.62013, 0.45456}, {0.207407, 1.46681, 1.26819}, {1.90211, 1.1321, 1.13397}, {0.754571, 0.97847, 1.94709}, {0.402976, 1.35465, 0.294257}, {1.14636, 0.753787, 0.}, {0., 0.601377, 0.812632}, {1.49801, 0.378059, 1.65271}, {1.69473, 0.265372, 0.679336}, {0.743607, 0., 0.518764}, {0.547534, 0.112681, 1.4927}};

  return Polyhedron(example4(),points,4);
}

Polyhedron Example5()
{
  // Extracted from Mathematica's GraphPlot3D
  vector<coord3d> points{{1.9765, 1.83722, 2.32688}, {2.61587, 1.81993, 1.50014}, {2.01955, 0.763659, 2.37434}, {1.0369, 1.8017, 2.47131}, {1.39925, 2.45161, 1.91728}, {2.65633, 0.751793, 1.5462}, {2.04614, 2.44157, 1.09983}, {2.46377, 1.76721, 0.593609}, {2.50182, 0.706189, 0.633354}, {2.14321, 0.0537827, 1.20207}, {1.49444, 0.0742579, 2.02037}, {1.07943, 0.739422, 2.51824}, {0.151892, 1.73832,  1.88159}, {1.15849, 2.37818, 0.498908}, {1.6197, 0.651043,  0.0387256}, {0.609667, 0.00590933, 1.4228}, {0.510259, 2.39391,  1.31902}, {1.5823, 1.71314, 0.}, {1.25732, 0., 0.599416}, {0.19506,  0.679217, 1.92928}, {0., 1.69371, 0.970814}, {0.642485, 1.67977,  0.143231}, {0.682721, 0.613335, 0.186023}, {0.0418463, 0.628505,  1.02086}};

  return Polyhedron(example5(),points,4);
}

Polyhedron Example6()
{
    // Extracted from Mathematica's GraphPlot3D
  vector<coord3d> points{{1.63916, 2.60875, 0.422664}, {1.04952, 1.85887, 0.0472124}, {2.35042, 2.7244, 1.00191}, {2.42996, 2.1053, 0.322226}, {0.891956, 2.76681, 1.16569}, {0.358389, 2.0311, 0.647543}, {0., 1.70183, 1.54335}, {0.485594, 1.18104, 0.320745}, {0.433313, 2.37944, 2.02187}, {1.32144, 2.55026, 2.49941}, {0.29854, 1.55519, 2.40983}, {1.71702, 2.9711, 1.76918}, {2.24393, 2.5084, 2.38789}, {1.91541, 1.24369, 0.}, {0.224142, 0.735448, 1.22728}, {1.23237, 1.58308, 2.90446}, {2.91654, 2.14707, 1.66936}, {3.06146, 1.48644, 0.86622}, {1.24728, 0.471266, 0.248902}, {0.497575, 0.565468, 2.22213}, {2.29848, 1.5736, 2.83708}, {2.65184, 0.711364, 0.540582}, {0.982172, 0.0272659, 1.01124}, {1.28117, 0.556162, 2.7255}, {2.93818, 1.2669, 2.25115}, {3.01272, 0.721563, 1.39468}, {1.88979, 0.044802, 0.771203}, {1.19857, 0., 2.00683}, {2.312, 0.671292, 2.61734}, {2.23247, 0.0965759, 1.75632}};

  return Polyhedron(example6(),points,4);
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
  default:
    break;
  }
  cerr << "invalid graph chosen, aborting ..." << endl;
  abort();    
}

Graph      examples[6] = {example1(), example2(), example3(), example4(), example5(), example6()};

int main(int ac, char **av)
{
  int Nex = ac>=2? strtol(av[1],0,0) : 1;


  string basename("gs-ex-"+to_string(Nex));

  ofstream output(("output/"+basename+".m").c_str());

  PlanarGraph g = ExampleGraph(Nex);
  Polyhedron  P = ExamplePolyhedron(Nex);

  cout << "g = " << g << ";\n";
  
  {
    ofstream mol2(("output/"+basename+"-P.mol2").c_str());
    mol2 << P.to_mol2();
    mol2.close();
  }
  
  
  return 0;
}
