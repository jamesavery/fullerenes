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
  vector<coord3d> points{{1.34874, 1.14808, 0.}, {1.63636, 0.296786, 0.439735}, {0.682595, 0., 0.489853}, {0.395479, 0.851524, 0.0504444}, {1.90795, 1.17722, 0.828152}, {1.2425, 0.029786, 1.3178}, {0., 0.583842, 0.929148}, {0.666174, 1.73162, 0.440428}, {1.22508, 1.76097, 1.26813}, {1.51408, 0.909684, 1.70714}, {0.559471, 0.613215, 1.75766}, {0.270617, 1.46379, 1.31917}};

  return Polyhedron(example4(),points,4);
}

Polyhedron Example5()
{
  // Extracted from Mathematica's GraphPlot3D
  vector<coord3d> points{{1.78772, 2.61595, 1.75696}, {2.56195, 2.13131, 1.2033}, {1.94714, 2.16539, 0.360628}, {1.17115, 2.65626, 0.916014}, {1.03699, 2.33375, 2.25703}, {2.70505, 1.51885, 1.88838}, {2.15292, 1.27314, 0.186658}, {0.51311, 2.09111, 0.544246}, {1.92402, 2.00134, 2.43788}, {2.79337, 1.24354, 1.01856}, {1.28856, 1.59919, 0.}, {0.423279, 2.36294, 1.41249}, {0.637782, 1.38124, 2.49494}, {2.27429, 0.57175, 2.14273}, {1.74467, 0.328324, 0.425834}, {0.0871675, 1.13635, 0.793483}, {1.49962, 1.06007, 2.6816}, {2.36759, 0.296375, 1.27945}, {0.862859, 0.660326, 0.238707}, {0., 1.41246, 1.66546}, {0.850044, 0.485334, 2.31301}, {1.6185, 0., 1.76505}, {1.00249, 0.044527, 0.915374}, {0.226174, 0.519761, 1.46912}};

  return Polyhedron(example5(),points,4);
}

Polyhedron Example6()
{
    // Extracted from Mathematica's GraphPlot3D
  vector<coord3d> points{{1.62467, 2.5837, 2.56598}, {2.79414, 2.39953, 1.73511}, {2.95009, 1.08749, 2.2223}, {1.78637, 1.25448, 3.03939}, {2.01964, 2.92525, 1.68637}, {3.17219, 1.48559, 1.41778}, {2.19241, 0.529478, 2.6359}, {1.02844, 1.94157, 2.88731}, {0.739802, 2.54613, 2.26015}, {2.42269, 2.60404, 0.919172}, {2.9296, 0.605446, 1.42089}, {1.25927, 0.540849, 2.76083}, {1.02759, 2.82887, 1.30509}, {2.69123, 1.67378, 0.506652}, {2.09066, 0.0349076, 1.71842}, {0.43165, 1.15306, 2.49559}, {0.0584776, 1.88739, 1.84951}, {1.51746, 2.56377, 0.41821}, {2.48601, 0.652223, 0.49452}, {1.04894, 0., 1.91624}, {0.245225, 2.17429, 0.997127}, {1.75773, 1.77649, 0.}, {1.78427, 0.0934925, 0.718707}, {0.278149, 0.470846, 1.70868}, {0., 1.28813, 1.14341}, {0.850026, 1.94977, 0.196111}, {1.64975, 0.762152, 0.0927389}, {0.807268, 0.12372, 1.03965}, {2.46007, 1.93508, 2.60695}, {0.653625, 0.938786, 0.401366}};

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

  PlanarGraph g = ExampleGraph(Nex);
  Polyhedron  P = ExamplePolyhedron(Nex);

  ofstream output(("output/"+basename+".m").c_str());
  output << "g = " << g << ";\n";
  output.close();
  
  {
    ofstream mol2(("output/"+basename+"-P.mol2").c_str());
    mol2 << P.to_mol2();
    mol2.close();
  }
  
  
  return 0;
}
