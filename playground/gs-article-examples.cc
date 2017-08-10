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

  return Graph(neighbours);
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
  return Graph(neighbours);
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
  return Graph(neighbours);
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
  return Graph(neighbours);
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
  return Graph(neighbours);
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
  return Graph(neighbours);
}


Graph examples[6] = {example1(), example2(), example3(), example4(), example5(), example6()};

int main(int ac, char **av)
{
  int Nex = ac>=2? strtol(av[1],0,0) : 1;

  if(Nex < 1 || Nex > 6){cerr << "invalid graph chosen, aborting ..." << endl; abort();}

  string basename("gs-ex-"+to_string(Nex));

  ofstream output(("output/"+basename+".m").c_str());

  PlanarGraph g(examples[Nex-1]);
  output << "g = " << g << ";\n";
  output.flush();

  return 0;
}
