#include "libgraph/fullerenegraph.hh"
#include "libgraph/polyhedron.hh"
#include "libgraph/isomerdb.hh"

int testN = 80;
int testRSPI[12] = {1,2,3,4,5,6,37,38,39,40,41,42};

// Smallest planar cubic graph with no face spiral
Graph example1()
{
  int M = 6, N = 18;
  neighbours_t neighbours(N,vector<node_t>(3));

  for(int i=0;i<3;i++){
    // Outer triangle
    neighbours[i*M+0][0] = i*M+1;
    neighbours[i*M+0][1] = i*M+2;
    neighbours[i*M+0][2] = ((i-1)*M+1+N)%N;
    
    neighbours[i*M+1][0] = i*M+0;
    neighbours[i*M+1][1] = i*M+2;
    neighbours[i*M+1][2] = ((i+1)*M+0)%N;

    neighbours[i*M+2][0] = i*M+0;
    neighbours[i*M+2][1] = i*M+1;
    neighbours[i*M+2][2] = i*M+3;

    // Inner triangle
    neighbours[i*M+3][0] = i*M+4;
    neighbours[i*M+3][1] = i*M+5;
    neighbours[i*M+3][2] = i*M+2;
    
    neighbours[i*M+4][0] = i*M+3;
    neighbours[i*M+4][1] = i*M+5;
    neighbours[i*M+4][2] = ((i+1)*M+5)%N;

    neighbours[i*M+5][0] = i*M+3;
    neighbours[i*M+5][1] = i*M+4;
    neighbours[i*M+5][2] = ((i-1)*M+4+N)%N;
  }
  return Graph(neighbours);
}


// Tutte graph
Graph example2(){
  const int m=3, n = 46;
  neighbours_t neighbours(n,vector<node_t>(3));

  for(int i=0; i!=3; ++i){
    neighbours[i][0] = 45;
    neighbours[i][1] = m+i;
    neighbours[i][2] = 2*m+i;
    
    neighbours[3+i][0] = 0+i;
    neighbours[3+i][1] = 3*m+i;

    neighbours[3+i][2] = 4*m+i;
    
    neighbours[6+i][0] = 0+i;
    neighbours[6+i][1] = 9+i;
    neighbours[6+i][2] = 18+i;
    
    neighbours[9+i][0] = 3+i;
    neighbours[9+i][1] = 6+i;
    neighbours[9+i][2] = 15+i;
    
    neighbours[12+i][0] = 3+i;
    neighbours[12+i][1] = 21+i;
    neighbours[12+i][2] = 24+i;
    
    neighbours[15+i][0] = 9+i;
    neighbours[15+i][1] = 24+i;
    neighbours[15+i][2] = 27+i;
    
    neighbours[6*m+i][0] = 2*m+i;
    neighbours[6*m+i][1] = 9*m+i;
    neighbours[6*m+i][2] = 10*m+i;
    
    neighbours[7*m+i][0] = 4*m+i;
    neighbours[7*m+i][1] = 11*m+i;
    neighbours[7*m+i][2] = 12*m+i;
    
    neighbours[8*m+i][0] = 4*m+i;
    neighbours[8*m+i][1] = 11*m+i;
    neighbours[8*m+i][2] = 12*m+i;
    
    neighbours[9*m+i][0] = 5*m+i;
    neighbours[9*m+i][1] = 6*m+i;
    neighbours[9*m+i][2] = 13*m+i;
    
    neighbours[10*m+i][0] = 6*m+i;
    neighbours[10*m+i][1] = 11*m+((i+1+3)%3);
    neighbours[10*m+i][2] = 14*m+i;
    
    neighbours[11*m+i][0] = 7*m+i;
    neighbours[11*m+i][1] = 10*m+((i-1+3)%3);
    neighbours[11*m+i][2] = 14*m+i;
    
    neighbours[12*m+i][0] = 7*m+i;
    neighbours[12*m+i][1] = 8*m+i;
    neighbours[12*m+i][2] = 13*m+i;
    
    neighbours[13*m+i][0] = 9*m+i;
    neighbours[13*m+i][1] = 12*m+i;
    neighbours[13*m+i][2] = 14*m+i;
    
    neighbours[14*m+i][0] = 10*m+i;
    neighbours[14*m+i][1] = 11*m+i;
    neighbours[14*m+i][2] = 13*m+i;
  }

  neighbours[45][0] = 0;
  neighbours[45][1] = 1;
  neighbours[45][2] = 2;

  return Graph(neighbours);
}

Graph examples[2] = {example1(), example2()};

int main(int ac, char **av)
{
  int Nex = ac>=2? strtol(av[1],0,0) : 1;
  string basename("gs-ex-"+to_string(Nex));

  ofstream output(("output/"+basename+".m").c_str());

  PlanarGraph g(examples[Nex-1]);
  output << "g = " << g << ";\n";
  output.flush();

  printf("Computing planar layout\n");
  g.layout2d = g.tutte_layout();

  printf("Computing dual\n");
  PlanarGraph dg(g.dual_graph(8,true));

  printf("Computing planar layour of dual\n");
  dg.layout2d = dg.tutte_layout();

  output << "g = " << g << ";\n";
  output << "dg = " << dg << ";\n";
  
  output.flush();

  Polyhedron P0(g,g.zero_order_geometry(),8);
  Polyhedron dP0(dg,dg.zero_order_geometry(),8);

  output << "P0coordinates = " << P0.points << ";\n";

  {
    ofstream mol2(("output/"+basename+"-P0.mol2").c_str());
    mol2 << P0.to_mol2();
    mol2.close();
  }

  {
    ofstream mol2(("output/"+basename+"-dP0.mol2").c_str());
    mol2 << dP0.to_mol2();
    mol2.close();
  }


  
  Polyhedron P(P0);

  printf("Optimizing polyhedron.\n");
  P.optimize();

  output << "Pcoordinates = "  << P.points << ";\n";
  output << "P = " << P << ";\n";
  
  Polyhedron D(P.dual(8,true));

  output << "PD = " << D << ";\n";
  
  output.close();


  {
    ofstream mol2(("output/"+basename+".mol2").c_str());
    mol2 << P.to_mol2();
    mol2.close();
  }

  {
    ofstream mol2(("output/"+basename+"-dual.mol2").c_str());
    mol2 << D.to_mol2();
    mol2.close();
  }

  return 0;
}
