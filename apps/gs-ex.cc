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


int main(int ac, char **av)
{
  string basename("gs-ex-"+to_string(1));

  PlanarGraph g(example1());
  g.layout2d = g.tutte_layout();
  PlanarGraph dg(g.dual_graph(8,true));
  dg.layout2d = dg.tutte_layout();


  Polyhedron P0(g,g.zero_order_geometry(),8);
  Polyhedron dP0(dg,dg.zero_order_geometry(),8);

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


  ofstream output(("output/"+basename+".m").c_str());

  output << "g = " << g << ";\n";
  output << "dg = " << dg << ";\n";
  output << "coordinates0 = " << P0.points << ";\n";

  
  output.close(); return 0;

  Polyhedron P(P0);
  P.optimize();

  output << "coordinates = "  << P.points << ";\n";
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
