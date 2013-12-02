#include "libgraph/fullerenegraph.hh"
#include "libgraph/polyhedron.hh"
#include "libgraph/isomerdb.hh"

int testN = 380;
int testjump[2] = {110,2};
int testRSPI[12] = {45, 70, 71, 82, 83, 110, 119, 120, 144, 184, 185, 192};

int main(int ac, char **av)
{
  int N;
  vector<int> RSPI(12);
  FullereneGraph::jumplist_t  jumps;

  if(ac>=14) {
    N = strtol(av[1],0,0);
    for(int i=0;i<12;i++) RSPI[i] = strtol(av[i+2],0,0)-1;
  } else {
    N = testN;
    for(int i=0;i<12;i++) RSPI[i] = testRSPI[i]-1;
  }

  if(ac>=16) {
    for(int i=0;i<ac-14;i++)
      jumps.push_back(make_pair(strtol(av[i+15],0,0),strtol(av[i+16],0,0)));
  } else 
    jumps.push_back(make_pair(testjump[0]-1,testjump[1]));

  string basename("polyhedron-"+to_string(N));

  FullereneGraph g(N,RSPI,jumps);
  g.layout2d = g.tutte_layout();

  g = g.halma_fullerene(2 ,true);
  g.layout2d = g.tutte_layout();

  Polyhedron P0(g,g.zero_order_geometry(),6);

  {
    ofstream mol2(("output/"+basename+"-P0.mol2").c_str());
    mol2 << P0.to_mol2();
    mol2.close();
  }
  
  Polyhedron P(P0);
  P.optimize();

  ofstream output(("output/"+basename+".m").c_str());

  output << "g = " << g << ";\n";
  output << "coordinates0 = " << P0.points << ";\n";
  output << "coordinates = "  << P.points << ";\n";

  output << "P = " << P << ";\n";
  
  Polyhedron D(P.dual(6,true));

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
