#include <string>

#include "fullerenes/triangulation.hh"
#include "fullerenes/spiral.hh"
#include "fullerenes/polyhedron.hh"

using namespace std;

int main(int ac, char **av)
{
  if(ac<=2){
    fprintf(stderr,"Syntax: %s <spiral> <filename.mol2>\n",av[0]);
    return -1;
  }
  string spiral_name = av[1];
  spiral_nomenclature fsn(spiral_name);
  FullereneGraph g(fsn);
  g.layout2d = g.tutte_layout();
  auto X = g.zero_order_geometry();
  Polyhedron P(g,X);
  P.optimize();
  Polyhedron::to_file(P,av[2]);  

  //  cout <<  P.to_povray() << "\n";

    
  return 0;
}
