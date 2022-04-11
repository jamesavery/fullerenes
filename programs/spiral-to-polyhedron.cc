#include <string>

#include "fullerenes/triangulation.hh"
#include "fullerenes/spiral.hh"
#include "fullerenes/polyhedron.hh"

using namespace std;

int main(int ac, char **av)
{
  string spiral_name = av[1];
  spiral_nomenclature fsn(spiral_name);
  PlanarGraph g(fsn);

  g.layout2d = g.tutte_layout();
  Polyhedron P(g,g.zero_order_geometry());
  P.optimize();
  Polyhedron::to_file(P,av[2]);  


    
  return 0;
}
