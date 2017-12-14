#include <string>

#include "libgraph/triangulation.hh"
#include "libgraph/spiral.hh"
#include "libgraph/polyhedron.hh"

using namespace std;

int main(int ac, char **av)
{
  string spiral_name = av[1];
  spiral_nomenclature fsn(spiral_name);
  PlanarGraph g(fsn);

  g.layout2d = g.tutte_layout();
  Polyhedron P(g,g.zero_order_geometry());
  P.optimize();
  
  cout << "spiralcode = \"" << spiral_name << "\";\n"
       << "fullname   = " << fsn << ";\n"
       << "g = " << g << ";\n"
       << "P = " << P << ";\n";
    
  return 0;
}
