#include <string>

#include "libgraph/triangulation.hh"
#include "libgraph/spiral.hh"
#include "libgraph/polyhedron.hh"

using namespace std;

int main(int ac, char **av)
{
  if(ac<=1 || string(av[1]) == "help"){
    fprintf(stderr,"Syntax: %s \"<spiral code>\" [output-file]\n"
	    "Graph-file formats: %s.\n",
	    av[0], to_string(PlanarGraph::output_formats).c_str());
    return -1;
  }

  string spiral_name = av[1];
  string filename    = ac>2? av[2] : "-";

  spiral_nomenclature sn(spiral_name);
  PlanarGraph g(sn);

  g.layout2d = g.tutte_layout();
  Polyhedron P(g,g.zero_order_geometry());
  P.optimize();

  // cout << "spiralcode = \"" << spiral_name << "\";\n"
  //      << "fullname   = " << sn << ";\n"
  //      << "g = " << g << ";\n"
  //      << "P = " << P << ";\n";

  FILE *file = filename=="-"? stdout : fopen(filename.c_str(),"wb");
  Polyhedron::to_mol2(P, file);
  fclose(file);

  return 0;
}
