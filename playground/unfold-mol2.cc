#include "libgraph/fullerenegraph.hh"
#include "libgraph/triangulation.hh"
#include "libgraph/unfold.hh"
#include "libgraph/polyhedron.hh"
#include "libgraph/auxiliary.hh"
#include <fstream>

int main(int ac, char **av)
{
  Polyhedron  P(av[1]);
  Polyhedron dP(P.dual());
  Triangulation dual(dP);
  
  dual = dual.sort_nodes();
  
  dual.layout2d = dual.tutte_layout();
  Unfolding uf(dual,true);

  ofstream output("output/C"+to_string(P.N)+"-unfold.m");
  output
    << "dual      = " << dual << ";\n"
    << "outline   = " << uf.outline << ";\n"
    << "triangles = " << dual.triangles << ";\n"
    << "P         = " << P << ";\n"
    << "PD        = " << dP << ";\n"
    << "arcs    = " << get_keys(uf.edgecoords) << ";\n"
    << "arcpos  = " << get_values(uf.edgecoords) << ";\n";
  output.close();

  ofstream latex_output("output/C"+to_string(P.N)+"-unfold.tex");

  latex_output << "\\newcommand{\\figunfolding}{"<<uf.to_latex(1,0,2) << "}\n\n";

  Unfolding UF = uf.straighten_lines();

  // output.open("output/C"+to_string(N)+"-star.m");

  // output << "dual    = " << dual << ";\n"
  // 	 << "outline = " << UF.outline << ";\n"
  // 	 << "arcs   = " << get_keys(UF.edgecoords) << ";\n"
  // 	 << "arcpos = " << get_values(UF.edgecoords) << ";\n";

  // output.close();

  // latex_output << "\\newcommand{\\figstar}{"<<UF.to_latex(1,0,1) << "}\n\n";
  
  latex_output.close();

  return 0;
}
