#include "libgraph/fullerenegraph.hh"
#include "libgraph/triangulation.hh"
#include "libgraph/unfold.hh"
#include "libgraph/auxiliary.hh"
#include <fstream>

int main(int ac, char **av)
{
  assert(ac >= 13);

  int N = strtol(av[1],0,0);
  vector<int> RSPI(12);
  for(int i=0;i<12;i++) RSPI[i] = strtol(av[i+2],0,0)-1;

  vector<int> spiral(N/2+2, 6);
  for (int i = 0; i < 12; i++){
    spiral[RSPI[i]] = 5;
  }

  Triangulation dual(spiral);
  dual = dual.sort_nodes();
  
  dual.layout2d = dual.tutte_layout();
  Unfolding uf(dual,true);

  ofstream output("output/C"+to_string(N)+"-unfold.m");
  output << "dual    = " << dual << ";\n"
	 << "outline = " << uf.outline << ";\n"
	 << "dedges   = " << get_keys(uf.edgecoords) << ";\n"
	 << "dedgepos = " << get_values(uf.edgecoords) << ";\n";
  output.close();

  ofstream latex_output("output/C"+to_string(N)+"-unfold.tex");

  latex_output << "\\newcommand{\\figunfolding}{"<<uf.to_latex(1,0,2) << "}\n\n";

  Unfolding UF = uf.straighten_lines();

  output.open("output/C"+to_string(N)+"-star.m");

  output << "dual    = " << dual << ";\n"
	 << "outline = " << UF.outline << ";\n"
	 << "dedges   = " << get_keys(UF.edgecoords) << ";\n"
	 << "dedgepos = " << get_values(UF.edgecoords) << ";\n";

  output.close();

  latex_output << "\\newcommand{\\figstar}{"<<UF.to_latex(1,0,1) << "}\n\n";
  
  latex_output.close();

  return 0;
}
