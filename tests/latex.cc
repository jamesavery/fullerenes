#include "libgraph/fullerenegraph.hh"
using namespace std;

int main(int ac, char **av)
{
  bool show_dual=false, number_vertices=false, include_latex_header=false;
  double width_cm=10, height_cm=10;
  // Parse command line options
  for(int i=1;i<ac;i++){
    if(string(av[i]) == "show_dual")            show_dual = true;
    if(string(av[i]) == "number_vertices")      number_vertices = true;
    if(string(av[i]) == "include_latex_header") include_latex_header = true;
    if(string(av[i]) == "size"){
      if(++i>=ac) abort();
      height_cm = strtod(av[i],0);
      width_cm = height_cm;
    }
  }

  CubicGraph g(stdin);

  if(g.layout2d.size() != g.N){
    g.layout2d = g.tutte_layout();
  }

  coord2d wh(g.width_height());
  g.scale(coord2d(1.0/wh.first,1.0/wh.second));

  cout << g.to_latex(width_cm, height_cm,
		     show_dual,number_vertices,include_latex_header) << endl;
  
  return 0;
}
