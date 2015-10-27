#include <stdlib.h>
#include <fstream>
#include "libgraph/cubicgraph.hh"
#include "libgraph/fullerenegraph.hh"

#include <vector>

using namespace std;


int main(int ac, char **av)
{
  vector<int> rspi(12);
  FullereneGraph::jumplist_t jumps;

  if(ac<13) return -1;

  int N = strtol(av[1],0,0);
  for(int i=0;i<12;i++) rspi[i] = strtol(av[2+i],0,0)-1;

  cout << "Attempting to create graph from spiral indices " << rspi << endl;
 
  FullereneGraph g(N, rspi, jumps);
  g.layout2d = g.tutte_layout();
  PlanarGraph  dg(g.dual_graph(6));
  dg.layout2d = dg.tutte_layout();

  cout << "g = " << g << ";\n";
  cout << "dg = " << dg << ";\n";

  ofstream g_latex("output/spiral-g.tex"), dg_latex("output/spiral-dg.tex");
  
  g_latex  << g.to_latex(20,20,true,true,false,0,0,0xffffff) << endl;
  dg_latex << dg.to_latex(20,20,false,true,false,0,0,0xffffff) << endl;

  g_latex.close(); 
  dg_latex.close();

  return 0;
}
