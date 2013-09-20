#include "libgraph/fullerenegraph.hh"
#include "libgraph/unfold.hh"


int main(int ac, char **av)
{
  assert(ac >= 13);

  int N = strtol(av[1],0,0);
  vector<int> spiral(12);
  for(int i=0;i<12;i++) spiral[i] = strtol(av[i+2],0,0)-1;

  FullereneGraph G(N,spiral);
  G.layout2d = G.tutte_layout();
  
  PlanarGraph dual(G.dual_graph());
  dual.layout2d = dual.tutte_layout();
  cout << "dual = " << dual << ";\n";

  Unfolding uf(dual,true);

  cout << "outline = " << uf.outline << ";\n";

  Unfolding UF = uf.straighten_lines();

  cout << "OUTLINE = " << UF.outline << ";\n";
  
  return 0;
}
