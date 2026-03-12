#include "fullerenes/fullerenegraph.hh"
#include "fullerenes/unfold.hh"


int main(int ac, char **av)
{
  assert(ac >= 13);

  int N = strtol(av[1],0,0);
  vector<int> spiral(12);
  for(int i=0;i<12;i++) spiral[i] = strtol(av[i+2],0,0)-1;

  FullereneGraph G(N,spiral);

  //  PlanarGraph dual(G.dual_graph());
  Triangulation dual(G.leapfrog_dual());
  cout << "dual = " << dual << ";\n";

  Unfolding uf(dual);

  cout << "outline = " << uf.outline << ";\n";

  Unfolding UF = uf.straighten_lines();

  cout << "OUTLINE = " << UF.outline << ";\n";

  return 0;
}
