#include "libgraph/fullerenegraph.hh"

int main(int ac, char **av)
{

  vector<int> rspi(12);
  FullereneGraph::jumplist_t jumps;

  if(ac<13) return -1;

  int N = strtol(av[1],0,0);
  for(int i=0;i<12;i++) rspi[i] = strtol(av[2+i],0,0)-1;

  FullereneGraph g(N,rspi);

  g.layout2d = g.tutte_layout();
  Graph dual(g.dual_graph(6));
 
  cout << "graph = " << g    << endl;
  cout << "dual  = " << dual << endl;

  return 0;
}
