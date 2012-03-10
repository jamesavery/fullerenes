#include "libgraph/fullerenegraph.hh"

int main()
{
  FullereneGraph g(stdin);
  g.layout2d = g.tutte_layout();
  cout << "graph = " << g << ";\n";
  //  cout << "dual  = " << g.dual_graph() << ";\n";

  return 0;
}
