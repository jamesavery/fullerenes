#include "libgraph/fullerenegraph.hh"
#include "libgraph/layout2d.hh"
using namespace std;

int main()
{
  FullereneGraph g(stdin);
  vector<coord2d> layout = g.tutte_layout();
  auto spherical_layout = layout2d::spherical_projection(g, layout);
  Graph triangles(g.triangulation(6));

  cout << "g = " << g << ";\n";
  cout << "gT = " << triangles << ";\n";

  return 0;
}
