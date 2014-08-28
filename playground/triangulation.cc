#include "libgraph/fullerenegraph.hh"
using namespace std;

int main()
{
  FullereneGraph g(stdin);
  g.layout2d = g.tutte_layout();
  g.spherical_layout = g.spherical_projection(g.layout2d);
  Graph triangles(g.triangulation(6));
  triangles.layout2d = g.layout2d;

  cout << "g = " << g << ";\n";
  cout << "gT = " << triangles << ";\n";

  return 0;
}
