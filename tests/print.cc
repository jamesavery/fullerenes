#include "libgraph/fullerenegraph.hh"
#include "libgraph/polyhedron.hh"

int main(int ac, char **av)
{
  CubicGraph g(stdin);

  if(ac>=2 && string(av[1]) == "do_layout")
    g.layout2d = g.tutte_layout();

  bool is_fullerene = g.this_is_a_fullerene();
  cerr << "Graph is "<<(is_fullerene?"":"not ")<<"a fullerene.\n";

  cout << "graph = " << g << ";\n";
  vector<coord3d> points(Polyhedron::polar_mapping(g.spherical_projection()));
  
  cout << "points = {"; for(int i=0;i<points.size();i++) cout << points[i] << (i+1<points.size()?",":"};\n");
  //  cout << "dual  = " << g.dual_graph() << ";\n";

  return 0;
}
