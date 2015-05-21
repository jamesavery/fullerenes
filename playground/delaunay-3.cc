#include "libgraph/triangulation.hh"
#include "libgraph/delaunay.hh"
#include "libgraph/polyhedron.hh"
#include "libgraph/fullerenegraph.hh"
#include "libgraph/debug.hh"

#include <fstream>

Triangulation regular_polygon(int n)
{
  vector<tri_t> triangles(n);
  Graph P(n+2,vector<int>((n+2)*(n+2),0));;
  node_t v = n, w = n+1;  
  for(node_t u=0;u<n;u++) P.insert_edge(edge_t(u,(u+1)%n)); // Connect cycle
  for(node_t u=0;u<n;u++){ P.insert_edge(edge_t(u,v)); P.insert_edge(edge_t(u,w)); }   // Connect centers
  
  return Triangulation(P);
}


int main(int ac, char **av) 
{
  const size_t n = 6;
  Triangulation pT(regular_polygon(n));

  Debug::channel_verbosity["Delaunay"] = Debug::INFO3;

  FulleroidDelaunay DY(pT);

  vector<double> angles(n);
  for(node_t u=0;u<n;u++) angles[u] = DY.angle_d6y(u,n,(u+1)%n);

  vector<node_t> hole(n);
  for(node_t u=0;u<n;u++) hole[u] = u;

  vector<double> distances = DY.new_distances_d6y(n,hole);

  cout << "angles    = " << angles << ";\n";
  cout << "distances = " << distances << ";\n";
  

  return 0;
}
