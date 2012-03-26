#include "libgraph/fullerenegraph.hh"
#include "libgraph/polyhedron.hh"

coord3d surface(const coord2d& polar_angle)
{
  double theta = polar_angle.first, phi = polar_angle.second;
  double r = 1.0;// + .2*cos(2*theta)*sin(2*phi);
  return coord3d(cos(theta)*sin(phi), sin(theta)*sin(phi), cos(phi))*r;
}

vector<coord3d> map_surface(const vector<coord2d>& angles)
{
  vector<coord3d> S(angles.size());
  for(unsigned int i=0;i<angles.size();i++) S[i] = surface(angles[i]);
  return S;
}

int main()
{
  FullereneGraph g(FullereneGraph::C20().leapfrog_fullerene());

  g.layout2d = g.tutte_layout();
  vector<face_t> faces(g.compute_faces_flat(6, g.layout2d));
  
  vector<coord2d> angles(g.spherical_projection(g.layout2d));
  vector<coord3d> points(map_surface(angles)); // Convex hull is {0,...,59} \ {0,2,14,20,50}

  Polyhedron P(g,points,6);

  cout << "P = " << P << endl;

  Polyhedron hull(P.convex_hull());

  cout << "hull = " << hull << endl;
  cout << "volume = " << P.volume() << endl;
  cout << "area = " << P.surface_area() << endl;

  return 0;
}
