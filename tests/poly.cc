#include "libgraph/fullerenegraph.hh"
#include "libgraph/polyhedron.hh"

coord3d surface(const coord2d& polar_angle)
{
  double theta = polar_angle.first, phi = polar_angle.second;
  double r = 1.0;// + .5*cos(2*theta)*sin(2*phi);
  return coord3d(cos(theta)*sin(phi), sin(theta)*sin(phi), cos(phi))*r;
}

vector<coord3d> map_surface(const vector<coord2d>& angles)
{
  vector<coord3d> S(angles.size());
  for(unsigned int i=0;i<angles.size();i++) S[i] = surface(angles[i]);
  return S;
}

int main(int ac, char **av)
{
  int name(ac<=1? 1 : atoi(av[1]));
  
  fprintf(stderr,"Constructing high order Halma fullerene.\n");
  FullereneGraph g(FullereneGraph::C20().halma_fullerene(name));
  fprintf(stderr,"Constructing Tutte layout\n");
  g.layout2d = g.tutte_layout();

  fprintf(stderr,"Computing faces\n");
  vector<face_t> faces(g.compute_faces_flat(6, g.layout2d));
  
  fprintf(stderr,"Spherical projection\n");
  vector<coord2d> angles(g.spherical_projection(g.layout2d));
  vector<coord3d> points(map_surface(angles)); // Convex hull is {0,...,59} \ {0,2,14,20,50}

  Polyhedron P(g,points,6);

  cout << "sph" << name << "P = " << P << endl;
  cout << "sph" << name <<"nodes = sph"<<name<<"P[[1]];\n";
  cout << "sph" << name <<"xs    = sph"<<name<<"P[[2]];\n";
  cout << "sph" << name <<"graph = sph"<<name<<"P[[3]];\n";
  cout << "sph" << name << "volume = " << P.volume() << endl;
  cout << "sph" << name << "volume2 = " << P.volume2() << endl;
  cout << "sph" << name << "area = " << P.surface_area() << endl;

  return 0;
}
