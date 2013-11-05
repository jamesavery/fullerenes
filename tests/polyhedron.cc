#include "libgraph/fullerenegraph.hh"
#include "libgraph/polyhedron.hh"
#include "libgraph/isomerdb.hh"


extern "C" void optff_(const FullereneGraph **graph, const int *N, const int *ihessian, const int *iprinthessian,
		       const int *iopt,double *Dist,double *ftol,double *force);
extern "C" void default_force_parameters_(const int *iopt, double *parameters);

vector<coord3d> &operator-=(vector<coord3d>& xs, const coord3d& y)
{
  for(int i=0;i<xs.size();i++) xs[i] -= y;
  return xs;
}

vector<coord3d> &operator*=(vector<coord3d>& xs, const double& y)
{
  for(int i=0;i<xs.size();i++) xs[i] *= y;
  return xs;
}

vector<coord3d> zero_order_geometry(const FullereneGraph& g, double scalerad)
{
  assert(g.layout2d.size() == g.N);
  vector<coord2d> angles(g.spherical_projection());

  // Spherical projection
  vector<coord3d> coordinates(g.N);
  for(int i=0;i<g.N;i++){
    double theta = angles[i].first, phi = angles[i].second;
    double x = cos(theta)*sin(phi), y = sin(theta)*sin(phi), z = cos(phi);
    coordinates[i] = coord3d(x,y,z);
  }

  // Move to centroid
  coord3d cm;
  for(node_t u=0;u<g.N;u++) cm += coordinates[u];
  cm /= double(g.N);
  coordinates -= cm;

  // Scale spherical projection
  double Ravg = 0;
  for(node_t u=0;u<g.N;u++)
    for(int i=0;i<3;i++) Ravg += (coordinates[u]-coordinates[g.neighbours[u][i]]).norm();
  Ravg /= (3.0*g.N);
  
  coordinates *= scalerad*1.5/Ravg;

  return coordinates;
}

vector<coord3d> optimize_fullerene(const FullereneGraph& g, int opt_method = 3, double ftol = 1e-10, double scalerad = 4)
{
  vector<coord3d> coordinates(zero_order_geometry(g,scalerad));

  vector<double> force_parameters(19);
  default_force_parameters_(&opt_method,&force_parameters[0]);

  cout << "force parameters: " << force_parameters << endl;

  int zero = 0, one = 1;
  const FullereneGraph *g_f = &g;
  optff_(&g_f,&g.N,&one,&zero,&opt_method,(double*)&coordinates[0],&ftol,&force_parameters[0]);

  return coordinates;
}

int main(int ac, char **av)
{
  int N = ac>=2? strtol(av[1],0,0) : 20;
  int isomer_number = ac>=3? strtol(av[2],0,0) : 1;
  bool IPR = ac>=4? strtol(av[3],0,0) : false;

  IsomerDB DB(N,IPR);

  FullereneGraph g(IsomerDB::makeIsomer(N,DB.getIsomer(N,isomer_number,IPR)));
  g.layout2d = g.tutte_layout();
  
  vector<coord3d> coordinates = optimize_fullerene(g,3,1e-12,4);

  ofstream output(("output/polyhedron-"+to_string(N)+"-"+to_string(isomer_number)+".m").c_str());

  output << "g = " << g << ";\n";
  output << "coordinates = " << coordinates << ";\n";

  Polyhedron P(g,coordinates,6);
  
  output << "P = " << P << ";\n";
  
  Polyhedron D(P.dual(6,true));

  output << "PD = " << D << ";\n";
  
  output.close();

  return 0;
}
