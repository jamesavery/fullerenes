#include <fstream>
#include <vector>
#include <list>
#include <vector>
#include <utility> //required for pair

#include "fullerenes/fullerenegraph.hh"
#include "fullerenes/triangulation.hh"
#include "fullerenes/layout2d.hh"

// Creates the m-point halma-fullerene from the current fullerene C_n with n(1+m)^2 vertices. (I.e. 4,9,16,25,36,... times)
FullereneGraph FullereneGraph::halma_fullerene(const int m, const bool) const {
  if(m<0) return FullereneGraph(*this);
  Triangulation dual(dual_graph(6));
  Triangulation halma = dual.halma_transform(m);
  return FullereneGraph(halma.dual_graph());
}

unsigned int gcd(unsigned int a, unsigned int b)
{
  unsigned int r = a % b;
  if(r == 0) return b; else return gcd(b,r);
}


// Actually works for all cubic graphs -- perhaps stick it there instead
// works for all CG, but here we know, the maximum ring size is 6
FullereneGraph FullereneGraph::GCtransform(unsigned k, unsigned l) const
{
  Triangulation t(dual_graph(6));
  Triangulation t_inflated(t.GCtransform(k,l));
  FullereneGraph fg(t_inflated.dual_graph());
  return fg;
}

// Creates the next leapfrog fullerene C_{3n} from the current fullerene C_n
FullereneGraph FullereneGraph::leapfrog_fullerene(const bool) const {
  // The leapfrog is equivalent to GC(1,1) = halma(1) on the dual triangulation
  Triangulation dual(dual_graph(6));
  Triangulation lf = dual.GCtransform(1,1);
  return FullereneGraph(lf.dual_graph());
}



// both the pentagon indices and the jumps start at 0
// n is the number of vertices
FullereneGraph::FullereneGraph(const int n, const vector<int>& spiral_indices, const jumplist_t& jumps) : CubicGraph() {
  assert(spiral_indices.size() == 12);

  const int n_faces = n/2 + 2;
  vector<int> spiral_string(n_faces,6);
  for(int i=0;i<spiral_indices.size();i++) spiral_string[spiral_indices[i]] = 5;

  Triangulation dual(spiral_string,jumps);
  Graph G(dual.dual_graph());

  *this = G;
}


// pentagon indices and jumps start to count at 0
// perform a general general spiral search and return 12 pentagon indices and the jump positions + their length
bool FullereneGraph::get_rspi_from_fg(const node_t f1, const node_t f2, const node_t f3, vector<int> &rspi, jumplist_t &jumps, const bool general) const
{
  rspi.clear();
  jumps.clear();

  FullereneDual FDual = Triangulation(this->dual_graph(6));

  if(!FDual.get_rspi(f1, f2, f3, rspi, jumps, general)) return false;
  assert(rspi.size()==12);
  return true;
}

// pentagon indices and jumps start to count at 0
// perform the canonical general general spiral search and return 12 pentagon indices and the jump positions + their length
bool FullereneGraph::get_rspi_from_fg(vector<int> &rspi, jumplist_t &jumps, const bool general, const bool pentagon_start) const
{
  rspi.clear();
  jumps.clear();

  FullereneDual FDual = Triangulation(this->dual_graph(6));

  if(!FDual.get_rspi(rspi, jumps, general, pentagon_start)) return false;
  assert(rspi.size()==12);
  return true;
}


// create a matrix that holds the topological distances between all pentagons
matrix<int> FullereneGraph::pentagon_distance_mtx() const {
  return Triangulation(dual_graph()).pentagon_distance_mtx();
}

vector<coord3d> FullereneGraph::zero_order_geometry(double scalerad) const
{
  vector<coord2d> flat_layout = tutte_layout();
  vector<coord2d> angles(layout2d::spherical_projection(*this, flat_layout));

  // Spherical projection
  vector<coord3d> coordinates(N);
  for(int i=0;i<N;i++){
    double theta = angles[i].first, phi = angles[i].second;
    double x = cos(theta)*sin(phi), y = sin(theta)*sin(phi), z = cos(phi);
    coordinates[i] = coord3d(x,y,z);
  }

  // Move to centroid
  coord3d cm;
  for(node_t u=0;u<N;u++) cm += coordinates[u];
  cm /= double(N);
  coordinates -= cm;

  // Scale spherical projection
  double Ravg = 0;
  for(node_t u=0;u<N;u++)
    for(int i=0;i<3;i++) Ravg += (coordinates[u]-coordinates[(*this)[u][i]]).norm();
  Ravg /= (3.0*N);

  coordinates *= scalerad*1.5/Ravg;

  return coordinates;
}

extern "C" void sa_optff_(const FullereneGraph **graph, const int *N, const int *ihessian, const int *iprinthessian,
                       const int *iopt,double *Dist,double *ftol,double *force);
extern "C" void default_force_parameters_(const int *iopt, double *parameters);

vector<coord3d> FullereneGraph::optimized_geometry(std::span<const coord3d> points, int opt_method, double ftol) const
{
  //  assert(layout2d.size() == N);
  vector<coord3d> coordinates(points.begin(),points.end());
  vector<double> force_parameters(19);

  default_force_parameters_(&opt_method,&force_parameters[0]);

  // cout << "force parameters: " << force_parameters << endl;
  // cout << "optimization method: " << opt_method << endl;

  // cout << "g = " << *this << ";\n";

  int zero = 0;
  int one = 1;
  const FullereneGraph *g = this;
  sa_optff_(&g,&N,&zero,&zero,&opt_method,(double*)&coordinates[0],&ftol,&force_parameters[0]);

  return coordinates;
}




