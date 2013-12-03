#include "planargraph.hh"
#include "polyhedron.hh"
#include "geometry.hh"

#include "math.h"

#ifdef HAS_GSL
#include <gsl/gsl_vector.h>
#include <gsl/gsl_multimin.h>
#endif

using namespace std;


#ifdef HAS_GSL

struct params_t {
  PlanarGraph *graph;
  vector<double> *zero_values_dist;
  vector<double> *force_constants_dist;
  vector<double> *force_constants_angle;
};

double pot(const gsl_vector* coordinates, void* parameters){
  // parameters is a pair of a vector of force constants and a graph

  params_t &params = *static_cast<params_t*>(parameters);
  PlanarGraph &graph = *params.graph;
  vector<double> &zero_values_dist = *params.zero_values_dist;
  vector<double> &force_constants_dist = *params.force_constants_dist;
  vector<double> &force_constants_angle = *params.force_constants_angle;

  assert(zero_values_dist.size() == graph.edge_set.size());
  assert(force_constants_dist.size() == graph.edge_set.size());
  assert(force_constants_angle.size() == 2*graph.edge_set.size()); // only for cubic graphs

  // iterate over edges
  //  V = k (r - r0)**2
  double potential_energy = 0.0;
  int i=0;
  set<edge_t>::const_iterator e=graph.edge_set.begin(), ee=graph.edge_set.end();
  for(; e!=ee; e++, i++){
    const double ax = gsl_vector_get(coordinates,3 * e->first);
    const double ay = gsl_vector_get(coordinates,3 * e->first +1);
    const double az = gsl_vector_get(coordinates,3 * e->first +2);
    const double bx = gsl_vector_get(coordinates,3 * e->second);
    const double by = gsl_vector_get(coordinates,3 * e->second +1);
    const double bz = gsl_vector_get(coordinates,3 * e->second +2);
    potential_energy += 0.5 * force_constants_dist[i] * pow(coord3d::dist(coord3d(ax,ay,az), coord3d(bx,by,bz)) - zero_values_dist[i], 2); //FIXME why cannot I just say dist?
  }

  //  determine size
  //  iterate over all angles in face
  //    V = k (a - a0)**2
  // get map<face size, face>
  facemap_t faces(graph.compute_faces_oriented());
  //iterate over all faces
  cout << "number of faces: " << faces.size() << endl;
  for(facemap_t::const_iterator it=faces.begin(); it!=faces.end(); ++it){
    // iterate over vertices in face
    vector<int> f(it->first,0); // it->first is the face size
    set<face_t>::const_iterator jt=it->second.begin();
    // cout << "it first: " << it->first << endl;
    // cout << "*jt: " << *jt << endl;
    for(int i=0; i!=it->first; ++i){f[i] = (*jt)[i];} // only because sets cannot be subscripted
    // cout << "f: " << f << endl;
    for(int i=0; i!=it->first; ++i){
      const double ax = gsl_vector_get(coordinates, 3*f[(i+ it->first -1) % it->first]);
      const double ay = gsl_vector_get(coordinates, 3*f[(i+ it->first -1) % it->first] +1);
      const double az = gsl_vector_get(coordinates, 3*f[(i+ it->first -1) % it->first] +2);
      const double bx = gsl_vector_get(coordinates, 3*f[i]);
      const double by = gsl_vector_get(coordinates, 3*f[i] +1);
      const double bz = gsl_vector_get(coordinates, 3*f[i] +2);
      const double cx = gsl_vector_get(coordinates, 3*f[(i+1) % it->first]);
      const double cy = gsl_vector_get(coordinates, 3*f[(i+1) % it->first] +1);
      const double cz = gsl_vector_get(coordinates, 3*f[(i+1) % it->first] +2);
      potential_energy += 0.5 * force_constants_angle[1] * pow(coord3d::angle(coord3d(ax,ay,az) - coord3d(bx,by,bz), coord3d(ax,ay,az) - coord3d(cx,cy,cz)) - M_PI*(1-2/it->first),2); // FIXME an angle function with three arguments would be more convenient
    }
  }

  return potential_energy;
}

void grad(const gsl_vector* coordinates, void* parameters, gsl_vector* gradient){


  params_t &params = *static_cast<params_t*>(parameters);
  PlanarGraph &graph = *params.graph;
  vector<double> &zero_values_dist = *params.zero_values_dist;
  vector<double> &force_constants_dist = *params.force_constants_dist;
  vector<double> &force_constants_angle = *params.force_constants_angle;

  assert(zero_values_dist.size() == graph.edge_set.size());
  assert(force_constants_dist.size() == graph.edge_set.size());
  assert(force_constants_angle.size() == 2*graph.edge_set.size()); // only for cubic graphs

  vector<coord3d> derivatives(graph.N);  

  int i=0;
  set<edge_t>::const_iterator e=graph.edge_set.begin(), ee=graph.edge_set.end();
  for(; e!=ee; e++, i++){
    const double ax = gsl_vector_get(coordinates,3 * e->first);
    const double ay = gsl_vector_get(coordinates,3 * e->first +1);
    const double az = gsl_vector_get(coordinates,3 * e->first +2);
    const double bx = gsl_vector_get(coordinates,3 * e->second);
    const double by = gsl_vector_get(coordinates,3 * e->second +1);
    const double bz = gsl_vector_get(coordinates,3 * e->second +2);
    derivatives[e->first] += coord3d::dnorm(coord3d(ax,ay,az) - coord3d(bx,by,bz)) * force_constants_dist[i] * (coord3d::dist(coord3d(ax,ay,az), coord3d(bx,by,bz)) - zero_values_dist[i]);
    derivatives[e->first] -= coord3d::dnorm(coord3d(ax,ay,az) - coord3d(bx,by,bz)) * force_constants_dist[i] * (coord3d::dist(coord3d(ax,ay,az), coord3d(bx,by,bz)) - zero_values_dist[i]);
//    derivatives[e->first] -= force_constants_dist[i] * (coord3d::dist(coord3d(ax,ay,az), coord3d(bx,by,bz)) - zero_values_dist[i]) * coord3d::dnorm(coord3d(ax,ay,az) - coord3d(bx,by,bz));
  }
  
  //  iterate over all angles in face
  //    V = k (a - a0)**2
  // get map<face size, face>
  facemap_t faces(graph.compute_faces_oriented());
  //iterate over all faces
  for(facemap_t::const_iterator it=faces.begin(); it!=faces.end(); ++it){
    // iterate over vertices in face
    vector<int> f(it->first,0); // it->first is the face size
    set<face_t>::const_iterator jt=it->second.begin();
    for(int i=0; i!=it->first; ++i, ++jt){f[i] = (*jt)[i];} // only because sets cannot be subscripted
    for(int i=0; i!=it->first; ++i){
      const double ax = gsl_vector_get(coordinates, 3*f[(i+ it->first -1) % it->first]);
      const double ay = gsl_vector_get(coordinates, 3*f[(i+ it->first -1) % it->first] +1);
      const double az = gsl_vector_get(coordinates, 3*f[(i+ it->first -1) % it->first] +2);
      const double bx = gsl_vector_get(coordinates, 3*f[i]);
      const double by = gsl_vector_get(coordinates, 3*f[i] +1);
      const double bz = gsl_vector_get(coordinates, 3*f[i] +2);
      const double cx = gsl_vector_get(coordinates, 3*f[(i+1) % it->first]);
      const double cy = gsl_vector_get(coordinates, 3*f[(i+1) % it->first] +1);
      const double cz = gsl_vector_get(coordinates, 3*f[(i+1) % it->first] +2);

      coord3d b(coord3d(ax,ay,az) - coord3d(bx,by,bz)), c(coord3d(ax,ay,az) - coord3d(cx,cy,cz)), db, dc;
      coord3d::dangle(b, c, db, dc);

      derivatives[f[(i+ it->first -1) % it->first]] += db * (coord3d::angle(coord3d(ax,ay,az) - coord3d(bx,by,bz), coord3d(ax,ay,az) - coord3d(cx,cy,cz)) - M_PI*(1-2/it->first)) * force_constants_angle[f[i]];
      derivatives[f[i]] -= (db+dc) * (coord3d::angle(coord3d(ax,ay,az) - coord3d(bx,by,bz), coord3d(ax,ay,az) - coord3d(cx,cy,cz)) - M_PI*(1-2/it->first)) * force_constants_angle[f[i]];
      derivatives[f[(i+1) % it->first]] += dc * (coord3d::angle(coord3d(ax,ay,az) - coord3d(bx,by,bz), coord3d(ax,ay,az) - coord3d(cx,cy,cz)) - M_PI*(1-2/it->first)) * force_constants_angle[f[i]];
    }
  }

  // return gradient
  for(int i = 0; i < graph.N; ++i) {
    gsl_vector_set(gradient, i, derivatives[i][0]);
    gsl_vector_set(gradient, i+1, derivatives[i][1]);
    gsl_vector_set(gradient, i+2, derivatives[i][2]);
  }
}


//FIXME redo
void pot_grad(const gsl_vector* coordinates, void* parameters, double* potential, gsl_vector* gradient) {
  *potential = pot(coordinates, parameters);
  grad(coordinates, parameters, gradient);
}


bool Polyhedron::optimize_other(){
  // cout << "entering opt other" << endl;

  // setings for the optimizations
  const double stepsize = 1e-3;// FIXME final value
  const double terminate_gradient = 1e-5;// FIXME final value
  const double tol = 1e-1; // accuracy of line minimization, the manual suggests 0.1
  const int max_iterations = 1000;// FIXME final value
  
  // init values
  vector<double> zero_values_dist(this->edge_set.size());
  vector<double> force_constants_dist(this->edge_set.size());
  for(int i=0; i!=this->edge_set.size(); i++){
    zero_values_dist[i] = 1.0; // FIXME possibly change later
    force_constants_dist[i] = 500.0; // FIXME possibly change later
  }
  vector<double> force_constants_angle(2 * this->edge_set.size(), 200.0); // FIXME possibly change later

  params_t params;
  params.graph = this; //FIXME is this legal?
  params.zero_values_dist = &zero_values_dist;
  params.force_constants_dist = &force_constants_dist;
  params.force_constants_angle = &force_constants_angle;
  
  // Create the minimisation function block, define the different functions and parameters
  gsl_multimin_function_fdf potential_function;
  potential_function.n = 3*N;
  potential_function.f = &pot;
  potential_function.df = &grad;
  potential_function.fdf = &pot_grad;
  potential_function.params = static_cast<void*>(&params);
  
  // Starting point
  gsl_vector* coordinates = gsl_vector_alloc(potential_function.n);

  for(int i=0; i < N; ++i){
    gsl_vector_set(coordinates, 3*i, points[i][0]);
    gsl_vector_set(coordinates, 3*i+1, points[i][1]);
    gsl_vector_set(coordinates, 3*i+2, points[i][2]);
  }

  const gsl_multimin_fdfminimizer_type *fT = gsl_multimin_fdfminimizer_conjugate_fr;
  gsl_multimin_fdfminimizer *s = gsl_multimin_fdfminimizer_alloc(fT, potential_function.n);
  gsl_multimin_fdfminimizer_set(s, &potential_function, coordinates, stepsize, tol); // runs fdf once
  size_t iter = 0;
  int status;
  do {
    ++iter;
    status = gsl_multimin_fdfminimizer_iterate(s);
  
    if(status) {
      break;
    }
  
    status = gsl_multimin_test_gradient(s->gradient, terminate_gradient);
    //cout << "Status 2: " << status << endl;
  
  } while (status == GSL_CONTINUE && iter < max_iterations);
  
  
  gsl_multimin_fdfminimizer_free(s);
  gsl_vector_free(coordinates);
  

  
  return status==0 ? true : false;
}

#else
bool Polyhedron::optimize_other(){
  cout << "Optimizing other polyhedra than fullerenes is only available through GSL." << endl;
  return 0;
}
#endif

