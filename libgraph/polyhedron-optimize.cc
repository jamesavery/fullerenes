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

struct params_t
{
  PlanarGraph *graph;
  vector<double> *zero_values_dist;
  vector<double> *force_constants_dist;
  vector<double> *force_constants_angle;
};

double polyhedron_pot(const gsl_vector* coordinates, void* parameters)
{

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
    potential_energy += 0.5 * force_constants_dist[i] * pow(coord3d::dist(coord3d(ax,ay,az), coord3d(bx,by,bz)) - zero_values_dist[i], 2);
  }

  facemap_t faces(graph.compute_faces_oriented());
  //iterate over all faces
  for(facemap_t::const_iterator it=faces.begin(); it!=faces.end(); ++it){
    // it->first is the face size
    // iterate over faces of equal size
    for(set<face_t>::const_iterator jt=it->second.begin(); jt!=it->second.end(); ++jt){
//      cout << " face of size-" << it->first << ": " << *jt << endl;
      // iterate over nodes in face
      for (int i=0; i!=it->first; ++i){
//        cout << " 3 nodes: " << (*jt)[(i+ it->first -1) % it->first] << ", " << (*jt)[i] <<", " <<  (*jt)[(i+1) % it->first] << endl;
        const double ax = gsl_vector_get(coordinates, 3* (*jt)[(i+ it->first -1) % it->first]);
        const double ay = gsl_vector_get(coordinates, 3* (*jt)[(i+ it->first -1) % it->first] +1);
        const double az = gsl_vector_get(coordinates, 3* (*jt)[(i+ it->first -1) % it->first] +2);
        const double bx = gsl_vector_get(coordinates, 3* (*jt)[i]);
        const double by = gsl_vector_get(coordinates, 3* (*jt)[i] +1);
        const double bz = gsl_vector_get(coordinates, 3* (*jt)[i] +2);
        const double cx = gsl_vector_get(coordinates, 3* (*jt)[(i+1) % it->first]);
        const double cy = gsl_vector_get(coordinates, 3* (*jt)[(i+1) % it->first] +1);
        const double cz = gsl_vector_get(coordinates, 3* (*jt)[(i+1) % it->first] +2);

        const double angle_beta = coord3d::angle(coord3d(ax,ay,az) - coord3d(bx,by,bz), coord3d(cx,cy,cz) - coord3d(bx,by,bz));
        potential_energy += 0.5 * force_constants_angle[1] * pow(angle_beta - M_PI*(1.0-2.0/it->first),2);
//        cout << angle_beta << ", " << M_PI*(1.0-2.0/it->first) << ", " << it->first << endl;
      }
    }
  }
  
//  cout << "pot_E: " << potential_energy << endl;
  return potential_energy;
}


void polyhedron_grad(const gsl_vector* coordinates, void* parameters, gsl_vector* gradient)
{

  params_t &params = *static_cast<params_t*>(parameters);
  PlanarGraph &graph = *params.graph;
  vector<double> &zero_values_dist = *params.zero_values_dist;
  vector<double> &force_constants_dist = *params.force_constants_dist;
  vector<double> &force_constants_angle = *params.force_constants_angle;

  assert(zero_values_dist.size() == graph.edge_set.size());
  assert(force_constants_dist.size() == graph.edge_set.size());
  assert(force_constants_angle.size() == 2*graph.edge_set.size()); // only for cubic graphs

  vector<coord3d> derivatives(graph.N);  

  set<edge_t>::const_iterator e=graph.edge_set.begin(), ee=graph.edge_set.end();
  for(int i=0; e!=ee; e++, i++){
    const double ax = gsl_vector_get(coordinates,3 * e->first);
    const double ay = gsl_vector_get(coordinates,3 * e->first +1);
    const double az = gsl_vector_get(coordinates,3 * e->first +2);
    const double bx = gsl_vector_get(coordinates,3 * e->second);
    const double by = gsl_vector_get(coordinates,3 * e->second +1);
    const double bz = gsl_vector_get(coordinates,3 * e->second +2);
//    cout << "ax " << ax << " ay " << ay << " az " << az << " bx " << bx << " by " << by << " bz " << bz << endl;
    derivatives[e->first]  += coord3d::dnorm(coord3d(ax,ay,az) - coord3d(bx,by,bz)) * force_constants_dist[i] * (coord3d::dist(coord3d(ax,ay,az), coord3d(bx,by,bz)) - zero_values_dist[i]);
    derivatives[e->second] -= coord3d::dnorm(coord3d(ax,ay,az) - coord3d(bx,by,bz)) * force_constants_dist[i] * (coord3d::dist(coord3d(ax,ay,az), coord3d(bx,by,bz)) - zero_values_dist[i]);
//    cout << "dist(" << i << "): " << coord3d::dist(coord3d(ax,ay,az), coord3d(bx,by,bz)) << endl;
  }
  
  facemap_t faces(graph.compute_faces_oriented());
  //iterate over all faces
  for(facemap_t::const_iterator it=faces.begin(); it!=faces.end(); ++it){
    // it->first is the face size
    // iterate over faces of equal size
    for(set<face_t>::const_iterator jt=it->second.begin(); jt!=it->second.end(); ++jt){
//      cout << " face of size-" << it->first << ": " << *jt << endl;
      // iterate over nodes in face
      for (int i=0; i!=it->first; ++i){
//        cout << " 3 nodes: " << (*jt)[(i+ it->first -1) % it->first] << ", " << (*jt)[i] <<", " <<  (*jt)[(i+1) % it->first] << endl;
        const double ax = gsl_vector_get(coordinates, 3* (*jt)[(i+ it->first -1) % it->first]);
        const double ay = gsl_vector_get(coordinates, 3* (*jt)[(i+ it->first -1) % it->first] +1);
        const double az = gsl_vector_get(coordinates, 3* (*jt)[(i+ it->first -1) % it->first] +2);
        const double bx = gsl_vector_get(coordinates, 3* (*jt)[i]);
        const double by = gsl_vector_get(coordinates, 3* (*jt)[i] +1);
        const double bz = gsl_vector_get(coordinates, 3* (*jt)[i] +2);
        const double cx = gsl_vector_get(coordinates, 3* (*jt)[(i+1) % it->first]);
        const double cy = gsl_vector_get(coordinates, 3* (*jt)[(i+1) % it->first] +1);
        const double cz = gsl_vector_get(coordinates, 3* (*jt)[(i+1) % it->first] +2);
  
        coord3d a(coord3d(ax,ay,az) - coord3d(bx,by,bz)), c(coord3d(cx,cy,cz) - coord3d(bx,by,bz)), da, dc;
        coord3d::dangle(a, c, da, dc);
  
        const double angle_beta = coord3d::angle(coord3d(ax,ay,az) - coord3d(bx,by,bz), coord3d(cx,cy,cz) - coord3d(bx,by,bz));
//        cout << angle_beta << ", " << M_PI*(1.0-2.0/it->first) << ", " << it->first << endl;
//        cout << "da: " << da << endl;

        derivatives[(*jt)[(i+ it->first -1) % it->first]] += da *       (angle_beta - M_PI*(1.0-2.0/it->first)) * force_constants_angle[(*jt)[i]];
        derivatives[(*jt)[i]]                             += -(da+dc) * (angle_beta - M_PI*(1.0-2.0/it->first)) * force_constants_angle[(*jt)[i]];
        derivatives[(*jt)[(i+1) % it->first]]             += dc *       (angle_beta - M_PI*(1.0-2.0/it->first)) * force_constants_angle[(*jt)[i]];
      }
    }
  }

  // return gradient
  for(int i = 0; i < graph.N; ++i) {
    gsl_vector_set(gradient, 3*i, derivatives[i][0]);
    gsl_vector_set(gradient, 3*i+1, derivatives[i][1]);
    gsl_vector_set(gradient, 3*i+2, derivatives[i][2]);
  }

//   double grad_debug = 0.0;
//   for(int i = 0; i < graph.N; ++i) {
//     grad_debug += sqrt(pow(gsl_vector_get(gradient, 3*i),2) + pow(gsl_vector_get(gradient, 3*i+1),2) + pow(gsl_vector_get(gradient, 3*i+2),2));
//   }
//   cout << "grad: " << grad_debug << endl;

}


//FIXME redo for performance reasons
void polyhedron_pot_grad(const gsl_vector* coordinates, void* parameters, double* potential, gsl_vector* gradient)
{
  *potential = polyhedron_pot(coordinates, parameters);
  polyhedron_grad(coordinates, parameters, gradient);
}


bool Polyhedron::optimize_other(vector<double> zero_values_dist)
{
//  cout << "entering opt other" << endl;

  // settings for the optimizations
  const double stepsize = 1e-3;// FIXME final value
  const double terminate_gradient = 1e-5;// FIXME final value
  const double tol = 1e-1; // accuracy of line minimization, the manual suggests 0.1
  const int max_iterations = 500;// FIXME final value
  
  // init values
  vector<double> force_constants_dist(edge_set.size(), 500.0);
  if(zero_values_dist.size() != edge_set.size())
  {
    zero_values_dist.resize(edge_set.size(), 1.4);
  }
  vector<double> force_constants_angle(2 * edge_set.size(), 200.0); // FIXME possibly change later

  params_t params;
  params.graph = this; //FIXME is this legal?
  params.zero_values_dist = &zero_values_dist;
  params.force_constants_dist = &force_constants_dist;
  params.force_constants_angle = &force_constants_angle;
  
  // Create the minimisation function block, define the different functions and parameters
  gsl_multimin_function_fdf potential_function;
  potential_function.n = 3*N;
  potential_function.f = &polyhedron_pot;
  potential_function.df = &polyhedron_grad;
  potential_function.fdf = &polyhedron_pot_grad;
  potential_function.params = static_cast<void*>(&params);
  
  // Starting point
  gsl_vector* coordinates = gsl_vector_alloc(potential_function.n);

  for(int i=0; i<N; ++i){
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
  
//    cout << "it: " << iter << " stat: " << status << endl;
//    printf ("error: %s\n", gsl_strerror (status));

    if(status) break;
  
    status = gsl_multimin_test_gradient(s->gradient, terminate_gradient);
//    cout << "Status 2: " << status << endl;
  
  } while (status == GSL_CONTINUE && iter < max_iterations);
  
  for(int i=0; i<N; ++i){
    points[i][0] = gsl_vector_get(s->x, 3*i);
    points[i][1] = gsl_vector_get(s->x, 3*i+1);
    points[i][2] = gsl_vector_get(s->x, 3*i+2);
  }

  gsl_multimin_fdfminimizer_free(s);
  gsl_vector_free(coordinates);
  
  return status==0 ? true : false;
}

#else
bool Polyhedron::optimize_other(vector<double> unused)
{
  cerr << "Optimizing other polyhedra than fullerenes is only available through GSL." << endl;
  return 0;
}
#endif

