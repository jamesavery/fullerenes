#include "planargraph.hh"
#include "polyhedron.hh"

#include "math.h"

#include <gsl/gsl_vector.h>
#include <gsl/gsl_multimin.h>

using namespace std;

inline double dist(const double ax, const double ay, const double az, const double bx, const double by, const double bz){
  return sqrt(pow(ax-bx,2) + pow(ay-by,2) + pow(az-bz,2));
}

// subroutine angle takes 9 reals (=3 coordinates) and yields an angel between 0 and +\pi (in radians)
// via law of cosines
double angle(const double ax, const double ay, const double az, const double bx, const double by, const double bz, const double cx, const double cy, const double cz){
      const double r2L=dist(ax,ay,az,bx,by,bz);
      const double r2M=dist(ax,ay,az,cx,cy,cz);
      const double r2R=dist(bx,by,bz,cx,cy,cz);
      const double den=2.0 * sqrt(r2L * r2R);
      double arg=(r2L+r2R-r2M)/den;
// the following two exceptions may be called in case of rounding errors
      if(arg > 1.0) arg=1.0;
      if(arg < -1.0) arg=-1.0;
      return acos(arg);
}



// double my_f (const gsl_vector *v, void *params)
// {
// double x, y;
// double *p = (double *)params;
// x = gsl_vector_get(v, 0);
// y = gsl_vector_get(v, 1);
// return p[2] * (x - p[0]) * (x - p[0]) +
// p[3] * (y - p[1]) * (y - p[1]) + p[4];
// }

double pot(const gsl_vector* coordinates, void* parameters){
  // parameters is a pair of a vector of force constants and a graph

  pair<vector<double>,PlanarGraph> *p = static_cast<pair<vector<double>,PlanarGraph> *>(parameters);
  
  // iterate over edges
  //  V = k (r - r0)**2
  double potential_energy = 0.0;
  double edge_zero_value = 1.0; // arbitrary number
  for(set<edge_t>::const_iterator e=p->second.edge_set.begin(); e!=p->second.edge_set.begin();e++){
    const double ax = gsl_vector_get(coordinates,3 * e->first);
    const double ay = gsl_vector_get(coordinates,3 * e->first +1);
    const double az = gsl_vector_get(coordinates,3 * e->first +2);
    const double bx = gsl_vector_get(coordinates,3 * e->second);
    const double by = gsl_vector_get(coordinates,3 * e->second +1);
    const double bz = gsl_vector_get(coordinates,3 * e->second +2);
    potential_energy += p->first[0] * pow(dist(ax,ay,az,bx,by,bz) - edge_zero_value,2);
  }

  //  determine size
  //  iterate over all angles in face
  //    V = k (a - a0)**2
  // get map<face size, face>
  facemap_t faces(p->second.compute_faces_oriented());
  //iterate over all faces
  for(facemap_t::const_iterator it = faces.begin(); it != faces.end(); ++it) {
    // iterate over vertices in face
    vector<int> f(it->first,0); set<face_t>::const_iterator jt=it->second.begin(); for(int i=0; i!=it->first; ++i, ++jt){f[i] = (*jt)[i];} // only because sets cannot be subscripted
    for(unsigned int i=0; i!=it->first; ++i){
      const double ax = gsl_vector_get(coordinates, 3*f[(i+ it->first -1) % it->first]);
      const double ay = gsl_vector_get(coordinates, 3*f[(i+ it->first -1) % it->first] +1);
      const double az = gsl_vector_get(coordinates, 3*f[(i+ it->first -1) % it->first] +2);
      const double bx = gsl_vector_get(coordinates, 3*f[i]);
      const double by = gsl_vector_get(coordinates, 3*f[i] +1);
      const double bz = gsl_vector_get(coordinates, 3*f[i] +2);
      const double cx = gsl_vector_get(coordinates, 3*f[(i+1) % it->first]);
      const double cy = gsl_vector_get(coordinates, 3*f[(i+1) % it->first] +1);
      const double cz = gsl_vector_get(coordinates, 3*f[(i+1) % it->first] +2);
      potential_energy += p->first[1] * pow(angle(ax,ay,az,bx,by,bz,cx,cy,cz) - M_PI*(1-2/it->first),2);
    }
  }

  return potential_energy;
}

// void my_df (const gsl_vector *v, void *params, gsl_vector *df)
// {
// double x, y;
// double *p = (double *)params;
// x = gsl_vector_get(v, 0);
// y = gsl_vector_get(v, 1);
// gsl_vector_set(df, 0, 2.0 * p[2] * (x - p[0]));
// gsl_vector_set(df, 1, 2.0 * p[3] * (y - p[1]));
// }

void grad(const gsl_vector* coordinates, void* parameters, gsl_vector* gradient){


}


//FIXME redo
void pot_grad(const gsl_vector* coordinates, void* parameters, double* potential, gsl_vector* gradient) {
  *potential = pot(coordinates, parameters);
  grad(coordinates, parameters, gradient);
}


bool Polyhedron::optimize_other(){

  const double stepsize = 1e-3;// FIXME final value
  const double terminate_gradient = 1e-5;// FIXME final value
  const double tol = 1e-1; // accuracy of line minimization, the manual suggests 0.1
  const int max_iterations = 1000;// FIXME final value
  
//  const int number_magnetic_field_vectors = conf.special().rdc.nmf;
//  const int number_magnetic_field_vector_components = 2*number_magnetic_field_vectors;
  // Create the parameters-struct, fill in all needed variabled for the potential function and the derivative
//  fit_param<B> parameters;
//  parameters.topo = &topo;
//  parameters.conf = &conf;
//  parameters.sim = &sim;

// gsl_multimin_function_fdf my_func;
// /* Paraboloid center at (1,2), scale factors (10, 20),
// minimum value 30 */
// double p[5] = { 1.0, 2.0, 10.0, 20.0, 30.0 };
// my_func.n = 2; /* number of function components */
// my_func.f = &my_f;
// my_func.df = &my_df;
// my_func.fdf = &my_fdf;
// my_func.params = (void *)p;


  double fc_array[] = {200.0, 500.0};
  vector<double> force_constants(fc_array, fc_array+2);

  pair<vector<double>,PlanarGraph> parameters = make_pair(force_constants, *this);

  // Create the minimisation function block, define the different functions and parameters
  gsl_multimin_function_fdf potential_function;
  potential_function.n = 3*N;
  potential_function.f = &pot;
  potential_function.df = &grad;
  potential_function.fdf = &pot_grad;
  potential_function.params = (void *) &parameters;
  
  // Starting point
  gsl_vector* coordinates = gsl_vector_alloc(3*N);

  for(int i=0; i < N; ++i){
    gsl_vector_set(coordinates, 3*i, points[i][0]);
    gsl_vector_set(coordinates, 3*i+1, points[i][1]);
    gsl_vector_set(coordinates, 3*i+2, points[i][2]);
  }


// int main (void)
// {
// size_t iter = 0;
// int status;
// const gsl_multimin_fdfminimizer_type *T;
// gsl_multimin_fdfminimizer *s;
// /* Position of the minimum (1,2), scale factors
// 10,20, height 30. */
// double par[5] = { 1.0, 2.0, 10.0, 20.0, 30.0 };
// gsl_vector *x;
// gsl_multimin_function_fdf my_func;
// my_func.n = 2;
// my_func.f = my_f;
// my_func.df = my_df;
// my_func.fdf = my_fdf;
// my_func.params = par;
// /* Starting point, x = (5,7) */
// x = gsl_vector_alloc (2);
// gsl_vector_set (x, 0, 5.0);
// gsl_vector_set (x, 1, 7.0);
// T = gsl_multimin_fdfminimizer_conjugate_fr;
// s = gsl_multimin_fdfminimizer_alloc (T, 2);
// gsl_multimin_fdfminimizer_set (s, &my_func, x, 0.01, 1e-4);
// do
// {
// iter++;
// status = gsl_multimin_fdfminimizer_iterate (s);
// if (status)
// break;
// status = gsl_multimin_test_gradient (s->gradient, 1e-3);
// if (status == GSL_SUCCESS)
// printf ("Minimum found at:\n");
// printf ("%5d %.5f %.5f %10.5f\n", iter,
// gsl_vector_get (s->x, 0),
// gsl_vector_get (s->x, 1),
// s->f);
// }
// while (status == GSL_CONTINUE && iter < 100);
// gsl_multimin_fdfminimizer_free (s);
// gsl_vector_free (x);
// return 0;
// }

  
  const gsl_multimin_fdfminimizer_type *fT = gsl_multimin_fdfminimizer_conjugate_fr;
  gsl_multimin_fdfminimizer *s = gsl_multimin_fdfminimizer_alloc(fT, 3*N);
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
