#include "planargraph.hh"
#include "cubicgraph.hh"
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
  vector<double> *k_dist;
  vector<double> *k_angle;
  vector<double> *k_area;
};


double layout_pot(const gsl_vector* coordinates, void* parameters)
{
//  cout << "entering layout_pot" << endl;

  params_t &params = *static_cast<params_t*>(parameters);
  PlanarGraph &graph = *params.graph;
  vector<double> &zero_values_dist = *params.zero_values_dist;
  vector<double> &k_dist = *params.k_dist;
  vector<double> &k_angle = *params.k_angle;
  vector<double> &k_area = *params.k_area;

  set<edge_t> edge_set = graph.undirected_edges();
  const int n_faces = 2 + edge_set.size() - graph.N;

  assert(zero_values_dist.size() == edge_set.size());
  assert(k_dist.size() == edge_set.size());
  assert(k_angle.size() == 2*edge_set.size());
  assert(k_area.size() == n_faces);

//
// DISTANCE TERM
//
  //get average edge_length
  double log_sum_edge_length=0;

  set<edge_t>::const_iterator e=edge_set.begin(), ee=edge_set.end();
  for(; e!=ee; e++){
    //cout << *e << endl;
    const double ax = gsl_vector_get(coordinates, 2 * e->first);
    const double ay = gsl_vector_get(coordinates, 2 * e->first +1);
    const double bx = gsl_vector_get(coordinates, 2 * e->second);
    const double by = gsl_vector_get(coordinates, 2 * e->second +1);
    log_sum_edge_length += log(coord2d(ax-bx,ay-by).norm());
//    cout << "l: " << coord2d(ax-bx,ay-by).norm() << " ll: " << log(coord2d(ax-bx,ay-by).norm()) << endl;
  }
  log_sum_edge_length /= edge_set.size();
  //const double zero_value = exp(log_sum_edge_length); // selfreferentiallity makes the optimization unstable
//  double zero_value = 0.25;
//  cout << "log average length: " << exp(log_sum_edge_length) << endl;

//  //  V = k (r - r0)**2
  double potential_energy = 0.0;
  e=edge_set.begin();
  for(int i=0; e!=ee; ++e, ++i){
    vector<node_t>::const_iterator it1, it2;
    it1 = find (graph.outer_face.begin(), graph.outer_face.end(), e->first);
    it2 = find (graph.outer_face.begin(), graph.outer_face.end(), e->second);
    if (it1 != graph.outer_face.end() && it2 != graph.outer_face.end() && ( it1 == it2+1 || it1 == it2-1 || (it1 == graph.outer_face.begin() && it2 == graph.outer_face.end()-1) || (it1 == graph.outer_face.end()-1 && it2 == graph.outer_face.begin()))){
//      cout << "omitting " << *it1 << "-" << *it2 << endl;
      continue; // edge is part of outer face
    }

    const double ax = gsl_vector_get(coordinates, 2 * e->first);
    const double ay = gsl_vector_get(coordinates, 2 * e->first +1);
    const double bx = gsl_vector_get(coordinates, 2 * e->second);
    const double by = gsl_vector_get(coordinates, 2 * e->second +1);
//    cout << "r: " << coord2d(ax-bx,ay-by).norm() << " r_0: " << zero_values_dist[i] << endl;
    potential_energy += 0.5 * k_dist[i] * pow(coord2d(ax-bx,ay-by).norm() - zero_values_dist[i], 2);
  }



//
// ANGLE TERM
//
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
        const double ax = gsl_vector_get(coordinates, 2* (*jt)[(i+ it->first -1) % it->first]);
        const double ay = gsl_vector_get(coordinates, 2* (*jt)[(i+ it->first -1) % it->first] +1);
        const double bx = gsl_vector_get(coordinates, 2* (*jt)[i]);
        const double by = gsl_vector_get(coordinates, 2* (*jt)[i] +1);
        const double cx = gsl_vector_get(coordinates, 2* (*jt)[(i+1) % it->first]);
        const double cy = gsl_vector_get(coordinates, 2* (*jt)[(i+1) % it->first] +1);

        const double angle_beta = coord3d::angle(coord3d(ax,ay,0) - coord3d(bx,by,0), coord3d(cx,cy,0) - coord3d(bx,by,0));
        potential_energy += 0.5 * k_angle[(*jt)[i]] * pow(angle_beta - M_PI*(1.0-2.0/it->first),2);
        //cout << angle_beta << ", " << M_PI*(1.0-2.0/it->first) << ", " << it->first << endl;
      }
    }
  }

//
// AREA TERM
//
  // find area of outer face:
  double A_tot=0;
  const double bx = gsl_vector_get(coordinates, 2* graph.outer_face[0]);
  const double by = gsl_vector_get(coordinates, 2* graph.outer_face[0] +1);
  // iterate over nodes in face
  for (int i=1; i!=graph.outer_face.size() -1; ++i){
//    cout << " 3 nodes: " << (*jt)[(i+ it->first -1) % it->first] << ", " << (*jt)[i] <<", " <<  (*jt)[(i+1) % it->first] << endl;
    const double ax = gsl_vector_get(coordinates, 2* graph.outer_face[(i)]);
    const double ay = gsl_vector_get(coordinates, 2* graph.outer_face[(i)] +1);
    const double cx = gsl_vector_get(coordinates, 2* graph.outer_face[(i+1)]);
    const double cy = gsl_vector_get(coordinates, 2* graph.outer_face[(i+1)] +1);

    A_tot += ((ax-bx)*(cy-by) - (ay-by)*(cx-bx))/2;
//    cout << "area of one triangle: " << ((ax-bx)*(cy-by) - (ay-by)*(cx-bx))/2 << endl;
  }
//  cout << "area of outer polygon: " << A_tot << endl;

  const double A_av = abs(A_tot)/(n_faces-1); // (excluding the outer face)
//  cout << "average area of polygon: " << A_av << endl;


  //iterate over all faces
  int i=0;
  for(facemap_t::const_iterator it=faces.begin(); it!=faces.end(); ++it){
    // it->first is the face size
    // iterate over faces of equal size
    for(set<face_t>::const_iterator jt=it->second.begin(); jt!=it->second.end(); ++jt){
//      cout << " face of size-" << it->first << ": " << *jt << endl;
      double A=0;
      const double bx = gsl_vector_get(coordinates, 2* (*jt)[0]);
      const double by = gsl_vector_get(coordinates, 2* (*jt)[0] +1);
      // iterate over nodes in face
      for (int i=1; i!=it->first -1; ++i){
//        cout << " 3 nodes: " << (*jt)[(i+ it->first -1) % it->first] << ", " << (*jt)[i] <<", " <<  (*jt)[(i+1) % it->first] << endl;
        const double ax = gsl_vector_get(coordinates, 2* (*jt)[(i)]);
        const double ay = gsl_vector_get(coordinates, 2* (*jt)[(i)] +1);
        const double cx = gsl_vector_get(coordinates, 2* (*jt)[(i+1)]);
        const double cy = gsl_vector_get(coordinates, 2* (*jt)[(i+1)] +1);

        A += ((ax-bx)*(cy-by) - (ay-by)*(cx-bx))/2;
//        cout << "area of one triangle: " << ((ax-bx)*(cy-by) - (ay-by)*(cx-bx))/2 << endl;
      }
//      cout << "area of one polygon: " << A << endl;
      potential_energy += 0.5 * k_area[i] * pow(abs(A) - A_av,2);
      ++i;
    }
    ++i;
  }

//  cout << "pot_E: " << potential_energy << endl;
//  cout << "leaving layout_pot" << endl;
  
  return potential_energy;
}


void layout_grad(const gsl_vector* coordinates, void* parameters, gsl_vector* gradient)
{
//  cout << "entering layout_grad" << endl;

//
// DIST TERM
//
  params_t &params = *static_cast<params_t*>(parameters);
  PlanarGraph &graph = *params.graph;
  vector<double> &zero_values_dist = *params.zero_values_dist;
  vector<double> &k_dist = *params.k_dist;
  vector<double> &k_angle = *params.k_angle;
  vector<double> &k_area = *params.k_area;

  set<edge_t> edge_set = graph.undirected_edges();
  const int n_faces = 2+edge_set.size() - graph.N;

  assert(zero_values_dist.size() == edge_set.size());
  assert(k_dist.size() == edge_set.size());
  assert(k_angle.size() == 2*edge_set.size());
  assert(k_area.size() == n_faces);

  vector<coord2d> derivatives(graph.N);
//  cout << "(ought to be all zero): " << derivatives << endl;

  //get average edge_length

  double log_sum_edge_length=0;
  set<edge_t>::const_iterator e=edge_set.begin(), ee=edge_set.end();
  for(; e!=ee; e++){
    //cout << *e << endl;
    const double ax = gsl_vector_get(coordinates, 2 * e->first);
    const double ay = gsl_vector_get(coordinates, 2 * e->first +1);
    const double bx = gsl_vector_get(coordinates, 2 * e->second);
    const double by = gsl_vector_get(coordinates, 2 * e->second +1);
    log_sum_edge_length += log(coord2d(ax-bx,ay-by).norm());
//    cout << "l: " << coord2d(ax-bx,ay-by).norm() << " ll: " << log(coord2d(ax-bx,ay-by).norm()) << endl;
  }
  log_sum_edge_length /= edge_set.size();
  //const double zero_value = exp(log_sum_edge_length); // selfreferentiallity makes the optimization unstable
//  const double zero_value = 0.25;
//  cout << "log average length: " << exp(log_sum_edge_length) << endl;

  e=edge_set.begin();
  for(int i=0; e!=ee; e++, i++){
    //cout << *e << endl;
    const double ax = gsl_vector_get(coordinates, 2 * e->first);
    const double ay = gsl_vector_get(coordinates, 2 * e->first +1);
    const double bx = gsl_vector_get(coordinates, 2 * e->second);
    const double by = gsl_vector_get(coordinates, 2 * e->second +1);
//    cout << "ax " << ax << " ay " << ay << " az " << az << " bx " << bx << " by " << by << " bz " << bz << endl;
    derivatives[e->first]  += coord2d::dnorm(coord2d(ax-bx,ay-by)) * k_dist[i] * (coord2d(ax-bx,ay-by).norm() - zero_values_dist[i]);
    derivatives[e->second] -= coord2d::dnorm(coord2d(ax-bx,ay-by)) * k_dist[i] * (coord2d(ax-bx,ay-by).norm() - zero_values_dist[i]);
//    cout << "dist(" << i << "): " << coord2d(ax-bx,ay-by).norm() << endl;
//    cout << "gradient: " << derivatives << endl;
  }
  

 //
 // ANGLE TERM
 //
   facemap_t faces(graph.compute_faces_oriented());
   //iterate over all faces
   for(facemap_t::const_iterator it=faces.begin(); it!=faces.end(); ++it){
     // it->first is the face size
     // iterate over faces of equal size
     for(set<face_t>::const_iterator jt=it->second.begin(); jt!=it->second.end(); ++jt){
       //cout << " face of size: " << it->first << ": " << *jt << endl;
       // iterate over nodes in face
       for (int i=0; i!=it->first; ++i){
 //        cout << " 3 nodes: " << (*jt)[(i+ it->first -1) % it->first] << ", " << (*jt)[i] <<", " <<  (*jt)[(i+1) % it->first] << endl;
         const double ax = gsl_vector_get(coordinates, 2* (*jt)[(i+ it->first -1) % it->first]);
         const double ay = gsl_vector_get(coordinates, 2* (*jt)[(i+ it->first -1) % it->first] +1);
         const double bx = gsl_vector_get(coordinates, 2* (*jt)[i]);
         const double by = gsl_vector_get(coordinates, 2* (*jt)[i] +1);
         const double cx = gsl_vector_get(coordinates, 2* (*jt)[(i+1) % it->first]);
         const double cy = gsl_vector_get(coordinates, 2* (*jt)[(i+1) % it->first] +1);
   
         coord3d a(coord3d(ax,ay,0) - coord3d(bx,by,0)), c(coord3d(cx,cy,0) - coord3d(bx,by,0)), da, dc;
         coord3d::dangle(a, c, da, dc);
 //        cout << da << dc << endl;
   
         const double angle_beta = coord3d::angle(coord3d(ax,ay,0) - coord3d(bx,by,0), coord3d(cx,cy,0) - coord3d(bx,by,0));
         //cout << angle_beta << ", " << M_PI*(1.0-2.0/it->first) << ", " << it->first << endl;
 //        cout << "da: " << da << endl;
 
         derivatives[(*jt)[(i+ it->first -1) % it->first]] += coord2d(da[0],da[1]) *       (angle_beta - M_PI*(1.0-2.0/it->first)) * k_angle[(*jt)[i]];
         derivatives[(*jt)[i]]                             += -coord2d((da+dc)[0],(da+dc)[1]) * (angle_beta - M_PI*(1.0-2.0/it->first)) * k_angle[(*jt)[i]];
         derivatives[(*jt)[(i+1) % it->first]]             += coord2d(dc[0],dc[1]) *       (angle_beta - M_PI*(1.0-2.0/it->first)) * k_angle[(*jt)[i]];
       }
     }
   }


 //
 // AREA TERM
 //
   // find area of outer face:
   double A_tot=0;
   const double bx = gsl_vector_get(coordinates, 2* graph.outer_face[0]);
   const double by = gsl_vector_get(coordinates, 2* graph.outer_face[0] +1);
   // iterate over nodes in face
   for (int i=1; i!=graph.outer_face.size() -1; ++i){
     const double ax = gsl_vector_get(coordinates, 2* graph.outer_face[(i)]);
     const double ay = gsl_vector_get(coordinates, 2* graph.outer_face[(i)] +1);
     const double cx = gsl_vector_get(coordinates, 2* graph.outer_face[(i+1)]);
     const double cy = gsl_vector_get(coordinates, 2* graph.outer_face[(i+1)] +1);
 
     A_tot += ((ax-bx)*(cy-by) - (ay-by)*(cx-bx))/2;
     //cout << "area of one triangle: " << ((ax-bx)*(cy-by) - (ay-by)*(cx-bx))/2 << endl;
   }
   //cout << "area of outer polygon: " << A_tot << endl;
   const double A_av = abs(A_tot)/(n_faces-1); // (excluding the outer face)
   //cout << "average area of polygon: " << A_av << endl;
 
   int i=0;
   //iterate over all faces
   for(facemap_t::const_iterator it=faces.begin(); it!=faces.end(); ++it){
     // it->first is the face size
     // iterate over faces of equal size
     for(set<face_t>::const_iterator jt=it->second.begin(); jt!=it->second.end(); ++jt){
 //      cout << " face of size-" << it->first << ": " << *jt << endl;
       double A=0;
       const double bx = gsl_vector_get(coordinates, 2* (*jt)[0]);
       const double by = gsl_vector_get(coordinates, 2* (*jt)[0] +1);
       // iterate over nodes in face
       for (int i=1; i!=it->first -1; ++i){
 //        cout << " 3 nodes: " << (*jt)[(i+ it->first -1) % it->first] << ", " << (*jt)[i] <<", " <<  (*jt)[(i+1) % it->first] << endl;
         const double ax = gsl_vector_get(coordinates, 2* (*jt)[(i)]);
         const double ay = gsl_vector_get(coordinates, 2* (*jt)[(i)] +1);
         const double cx = gsl_vector_get(coordinates, 2* (*jt)[(i+1)]);
         const double cy = gsl_vector_get(coordinates, 2* (*jt)[(i+1)] +1);
 
         A += ((ax-bx)*(cy-by) - (ay-by)*(cx-bx))/2;
         //cout << "area of one triangle: " << ((ax-bx)*(cy-by) - (ay-by)*(cx-bx))/2 << endl;
       }
 //      cout << "area of one polygon: " << A << endl;
       const double sign = A/abs(A);
 //      potential_energy += 0.5 * k_area * pow(sign*A - A_av,2);
 
 //      cout << "face: " << *jt << endl;
 //      cout << "sign: " << sign << endl;
       // iterate over nodes in face
       for (int i=0; i!=it->first; ++i){
         const double ax = gsl_vector_get(coordinates, 2* (*jt)[(i+ it->first -1) % it->first]);
         const double ay = gsl_vector_get(coordinates, 2* (*jt)[(i+ it->first -1) % it->first] +1);
         const double cx = gsl_vector_get(coordinates, 2* (*jt)[(i+1) % it->first]);
         const double cy = gsl_vector_get(coordinates, 2* (*jt)[(i+1) % it->first] +1);
 
 //        cout << "dA/dr(" << i << "): " << coord2d(cy-ay, -(cx-ax))/2 << endl;
 //        cout << "dE/dr(" << i << "): " << coord2d(cy-ay, -(cx-ax))/2 * k_area[i] * (abs(A) - A_av) * sign << endl;
         derivatives[(*jt)[i]] += coord2d(cy-ay, ax-cx)/2 * k_area[i] * (abs(A) - A_av) * sign;
       }
       ++i;
     }
     ++i;
   }

  // fix outer face
//  cout << "d: " << derivatives << endl;
//  cout << "outer face: " << graph.outer_face << endl;
  for(vector<node_t>::iterator it = graph.outer_face.begin(); it != graph.outer_face.end(); ++it){
    derivatives[*it].first = 0;
    derivatives[*it].second = 0;
  }
//  cout << "d: " << derivatives << endl;
  // return gradient
  for(int i=0; i < graph.N; ++i) {
    gsl_vector_set(gradient, 2*i, derivatives[i].first);
    gsl_vector_set(gradient, 2*i+1, derivatives[i].second);
  }

  double grad_debug = 0.0;
  for(int i = 0; i < graph.N; ++i) {
    grad_debug += sqrt(pow(gsl_vector_get(gradient, 2*i),2) + pow(gsl_vector_get(gradient, 2*i+1),2));
  }
//  cout << "norm grad: " << grad_debug << endl;

//  cout << "leaving layout_grad" << endl;
}


void layout_pot_grad(const gsl_vector* coordinates, void* parameters, double* potential, gsl_vector* gradient)
{
  *potential = layout_pot(coordinates, parameters);
  layout_grad(coordinates, parameters, gradient);
}



bool PlanarGraph::optimize_layout(const double zv_dist_inp, const double k_dist_inp, const double k_angle_inp, const double k_area_inp)
{
//  cout << "entering opt layout" << endl;

  if(layout2d.size()!=N){
    layout2d = this->tutte_layout();
  }

  // settings for the optimizations
  const double stepsize = 1e-2;// FIXME final value
  const double terminate_gradient = 1e-2;// FIXME final value
  const double tol = 1e-1; // accuracy of line minimization, the manual suggests 0.1
  const int max_iterations = 5000;// FIXME final value

  set<edge_t> edge_set = undirected_edges();
  const int n_faces = 2 + edge_set.size() - N;
  
  // init values
  vector<double> zero_values_dist(edge_set.size(), zv_dist_inp); // FIXME choose helpful values
  vector<double> k_dist(edge_set.size(), k_dist_inp);
  vector<double> k_angle(3*N, k_angle_inp); // FIXME possibly change later
  vector<double> k_area(n_faces, k_area_inp); // FIXME possibly change later

  params_t params;
  params.graph = this;
  params.zero_values_dist = &zero_values_dist;
  params.k_dist = &k_dist;
  params.k_angle = &k_angle;
  params.k_area = &k_area;
  
  // Create the minimisation function block, define the different functions and parameters
  gsl_multimin_function_fdf potential_function;
  potential_function.n = 2*N;
  potential_function.f = &layout_pot;
  potential_function.df = &layout_grad;
  potential_function.fdf = &layout_pot_grad;
  potential_function.params = static_cast<void*>(&params);
  
  // Starting point
  gsl_vector* coordinates = gsl_vector_alloc(potential_function.n);
  for(int i=0; i<N; ++i){
    gsl_vector_set(coordinates, 2*i, layout2d[i].first);
    gsl_vector_set(coordinates, 2*i+1, layout2d[i].second);
  }

  const gsl_multimin_fdfminimizer_type *fT = gsl_multimin_fdfminimizer_conjugate_fr;
  gsl_multimin_fdfminimizer *s = gsl_multimin_fdfminimizer_alloc(fT, potential_function.n);
  gsl_multimin_fdfminimizer_set(s, &potential_function, coordinates, stepsize, tol); // runs fdf once
  size_t iter = 0;
  int status=0;
  do {
    ++iter;
    status = gsl_multimin_fdfminimizer_iterate(s);
  
//    cout << "it: " << iter << " Status 1: " << status << endl;
//    printf ("error: %s\n", gsl_strerror (status));

//     for(int i=0; i<N; ++i){
//       cout << "intemediate layout ..." << endl;
//       cout << gsl_vector_get(s->x, 2*i) << ", " <<  gsl_vector_get(s->x, 2*i+1) << endl;
//     }
  

    if(status) break;
  
    status = gsl_multimin_test_gradient(s->gradient, terminate_gradient);
//    printf ("error: %s\n", gsl_strerror (status));
    //cout << "Status 2: " << status << endl;
  } while (status == GSL_CONTINUE && iter < max_iterations);
  
  for(int i=0; i<N; ++i){
//    cout << "returning layout ..." << endl;
//    cout << gsl_vector_get(s->x, 2*i) << "(was: " << layout2d[i].first << ")" << endl;
//    cout << gsl_vector_get(s->x, 2*i+1) << "(was: " << layout2d[i].second << ")" << endl;
    layout2d[i].first = gsl_vector_get(s->x, 2*i);
    layout2d[i].second = gsl_vector_get(s->x, 2*i+1);
  }

  gsl_multimin_fdfminimizer_free(s);
  gsl_vector_free(coordinates);
  
  return status==0 ? true : false;
}

#else
bool PlanarGraph::optimize_layout(const double a, const double b, const double c, const double d){
  cerr << "Optimizing layouts is only available through GSL." << endl;
  return 0;
}
#endif

