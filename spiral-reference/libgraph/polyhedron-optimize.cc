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
  Polyhedron *P;
  map<edge_t, double> *zero_values_dist;
  vector<double> *zero_values_dihedral;
  vector<double> *force_constants_dist;
  vector<double> *force_constants_angle;
  vector<double> *force_constants_dihedral;
  set<edge_t> *edge_set;
  bool optimize_angles;
};

double polyhedron_pot(const gsl_vector* coordinates, void* parameters)
{
  //cout << "entering polyhedron pot" << endl;

  params_t &params = *static_cast<params_t*>(parameters);
  Polyhedron &P = *params.P;
  map<edge_t, double> &zero_values_dist = *params.zero_values_dist;
  vector<double> &zero_values_dihedral = *params.zero_values_dihedral;
  vector<double> &force_constants_dist = *params.force_constants_dist;
  vector<double> &force_constants_angle = *params.force_constants_angle;
  vector<double> &force_constants_dihedral = *params.force_constants_dihedral;
  set<edge_t> &edge_set = *params.edge_set;

  assert(zero_values_dist.size() == edge_set.size());
  assert(force_constants_dist.size() == edge_set.size());
  assert(force_constants_angle.size() == 2*edge_set.size()); // only for cubic graphs
  assert(force_constants_dihedral.size() == P.N);


  // iterate over edges
  //  V = k (r - r0)**2
  double potential_energy = 0.0;
  int i=0;
  set<edge_t>::const_iterator e=edge_set.begin(), ee=edge_set.end();
  for(; e!=ee; e++, i++){
    const double ax = gsl_vector_get(coordinates,3 * e->first);
    const double ay = gsl_vector_get(coordinates,3 * e->first +1);
    const double az = gsl_vector_get(coordinates,3 * e->first +2);
    const double bx = gsl_vector_get(coordinates,3 * e->second);
    const double by = gsl_vector_get(coordinates,3 * e->second +1);
    const double bz = gsl_vector_get(coordinates,3 * e->second +2);
    potential_energy += 0.5 * force_constants_dist[i] * pow(coord3d::dist(coord3d(ax,ay,az), coord3d(bx,by,bz)) - zero_values_dist[*e], 2);
  }
  
  
  // it->first is the face size
  // iterate over faces of equal size
  if(params.optimize_angles){
    for(int i=0;i<P.faces.size();i++){
      const face_t &f(P.faces[i]);
      const int face_size = f.size();
    
//      cout << " face of size-" << it->first << ": " << *jt << endl;
      // iterate over nodes in face
      for (int i=0; i<face_size; ++i)
      {
        int f0 = f[(i-1+face_size)%face_size], f1 = f[i], f2 = f[(i+1)%face_size];
        // cout << " 3 nodes: " << f[(i+ it->first -1) % it->first] << ", " << f[i] <<", " <<  f[(i+1) % it->first] << endl;
        const double ax = gsl_vector_get(coordinates, 3*f0    );
        const double ay = gsl_vector_get(coordinates, 3*f0 + 1);
        const double az = gsl_vector_get(coordinates, 3*f0 + 2);
        const double bx = gsl_vector_get(coordinates, 3*f1    );
        const double by = gsl_vector_get(coordinates, 3*f1 + 1);
        const double bz = gsl_vector_get(coordinates, 3*f1 + 2);
        const double cx = gsl_vector_get(coordinates, 3*f2    );
        const double cy = gsl_vector_get(coordinates, 3*f2 + 1);
        const double cz = gsl_vector_get(coordinates, 3*f2 + 2);

        const double angle_beta = coord3d::angle(coord3d(ax,ay,az) - coord3d(bx,by,bz), coord3d(cx,cy,cz) - coord3d(bx,by,bz));
        potential_energy += 0.5 * force_constants_angle[1] * pow(angle_beta - M_PI*(1.0-2.0/double(face_size)),2);
        // cout << angle_beta << ", " << M_PI*(1.0-2.0/it->first) << ", " << it->first << endl;
      }
    }
  }

  // iterate over all dihedrals, i.e., all points
  if(P.is_cubic()){
    for(int u=0; u<P.N; u++)
    {
//        t   B   s
//          \   /
//        C   u   A
//            |
//            r
      node_t r = P.neighbours[u][0];
      node_t s = P.neighbours[u][1];
      node_t t = P.neighbours[u][2];

      const double ax = gsl_vector_get(coordinates, 3*u    );
      const double ay = gsl_vector_get(coordinates, 3*u + 1);
      const double az = gsl_vector_get(coordinates, 3*u + 2);
      const double bx = gsl_vector_get(coordinates, 3*r    );
      const double by = gsl_vector_get(coordinates, 3*r + 1);
      const double bz = gsl_vector_get(coordinates, 3*r + 2);
      const double cx = gsl_vector_get(coordinates, 3*s    );
      const double cy = gsl_vector_get(coordinates, 3*s + 1);
      const double cz = gsl_vector_get(coordinates, 3*s + 2);
      const double dx = gsl_vector_get(coordinates, 3*t    );
      const double dy = gsl_vector_get(coordinates, 3*t + 1);
      const double dz = gsl_vector_get(coordinates, 3*t + 2);

      const double dihedral_abcd = coord3d::dihedral(coord3d(bx,by,bz) - coord3d(ax,ay,az), coord3d(cx,cy,cz) - coord3d(ax,ay,az), coord3d(dx,dy,dz) - coord3d(ax,ay,az));
      potential_energy += 0.5 * force_constants_dihedral[u] * pow(dihedral_abcd - zero_values_dihedral[u],2);
    
//      cout << "value: " << dihedral_abcd << ", " << zero_values_dihedral[u] << ", " << u << endl;
    }
  }

//  // NB: Try out, then remove or rewrite
//  for(int i=0;i<P.points.size();i++){
//    const coord3d x(gsl_vector_get(coordinates,3*i),gsl_vector_get(coordinates,3*i+1),gsl_vector_get(coordinates,3*i+2));
//    potential_energy += 2000/x.dot(x);
//  }

  //  cout << "pot_E: " << potential_energy << endl;
  return potential_energy;
}


void polyhedron_grad(const gsl_vector* coordinates, void* parameters, gsl_vector* gradient)
{

  params_t &params = *static_cast<params_t*>(parameters);
  Polyhedron &P = *params.P;
  map<edge_t, double> &zero_values_dist = *params.zero_values_dist;
  vector<double> &zero_values_dihedral = *params.zero_values_dihedral;
  vector<double> &force_constants_dist = *params.force_constants_dist;
  vector<double> &force_constants_angle = *params.force_constants_angle;
  vector<double> &force_constants_dihedral = *params.force_constants_dihedral;
  set<edge_t> &edge_set = *params.edge_set;

  assert(zero_values_dist.size() == edge_set.size());
  assert(force_constants_dist.size() == edge_set.size());
  assert(force_constants_angle.size() == 2*edge_set.size()); // only for cubic graphs
  assert(force_constants_dihedral.size() == P.points.size());

  vector<coord3d> derivatives(P.N, coord3d(0.0,0.0,0.0));

  set<edge_t>::const_iterator e=edge_set.begin(), ee=edge_set.end();
  for(int i=0; e!=ee; e++, i++){
    const double ax = gsl_vector_get(coordinates,3 * e->first);
    const double ay = gsl_vector_get(coordinates,3 * e->first +1);
    const double az = gsl_vector_get(coordinates,3 * e->first +2);
    const double bx = gsl_vector_get(coordinates,3 * e->second);
    const double by = gsl_vector_get(coordinates,3 * e->second +1);
    const double bz = gsl_vector_get(coordinates,3 * e->second +2);
    // cout << "edge: " << *e << endl;
    // cout << "ax " << ax << " ay " << ay << " az " << az << " bx " << bx << " by " << by << " bz " << bz << endl;
    // cout << "norm(" << i << "): " << coord3d::dnorm(coord3d(ax,ay,az) - coord3d(bx,by,bz)) << endl;
    derivatives[e->first]  += coord3d::dnorm(coord3d(ax,ay,az) - coord3d(bx,by,bz)) * force_constants_dist[i] * (coord3d::dist(coord3d(ax,ay,az), coord3d(bx,by,bz)) - zero_values_dist[*e]);
    derivatives[e->second] -= coord3d::dnorm(coord3d(ax,ay,az) - coord3d(bx,by,bz)) * force_constants_dist[i] * (coord3d::dist(coord3d(ax,ay,az), coord3d(bx,by,bz)) - zero_values_dist[*e]);
    // cout << "dist(" << i << "): " << coord3d::dist(coord3d(ax,ay,az), coord3d(bx,by,bz)) << endl;
  }

  //iterate over all faces
  if(params.optimize_angles){
    for(int i=0;i<P.faces.size();i++){
      const face_t &f(P.faces[i]);
      const int face_size = f.size();
    
      // cout << " face of size-" << it->first << ": " << *jt << endl;
      // iterate over nodes in face
      for (int i=0; i<face_size; ++i){
        int f0 = f[(i-1+face_size)%face_size], f1 = f[i], f2 = f[(i+1)%face_size];
        //        cout << " 3 nodes: " << (*jt)[(i+ it->first -1) % it->first] << ", " << (*jt)[i] <<", " <<  (*jt)[(i+1) % it->first] << endl;
        const double ax = gsl_vector_get(coordinates,   3*f0);
        const double ay = gsl_vector_get(coordinates,   3*f0 + 1);
        const double az = gsl_vector_get(coordinates,   3*f0 + 2);
        const double bx = gsl_vector_get(coordinates,   3*f1   );
        const double by = gsl_vector_get(coordinates,   3*f1 + 1);
        const double bz = gsl_vector_get(coordinates,   3*f1 + 2);
        const double cx = gsl_vector_get(coordinates,   3*f2   );
        const double cy = gsl_vector_get(coordinates,   3*f2 + 1);
        const double cz = gsl_vector_get(coordinates,   3*f2 + 2);
        
        coord3d a(coord3d(ax,ay,az) - coord3d(bx,by,bz)), c(coord3d(cx,cy,cz) - coord3d(bx,by,bz)), da, dc;
        coord3d::dangle(a, c, da, dc);
        
        const double angle_beta = coord3d::angle(coord3d(ax,ay,az) - coord3d(bx,by,bz), coord3d(cx,cy,cz) - coord3d(bx,by,bz));
        //        cout << angle_beta << ", " << M_PI*(1.0-2.0/it->first) << ", " << it->first << endl;
        //        cout << "da: " << da << endl;
        
        derivatives[f0] += da *       (angle_beta - M_PI*(1.0-2.0/face_size)) * force_constants_angle[f1];
        derivatives[f1] += -(da+dc) * (angle_beta - M_PI*(1.0-2.0/face_size)) * force_constants_angle[f1];
        derivatives[f2] += dc *       (angle_beta - M_PI*(1.0-2.0/face_size)) * force_constants_angle[f1];
      }
    }
  }

  // iterate over all dihedrals, i.e., all points
  if(P.is_cubic()){
    for(int u=0; u<P.N; u++)
    {
//        t   B   s
//          \   /
//        C   u   A
//            |
//            r

      node_t r = P.neighbours[u][0];
      node_t s = P.neighbours[u][1];
      node_t t = P.neighbours[u][2];

      const double ax = gsl_vector_get(coordinates, 3*u    );
      const double ay = gsl_vector_get(coordinates, 3*u + 1);
      const double az = gsl_vector_get(coordinates, 3*u + 2);
      const double bx = gsl_vector_get(coordinates, 3*r    );
      const double by = gsl_vector_get(coordinates, 3*r + 1);
      const double bz = gsl_vector_get(coordinates, 3*r + 2);
      const double cx = gsl_vector_get(coordinates, 3*s    );
      const double cy = gsl_vector_get(coordinates, 3*s + 1);
      const double cz = gsl_vector_get(coordinates, 3*s + 2);
      const double dx = gsl_vector_get(coordinates, 3*t    );
      const double dy = gsl_vector_get(coordinates, 3*t + 1);
      const double dz = gsl_vector_get(coordinates, 3*t + 2);

      coord3d b(coord3d(bx,by,bz) - coord3d(ax,ay,az)), c(coord3d(cx,cy,cz) - coord3d(ax,ay,az)), d(coord3d(dx,dy,dz) - coord3d(ax,ay,az)), db, dc, dd;
      const double dihedral_abcd = coord3d::dihedral(b, c, d);
      coord3d::ddihedral(b, c, d, db, dc, dd);
        
      derivatives[u] += -(db+dc+dd) * (dihedral_abcd - zero_values_dihedral[u]) * force_constants_dihedral[u];
      derivatives[r] += db *          (dihedral_abcd - zero_values_dihedral[u]) * force_constants_dihedral[u];
      derivatives[s] += dc *          (dihedral_abcd - zero_values_dihedral[u]) * force_constants_dihedral[u];
      derivatives[t] += dd *          (dihedral_abcd - zero_values_dihedral[u]) * force_constants_dihedral[u];

//      cout << "derivarive: " << dihedral_abcd << ", " << zero_values_dihedral[u] << ", " << u << endl;
    }
  }

//   // NB: Try out, then remove or rewrite
//   for(int i=0;i<P.points.size();i++){
//     const coord3d x(gsl_vector_get(coordinates,3*i),gsl_vector_get(coordinates,3*i+1),gsl_vector_get(coordinates,3*i+2));
//     const double norm2 = x.dot(x);
//     const double dcoul = -2000.0*2/(norm2*norm2);
//     for(int j=0;j<3;j++) derivatives[i][j] += x[j]*dcoul;
//   }

  // return gradient
  for(int i = 0; i < P.N; ++i) {
    gsl_vector_set(gradient, 3*i, derivatives[i][0]);
    gsl_vector_set(gradient, 3*i+1, derivatives[i][1]);
    gsl_vector_set(gradient, 3*i+2, derivatives[i][2]);
  }

//   double grad_debug = 0.0;
//   for(int i = 0; i < P.N; ++i) {
//     grad_debug += sqrt(pow(gsl_vector_get(gradient, 3*i),2) + pow(gsl_vector_get(gradient, 3*i+1),2) + pow(gsl_vector_get(gradient, 3*i+2),2));
//   }
//   cout << "grad: " << grad_debug << endl;

}


void polyhedron_pot_grad(const gsl_vector* coordinates, void* parameters, double* potential, gsl_vector* gradient)
{
  *potential = polyhedron_pot(coordinates, parameters);
  polyhedron_grad(coordinates, parameters, gradient);
}


bool Polyhedron::optimize_other(bool optimize_angles, map<edge_t, double> zero_values_dist)
{
  //cout << "entering opt other" << endl;

  assert(points.size() == N);

  const double eps=1e-4;

  // settings for the optimizations
  const double stepsize = 1e-3;// FIXME final value
  const double terminate_gradient = 1e-10;// FIXME final value
  const double tol = 1e-1; // accuracy of line minimization, the manual suggests 0.1
  const int max_iterations = 10000;// FIXME final value
  
  // init values
  set<edge_t> edge_set = undirected_edges();

  const double default_edge_length = 1.4;
  if(zero_values_dist.size() != edge_set.size()){
    for(set<edge_t>::iterator it=edge_set.begin(), to=edge_set.end(); it!=to; it++){
      zero_values_dist.insert(make_pair(*it, default_edge_length));
    }
  }

  vector<double> zero_values_dihedral(N);
  if(is_cubic())
  {
    orient_neighbours(); // CCW
    const int fmax = 10; // maximum face size
//          t   B   s
//            \   /
//          C   u   A
//              |
//              r
    for(node_t u=0; u<N; ++u){
      node_t r = neighbours[u][0];
      node_t s = neighbours[u][1];
      node_t t = neighbours[u][2];
      // face sizes
      int lA = shortest_cycle(r,u,s,fmax).size();
      int lB = shortest_cycle(s,u,t,fmax).size();
      int lC = shortest_cycle(t,u,r,fmax).size();
      // bond lengths
      double ur = zero_values_dist.find(edge_t(u,r))->second;
      double us = zero_values_dist.find(edge_t(u,s))->second;
      double ut = zero_values_dist.find(edge_t(u,t))->second;

      // we may have to rearrange the neighbour list to preserve symmetry
      auto right_shift = [&](){
        neighbours[u][0] = t; neighbours[u][1] = r; neighbours[u][2] = s;
        int tmp_f=lA; lA=lC; lC=lB; lB=tmp_f;
        double tmp_l=ur; ur=ut; ut=us; us=tmp_l;
      };
      auto left_shift = [&](){
        neighbours[u][0] = s; neighbours[u][1] = t; neighbours[u][2] = r;
        int tmp_f=lA; lA=lB; lB=lC; lC=tmp_f;
        double tmp_l=ur; ur=us; us=ut; ut=tmp_l;
      };
      // check for face sizes
      if(lA==lB && lA!=lC){ // AAB --> ABA
        left_shift();
      }
      else if(lA!=lB && lB==lC){ // BAA --> ABA
        right_shift();
      }
      else if(lA!=lB && lB!=lC && lA!=lC){ // rotate the smallest to the front
        if(lB<lA && lB<lC){
          left_shift();
        }
        else if(lC<lA && lC<lA){
          right_shift();
        }
      }
      // and if the face sizes are all the same, check for edge length
      else if(lA==lB && lA==lC){
        if(ur-us<eps && abs(ur-ut)>eps){ // aab --> baa
          right_shift();
        }
        else if(abs(ur-us)>eps && us-ut<eps){ // aba --> baa
          left_shift();
        }
        else if(abs(ur-us)>eps && abs(us-ut)>eps && abs(ur-ut)>eps){ // rotate the shortest to the front
          if(us<ur && us<ut){
            left_shift();
          }
          else if(ut<ur && ut<ur){
            right_shift();
          }
        }
      }

      zero_values_dihedral[u] = coord3d::ideal_dihedral(lA, lB, lC, ur, us, ut);
//      cout << "opt-main, r,s,t, a,lb,lc, th_0: " <<r<<" " <<s<<" " <<t<<" " << lA<<" " <<lB<<" "<<lC<<" "<<zero_values_dihedral[u] << endl;
    }
  }
  vector<double> force_constants_dist(edge_set.size(), 500.0);
  vector<double> force_constants_angle(2 * edge_set.size(), 200.0); // FIXME possibly change later // only for cubic graphs
  vector<double> force_constants_dihedral(N, 50.0);

  params_t params;
  params.P = this; 
  params.zero_values_dist = &zero_values_dist;
  params.zero_values_dihedral = &zero_values_dihedral;
  params.force_constants_dist = &force_constants_dist;
  params.force_constants_angle = &force_constants_angle;
  params.force_constants_dihedral = &force_constants_dihedral;
  params.edge_set = &edge_set;
  params.optimize_angles = optimize_angles;
  
  // Create the minimisation function block, define the different functions and parameters
  gsl_multimin_function_fdf potential_function;
  potential_function.n = 3*N;
  potential_function.f = &polyhedron_pot;
  potential_function.df = &polyhedron_grad;
  potential_function.fdf = &polyhedron_pot_grad;
  potential_function.params = static_cast<void*>(&params);
  
  // our zero-th order geometry sometimes places two vertices at exactly the same coordinates
  // solution:  displace randomly and let the FF do the rest
  set<edge_t>::const_iterator e=edge_set.begin(), ee=edge_set.end();
  for(int i=0; e!=ee; e++, i++){
    const double ax = points[e->first][0];
    const double ay = points[e->first][1];
    const double az = points[e->first][2];
    const double bx = points[e->second][0];
    const double by = points[e->second][1];
    const double bz = points[e->second][2];
    // cout << "edge: " << *e << endl;
    // cout << "ax " << ax << " ay " << ay << " az " << az << " bx " << bx << " by " << by << " bz " << bz << endl;
    // cout << "dist(" << i << "): " << coord3d::dist(coord3d(ax,ay,az), coord3d(bx,by,bz)) << endl;
    if(coord3d::dist(coord3d(ax,ay,az), coord3d(bx,by,bz)) < 1e-1 ){
      // cout << "short bond found ... displacing." << endl;
      const double displacement = 0.5;
      points[e->first][0]  += displacement;
      points[e->first][1]  += displacement;
      points[e->first][2]  += displacement;
      points[e->second][0] -= displacement;
      points[e->second][1] -= displacement;
      points[e->second][2] -= displacement;
    }
  }

  // Starting point
  gsl_vector* coordinates = gsl_vector_alloc(potential_function.n);

  for(int i=0; i<N; ++i){
    gsl_vector_set(coordinates, 3*i, points[i][0]);
    gsl_vector_set(coordinates, 3*i+1, points[i][1]);
    gsl_vector_set(coordinates, 3*i+2, points[i][2]);
  }

  const gsl_multimin_fdfminimizer_type *fT = gsl_multimin_fdfminimizer_conjugate_pr;
  gsl_multimin_fdfminimizer *s = gsl_multimin_fdfminimizer_alloc(fT, potential_function.n);
  gsl_multimin_fdfminimizer_set(s, &potential_function, coordinates, stepsize, tol); // runs fdf once
  size_t iter = 0;
  int status;
  do {
    ++iter;
    status = gsl_multimin_fdfminimizer_iterate(s);
    
    if(iter % 1000 == 0){
      cout << "it: " << iter << " stat: " << status << endl;
      printf ("error: %s\n", gsl_strerror (status));
    }
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
bool Polyhedron::optimize_other(bool, map<edge_t, double>)
{
  cerr << "Optimizing other polyhedra than fullerenes is only available through GSL." << endl;
  return 0;
}
#endif

