#include <chrono>
#include <assert.h>
#include <inttypes.h>
#include <stdio.h>
#include <math.h>
#include <iostream>
#include <vector>
using namespace std;
using namespace std::literals;

#include <array>
#include <tuple>


#include <fullerenes/spiral.hh>
#include <fullerenes/triangulation.hh>
#include <fullerenes/fullerenegraph.hh>
#include <fullerenes/polyhedron.hh>

#include "coord3d.cc"

// TODO: Place somewhere sensible
template <typename T, size_t N> 
inline T sum(const array<T,N>& a){
  T result = 0;
  for(size_t i=0;i<N;i++) result += a[i];
  return result;
}

inline real_t array_dot(const vector<Coord3d>& a, const vector<Coord3d>& b){
  real_t result = 0;
  for(size_t i=0;i<a.size();i++) result += dot(a[i],b[i]);
  return result;
}


//Parallelizable euclidean norm function.
inline real_t norm(const vector<Coord3d>& array){
  return sqrt(array_dot(array,array));
}

class FullereneForcefield
{
public:
  const int N;
  const node_t *neighbours; // Sparse fullerene bond graph Nx3 array of neighbours b,c,d
  Coord3d *X;               // Current nucleus positions, N x Coord3d (Nx3 reals)
  
  const uint8_t *face_right;                   // Face-size to the right of the edges ab, ac, ad; Nx3 array [TODO: 1 bit/arc]
  const node_t  *next_on_face,  *prev_on_face; // Next node on the face to the left of the edges ab, ac, ad; Nx3 array
  
  // TODO: Separate data from methods
  //All parameter arrays are indexed by the total number of hexagons
  //Pentagon = 0
  //Hexagon  = 1
  //PPP = 0, {HPP, PHP, PPH} = 1, {PHH, HPH, HHP} = 2, {HHH} = 3
  const array<real_t, 2> optimal_corner_cos_angles = {cos(M_PI * 108 / 180), cos(M_PI * 120 / 180)}; 
  const array<real_t, 3> optimal_bond_lengths      = {1.479, 1.458, 1.401}; 
  const array<real_t, 4> optimal_dih_cos_angles    = {cos(0.652358), cos(0.509674), cos(0.417884), cos(0)}; 
  
  const array<real_t, 2> angle_forces = {100.0, 100.0}; 
  const array<real_t, 3> bond_forces  = {260.0, 390.0, 450.0}; 
  const array<real_t, 4> dih_forces   = {35.0, 65.0, 85.0, 270.0}; 

  struct ArcData
  {
    static inline real_t harmonic_energy(const real_t p0, const real_t p){
      return 0.5*(p-p0)*(p-p0);
    }
    static inline Coord3d  harmonic_energy_gradient(const real_t p0, const real_t p, const Coord3d gradp){
      return (p-p0)*gradp;     
    }
    //Helper Function for dihedral calculation.
    inline tuple<Coord3d, real_t, real_t> dihedral_terms(const Coord3d& ba_hat, const Coord3d& bc_hat) const{
      real_t cos_b   = dot(ba_hat,bc_hat);
      real_t r_sin_b = 1.0/sqrt(1 - cos_b*cos_b);
      Coord3d nabc = cross(ba_hat, bc_hat) * r_sin_b;
      return {nabc,cos_b,r_sin_b};
    }

    inline real_t bond_length() const {return 1/r_rab;}
    inline real_t angle()       const {return dot(ab_hat,ac_hat);}
    //Returns the inner dihedral angle for the current arc. Used here only for energy calculation, 
    //otherwise embedded in dihedral computation because the planes and angles that make up the dihedral angle computation are required for derivative computation.
    inline real_t dihedral()    const 
    { 
      Coord3d nabc, nbcd; real_t cos_b, cos_c, r_sin_b, r_sin_c;
      tie(nabc, cos_b, r_sin_b) = dihedral_terms(ba_hat, bc_hat);
      tie(nbcd, cos_c, r_sin_c) = dihedral_terms(-bc_hat, cd_hat);

      return dot(nabc, nbcd);
    }

    // Chain rule terms for angle calculation
    //Computes gradient related to bending term. ~24 FLOPs
    inline Coord3d inner_angle_gradient() const
    {
      real_t cos_angle = angle(); //Inner angle of arcs ab,ac.
      Coord3d grad = cos_angle * (ab_hat * r_rab + ac_hat * r_rac) - ab_hat * r_rac - ac_hat* r_rab; //Derivative of inner angle: Eq. 21. 
      return f_inner_angle * harmonic_energy_gradient(angle0, cos_angle, grad); //Harmonic Energy Gradient: Eq. 21. multiplied by harmonic term.
    }
    //Computes gradient related to bending of outer angles. ~20 FLOPs
    inline Coord3d outer_angle_gradient_m() const
    {
      real_t cos_angle = -dot(ab_hat, bm_hat); //Compute outer angle. ab,bm
      Coord3d grad = (bm_hat + ab_hat * cos_angle) * r_rab; //Derivative of outer angles Eq. 30. Buster Thesis
      return f_outer_angle_m * harmonic_energy_gradient(outer_angle_m0,cos_angle,grad); //Harmonic Energy Gradient: Eq. 30 multiplied by harmonic term.
    }
    inline Coord3d outer_angle_gradient_p() const
    {
      real_t cos_angle = -dot(ab_hat, bp_hat);              //Compute outer angle. ab,bp
      Coord3d grad = (bp_hat + ab_hat * cos_angle) * r_rab; //Derivative of outer angles Eq. 28. Buster Thesis
      return f_outer_angle_p * harmonic_energy_gradient(outer_angle_p0,cos_angle,grad); //Harmonic Energy Gradient: Eq. 28 multiplied by harmonic term.
    }
    // Chain rule terms for dihedral calculation
    //Computes gradient related to dihedral/out-of-plane term. ~75 FLOPs
    inline Coord3d inner_dihedral_gradient() const
    {
      Coord3d nabc, nbcd; real_t cos_b, cos_c, r_sin_b, r_sin_c;
      tie(nabc, cos_b, r_sin_b) = dihedral_terms(ba_hat, bc_hat);
      tie(nbcd, cos_c, r_sin_c) = dihedral_terms(-bc_hat, cd_hat);

      real_t cos_beta = dot(nabc, nbcd); //Inner dihedral angle from planes abc,bcd.
            
      real_t cot_b = cos_b * r_sin_b * r_sin_b; //cos(b)/sin(b)^2

      //Derivative w.r.t. inner dihedral angle F and G in Eq. 26
      Coord3d grad = cross(bc_hat, nbcd) * r_sin_b * r_rab - ba_hat * cos_beta * r_rab + (cot_b * cos_beta * r_rab) * (bc_hat - ba_hat * cos_b);

      return f_inner_dihedral * harmonic_energy_gradient(inner_dih0, cos_beta, grad); //Eq. 26.
    }

    //Computes gradient from dihedral angles constituted by the planes bam, amp ~162 FLOPs
    inline Coord3d outer_a_dihedral_gradient() const
    {
      Coord3d nbam_hat, namp_hat; real_t cos_a, cos_m, r_sin_a, r_sin_m;
      tie(nbam_hat, cos_a, r_sin_a) = dihedral_terms(ab_hat, am_hat);
      tie(namp_hat, cos_m, r_sin_m) = dihedral_terms(-am_hat, mp_hat);
            
      real_t cos_beta = dot(nbam_hat, namp_hat); //Outer Dihedral angle bam, amp
      real_t cot_a = cos_a * r_sin_a * r_sin_a;
      real_t cot_m = cos_m * r_sin_m * r_sin_m;

      //Derivative w.r.t. outer dihedral angle, factorized version of Eq. 31.
      Coord3d grad = cross(mp_hat,nbam_hat)*r_ram*r_sin_m - (cross(namp_hat,ab_hat)*r_ram + cross(am_hat,namp_hat)*r_rab)*r_sin_a +
	cos_beta*(ab_hat*r_rab + r_ram * (2*am_hat + cot_m*(mp_hat+cos_m*am_hat)) - cot_a*(r_ram*(ab_hat - am_hat*cos_a) + r_rab*(am_hat-ab_hat*cos_a)));
            
      //Eq. 31 multiplied by harmonic term.
      return f_outer_dihedral * harmonic_energy_gradient(outer_dih0, cos_beta, grad);
    }

    //Computes gradient from dihedral angles constituted by the planes nbmp, nmpa ~92 FLOPs
    inline Coord3d outer_m_dihedral_gradient() const
    {
      Coord3d nbmp_hat, nmpa_hat; real_t cos_m, cos_p, r_sin_m, r_sin_p;
      tie(nbmp_hat, cos_m, r_sin_m) = dihedral_terms(mb_hat, mp_hat);
      tie(nmpa_hat, cos_p, r_sin_p) = dihedral_terms(-mp_hat, pa_hat);
            
      //Cosine to the outer dihedral angle constituted by the planes bmp and mpa
      real_t cos_beta = dot(nbmp_hat, nmpa_hat); //Outer dihedral angle bmp,mpa.
      real_t cot_p = cos_p * r_sin_p * r_sin_p;
            
      //Derivative w.r.t. outer dihedral angle, factorized version of Eq. 32.
      Coord3d grad = r_rap * (cot_p*cos_beta * (-mp_hat - pa_hat*cos_p) - cross(nbmp_hat, mp_hat)*r_sin_p - pa_hat*cos_beta );

      //Eq. 32 multiplied by harmonic term.
      return f_outer_dihedral * harmonic_energy_gradient(outer_dih0, cos_beta, grad);
    }

    //Computes gradient from dihedral angles constituted by the planes bpa, pam ~162 FLOPs
    inline Coord3d outer_p_dihedral_gradient() const
    {
      Coord3d nbpa_hat, npam_hat; real_t cos_p, cos_a, r_sin_p, r_sin_a;
      tie(npam_hat, cos_a, r_sin_a) = dihedral_terms(ap_hat, am_hat);
      tie(nbpa_hat, cos_p, r_sin_p) = dihedral_terms(pb_hat, -ap_hat);
      real_t cos_beta = dot(nbpa_hat, npam_hat); //Outer dihedral angle bpa, pam.
      real_t cot_p = cos_p * r_sin_p * r_sin_p;
      real_t cot_a = cos_a * r_sin_a * r_sin_a;

      //Derivative w.r.t. outer dihedral angle, factorized version of Eq. 33.
      Coord3d grad = cross(npam_hat,pb_hat)*r_rap*r_sin_p - (cross(am_hat,nbpa_hat)*r_rap + cross(nbpa_hat,ap_hat)*r_ram)*r_sin_a +
	cos_beta*(am_hat*r_ram + r_rap * (2*ap_hat + cot_p*(pb_hat+cos_p*ap_hat)) - cot_a*(r_rap*(am_hat - ap_hat*cos_a) + r_ram*(ap_hat-am_hat*cos_a)));
            
      //Eq. 33 multiplied by harmonic term.
      return f_outer_dihedral * harmonic_energy_gradient(outer_dih0, cos_beta, grad);
    }

    // Internal coordinate gradients
    // [JA: TODO: Nu er det ikke længere gradienter af de interne koordinater: bond length, angle, osv., men af de harmoniske energi-udtryk => skal bruge nye navne]
    inline Coord3d bond_length_gradient() const { return -f_bond * harmonic_energy_gradient(r0,bond_length(),ab_hat);}
    //Sum of angular gradient components.
    inline Coord3d angle_gradient()       const { return inner_angle_gradient() + outer_angle_gradient_p() + outer_angle_gradient_m();}
    //Sum of inner and outer dihedral gradient components.
    inline Coord3d dihedral_gradient()    const { return inner_dihedral_gradient() + outer_a_dihedral_gradient() + outer_m_dihedral_gradient() + outer_p_dihedral_gradient();}
    //inline Coord3d flatness()             const { return ;  }   
        

    //Harmonic energy contribution from bond stretching, angular bending and dihedral angle bending.
    inline real_t energy() const {return 0.5*f_bond *harmonic_energy(bond_length(),r0)+f_inner_angle* harmonic_energy(angle(),angle0)+f_inner_dihedral* harmonic_energy(dihedral(),inner_dih0);}
    //Sum of bond, angular and dihedral gradient components.
    inline Coord3d gradient() const{return bond_length_gradient() + angle_gradient() + dihedral_gradient();}


    //Force constants for all paremeters.
    real_t  
    f_outer_dihedral,
      f_inner_dihedral,
      f_inner_angle,
      f_outer_angle_p,
      f_outer_angle_m,
      f_bond;
        
    //Residual lengths of arcs ab, ac, am, ap.
    real_t
    r_rab,
      r_rac,
      r_ram,
      r_rap;

    //Equillibrium parameters.
    real_t
    r0,
      angle0,
      outer_angle_m0,
      outer_angle_p0,
      inner_dih0,
      outer_dih0;

    /*
      All normalized arcs required to perform energy & gradient calculations.
      Note that all these arcs are cyclical the arc ab becomes: ab->ac->ad,  the arc ac becomes: ac->ad->ab , the arc bc becomes: bc->cd->db (For iterations 0, 1, 2)
      As such the naming convention here is related to the arcs as they are used in the 0th iteration. */
    Coord3d 
      ab_hat,			// Direction vector of current arc
      ac_hat,			// Direction vector of next arc around node a
      ad_hat,			// Direction vector of previous arc around node a
      am_hat,			// See Figure...
      ap_hat,			// See Figure...
      ba_hat,			// No need to store this
      bc_hat,			// Direction 
      bp_hat,			// 
      bm_hat,			// Many of these should probably be calculated on the fly...
      cd_hat,			// 
      mp_hat,
      mb_hat,
      pa_hat,
      pb_hat;
  };


  FullereneForcefield(size_t N,
		      const node_t *neighbours,
		      Coord3d *X,
		      const uint8_t *face_right,
		      const node_t  *next_on_face,
		      const node_t  *prev_on_face) : 
    N(N), neighbours(neighbours), X(X), face_right(face_right), next_on_face(next_on_face), prev_on_face(prev_on_face) {}


  //~393*N FLOPs
  real_t energy(const Coord3d *X){
    real_t total_energy = 0;
    
    for (node_t a = 0; a < N; a++)
      {   
	real_t node_energy = 0.0;

	//Fetching neighbour indices.
	array<node_t, 3> na = {neighbours[a * 3], neighbours[a * 3 + 1], neighbours[a * 3 + 2]}; // Neighbours b,c,d to a

	//Fetching current node coordinate.
	Coord3d Xa = X[a];

	//Fetching face information to the left and right of the edges: ab, ac, ad.
	array<uint8_t, 3> f_r = {face_right[a * 3], face_right[a * 3 + 1], face_right[a * 3 + 2]};
	array<uint8_t, 3> f_l = {face_right[a * 3 + 2], face_right[a * 3], face_right[a * 3 + 1]};

	//Fetching coordinates for neighbour nodes b,c,d.
	array<Coord3d, 3> Xbcd = {X[na[0]], X[na[1]], X[na[2]]};
	    
	array<real_t,3> r_rab;
	array<Coord3d,3> abs;
	array<Coord3d,3> ab_hats;

	for (size_t i = 0; i < 3; i++)// 30 FLOPs
	  {
	    tie(r_rab[i],abs[i],ab_hats[i]) = split_norm(Xbcd[i] - Xa); //10 FLOPs
	  }
	for (size_t j = 0; j < 3; j++)// 121*3 = 363 FLOPs
	  {   
	    ArcData arc;
	    arc.r0 = optimal_bond_lengths[ f_l[j] + f_r[j] ];
	    arc.angle0 = optimal_corner_cos_angles[ f_r[j] ];
	    arc.inner_dih0 = optimal_dih_cos_angles[ f_l[0] + f_l[1] + f_l[2]];
	    arc.f_bond = bond_forces[ f_l[j] + f_r[j] ];
	    arc.f_inner_angle = angle_forces[ f_l[j] ];
	    arc.f_inner_dihedral = dih_forces[ f_l[0] + f_l[1] + f_l[2] ];
                
	    arc.ab_hat = ab_hats[j];
	    arc.ba_hat = -ab_hats[j];
	    arc.ac_hat = ab_hats[(j+1)%3];
	    arc.bc_hat = unit_vector(abs[(j+1)%3] - abs[j]);
	    arc.cd_hat = unit_vector(abs[(j+2)%3] - abs[(j+1)%3]);
	    arc.r_rab = r_rab[j];
                
	    node_energy += arc.energy();
	  }
	//All bond energies are computed twice, (node a computes ab and node b computes ba) , hence the factor 0.5;
	total_energy += node_energy;
      }

    return total_energy;
  }

  //~1914 * N -FLOPs
  void gradient(const Coord3d *X, Coord3d *gradient)
  {   
    for (node_t a = 0; a < N; a++) 
      {   
	//Fetching neighbour indices.
	array<node_t, 3> na = {neighbours[a * 3], neighbours[a * 3 + 1], neighbours[a * 3 + 2]}; // Neighbours b,c,d to a

	//Fetching current node coordinate.
	Coord3d Xa = X[a];


	//Fetching face information to the left and right of the edges: ab, ac, ad.
	array<uint8_t, 3> f_r = {face_right[a * 3] , face_right[a * 3 + 1] , face_right[a * 3 + 2] };
	array<uint8_t, 3> f_l = {face_right[a * 3 + 2] , face_right[a * 3] , face_right[a * 3 + 1] };

	//Fetching face information to the right of the edges {ba, bm, bp}, {ca, cm, cp}, {da, dm, dp}
	array<array<uint8_t,3>, 3> f_r_bcd = {face_right[na[0] * 3], face_right[na[0] * 3 + 1], face_right[na[0] * 3 + 2], 
					      face_right[na[1] * 3], face_right[na[1] * 3 + 1], face_right[na[1] * 3 + 2],
					      face_right[na[2] * 3], face_right[na[2] * 3 + 1], face_right[na[2] * 3 + 2]};
           
	//Pre calculate the sum of the face information for each neighbour. Used to index dihedral constants.
	array<uint8_t,3> dihedral_face_sum_bcd = {sum(f_r_bcd[0]), sum(f_r_bcd[1]), sum(f_r_bcd[2])};

	//Fetching coordinates for neighbour nodes b,c,d.
	array<Coord3d, 3> Xbcd = {X[na[0]], X[na[1]], X[na[2]]};

	//Fetching coordinates for outer nodes.
	array<Coord3d, 3> Xbp = {X[next_on_face[a*3]], X[next_on_face[a*3 + 1]], X[next_on_face[a*3 + 2]]};
	array<Coord3d, 3> Xbm = {X[prev_on_face[a*3]], X[prev_on_face[a*3 + 1]], X[prev_on_face[a*3 + 2]]};

	//Vector reciprocal norms.
	array<real_t, 3> r_rab, r_rbp, r_rbm;

	//Common vectors.
	array<Coord3d, 3> abs, bps, bms, ab_hats, bp_hats, bm_hats;

	Coord3d am, mp, ap, am_hat, mp_hat, ap_hat;

	real_t r_ram, r_rmp, r_rap;

	Coord3d node_gradient = {0,0,0};
            
	for (size_t i = 0; i < 3; i++) // 90 FLOPs
	  {
	    tie(r_rab[i], abs[i], ab_hats[i]) = split_norm(Xbcd[i] - Xa); // 10 FLOPs
	    tie(r_rbm[i], bms[i], bm_hats[i]) = split_norm(Xbm[i] - Xbcd[i]);
	    tie(r_rbp[i], bps[i], bp_hats[i]) = split_norm(Xbp[i] - Xbcd[i]);
	  }
	for (size_t j = 0; j < 3; j++) // 608*3 = 1824 FLOPs
	  {       
	    ArcData arc;
	    //Indexing of 'Equillibrium' parameters and related force constants, through face information.
	    arc.r0 = optimal_bond_lengths[ f_l[j] + f_r[j] ];
	    arc.angle0 = optimal_corner_cos_angles[ f_r[j] ];
	    arc.inner_dih0 = optimal_dih_cos_angles[ f_l[0] + f_l[1] + f_l[2]];
	    arc.outer_angle_m0 = optimal_corner_cos_angles[ f_l[j] ];
	    arc.outer_angle_p0 = optimal_corner_cos_angles[ f_r[j] ];
	    arc.outer_dih0 = optimal_dih_cos_angles[ dihedral_face_sum_bcd[j] ];

	    arc.f_bond = bond_forces[ f_l[j] + f_r[j] ];
	    arc.f_inner_angle = angle_forces[ f_l[j] ];
	    arc.f_inner_dihedral = dih_forces[ f_l[0] + f_l[1] + f_l[2] ];
	    arc.f_outer_angle_m = angle_forces[ f_r[j] ];
	    arc.f_outer_angle_p = angle_forces[ f_l[j] ];
	    arc.f_outer_dihedral = dih_forces[ dihedral_face_sum_bcd[j] ];

	    //Compute relevant dihedral vectors.
	    tie(r_ram, am, am_hat) = split_norm(bms[j] + abs[j]); //10 FLOPs
	    tie(r_rmp, mp, mp_hat) = split_norm(bps[j] - bms[j]); //10 FLOPs
	    tie(r_rap, ap, ap_hat) = split_norm(bps[j] + abs[j]); //10 FLOPs

	    arc.ab_hat = ab_hats[j];
	    arc.ac_hat = ab_hats[(j+1)%3];
	    arc.ad_hat = ab_hats[(j+2)%3];
	    arc.am_hat = am_hat;
	    arc.ap_hat = ap_hat;
	    arc.ba_hat = -ab_hats[j];
	    arc.bc_hat = unit_vector(abs[(j+1)%3] - abs[j]);
	    arc.bp_hat = bp_hats[j];
	    arc.bm_hat = bm_hats[j];
	    arc.cd_hat = unit_vector(abs[(j+2)%3] - abs[(j+1)%3]);
	    arc.mp_hat = mp_hat;
	    arc.mb_hat = -bm_hats[j];
	    arc.pa_hat = -ap_hat;
	    arc.pb_hat = -bp_hats[j];

	    arc.r_rab = r_rab[j];
	    arc.r_rac = r_rab[(j+1)%3];
	    arc.r_ram = r_ram;
	    arc.r_rap = r_rap;
                
	    node_gradient = node_gradient + arc.gradient(); 

	  }

	gradient[a] = node_gradient;
      }
  }

  size_t golden_section_search(Coord3d *X, const Coord3d *direction,
			       Coord3d *new_direction,
			       real_t a, real_t b, real_t tol)
  {
    real_t tau = (sqrt(5) - 1) / 2;
        
    //Actual coordinates resulting from each traversal 
    vector<Coord3d> X1(N), X2(N);
    
    //Line search x - values;
    real_t x1,  x2;
    x1 = (a + (1 - tau) * (b - a));
    x2 = (a + tau * (b - a));
    for (node_t a = 0; a < N; a++)
      {
	X1[a] = X[a] + x1 * direction[a];
	X2[a] = X[a] + x2 * direction[a];
      }
    real_t f1 = energy(&X1[0]);
    real_t f2 = energy(&X2[0]);
    size_t evals = 2;

    while (abs(b - a) > tol){
      evals++;
      if (f1 > f2){
	a = x1;
	x1 = x2;
	f1 = f2;
	x2 = a + tau * (b - a);
	for (node_t a = 0; a < N; a++)
	  {
	    X2[a] = X[a] + x2 * direction[a];
	  }
	f2 = energy(&X2[0]);
      }else
	{
	  b = x2;
	  x2 = x1;
	  f2 = f1;
	  x1 = a + (1 - tau) * (b - a);
	  for (node_t a = 0; a < N; a++){
	    X1[a] = X[a] + x1 * direction[a];
	  }
	  f1 = energy(&X1[0]);
	}
    }
    //Line search coefficient
    real_t alfa = (a+b)/2;

    // TODO: Why is the input being updated inside a line-search? Should be const.
    for (node_t a = 0; a < N; a++)
      {
	X[a] = X[a] + alfa*direction[a];
      }
    gradient(X,new_direction);
    for (node_t a = 0; a < N; a++){ new_direction[a] = -new_direction[a];}
    return evals;
  }

  // size_t bisection_search(const Coord3d *X, const Coord3d *direction,
  // 			  Coord3d *fc, real_t a, real_t b, real_t eps,
  // 			  size_t max_iter)
  // {
  //   size_t count = 0;
  //   size_t nBrack = 0;
  //   real_t sum = -1;
  //   real_t dfc = 0.0;
  //   real_t coefficient = 0.0;
  //   size_t gradient_evals = 0;

  //   vector<Coord3d>  temp_coords(N);
  //   vector<Coord3d> temp_direction(N);

  //   //Copy coordinates into temporary coordinates. TODO: Is this necessary?
  //   for(node_t a = 0; a < N; a++)
  //     temp_coords[a]= X[a];

  //   for (node_t a = 0; a < N; a++)
  //     {
  // 	dfc -= dot(direction[a],direction[a]);
  //     }

  //   // TODO: Bedre navn end 'sum' (hvad er det, der tjekkes her?)
  //   while (sum < 0)
  //     {   
  // 	if (nBrack == 0)
  // 	  {
  // 	    coefficient = 0.0;
  // 	  } else
  // 	  {
  // 	    coefficient = b;
  // 	    b *= 1.5;
  // 	  }
  // 	coefficient = b - coefficient;
            
           
  // 	nBrack ++;
  // 	for (node_t a = 0; a < N; a++)
  // 	  {
  // 	    temp_coords[a] = temp_coords[a] + coefficient*direction[a];
  // 	  }
  // 	gradient(&temp_coords[0],&temp_direction[0]);
  // 	gradient_evals++;
  // 	//	sum = multiply_reduce(temp_direction,direction);
	
  //     }
        
  //   while (abs(dfc) > eps)
  //     {
  // 	//Copy temporary coordinates back to coordinates. TODO: Is this necessary?
  // 	for(node_t a = 0; a < N; a++)
  // 	  temp_coords[a]= X[a];
	
  // 	count++;
  // 	real_t c = (a+b)/2;

  // 	for (node_t i = 0; i < N; i++)
  // 	  {   
  // 	    temp_coords[i] = temp_coords[i] + c*direction[i];
  // 	  }
  // 	gradient(temp_coords,fc);
  // 	gradient_evals++;
  // 	// multiply_reduce(fc,direction);
  // 	dfc = 0;
  // 	for(node_t a=0;a<N;a++)
  // 	  dfc += dot(fc[a],direction[a]);
            
  // 	if (count > max_iter){ break;}
  // 	if (dfc < 0){a = c;} else{b = c;}
  //     }
  //   parallel_copy(X, temp_coords);
  //   for (node_t i = 0; i < N; i++)
  //     {
  // 	fc[i]= -fc[i];
  //     }
  //   return gradient_evals;
  // }


  void conjugate_gradient(){
    size_t iter_count = 0;
    size_t max_iter = N*10;
    real_t beta = 0.0;
    real_t dnorm = 1.0;
    size_t gradient_evals = 0;
    size_t energy_evals = 0;

    vector<Coord3d> X_temp(N);	// TODO: Nødvendig?
    vector<Coord3d> delta_x0(N);
    vector<Coord3d> delta_x1(N);
    vector<Coord3d> direction(N);

    gradient(X, &direction[0]);
    gradient_evals ++;
    dnorm = norm(direction);
    for(node_t a = 0; a < N; a++)
      direction[a] = -direction[a]/dnorm;

    //    parallel_copy(X_temp,X);
    //    parallel_copy(delta_x0, direction);
    for(node_t a=0; a<N ; a++){
      X_temp[a]   = X[a];
      delta_x0[a] = direction[a];
    }
    
    while (dnorm > 1e-5){   
	beta = 0.0;
	energy_evals += golden_section_search(&X_temp[0], &(direction[0]),
					      &(delta_x1[0]), 0, 1, 1e-10);
	gradient_evals++;
	//gradient_evals += bisection_search(X_temp, direction, delta_x1,0,1e-5,1e-10,N);
	//Polak Ribiere method
	for (node_t a = 0; a < N; a++)
	  {
	    beta += dot(delta_x1[a], (delta_x1[a] - delta_x0[a]));
	  }
	beta /= array_dot(delta_x0,delta_x0);
        
	if (energy(&(X_temp[0])) > energy(X))
	  {
	    //parallel_copy(X_temp, X);
	    //parallel_copy(delta_x1, delta_x0);
	    for(node_t a=0;a<N;a++){
	      X_temp[a]   = X[a];
	      delta_x1[a] = delta_x0[a];
	    }
	    beta = 0.0;
	  }
	else
	  {   
	    //parallel_copy(X, X_temp);
	    //parallel_copy(delta_x0,delta_x1);
	    for(node_t a=0;a<N;a++){
	      X[a]        = X_temp[a];
	      delta_x0[a] = delta_x1[a];
	    }
	  }
	for(node_t a = 0; a < N; a++)
	    direction[a] = delta_x1[a] + beta*direction[a];

	//	normalize(direction);
	//	dnorm = norm(delta_x1);
	//print_real(dnorm);
	real_t direction_norm2 = 0, dnorm2 = 0;
	for(node_t a=0; a<N; a++){
	  direction_norm2 += dot(direction[a],direction[a]);
	  dnorm2          += dot(delta_x1[a],delta_x1[a]);
	}
	
	if (iter_count > max_iter)
	  {
	    cout << "Conjugated Gradient Terminated Due to Max Iterations :" << N*10 << "\n";
	    cout << "Gradient Evaluations: " << gradient_evals << "\n";
	    cout << "Energy Evaluations: " << energy_evals << "\n";
	    return;
	  }
	iter_count++;
      }
    printf("Conjugated Gradient Finished in %ld iterations\n"
	   "Gradient evaluations: %ld\n"
	   "Energy   evaluations: %ld\n", iter_count, gradient_evals, energy_evals);
  }
};


int face_size(const Graph &g, node_t u, node_t v)
{
  int d = 1;
  node_t u0 = u;
  while(v != u0){
    node_t w = v;
    v = g.next_on_face(u,v);
    u = w;
    d++;
  }
  return d;
}


FullereneForcefield new_fullereneforcefield(string name)
{
  spiral_nomenclature spiral_code(name);
  
  // Generate dual bond graph from spiral code and dualize to get fullerene graph
  Triangulation dual(spiral_code);
  FullereneGraph g = dual.dual_graph();
  
  // Compute initial geometry: First get 2D embedding of graph, then wrap around sphere
  g.layout2d = g.tutte_layout();
  Polyhedron P0(g, g.zero_order_geometry());

  int N = g.N;
  Coord3d *X            = (Coord3d*) calloc(N,sizeof(Coord3d));
  uint8_t *face_right   = (uint8_t*) calloc(3*N,sizeof(uint8_t));
  node_t 
    *neighbours   = (node_t*) calloc(3*N,sizeof(node_t)),
    *next_on_face = (node_t*) calloc(3*N,sizeof(node_t)),
    *prev_on_face = (node_t*) calloc(3*N,sizeof(node_t));
    

  for(node_t a=0;a<N;a++){
    for(int j=0;j<3;j++){
      node_t b = g.neighbours[a][j];

      X[a][j] = P0.points[a][j];
      neighbours  [a*3+j] = g.neighbours[a][j];
      next_on_face[a*3+j] = g.next_on_face(a,b);
      prev_on_face[a*3+j] = g.prev_on_face(a,b);
      face_right  [a*3+j] = face_size(g,a,b);
    }
  }

  FullereneForcefield *F = new   FullereneForcefield(N,neighbours,X,face_right,next_on_face,prev_on_face);

  return F;
}

void destroy_fullereneforcefield(FullereneForcefield *F){
  free((void*)F->neighbours);  
  free((void*)F->X);
  free((void*)F->face_right);
  free((void*)F->next_on_face);
  free((void*)F->prev_on_face);
  free((void*)F);
}
