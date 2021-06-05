#include <chrono>
#include <assert.h>
#include <inttypes.h>
#include <stdio.h>
#include <math.h>
using namespace std;

#include <array>
#include <tuple>
typedef double real_t;

/*********** DATA TYPES ***********/
#include "coord3d.cc"
typedef uint16_t              node_t;  // NB: Increase to int32 to do more than 65k atoms
template <int N> using CubicArcs = array<node_t, N*3>;



template <int N> class FullereneForcefieldEnergy {
public:

  /************** DATA *(************/
  const CubicArcs<N> neighbours;
  const array<uint8_t, N*3>  face_left, face_right;
  array<coord3d,N>      X;	   // Current nucleus positions  
   
  // Optimal parameters. TODO: Fill out with the correct constants as constexpr's.
  array<real_t,2>          optimal_corner_cos_angles;	    // indexed by face size: {p,h}  
  array<array<real_t,2>,2> optimal_bond_lengths; 	    // indexed by incident faces: {{pp,ph},{hp,pp}}
  array<array<array<real_t,2>,2>,2> optimal_dih_cos_angles; // indexed by ... faces: {{{ppp,pph},{php,phh}}, {{hpp,hph},{hhp,hhh}}} 

  // Force constants. TODO: Add them to the code, both here as constexpr and use them down in the energy function
  real_t bond_force, angle_force, dih_force;

  /************** CONSTRUCTORS **************/
  FullereneForcefieldEnergy(const CubicArcs<N> &neighbours, const array<coord3d,N> &X, const array<uint8_t,N*3> &face_left, const array<uint8_t,N*3> &face_right) :
    neighbours(neighbours), X(X), face_left(face_left), face_right(face_right) {}
  

  /************* IMPLEMENTATION **************/
  // The three parameters for the force field are bond-length stretching (Buster Eq. (15)),
  // bending of inner angles (Eq. (18)), and bending of dihedral angles (Eq. (22)).

  // The bond length is simply the length of the displacement vector ||b-a||
  inline real_t bond_length(const coord3d &ab){
    return sqrt(dot(ab,ab));
  }

  // The cosine of the angle between the vectors b-a and c-a is the dot product of the unit vectors
  inline real_t corner_cos_angle(const coord3d &ab_hat, const coord3d &ac_hat) {
    return dot(ab_hat, ac_hat);
  }

  // The calculation of the cosine to the dihedral angle is described in Buster, Eq. (46) (Appendix A)
  inline real_t dihedral_cos_angle(const coord3d& ab, const coord3d& ac, const coord3d &ad)
  {
    coord3d 
      bc = ac-ab,
      cd = ad-ac;

    // Unit vectors of arcs in BusterThesis Fig. 10
    coord3d 
      bc_hat = bc/sqrt(dot(bc,bc)),
      cd_hat = cd/sqrt(dot(cd,cd)),
      ba_hat = -ab/sqrt(dot(ab,ab)),
      cb_hat = -bc/sqrt(dot(bc,bc));

    // Find angles between vectors ba,bc and between vectors cb,cd
    real_t
      cos_abc = dot(ba_hat,bc_hat),
      cos_bcd = dot(cb_hat,cd_hat);

    // Sines for normalizing cross products. NB: Assumes angles < pi
    real_t sin_abc = sqrt(1 - cos_abc*cos_abc);
    real_t sin_bcd = sqrt(1 - cos_bcd*cos_bcd);
        
    // Normals to abc and bcd triangles
    coord3d n_abc = cross(ba_hat,bc_hat)/sin_abc;
    coord3d n_bcd = cross(cb_hat,cd_hat)/sin_bcd;

    real_t cos_beta = dot(n_abc,n_bcd);

    return cos_beta;
  }

  // Generic harmonic energy term, appropriate for any parameter.
  inline real_t harmonic_energy(real_t p, real_t p0) {  return 0.5*(p-p0)*(p-p0);  }    
  
  // Now we put it all together to calculate the total energy of a molecular geometry X
  real_t energy()
  {
    real_t energy = 0;
    
#pragma acc parallel loop reduction(+:energy) copyin(neighbours[:N*3], X[:N*3], face_left[:4], face_right[4:])
    for(node_t a=0;a<N;a++){	// For each node a
      // Data needed to process this node
      array<node_t,3>  na            = {neighbours[a*3],neighbours[a*3+1],neighbours[a*3+2]}; // Neighbours b,c,d to a 
      // Atomic positions
      coord3d          Xa           = X[a];
      array<coord3d,3> Xbcd         = {X[na[0]], X[na[1]], X[na[2]]};

      // Degree of faces incident to arcs a->b, a->c, a->d
      array<int,3>     f_r = {face_right[a*3+0]-5, face_right[a*3+1]-5, face_right[a*3+2]-5};
      array<int,3>     f_l = {face_left[a*3+0]-5,  face_left[a*3+1]-5,  face_left[a*3+2]-5};

      array<real_t,3>  rab;
      array<coord3d,3> ab_hats;
      
      for(int j=0;j<3;j++){	// Precalculate stuff needed by everyone
	tie(rab[j],ab_hats[j]) = split_norm(Xbcd[j]-Xa); // Length and direction from a to jth neighbour
      }

      for(int j=0;j<3;j++){
	// Energy contribution from bond length
	real_t bond_length0  = optimal_bond_lengths[ f_r[j] ][ f_l[j] ];
	real_t bond_length   = rab[j];
	energy += bond_force*harmonic_energy(bond_length,bond_length0); 

	// Energy contribution from angles
	real_t cos_angle0  = optimal_corner_cos_angles[ f_r[j] ];
	real_t cos_angle   = corner_cos_angle(ab_hats[j],ab_hats[(j+1)%3]);
	energy += angle_force*harmonic_energy(cos_angle, cos_angle0);
	
	// Energy contribution from dihedrals
	real_t dih_cos_angle0  = optimal_dih_cos_angles[ f_r[j] ][ f_r[(j+1)%3] ][ f_r[(j+2)%2] ];
	real_t dih_cos_angle   = dihedral_cos_angle(rab[j]*ab_hats[j], rab[(j+1)%3]*ab_hats[(j+1)%3], rab[(j+2)%3]*ab_hats[(j+2)%3]);
	energy += dih_force*harmonic_energy(dih_cos_angle,dih_cos_angle0);
      }
    }

    return energy;
  }

};


int main()
{
  // TODO: Read in.
  CubicArcs<60>        neighbours;
  array<coord3d,60>    X0;
  array<uint8_t,3*60>  face_left, face_right;
  
  FullereneForcefieldEnergy<60> F(neighbours,X0,face_left,face_right);

  printf("%g\n",F.energy());	// Just to generate the code
  
  return 0;
}
