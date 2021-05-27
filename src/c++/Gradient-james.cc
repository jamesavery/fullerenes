#include <chrono>
#include <assert.h>
#include <inttypes.h>
#include <stdio.h>
#include <math.h>
#include <iostream>
using namespace std;
using namespace std::literals;

#include <array>
#include <tuple>

#include "coord3d.cc"

typedef array<real_t, 3> coord3d;
typedef uint16_t node_t; // NB: Increase to int32 to do more than 65k atoms


// TODO: Move MxN array operations to somewhere sensible + make nicer
template <typename T, size_t N> inline T sum(const array<T,N>& A){
  T result=0;			// T default value must be 0
  for(size_t i=0;i<N;i++) result += A[i];
  return result;
}

template <typename T, size_t N> inline T dot(const array<T,N>& A, const array<T,N>& B){
  T result= 0;			// T default value must be 0
  for(size_t i=0;i<N;i++) result += A[i]*B[i];
  return result;
}


template <typename T, size_t N> inline T norm(const array<T,N>& A){
  return sqrt(dot(A,A));
}

template <typename T, size_t M, size_t N> T flatten_dot(const array<array<T,N>,M>& A, const array<array<T,N>,M>& B)
{
  T result=0;
  for(size_t i=0;i<M;i++)
    for(size_t j=0;j<N;j++) result += A[i][j] * B[i][j];
  return result;
}


template <typename T, size_t M, size_t N> 
array<array<T,N>,M>& operator*=(array<array<T,N>,M>& A, const real_t& s){
  for(size_t i=0;i<M;i++)
    for(size_t j=0;j<N;j++) A[i][j] *= s;
  return A;
}

template <typename T, size_t M, size_t N> 
array<array<T,N>,M>& operator/=(array<array<T,N>,M>& A, const real_t& s){
  for(size_t i=0;i<M;i++)
    for(size_t j=0;j<N;j++) A[i][j] /= s;
  return A;
}


template <int N>
using CubicArcs = array<node_t, N * 3>;

template <int N>
class FullereneForcefield
{
public:
    const CubicArcs<N> neighbours; //Bookkeeping array for storing the indices of neighbour nodes: b,c,d
    
    array<coord3d, N> X; // Current nucleus positions
    array<coord3d, N> X_temp; // Temporary positions to be evaluated
    array<real_t, N> energy_array; // Energy contributions from each node stored here.

    const array<uint8_t, N * 3> face_right; //Faces to the right of the edges ab, ac, ad
    const array<node_t, N * 3> next_on_face, prev_on_face; //Next node on the face to the left of the edges ab, ac, ad

    
    //All parameter arrays are indexed by a binary sum, 0,1,2,3,4,...
    //Pentagons = 0
    //Hexagons = 1
    //PPP = 0, {HPP, PHP, PPH} = 1, {PHH, HPH, HHP} = 2, {HHH} = 3
    const array<real_t, 2> optimal_corner_cos_angles = {cos(M_PI * 108 / 180), cos(M_PI * 120 / 180)}; 
    const array<real_t, 3> optimal_bond_lengths = {1.479, 1.458, 1.401}; 
    const array<real_t, 4> optimal_dih_cos_angles = {cos(0.652358), cos(0.509674), cos(0.417884), cos(0)}; 

    const array<real_t, 2> angle_forces = {100.0, 100.0}; 
    const array<real_t, 3> bond_forces = {260.0, 390.0, 450.0}; 
    const array<real_t, 4> dih_forces = {35.0, 65.0, 85.0, 270.0}; 

    FullereneForcefield(const CubicArcs<N> &neighbours, const array<coord3d, N> &X, const array<uint8_t, N * 3> &face_right, const array<node_t, N * 3> &next_on_face, const array<node_t, N * 3> &prev_on_face) : 
        neighbours(neighbours), X(X), face_right(face_right), next_on_face(next_on_face), prev_on_face(prev_on_face) {}


    //Parallelizable copying,  a = b;
    template <typename T>
    static inline void parallel_copy(array<T,N>& a, array<T,N>& b){
        for (node_t i = 0; i < N; i++)
        {
            a[i]= b[i];
        }
    }
    
    
    static inline real_t multiply_reduce(array<coord3d,N>& array1, array<coord3d,N>& array2){
        real_t total = 0.0;        
        for (node_t a = 0; a < N; a++)
        {
            total += dot(array1[a],array2[a]);
        }
        return total;
    }

    //Computes gradient related to stretch term.
    static inline coord3d bond_gradient(const real_t force_const, const real_t r0, const real_t rab, const coord3d ab_hat){
        return - force_const * (rab - r0) * ab_hat;
    }

    //Computes gradient related to bending term.
    static inline coord3d angle_gradient(const real_t force_const, const real_t ang0, const real_t r_rab, const real_t r_rac, const coord3d &ab_hat, const coord3d &ac_hat)
    {   
        real_t cos_angle = dot(ab_hat, ac_hat);
        coord3d grad = cos_angle * (ab_hat * r_rab + ac_hat * r_rac) - ab_hat * r_rac - ac_hat* r_rab;
        return force_const * (cos_angle - ang0) * grad;
    }

    //Computes gradient related to bending of angles.
    static inline coord3d outer_angle_gradient(const real_t force_const, const real_t ang0, const real_t r_rab, coord3d &ab_hat, const coord3d& bp_hat)
    {   
        real_t cos_angle = -dot(ab_hat, bp_hat);
        coord3d grad = (bp_hat + ab_hat * cos_angle) * r_rab;
        return force_const * (cos_angle - ang0) * grad;
    }
    
    //Helper Function for dihedral calculation.
    static inline tuple<coord3d, real_t, real_t> dihedral_terms(const coord3d& ba_hat, const coord3d& bc_hat){
        real_t cos_b = dot(ba_hat,bc_hat);
        real_t r_sin_b = 1.0/sqrt(1 - cos_b*cos_b);
        coord3d nabc = cross(ba_hat, bc_hat) * r_sin_b;
        return {nabc,cos_b,r_sin_b};
    }

    //Computes gradient related to dihedral/out-of-plane term.
    static inline pair<coord3d, real_t> dihedral_gradient(const real_t force_const, const real_t dih0, const real_t r_rab, const coord3d& ba_hat, const coord3d& bc_hat, const coord3d& cd_hat)
    {
        coord3d nabc, nbcd; real_t cos_b, cos_c, r_sin_b, r_sin_c;
        tie(nabc, cos_b, r_sin_b) = dihedral_terms(ba_hat, bc_hat);
        tie(nbcd, cos_c, r_sin_c) = dihedral_terms(-bc_hat, cd_hat);

        real_t cos_beta = dot(nabc, nbcd);
        real_t cot_b = cos_b * r_sin_b * r_sin_b;

        coord3d d_dih_a = cross(bc_hat, nbcd) * r_sin_b * r_rab - ba_hat * cos_beta * r_rab + (cot_b * cos_beta * r_rab) * (bc_hat - ba_hat * cos_b);

        return {force_const * (cos_beta - dih0) * d_dih_a, cos_beta};
    }
    
    //Computes gradient from dihedral angles constituted by the planes nbam, namp
    static inline coord3d outer_dih_a_gradient(const real_t force_const, const real_t dih0, const real_t r_rab, const real_t r_ram, const coord3d& ab_hat, const coord3d  am_hat, const coord3d mp_hat)//const real_t force_const, const real_t dih0, const real_t rab, const coord3d &ab, const coord3d &bm, const coord3d &bp, coord3d &ab_hat, coord3d &bp_hat, coord3d &bm_hat)
    {   
        coord3d nbam_hat, namp_hat; real_t cos_a, cos_m, r_sin_a, r_sin_m;
        tie(nbam_hat, cos_a, r_sin_a) = dihedral_terms(ab_hat, am_hat);
        tie(namp_hat, cos_m, r_sin_m) = dihedral_terms(-am_hat, mp_hat);
        
        real_t cos_beta = dot(nbam_hat, namp_hat);
        real_t cot_a = cos_a * r_sin_a * r_sin_a;
        real_t cot_m = cos_m * r_sin_m * r_sin_m;

        coord3d d_dih = cross(mp_hat,nbam_hat)*r_ram*r_sin_m - (cross(namp_hat,ab_hat)*r_ram + cross(am_hat,namp_hat)*r_rab)*r_sin_a +
                        cos_beta*(ab_hat*r_rab + r_ram * (2*am_hat + cot_m*(mp_hat+cos_m*am_hat)) - cot_a*(r_ram*(ab_hat - am_hat*cos_a) + r_rab*(am_hat-ab_hat*cos_a)));
        
        return force_const * (cos_beta - dih0) * d_dih;     
    }
    
    //Computes gradient from dihedral angles constituted by the planes nbmp, nmpa
    static inline coord3d outer_dih_m_gradient(const real_t force_const, const real_t dih0, const real_t r_rpa, const coord3d& mb_hat, const coord3d& mp_hat, const coord3d& pa_hat)
    {
        coord3d nbmp_hat, nmpa_hat; real_t cos_m, cos_p, r_sin_m, r_sin_p;
        tie(nbmp_hat, cos_m, r_sin_m) = dihedral_terms(mb_hat, mp_hat);
        tie(nmpa_hat, cos_p, r_sin_p) = dihedral_terms(-mp_hat, pa_hat);
        
        real_t cos_beta = dot(nbmp_hat, nmpa_hat);
        real_t cot_p = cos_p * r_sin_p * r_sin_p;
        
        coord3d d_dih_m = r_rpa * (cot_p*cos_beta * (-mp_hat - pa_hat*cos_p) - cross(nbmp_hat, mp_hat)*r_sin_p - pa_hat*cos_beta );
        return force_const * (cos_beta - dih0) * d_dih_m;
    }
    
    //Computes gradient from dihedral angles constituted by the planes nbpa, npam
    static inline coord3d outer_dih_p_gradient(const real_t force_const, const real_t dih0, const real_t r_ram, const real_t r_rpa, const coord3d& ap_hat, const coord3d& am_hat, const coord3d& pb_hat)
    {
        coord3d nbpa_hat, npam_hat; real_t cos_p, cos_a, r_sin_p, r_sin_a;
        tie(npam_hat, cos_a, r_sin_a) = dihedral_terms(ap_hat, am_hat);
        tie(nbpa_hat, cos_p, r_sin_p) = dihedral_terms(pb_hat, -ap_hat);
        real_t cos_beta = dot(nbpa_hat, npam_hat);
        real_t cot_p = cos_p * r_sin_p * r_sin_p;
        real_t cot_a = cos_a * r_sin_a * r_sin_a;

        coord3d d_dih = cross(npam_hat,pb_hat)*r_rpa*r_sin_p - (cross(am_hat,nbpa_hat)*r_rpa + cross(nbpa_hat,ap_hat)*r_ram)*r_sin_a +
                        cos_beta*(am_hat*r_ram + r_rpa * (2*ap_hat + cot_p*(pb_hat+cos_p*ap_hat)) - cot_a*(r_rpa*(am_hat - ap_hat*cos_a) + r_ram*(ap_hat-am_hat*cos_a)));
        
        return force_const * (cos_beta - dih0) * d_dih;
    }

    static inline real_t harmonic_energy(real_t param, real_t param0){
        return 0.5*(param-param0)*(param-param0);
    }

    real_t energy(array<coord3d,N>& X){
        for (node_t a = 0; a < N; a++)
        {   
            real_t node_bond_energy = 0.0;
            real_t node_bend_energy = 0.0;
            real_t node_dihedral_energy = 0.0;

             //Fetching neighbour indices.
            array<node_t, 3> na = {neighbours[a * 3], neighbours[a * 3 + 1], neighbours[a * 3 + 2]}; // Neighbours b,c,d to a

            //Fetching current node coordinate.
            coord3d Xa = X[a];

            //Fetching face information to the left and right of the edges: ab, ac, ad.
            array<uint8_t, 3> f_r = {face_right[a * 3],     face_right[a * 3 + 1] , face_right[a * 3 + 2] };
            array<uint8_t, 3> f_l = {face_right[a * 3 + 2], face_right[a * 3] ,     face_right[a * 3 + 1] };

            //Fetching coordinates for neighbour nodes b,c,d.
            array<coord3d, 3> Xbcd = {X[na[0]], X[na[1]], X[na[2]]};

            array<real_t,3> r_rab;
            array<coord3d,3> abs;
            array<coord3d,3> ab_hats;

            for (size_t i = 0; i < 3; i++)
            {
                tie(r_rab[i],abs[i],ab_hats[i]) = split_norm(Xbcd[i] - Xa);
            }
            for (size_t j = 0; j < 3; j++)
            {
                real_t r0 = optimal_bond_lengths[ f_l[j] + f_r[j] ];
                real_t ang0 = optimal_corner_cos_angles[ f_r[j] ];
                real_t dih0 = optimal_dih_cos_angles[ f_l[0] + f_l[1] + f_l[2]];
                real_t bond_fc = bond_forces[ f_l[j] + f_r[j] ];
                real_t ang_fc = angle_forces[ f_l[j] ];
                real_t dih_fc = dih_forces[ f_l[0] + f_l[1] + f_l[2] ];

                real_t cos_dih_angle, r_rbc, r_rcd;
                coord3d null, bc, bc_hat, cd, cd_hat;

                tie(null, cos_dih_angle) = dihedral_gradient(dih_fc, dih0, r_rab[j], -ab_hats[j], unit_vector(abs[(j+1)%3] - abs[j]), unit_vector(abs[(j+2)%3] - abs[(j+1)%3]));

                node_bond_energy += bond_fc * harmonic_energy(1.0/r_rab[j],r0);
                node_bend_energy += ang_fc * harmonic_energy(dot(ab_hats[j],ab_hats[(j+1)%3]), ang0);
                node_dihedral_energy += dih_fc * harmonic_energy(cos_dih_angle,dih0);
            }
            //All bond energies are computed twice, (node a computes ab and node b computes ba) , hence the factor 0.5;
            energy_array[a] = node_bond_energy*0.5 + node_bend_energy + node_dihedral_energy;
        }

        return sum(energy_array);
    }

    void gradient(array<coord3d,N>& X,array<coord3d,N>& gradient)
    {   
        #pragma acc parallel loop
        for (node_t a = 0; a < N; a++) 
        {   
            //Fetching neighbour indices.
            array<node_t, 3> na = {neighbours[a * 3], neighbours[a * 3 + 1], neighbours[a * 3 + 2]}; // Neighbours b,c,d to a

            //Fetching current node coordinate.
            coord3d Xa = X[a];

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
            array<coord3d, 3> Xbcd = {X[na[0]], X[na[1]], X[na[2]]};

            //Fetching coordinates for outer nodes.
            array<coord3d, 3> Xbp = {X[next_on_face[a*3]], X[next_on_face[a*3 + 1]], X[next_on_face[a*3 + 2]]};
            array<coord3d, 3> Xbm = {X[prev_on_face[a*3]], X[prev_on_face[a*3 + 1]], X[prev_on_face[a*3 + 2]]};

            //Vector reciprocal norms.
            array<real_t, 3> r_rab, r_rbp, r_rbm;

            //Common vectors.
            array<coord3d, 3> abs, bps, bms, ab_hats, bp_hats, bm_hats;

            coord3d am, mp, ap, bc, cd, am_hat, mp_hat, ap_hat, bc_hat, cd_hat;

            real_t r_ram, r_rmp, r_rap, r_rbc, r_rcd;

            coord3d node_gradient = {0,0,0};
            
            for (size_t i = 0; i < 3; i++)
            {
                tie(r_rab[i], abs[i], ab_hats[i]) = split_norm(Xbcd[i] - Xa);
                tie(r_rbm[i], bms[i], bm_hats[i]) = split_norm(Xbm[i] - Xbcd[i]);
                tie(r_rbp[i], bps[i], bp_hats[i]) = split_norm(Xbp[i] - Xbcd[i]);
            }
            for (size_t j = 0; j < 3; j++)
            {   
                //Indexing of 'Equillibrium' parameters and related force constants, through face information.
                real_t r0 = optimal_bond_lengths[ f_l[j] + f_r[j] ];
                real_t ang0 = optimal_corner_cos_angles[ f_r[j] ];
                real_t dih0 = optimal_dih_cos_angles[ f_l[0] + f_l[1] + f_l[2]];
                real_t outer_ang0_m = optimal_corner_cos_angles[ f_l[j] ];
                real_t outer_ang0_p = optimal_corner_cos_angles[ f_r[j] ];
                real_t outer_dih0 = optimal_dih_cos_angles[ dihedral_face_sum_bcd[j] ];
                real_t bond_fc = bond_forces[ f_l[j] + f_r[j] ];
                real_t ang_fc = angle_forces[ f_l[j] ];
                real_t dih_fc = dih_forces[ f_l[0] + f_l[1] + f_l[2] ];
                real_t out_ang_fc_m = angle_forces[ f_r[j] ];
                real_t out_ang_fc_p = angle_forces[ f_l[j] ];
                real_t out_dih_fc = dih_forces[ dihedral_face_sum_bcd[j] ];

                //Computation of gradient terms.
                coord3d bond_grad = bond_gradient(bond_fc, r0, 1.0/r_rab[j], ab_hats[j]);
                coord3d angle_grad = angle_gradient(ang_fc, ang0, r_rab[j], r_rab[(j+1)%3], ab_hats[j], ab_hats[(j+1)%3]);
                
                real_t cos_beta; coord3d dihedral_grad;

                tie(dihedral_grad, cos_beta) = dihedral_gradient(dih_fc, dih0, r_rab[j], -ab_hats[j], unit_vector(abs[(j+1)%3] - abs[j]), unit_vector(abs[(j+2)%3] - abs[(j+1)%3]));

                coord3d outer_ang_grad_m = outer_angle_gradient(out_ang_fc_m, outer_ang0_m, r_rab[j], ab_hats[j], bm_hats[j]);
                coord3d outer_ang_grad_p = outer_angle_gradient(out_ang_fc_p, outer_ang0_p, r_rab[j], ab_hats[j], bp_hats[j]);
                
                //Compute relevant dihedral vectors.
                tie(r_ram, am, am_hat) = split_norm(bms[j] + abs[j]);
                tie(r_rmp, mp, mp_hat) = split_norm(bps[j] - bms[j]);
                tie(r_rap, ap, ap_hat) = split_norm(bps[j] + abs[j]);

                coord3d outer_dihedral_a_grad = outer_dih_a_gradient(out_dih_fc, outer_dih0, r_rab[j], r_ram, ab_hats[j], am_hat, mp_hat);
                coord3d outer_dihedral_m_grad = outer_dih_m_gradient(out_dih_fc, outer_dih0, r_rap, -bm_hats[j], mp_hat, -ap_hat);
                coord3d outer_dihedral_p_grad = outer_dih_p_gradient(out_dih_fc, outer_dih0, r_ram, r_rap, ap_hat, am_hat, -bp_hats[j]);
                
                node_gradient = node_gradient + bond_grad + angle_grad + dihedral_grad + outer_ang_grad_m + outer_ang_grad_p + outer_dihedral_a_grad + outer_dihedral_m_grad + outer_dihedral_p_grad; 

            }
            gradient[a] = node_gradient;
        }
    }

    void bisection_search(array<coord3d,N>& X, array<coord3d,N>& direction, array<coord3d,N>& fc, real_t a, real_t b, real_t eps, size_t max_iter){
        size_t count = 0;
        size_t nBrack = 0;
        real_t sum = -1;
        real_t dfc = 0.0;
        real_t coefficient = 0.0;

        array<coord3d,N> temp_coords;
        array<coord3d,N> temp_direction;

        //Copy coordinates into temporary coordinates.
        parallel_copy(temp_coords, X);

        for (node_t a = 0; a < N; a++)
        {
            dfc -= dot(direction[a],direction[a]);
        }

        while (sum < 0)
        {   
            if (nBrack == 0)
            {
                coefficient = 0.0;
            } else
            {
                coefficient = b;
                b *= 1.5;
            }
            coefficient = b - coefficient;
            
           
            nBrack ++;
            for (node_t a = 0; a < N; a++)
            {
                temp_coords[a] = temp_coords[a] + coefficient*direction[a];
            }
            gradient(temp_coords,temp_direction);

            sum = multiply_reduce(temp_direction,direction);
        }
        
        while (abs(dfc) > eps)
        {   
            parallel_copy(temp_coords,X);
            count++;
            real_t c = (a+b)/2;

            for (node_t i = 0; i < N; i++)
            {   
                temp_coords[i] = temp_coords[i] + c*direction[i];
            }
            gradient(temp_coords,fc);
            dfc = multiply_reduce(fc,direction);
            
            if (count > max_iter){ break;}
            if (dfc < 0){a = c;} else{b = c;}
        }
        parallel_copy(X, temp_coords);
        for (node_t i = 0; i < N; i++)
        {
            fc[i]= -fc[i];
        }
        
    }



    //Parallelizable euclidean norm function.
    static inline real_t norm(array<coord3d,N>& array){
        return sqrt(array_dot(array,array));
    }

    static inline real_t array_dot(array<coord3d,N>& array1, array<coord3d,N>& array2){
        real_t result = 0.0;
        for (node_t a = 0; a < N; a++)
        {
            result += dot(array1[a],array2[a]);
        }
        return result;
        
    }
  
    void conjugate_gradient(){
        size_t iter_count = 0;
        size_t max_iter = N*10;
        real_t beta = 0.0;
        real_t dnorm = 1.0;

        array<coord3d, N> delta_x0;
        array<coord3d, N> delta_x1;
        array<coord3d, N> direction;

        gradient(X, direction);
        dnorm = norm(direction);
	direction *= -1/dnorm;

        parallel_copy(X_temp,X);
        parallel_copy(delta_x0, direction);
        while (dnorm > 1e-7)
        {   
            beta = 0.0;
            bisection_search(X_temp, direction,delta_x1,0,1e-5,1e-10,N);
            //Polak Ribiere method
            for (node_t a = 0; a < N; a++)
            {
                beta += dot(delta_x1[a], (delta_x1[a] - delta_x0[a]));
            }
            beta /= array_dot(delta_x0,delta_x0);
        
            if (energy(X_temp) > energy(X))
            {
                parallel_copy(X_temp, X);
                parallel_copy(delta_x1, delta_x0);
                beta = 0.0;
            }
            else
            {   
                parallel_copy(X, X_temp);
                parallel_copy(delta_x0,delta_x1);
            }
            for (node_t a = 0; a < N; a++)
            {
                direction[a] = delta_x1[a] + beta*direction[a];
            }
            dnorm = norm(direction);
	    direction *= 1/dnorm;

            dnorm = norm(delta_x1);
            if (iter_count > N*10){return;}
            iter_count++;
        }   
    }
};

int main()
{
    const int size = 60;

    //Data necessary to compute gradients for the C60ih Fullerene
    CubicArcs<size> cubic_neighbours  = {4,15,1,0,12,2,1,9,3,2,5,4,3,8,0,3,11,6,5,22,7,6,20,8,7,18,4,2,14,10,9,26,11,10,25,5,1,17,13,12,30,14,13,29,9,0,19,16,15,34,17,16,33,12,8,21,19,18,37,15,7,24,21,20,38,18,6,25,23,22,42,24,23,40,20,11,28,22,10,29,27,26,45,28,27,44,25,14,32,26,13,33,31,30,48,32,31,47,29,17,36,30,16,37,35,34,51,36,35,50,33,19,39,34,21,41,39,38,53,37,24,43,41,40,54,38,23,44,43,42,55,40,28,46,42,27,47,46,45,57,44,32,49,45,31,50,49,48,58,47,36,52,48,35,53,52,51,59,50,39,54,51,41,56,53,43,57,56,55,59,54,46,58,55,49,59,57,52,56,58};    
    array<node_t, size*3> next_on_face = {8,16,2,15,13,3,12,10,4,9,6,0,5,18,1,2,25,7,11,23,8,22,21,4,20,19,3,1,29,11,14,27,5,26,22,3,0,33,14,17,31,9,30,26,2,4,37,17,19,35,12,34,30,1,7,38,15,21,34,0,6,40,18,24,39,8,5,28,24,25,43,20,42,41,7,10,44,6,9,32,28,29,46,25,45,42,11,13,47,10,12,36,32,33,49,29,48,45,14,16,50,13,15,39,36,37,52,33,51,48,17,18,53,16,20,54,37,41,51,19,23,55,38,43,53,21,22,46,40,44,56,24,27,57,23,26,49,44,47,55,28,31,58,27,30,52,47,50,57,32,35,59,31,34,54,50,53,58,36,38,56,35,40,59,39,42,58,54,57,52,41,45,59,43,48,56,46,51,55,49};
    array<node_t, size*3> prev_on_face = {3,19,12,4,17,9,0,14,5,1,11,8,2,7,15,4,10,22,3,25,20,5,24,18,6,21,0,3,13,26,2,29,25,9,28,6,2,16,30,1,33,29,12,32,10,1,18,34,0,37,33,15,36,13,4,20,37,8,39,16,8,23,38,7,41,19,7,11,42,6,44,40,22,43,21,5,27,23,11,14,45,10,47,44,26,46,22,9,31,27,14,17,48,13,50,47,30,49,26,12,35,31,17,19,51,16,53,50,34,52,30,15,38,35,18,40,53,21,54,34,20,42,54,24,56,39,24,28,55,23,57,41,25,45,43,28,32,57,27,58,42,29,48,46,32,36,58,31,59,45,33,51,49,36,39,59,35,56,48,37,41,52,38,55,51,40,46,59,43,58,53,44,49,56,47,52,55,50,54,57};
    array<uint8_t, size*3> face_right   = {1,1,0,1,1,0,1,1,0,1,1,0,1,1,0,1,0,1,0,1,1,1,0,1,0,1,1,1,0,1,0,1,1,1,0,1,1,0,1,0,1,1,1,0,1,1,0,1,0,1,1,1,0,1,0,1,1,1,0,1,1,1,0,1,1,0,0,1,1,1,0,1,0,1,1,1,1,0,0,1,1,1,0,1,0,1,1,1,1,0,0,1,1,1,0,1,0,1,1,1,1,0,0,1,1,1,0,1,0,1,1,1,1,0,1,0,1,0,1,1,0,1,1,1,0,1,1,1,0,1,1,0,0,1,1,1,1,0,1,1,0,0,1,1,1,1,0,1,1,0,0,1,1,1,1,0,1,1,0,0,1,1,1,1,0,1,0,1,0,1,1,1,0,1,1,0,1,1,0,1}; // 0: pentagon, 1: hexagon
    array<coord3d, size> X = {2.82489,8.43261e-17,14.2017,8.0219,0.604118,12.0396,10.4266,6.01981,8.04461,4.53413,6.64511,12.0396,1.41245,2.44643,14.2017,3.29447,11.5801,8.04461,-3.29447,11.5801,8.04461,-4.53413,6.64511,12.0396,-1.41245,2.44643,14.2017,12.299,7.10085,2.82489,9.58483,10.4795,-2.82489,4.63059,13.4256,2.82489,11.6759,-2.93695,8.04461,13.9422,-2.70257,2.82489,13.8679,3.06098,-2.82489,1.41245,-2.44643,14.2017,3.48777,-7.24923,12.0396,8.38143,-8.64315,8.04461,-2.82489,5.87089e-16,14.2017,-1.41245,-2.44643,14.2017,-10.4266,6.01981,8.04461,-8.0219,0.604118,12.0396,-4.63059,13.4256,2.82489,-9.58483,10.4795,-2.82489,-12.299,7.10085,2.82489,2.39962e-16,14.2017,-2.82489,8.97979,8.0197,-8.04461,3.58009,7.20408,-12.0396,1.07573e-16,12.0396,-8.04461,11.4352,3.76688,-8.04461,12.299,-7.10085,-2.82489,10.4266,-6.01981,-8.04461,8.02896,-0.501593,-12.0396,9.31158,-10.723,2.82489,-2.56576e-15,-12.0396,8.04461,-2.91345e-15,-14.2017,2.82489,4.28306,-13.5404,-2.82489,-3.48777,-7.24923,12.0396,-11.6759,-2.93695,8.04461,-8.38143,-8.64315,8.04461,-13.8679,3.06098,-2.82489,-13.9422,-2.70257,2.82489,-8.97979,8.0197,-8.04461,-11.4352,3.76688,-8.04461,-3.58009,7.20408,-12.0396,1.33536,2.48934,-14.2017,-1.33536,2.48934,-14.2017,2.82351,-0.088213,-14.2017,4.44887,-6.70249,-12.0396,1.48815,-2.40113,-14.2017,2.45537,-11.7866,-8.04461,-4.28306,-13.5404,-2.82489,-2.45537,-11.7866,-8.04461,-9.31158,-10.723,2.82489,-12.299,-7.10085,-2.82489,-8.02896,-0.501593,-12.0396,-10.4266,-6.01981,-8.04461,-2.82351,-0.088213,-14.2017,-1.48815,-2.40113,-14.2017,-4.44887,-6.70249,-12.0396};

    //Gradient container
    array<coord3d,size> grad;

    //Test gradient computation
    FullereneForcefield<60> forcefield = FullereneForcefield<60>(cubic_neighbours, X, face_right, next_on_face, prev_on_face);
    
    
    forcefield.gradient(X,grad);



    for (size_t i = 0; i < size; i++)
    {
        //print_coord(grad[i]);
    }
    
    
    forcefield.X = X;

    auto start = chrono::system_clock::now();
    forcefield.conjugate_gradient();
    auto end = chrono::system_clock::now();
    cout << "Elapsed time: " << (end-start)/ 1ms << "ms\n" ;
    printf("%.16g\n",double(forcefield.energy(forcefield.X)));
    
    //write_to_file<size>(forcefield.X);

    for (size_t i = 0; i < size; i++)
    {
        print_coord(forcefield.X[i]);
    }
    
    
}
