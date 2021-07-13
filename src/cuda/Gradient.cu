#include "fullerenes/gpu/isomerspace_forcefield.hh"
#include <chrono>
#include <assert.h>
#include <inttypes.h>
#include <stdio.h>
#include <math.h>
#include <iostream>
using namespace std;
using namespace std::literals;

namespace IsomerspaceForcefield {

#include "coord3d.cu"
typedef uint16_t node_t;




// NB: Increase to int32 to do more than 65k atoms

//Convert neighbouring face information into a unique linear index.
inline uint8_t face_index(const uint8_t f1, const uint8_t f2, const uint8_t f3){
    return f1*4 + f2*2 + f3;
}


inline void print_real(real_t a)
{
    printf("%.16e \n",a);
}

class FullereneForcefield
{
public:
    const node_t* neighbours; //Bookkeeping array for storing the indices of neighbour nodes: b,c,d
    coord3d* X_temp; // Temporary positions to be evaluated
    coord3d* X; // Current nucleus positions
    const size_t N;

    const uint8_t* face_right; //Faces to the right of the edges ab, ac, ad
    const node_t* next_on_face; const node_t* prev_on_face; //Next node on the face to the left of the edges ab, ac, ad

    struct ArcData
    {   
        //All parameter arrays are indexed by a binary sum, 0,1,2,3,4,...
        //Pentagons = 0
        //Hexagons = 1
        //PPP = 0, {HPP, PHP, PPH} = 1, {PHH, HPH, HHP} = 2, {HHH} = 3
        const real_t optimal_corner_cos_angles[2] = {-0.3090169944, -0.5}; 
        const real_t optimal_bond_lengths[3] = {1.479, 1.458, 1.401}; 
        const real_t optimal_dih_cos_angles[8] = {0.79465455715, 0.87290360705, 0.87290360705, 0.9410338473, 0.816287936, 0.913965949, 0.913965949, 1}; 

        const real_t angle_forces[2] = {207.924,216.787}; 
        const real_t bond_forces[3] = {260.0, 353.377, 518.992}; 
        const real_t dih_forces[4] = {35.0,65.0,3.772,270.0}; 

        ArcData(const node_t a, const uint8_t j, const node_t* neighbours, const coord3d* X, const uint8_t* face_right, const node_t* next_on_face, const node_t* prev_on_face)
        {   
            real_t r_rmp;
            coord3d ap, am, ab, ac, ad, mp;
            //printf("Index: %d \n", a*3 + j);

            //Compute the arcs ab, ac, ad, bp, bm, ap, am, mp, bc and cd
            ab = (X[neighbours[(node_t)(a*3) + j]] - X[a]);  
            r_rab = bond_length(ab); 
            ab_hat = r_rab * ab;
            ac = (X[neighbours[(node_t)(a*3) + (j+1)%3]] - X[a]); r_rac = bond_length(ac); ac_hat = r_rac * ac;
            ad = (X[neighbours[(node_t)(a*3) + (j+2)%3]] - X[a]); r_rad = bond_length(ad); ad_hat = r_rad * ad;
            
            coord3d bp = (X[next_on_face[(node_t)a*3 + j]] - X[neighbours[(node_t)a*3 + j]]); bp_hat = unit_vector(bp);
            coord3d bm = (X[prev_on_face[(node_t)a*3 + j]] - X[neighbours[(node_t)a*3 + j]]); bm_hat = unit_vector(bm);

            ap = bp + ab; r_rap = bond_length(ap); ap_hat = r_rap * ap;
            am = bm + ab; r_ram = bond_length(am); am_hat = r_ram * am;
            mp = bp - bm; r_rmp = bond_length(mp); mp_hat = r_rmp * mp;

            bc_hat = unit_vector(ac - ab);
            cd_hat = unit_vector(ad - ac);

            //Compute inverses of some arcs, these are subject to be omitted if the equations are adapted appropriately with inversion of signs.
            ba_hat = -ab_hat;
            mb_hat = -bm_hat;
            pa_hat = -ap_hat;
            pb_hat = -bp_hat;
            
            uint8_t f_r = face_right[a * 3 + j] - 5;
            uint8_t f_l = face_right[a * 3 + (2 + j)%3] - 5;

            uint8_t face_sum = face_right[a * 3] - 5 + face_right[a * 3 + 1] - 5 + face_right[a * 3 + 2] - 5;
            uint8_t dihedral_face_sum = face_right[neighbours[a*3 + j] * 3]-5 + face_right[neighbours[a*3 + j] * 3 + 1]-5 +  face_right[neighbours[a*3 + j] * 3 + 2]-5;

            uint8_t dihedral_index_a = face_index(f_l,face_right[a * 3 + (1 + j)%3] - 5,f_r);
            uint8_t dihedral_index_m =  face_index(face_right[a * 3 + (1 + j)%3] - 5, f_r, f_l);
            uint8_t dihedral_index_p = face_index(f_r,f_l, face_right[a * 3 + (1 + j)%3] - 5);

            outer_dih0_a = optimal_dih_cos_angles[dihedral_index_a];
            outer_dih0_m = optimal_dih_cos_angles[dihedral_index_m];
            outer_dih0_p = optimal_dih_cos_angles[dihedral_index_p];

            //Load equillibirium distance, angles and dihedral angles from face information.
            r0 = optimal_bond_lengths[ f_l + f_r ];
            angle0 = optimal_corner_cos_angles[ f_r ];
            inner_dih0 = optimal_dih_cos_angles[ face_sum ];
            outer_angle_m0 = optimal_corner_cos_angles[ f_l ];
            outer_angle_p0 = optimal_corner_cos_angles[ f_r ];

            //Load force constants from neighbouring face information.
            f_bond = bond_forces[ f_l + f_r ];
            f_inner_angle = angle_forces[ f_l ];
            f_inner_dihedral = dih_forces[ face_sum];
            f_outer_angle_m = angle_forces[ f_r ];
            f_outer_angle_p = angle_forces[ f_l ];
            f_outer_dihedral = dih_forces[ dihedral_face_sum];
            
        }
        static inline real_t harmonic_energy(const real_t p0, const real_t p){
            return 0.5*(p-p0)*(p-p0);
        }
        static inline coord3d  harmonic_energy_gradient(const real_t p0, const real_t p, const coord3d gradp){
            return (p-p0)*gradp;     
        }

        inline real_t bond() const {return 1/r_rab;}
        inline real_t angle()       const {return dot(ab_hat,ac_hat);}
        //Returns the inner dihedral angle for the current arc. Used here only for energy calculation, 
        //otherwise embedded in dihedral computation because the planes and angles that make up the dihedral angle computation are required for derivative computation.
        inline real_t dihedral()    const 
        { 
            coord3d nabc, nbcd; real_t cos_b, cos_c, r_sin_b, r_sin_c;
            cos_b = dot(ba_hat,bc_hat); r_sin_b = rsqrt((real_t)1.0 - cos_b*cos_b); nabc = cross(ba_hat, bc_hat) * r_sin_b;
            cos_c = dot(-bc_hat,cd_hat); r_sin_c = rsqrt((real_t)1.0 - cos_c*cos_c); nbcd = cross(-bc_hat,cd_hat) * r_sin_c;

            return dot(nabc, nbcd);
        }
        
        // Chain rule terms for angle calculation
        //Computes gradient related to bending term. ~24 FLOPs
        inline coord3d inner_angle_gradient() const
        {   
            
            real_t cos_angle = angle(); //Inner angle of arcs ab,ac.
            coord3d grad = cos_angle * (ab_hat * r_rab + ac_hat * r_rac) - ab_hat * r_rac - ac_hat* r_rab; //Derivative of inner angle: Eq. 21. 
            return f_inner_angle * harmonic_energy_gradient(angle0, cos_angle, grad); //Harmonic Energy Gradient: Eq. 21. multiplied by harmonic term.
        }
        //Computes gradient related to bending of outer angles. ~20 FLOPs
        inline coord3d outer_angle_gradient_m() const
        {
            real_t cos_angle = -dot(ab_hat, bm_hat); //Compute outer angle. ab,bm
            coord3d grad = (bm_hat + ab_hat * cos_angle) * r_rab; //Derivative of outer angles Eq. 30. Buster Thesis
            return f_outer_angle_m * harmonic_energy_gradient(outer_angle_m0,cos_angle,grad); //Harmonic Energy Gradient: Eq. 30 multiplied by harmonic term.
        }
        inline coord3d outer_angle_gradient_p() const
        {
            real_t cos_angle = -dot(ab_hat, bp_hat); //Compute outer angle. ab,bp
            coord3d grad = (bp_hat + ab_hat * cos_angle) * r_rab; //Derivative of outer angles Eq. 28. Buster Thesis
            return f_outer_angle_p * harmonic_energy_gradient(outer_angle_p0,cos_angle,grad); //Harmonic Energy Gradient: Eq. 28 multiplied by harmonic term.
        }
        // Chain rule terms for dihedral calculation
        //Computes gradient related to dihedral/out-of-plane term. ~75 FLOPs
        inline coord3d inner_dihedral_gradient() const
        {
            coord3d nabc, nbcd; real_t cos_b, cos_c, r_sin_b, r_sin_c;
            cos_b = dot(ba_hat,bc_hat); r_sin_b = (real_t)1.0/sqrt((real_t)1.0 - cos_b*cos_b); nabc = cross(ba_hat, bc_hat) * r_sin_b;
            cos_c = dot(-bc_hat,cd_hat); r_sin_c = (real_t)1.0/sqrt((real_t)1.0 - cos_c*cos_c); nbcd = cross(-bc_hat,cd_hat) * r_sin_c;

            real_t cos_beta = dot(nabc, nbcd); //Inner dihedral angle from planes abc,bcd.
            real_t cot_b = cos_b * r_sin_b * r_sin_b; //cos(b)/sin(b)^2

            //Derivative w.r.t. inner dihedral angle F and G in Eq. 26
            coord3d grad = cross(bc_hat, nbcd) * r_sin_b * r_rab - ba_hat * cos_beta * r_rab + (cot_b * cos_beta * r_rab) * (bc_hat - ba_hat * cos_b);

            return f_inner_dihedral * harmonic_energy_gradient(inner_dih0, cos_beta, grad); //Eq. 26.
        }

        //Computes gradient from dihedral angles constituted by the planes bam, amp ~162 FLOPs
        inline coord3d outer_a_dihedral_gradient() const
        {
            coord3d nbam_hat, namp_hat; real_t cos_a, cos_m, r_sin_a, r_sin_m;

            cos_a = dot(ab_hat,am_hat); r_sin_a = (real_t)1.0/sqrtf((real_t)1.0 - cos_a*cos_a); nbam_hat = cross(ab_hat,am_hat) * r_sin_a;
            cos_m = dot(-am_hat,mp_hat); r_sin_m = (real_t)1.0/sqrtf((real_t)1.0 - cos_m*cos_m); namp_hat = cross(-am_hat,mp_hat) * r_sin_m;
        
            real_t cos_beta = dot(nbam_hat, namp_hat); //Outer Dihedral angle bam, amp
            real_t cot_a = cos_a * r_sin_a * r_sin_a;
            real_t cot_m = cos_m * r_sin_m * r_sin_m;

            //Derivative w.r.t. outer dihedral angle, factorized version of Eq. 31.
            coord3d grad = cross(mp_hat,nbam_hat)*r_ram*r_sin_m - (cross(namp_hat,ab_hat)*r_ram + cross(am_hat,namp_hat)*r_rab)*r_sin_a +
                            cos_beta*(ab_hat*r_rab + r_ram * (2*am_hat + cot_m*(mp_hat+cos_m*am_hat)) - cot_a*(r_ram*(ab_hat - am_hat*cos_a) + r_rab*(am_hat-ab_hat*cos_a)));
            
            //Eq. 31 multiplied by harmonic term.
            return f_outer_dihedral * harmonic_energy_gradient(outer_dih0_a, cos_beta, grad);
        }

        //Computes gradient from dihedral angles constituted by the planes nbmp, nmpa ~92 FLOPs
        inline coord3d outer_m_dihedral_gradient() const
        {
            coord3d nbmp_hat, nmpa_hat; real_t cos_m, cos_p, r_sin_m, r_sin_p;
            cos_m = dot(mb_hat,mp_hat);  r_sin_m = (real_t)1.0/sqrtf((real_t)1.0 - cos_m*cos_m); nbmp_hat = cross(mb_hat,mp_hat) * r_sin_m;
            cos_p = dot(-mp_hat,pa_hat); r_sin_p = (real_t)1.0/sqrtf((real_t)1.0 - cos_p*cos_p); nmpa_hat = cross(-mp_hat,pa_hat) * r_sin_p;
        
            //Cosine to the outer dihedral angle constituted by the planes bmp and mpa
            real_t cos_beta = dot(nbmp_hat, nmpa_hat); //Outer dihedral angle bmp,mpa.
            real_t cot_p = cos_p * r_sin_p * r_sin_p;
            
            //Derivative w.r.t. outer dihedral angle, factorized version of Eq. 32.
            coord3d grad = r_rap * (cot_p*cos_beta * (-mp_hat - pa_hat*cos_p) - cross(nbmp_hat, mp_hat)*r_sin_p - pa_hat*cos_beta );

            //Eq. 32 multiplied by harmonic term.
            return f_outer_dihedral * harmonic_energy_gradient(outer_dih0_m, cos_beta, grad);
        }

        //Computes gradient from dihedral angles constituted by the planes bpa, pam ~162 FLOPs
        inline coord3d outer_p_dihedral_gradient() const
        {
            coord3d nbpa_hat, npam_hat; real_t cos_p, cos_a, r_sin_p, r_sin_a;
            cos_a = dot(ap_hat,am_hat);  r_sin_a = (real_t)1.0/sqrtf((real_t)1.0 - cos_a*cos_a); npam_hat = cross(ap_hat,am_hat)  * r_sin_a;
            cos_p = dot(pb_hat,-ap_hat); r_sin_p = (real_t)1.0/sqrtf((real_t)1.0 - cos_p*cos_p); nbpa_hat = cross(pb_hat,-ap_hat) * r_sin_p;

            real_t cos_beta = dot(nbpa_hat, npam_hat); //Outer dihedral angle bpa, pam.
            real_t cot_p = cos_p * r_sin_p * r_sin_p;
            real_t cot_a = cos_a * r_sin_a * r_sin_a;

            //Derivative w.r.t. outer dihedral angle, factorized version of Eq. 33.
            coord3d grad = cross(npam_hat,pb_hat)*r_rap*r_sin_p - (cross(am_hat,nbpa_hat)*r_rap + cross(nbpa_hat,ap_hat)*r_ram)*r_sin_a +
                            cos_beta*(am_hat*r_ram + r_rap * (2*ap_hat + cot_p*(pb_hat+cos_p*ap_hat)) - cot_a*(r_rap*(am_hat - ap_hat*cos_a) + r_ram*(ap_hat-am_hat*cos_a)));
            
            //Eq. 33 multiplied by harmonic term.
            return f_outer_dihedral * harmonic_energy_gradient(outer_dih0_p, cos_beta, grad);
        }
        // Internal coordinate gradients
        inline coord3d bond_length_gradient() const { return - f_bond * harmonic_energy_gradient(r0,bond(),ab_hat);}
        //Sum of angular gradient components.
        inline coord3d angle_gradient()       const { return inner_angle_gradient() + outer_angle_gradient_p() + outer_angle_gradient_m();}
        //Sum of inner and outer dihedral gradient components.
        inline coord3d dihedral_gradient()    const { return inner_dihedral_gradient() + outer_a_dihedral_gradient() + outer_m_dihedral_gradient() + outer_p_dihedral_gradient();}
        //inline coord3d flatness()             const { return ;  }   
        

        //Harmonic energy contribution from bond stretching, angular bending and dihedral angle bending.
        inline real_t energy() const {return 0.5*f_bond *harmonic_energy(bond(),r0)+ f_inner_angle* harmonic_energy(angle(),angle0) +f_inner_dihedral* harmonic_energy(dihedral(),inner_dih0);}
        //Sum of bond, angular and dihedral gradient components.
        inline coord3d gradient() const{return bond_length_gradient() + angle_gradient() + dihedral_gradient();}


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
            r_rad,
            r_ram,
            r_rap;

        //Equillibrium parameters.
        real_t
            r0,
            angle0,
            outer_angle_m0,
            outer_angle_p0,
            inner_dih0,
            outer_dih0_a,
            outer_dih0_m,
            outer_dih0_p;

        //Base Arcs,
        coord3d
            ab,
            ac,
            ad;

        /*
        All normalized arcs required to perform energy & gradient calculations.
        Note that all these arcs are cyclical the arc ab becomes: ab->ac->ad,  the arc ac becomes: ac->ad->ab , the arc bc becomes: bc->cd->db (For iterations 0, 1, 2)
        As such the naming convention here is related to the arcs as they are used in the 0th iteration. */
        coord3d 
            ab_hat,
            ac_hat,
            ad_hat,
            bp_hat,
            bm_hat,
            am_hat,
            ap_hat,
            ba_hat,
            bc_hat,
            cd_hat,
            mp_hat,
            mb_hat,
            pa_hat,
            pb_hat;
    };

    FullereneForcefield(const node_t* neighbours, coord3d* X, const uint8_t* face_right, const node_t* next_on_face, const node_t* prev_on_face, size_t N) : 
        neighbours(neighbours), X(X), face_right(face_right), next_on_face(next_on_face), prev_on_face(prev_on_face), N(N) {
            X_temp = new coord3d[N];

        }
    ~FullereneForcefield(){free(X_temp);}

    //Parallelizable copying,  a = b;
    template <typename T>
    inline void parallel_copy(T* a, T* b){
        for (node_t i = 0; i < N; i++)
        {
            a[i]= b[i];
        }
    }
       
    
    //~393*N FLOPs
    real_t energy(coord3d* X){
        real_t energy = 0.0;
        for (node_t a = 0; a < N; a++)
        {   
            for (size_t j = 0; j < 3; j++)
            {   
                ArcData arc = ArcData(a,j,neighbours,X,face_right,next_on_face,prev_on_face);
                energy += arc.energy();
            }
        }
        return energy;
    }

    //~1914 * N -FLOPs
    void gradient(coord3d* X,coord3d* gradient)
    {   
        //#pragma acc parallel loop copyin(X[:N], neighbours[:3*N], face_right[:3*N], next_on_face[:3*N], prev_on_face[:3*N]) copy(gradient[:N])
        //#pragma omp parallel for
        for (node_t a = 0; a < N; a++) 
        {       
            coord3d node_gradient = {0,0,0};
            
            for (size_t j = 0; j < 3; j++)
            {       
                ArcData arc = ArcData(a,j,neighbours,X,face_right,next_on_face,prev_on_face);
                node_gradient += arc.gradient(); 
            }
            gradient[a] = node_gradient;
        }
        //#pragma acc wait
    }

    size_t golden_section_search(coord3d* X, coord3d* direction, coord3d* new_direction,coord3d* X1, coord3d* X2, real_t a, real_t b, real_t tol){
        real_t tau = (sqrt(5) - 1) / 2;
        
        //Actual coordinates resulting from each traversal 
        //Line search x - values;
        real_t x1,  x2, dfc;
        x1 = (a + (1 - tau) * (b - a));
        x2 = (a + tau * (b - a));
        for (node_t a = 0; a < N; a++)
        {
            X1[a] = X[a] + x1 * direction[a];
            X2[a] = X[a] + x2 * direction[a];
        }
        real_t f1 = energy(X1);
        real_t f2 = energy(X2);
        size_t evals = 2;

        for (uint8_t i = 0; i < 40; i++){
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
                f2 = energy(X2);
            }else
            {
                b = x2;
                x2 = x1;
                f2 = f1;
                x1 = a + (1 - tau) * (b - a);
                for (node_t a = 0; a < N; a++){
                    X1[a] = X[a] + x1 * direction[a];
                }
                f1 = energy(X1);
            }
        }
        //Line search coefficient
        real_t alfa = (a+b)/2;

        for (node_t a = 0; a < N; a++)
        {
            X[a] = X[a] + alfa*direction[a];
        }
        gradient(X,new_direction);
        for (node_t a = 0; a < N; a++){ new_direction[a] = -new_direction[a];}
        return evals;
    }


    void conjugate_gradient(){
        size_t iter_count = 0;
        size_t max_iter = N*3;
        real_t beta = 0.0;
        real_t dnorm = 0.0;
        real_t r0_norm = 0.0;
        real_t direction_norm = 0.0;
        size_t gradient_evals = 0;
        size_t energy_evals = 0;

        coord3d* delta_x0 = new coord3d[N];
        coord3d* delta_x1 = new coord3d[N];
        coord3d* direction = new coord3d[N];
        coord3d* X1 = new coord3d[N]; // Trial coordinates.
        coord3d* X2 = new coord3d[N]; // Trial coordinates.
        //if (X[0][0] > 0) { X[1][1] = 0;}
        gradient(X, direction);
        //if (direction[0][0] > 0) { direction[1][1] = 0;}
        
        gradient_evals ++;
        for (node_t a = 0; a < N; a++)
        {
            dnorm += dot(direction[a],direction[a]);
        }

        dnorm = sqrt(dnorm);
        for (node_t a = 0; a < N; a++)
        {   
            direction[a] = -direction[a]/dnorm;
        }
        parallel_copy(X_temp,X);
        parallel_copy(delta_x0, direction);

        
        for (node_t i = 0; i < max_iter; i++)
        {   
            beta = 0.0; direction_norm = 0.0; dnorm=0.0; r0_norm = 0.0;
            energy_evals += golden_section_search(X_temp,direction, delta_x1, X1, X2, 0, 1, 1e-10);
            gradient_evals++;
            //gradient_evals += bisection_search(X_temp, direction, delta_x1,0,1e-5,1e-10,N);
            //Polak Ribiere method
            for (node_t a = 0; a < N; a++)
            {
                beta += dot(delta_x1[a], (delta_x1[a] - delta_x0[a]));
                r0_norm += dot(delta_x0[a], delta_x0[a]);
            }
            beta /= r0_norm;
        
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

            //Calculate gradient and residual gradient norms..
            for (node_t a = 0; a < N; a++)
            {
                direction_norm += dot(direction[a],direction[a]);
                dnorm += dot(delta_x1[a],delta_x1[a]);
            }
            direction_norm = sqrt(direction_norm);
            dnorm = sqrt(dnorm);

            //Normalize gradient.
            for (node_t a = 0; a < N; a++)
            {
                direction[a] /= direction_norm;
            }
            iter_count++;
        }
        free(delta_x0); free(delta_x1); free(direction); free(X1); free(X2);
        
        
    }
};

size_t computeBatchSize(size_t N){
    return (size_t)1;
}

void OptimizeBatch(real_t* h_X, node_t* h_cubic_neighbours, node_t* h_next_on_face, node_t* h_prev_on_face, uint8_t* h_face_right, const size_t N, const size_t batch_size){  
    //Test gradient computation
    FullereneForcefield forcefield = FullereneForcefield(h_cubic_neighbours, reinterpret_cast<coord3d*>(h_X), h_face_right, h_next_on_face, h_prev_on_face, N);

    auto start = chrono::system_clock::now();
    forcefield.conjugate_gradient();
    auto end = chrono::system_clock::now();
    cout << "Elapsed time: " << (end-start)/ 1ms << "ms\n" ;
}

};