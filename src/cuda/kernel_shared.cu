
#include "cuda_runtime.h"
#include <cuda_runtime_api.h>
#include <cuda.h>
#include "device_launch_parameters.h"
#include <stdio.h>
#include <cooperative_groups.h>
#include <cooperative_groups/reduce.h>
#include <helper_cuda.h>
#include <iostream>
#include <fstream>
#include <chrono>
#include "C60ih.cu"
#include "C960ih.cu"
#include "coord3d.cu"

using namespace std::literals;
namespace cg = cooperative_groups;

typedef uint16_t node_t; 
typedef double3 coord3d;


//All parameter arrays are indexed by a binary sum, 0,1,2,3,4,...
//Pentagons = 0
//Hexagons = 1
//PPP = 0, {HPP, PHP, PPH} = 1, {PHH, HPH, HHP} = 2, {HHH} = 3
__constant__ real_t optimal_corner_cos_angles[2] = {-0.3090169944, -0.5}; 
__constant__ real_t optimal_bond_lengths[3] = {1.479, 1.458, 1.401}; 
__constant__ real_t optimal_dih_cos_angles[4] = {0.511161369, 0.6543772853, -0.3342353205, 1}; 

__constant__ real_t angle_forces[2] = {100.0, 100.0}; 
__constant__ real_t bond_forces[3] = {260.0, 390.0, 450.0}; 
__constant__ real_t dih_forces[4] = {35.0, 65.0, 85.0, 270.0}; 


__device__ __host__ struct BookkeepingData{
    const node_t* neighbours;
    const uint8_t* face_right;
    const node_t* next_on_face;
    const node_t* prev_on_face;
    __device__ __host__ BookkeepingData (const node_t* neighbours, const uint8_t* face_right, const node_t* next_on_face, const node_t* prev_on_face) : 
        neighbours(neighbours), face_right(face_right), next_on_face(next_on_face), prev_on_face(prev_on_face) {}
};

__device__ struct Constants{
    const real_t* f_bond;
    const real_t* f_inner_angle;
    const real_t* f_inner_dihedral;
    const real_t* f_outer_angle_m;
    const real_t* f_outer_angle_p;
    const real_t* f_outer_dihedral;

    const real_t* r0;
    const real_t* angle0;
    const real_t* outer_angle_m0;
    const real_t* outer_angle_p0;
    const real_t* inner_dih0;
    const real_t* outer_dih0;

    __device__ Constants(){}
    __device__ Constants(const real_t* f_bond, const real_t* f_inner_angle, const real_t* f_inner_dihedral, const real_t* f_outer_angle_m, const real_t* f_outer_angle_p, const real_t* f_outer_dihedral,
                            const real_t* r0, const real_t* angle0, const real_t* outer_angle_m0, const real_t* outer_angle_p0, const real_t* inner_dih0, const real_t* outer_dih0): f_bond(f_bond), f_inner_angle(f_inner_angle),
                            f_inner_dihedral(f_inner_dihedral), f_outer_angle_m(f_outer_angle_m), f_outer_angle_p(f_outer_angle_p), f_outer_dihedral(f_outer_dihedral), r0(r0), angle0(angle0), outer_angle_m0(outer_angle_m0), outer_angle_p0(outer_angle_p0),
                            inner_dih0(inner_dih0), outer_dih0(outer_dih0) {}

};

struct DevicePointers{
    coord3d* X;
    coord3d* X_temp;
    coord3d* X1;
    coord3d* X2;
    __device__ __host__ DevicePointers (coord3d* X, coord3d* X_temp, coord3d* X1, coord3d* X2) : 
        X(X), X_temp(X_temp), X1(X1), X2(X2) {}
};



__device__ struct ArcData{
    __device__ __forceinline__ ArcData(const node_t a, const uint8_t j, const coord3d* X, const BookkeepingData& bdat, const Constants& constants){   
        real_t r_rmp;
        coord3d ap, am, ab, ac, ad, mp;
        //printf("Index: %d \n", a*3 + j);
        this->j = j;
        this->c = constants;

        //Compute the arcs ab, ac, ad, bp, bm, ap, am, mp, bc and cd
        ab = (X[bdat.neighbours[a*3 + j]] - X[a]);  r_rab = bond_length(ab); ab_hat = r_rab * ab;
        ac = (X[bdat.neighbours[a*3 + (j+1)%3]] - X[a]); r_rac = bond_length(ac); ac_hat = r_rac * ac;
        ad = (X[bdat.neighbours[a*3 + (j+2)%3]] - X[a]); r_rad = bond_length(ad); ad_hat = r_rad * ad;
        
        coord3d bp = (X[bdat.next_on_face[a*3 + j]] - X[bdat.neighbours[a*3 + j]]); bp_hat = unit_vector(bp);
        coord3d bm = (X[bdat.prev_on_face[a*3 + j]] - X[bdat.neighbours[a*3 + j]]); bm_hat = unit_vector(bm);

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
    }
    __device__ __forceinline__ real_t harmonic_energy(const real_t p0, const real_t p) const{
        return 0.5*(p-p0)*(p-p0);
    }
    __device__ __forceinline__ coord3d  harmonic_energy_gradient(const real_t p0, const real_t p, const coord3d gradp) const{
        return (p-p0)*gradp;     
    }

    __device__ __forceinline__ real_t bond() const {return 1/r_rab;}
    __device__ __forceinline__ real_t angle() const {return dot(ab_hat,ac_hat);}
    //Returns the inner dihedral angle for the current arc. Used here only for energy calculation, 
    //otherwise embedded in dihedral computation because the planes and angles that make up the dihedral angle computation are required for derivative computation.
    __device__ __forceinline__ real_t dihedral() const 
    { 
        coord3d nabc, nbcd; real_t cos_b, cos_c, r_sin_b, r_sin_c;
        cos_b = dot(ba_hat,bc_hat); r_sin_b = rsqrt(1 - cos_b*cos_b); nabc = cross(ba_hat, bc_hat) * r_sin_b;
        cos_c = dot(-bc_hat,cd_hat); r_sin_c = rsqrt(1 - cos_c*cos_c); nbcd = cross(-bc_hat,cd_hat) * r_sin_c;
        return dot(nabc, nbcd);
    }
    
    // Chain rule terms for angle calculation
    //Computes gradient related to bending term. ~24 FLOPs
    __device__ coord3d inner_angle_gradient() const
    {
        real_t cos_angle = angle(); //Inner angle of arcs ab,ac.
        coord3d grad = cos_angle * (ab_hat * r_rab + ac_hat * r_rac) - ab_hat * r_rac - ac_hat* r_rab; //Derivative of inner angle: Eq. 21. 
        return c.f_inner_angle[j] * harmonic_energy_gradient(c.angle0[j], cos_angle, grad); //Harmonic Energy Gradient: Eq. 21. multiplied by harmonic term.
    }
    //Computes gradient related to bending of outer angles. ~20 FLOPs
    __device__ coord3d outer_angle_gradient_m() const
    {
        real_t cos_angle = -dot(ab_hat, bm_hat); //Compute outer angle. ab,bm
        coord3d grad = (bm_hat + ab_hat * cos_angle) * r_rab; //Derivative of outer angles Eq. 30. Buster Thesis
        return c.f_outer_angle_m[j] * harmonic_energy_gradient(c.outer_angle_m0[j],cos_angle,grad); //Harmonic Energy Gradient: Eq. 30 multiplied by harmonic term.
    }
    __device__ coord3d outer_angle_gradient_p() const
    {
        real_t cos_angle = -dot(ab_hat, bp_hat); //Compute outer angle. ab,bp
        coord3d grad = (bp_hat + ab_hat * cos_angle) * r_rab; //Derivative of outer angles Eq. 28. Buster Thesis
        return c.f_outer_angle_p[j] * harmonic_energy_gradient(c.outer_angle_p0[j],cos_angle,grad); //Harmonic Energy Gradient: Eq. 28 multiplied by harmonic term.
    }
    // Chain rule terms for dihedral calculation
    //Computes gradient related to dihedral/out-of-plane term. ~75 FLOPs
    __device__ coord3d inner_dihedral_gradient() const
    {
        coord3d nabc, nbcd; real_t cos_b, cos_c, r_sin_b, r_sin_c;
        cos_b = dot(ba_hat,bc_hat); r_sin_b = rsqrtf(1 - cos_b*cos_b); nabc = cross(ba_hat, bc_hat) * r_sin_b;
        cos_c = dot(-bc_hat,cd_hat); r_sin_c = rsqrtf(1 - cos_c*cos_c); nbcd = cross(-bc_hat,cd_hat) * r_sin_c;

        real_t cos_beta = dot(nabc, nbcd); //Inner dihedral angle from planes abc,bcd.
        
        real_t cot_b = cos_b * r_sin_b * r_sin_b; //cos(b)/sin(b)^2

        //Derivative w.r.t. inner dihedral angle F and G in Eq. 26
        coord3d grad = cross(bc_hat, nbcd) * r_sin_b * r_rab - ba_hat * cos_beta * r_rab + (cot_b * cos_beta * r_rab) * (bc_hat - ba_hat * cos_b);

        return c.f_inner_dihedral[j] * harmonic_energy_gradient(c.inner_dih0[j], cos_beta, grad); //Eq. 26.
    }

    //Computes gradient from dihedral angles constituted by the planes bam, amp ~162 FLOPs
    __device__ coord3d outer_a_dihedral_gradient() const
    {
        coord3d nbam_hat, namp_hat; real_t cos_a, cos_m, r_sin_a, r_sin_m;

        cos_a = dot(ab_hat,am_hat); r_sin_a = rsqrtf(1 - cos_a*cos_a); nbam_hat = cross(ab_hat,am_hat) * r_sin_a;
        cos_m = dot(-am_hat,mp_hat); r_sin_m = rsqrtf(1 - cos_m*cos_m); namp_hat = cross(-am_hat,mp_hat) * r_sin_m;
        
        real_t cos_beta = dot(nbam_hat, namp_hat); //Outer Dihedral angle bam, amp
        real_t cot_a = cos_a * r_sin_a * r_sin_a;
        real_t cot_m = cos_m * r_sin_m * r_sin_m;

        //Derivative w.r.t. outer dihedral angle, factorized version of Eq. 31.
        coord3d grad = cross(mp_hat,nbam_hat)*r_ram*r_sin_m - (cross(namp_hat,ab_hat)*r_ram + cross(am_hat,namp_hat)*r_rab)*r_sin_a +
                        cos_beta*(ab_hat*r_rab + r_ram * (2*am_hat + cot_m*(mp_hat+cos_m*am_hat)) - cot_a*(r_ram*(ab_hat - am_hat*cos_a) + r_rab*(am_hat-ab_hat*cos_a)));
        
        //Eq. 31 multiplied by harmonic term.
        return c.f_outer_dihedral[j] * harmonic_energy_gradient(c.outer_dih0[j], cos_beta, grad);
    }

    //Computes gradient from dihedral angles constituted by the planes nbmp, nmpa ~92 FLOPs
    __device__ coord3d outer_m_dihedral_gradient() const
    {
        coord3d nbmp_hat, nmpa_hat; real_t cos_m, cos_p, r_sin_m, r_sin_p;
        cos_m = dot(mb_hat,mp_hat); r_sin_m = rsqrtf(1 - cos_m*cos_m); nbmp_hat = cross(mb_hat,mp_hat) * r_sin_m;
        cos_p = dot(-mp_hat,pa_hat); r_sin_p = rsqrtf(1 - cos_p*cos_p); nmpa_hat = cross(-mp_hat,pa_hat) * r_sin_p;
        
        //Cosine to the outer dihedral angle constituted by the planes bmp and mpa
        real_t cos_beta = dot(nbmp_hat, nmpa_hat); //Outer dihedral angle bmp,mpa.
        real_t cot_p = cos_p * r_sin_p * r_sin_p;
        
        //Derivative w.r.t. outer dihedral angle, factorized version of Eq. 32.
        coord3d grad = r_rap * (cot_p*cos_beta * (-mp_hat - pa_hat*cos_p) - cross(nbmp_hat, mp_hat)*r_sin_p - pa_hat*cos_beta );

        //Eq. 32 multiplied by harmonic term.
        return c.f_outer_dihedral[j] * harmonic_energy_gradient(c.outer_dih0[j], cos_beta, grad);
    }

    //Computes gradient from dihedral angles constituted by the planes bpa, pam ~162 FLOPs
    __device__ coord3d outer_p_dihedral_gradient() const
    {
        coord3d nbpa_hat, npam_hat; real_t cos_p, cos_a, r_sin_p, r_sin_a;
        cos_a = dot(ap_hat,am_hat); r_sin_a = rsqrtf(1 - cos_a*cos_a); npam_hat = cross(ap_hat,am_hat) * r_sin_a;
        cos_p = dot(pb_hat,-ap_hat); r_sin_p = rsqrtf(1 - cos_p*cos_p); nbpa_hat = cross(pb_hat,-ap_hat) * r_sin_p;

        real_t cos_beta = dot(nbpa_hat, npam_hat); //Outer dihedral angle bpa, pam.
        real_t cot_p = cos_p * r_sin_p * r_sin_p;
        real_t cot_a = cos_a * r_sin_a * r_sin_a;

        //Derivative w.r.t. outer dihedral angle, factorized version of Eq. 33.
        coord3d grad = cross(npam_hat,pb_hat)*r_rap*r_sin_p - (cross(am_hat,nbpa_hat)*r_rap + cross(nbpa_hat,ap_hat)*r_ram)*r_sin_a +
                        cos_beta*(am_hat*r_ram + r_rap * (2*ap_hat + cot_p*(pb_hat+cos_p*ap_hat)) - cot_a*(r_rap*(am_hat - ap_hat*cos_a) + r_ram*(ap_hat-am_hat*cos_a)));
        
        //Eq. 33 multiplied by harmonic term.
        return c.f_outer_dihedral[j] * harmonic_energy_gradient(c.outer_dih0[j], cos_beta, grad);
    }
    // Internal coordinate gradients
    __device__ coord3d bond_length_gradient() const { return - c.f_bond[j] * harmonic_energy_gradient(c.r0[j],bond(),ab_hat);}
    //Sum of angular gradient components.
    __device__ coord3d angle_gradient() const { return inner_angle_gradient() + outer_angle_gradient_p() + outer_angle_gradient_m();}
    //Sum of inner and outer dihedral gradient components.
    __device__ coord3d dihedral_gradient() const { return inner_dihedral_gradient() + outer_a_dihedral_gradient() + outer_m_dihedral_gradient() + outer_p_dihedral_gradient();}
    //coord3d flatness()             const { return ;  }   
    

    //Harmonic energy contribution from bond stretching, angular bending and dihedral angle bending.
    __device__ __forceinline__ real_t energy() const {return 0.5*c.f_bond[j] *harmonic_energy(bond(),c.r0[j])+c.f_inner_angle[j]* harmonic_energy(angle(),c.angle0[j])+c.f_inner_dihedral[j]* harmonic_energy(dihedral(),c.inner_dih0[j]);}
    //Sum of bond, angular and dihedral gradient components.
    __device__ coord3d gradient() const{return bond_length_gradient() + angle_gradient() + dihedral_gradient();}
    
    uint8_t j; Constants c;

    //Residual lengths of arcs ab, ac, am, ap.
    real_t
        r_rab,
        r_rac,
        r_rad,
        r_ram,
        r_rap;

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

__device__ void reduction(real_t *sdata, const node_t N, const node_t node_id, const cg::thread_block& block){
    
    cg::sync(block);
    cg::thread_block_tile<32> tile32 = cg::tiled_partition<32>(block);
    sdata[node_id] = cg::reduce(tile32, sdata[node_id], cg::plus<real_t>());
    cg::sync(block);

    real_t beta = 0.0;
    if (block.thread_rank() == 0) {
        beta  = 0;
        for (uint16_t i = 0; i < N; i += tile32.size()) {
            beta  += sdata[i];
        }
        sdata[0] = beta;
    }
    cg::sync(block);
}



__device__ coord3d gradient(const coord3d* X, const node_t node_id, const BookkeepingData &dat, const Constants &constants) {
    coord3d grad = {0.0, 0.0, 0.0};
    #pragma unroll (3)
    for (uint8_t j = 0; j < 3; j++ ){
        ArcData arc = ArcData::ArcData(node_id, j, X, dat, constants);
        grad += arc.gradient();
    }
    return grad;
}

__device__ real_t energy(const coord3d* X, const node_t node_id, const BookkeepingData &dat, const Constants &constants, real_t* reduction_array, node_t N) {
    reduction_array[node_id] = 0.0;
    #pragma unroll (3)
    for (uint8_t j = 0; j < 3; j++ ){
        ArcData arc = ArcData::ArcData(node_id, j, X, dat, constants);
        reduction_array[node_id] += arc.energy();
    }
    reduction(reduction_array, N, node_id, cg::this_thread_block());
    return reduction_array[0];
}

__device__ void golden_section_search(coord3d* X, coord3d& direction, coord3d& new_direction,coord3d* X1, coord3d* X2, real_t* reduction_array, real_t a, real_t b, const node_t node_id, const node_t N, const BookkeepingData &dat, const Constants& constants){
    real_t tau = (sqrtf(5) - 1) / 2;
    cg::thread_block block = cg::this_thread_block();
    //Actual coordinates resulting from each traversal 
    //Line search x - values;
    real_t x1,  x2, dfc;
    x1 = (a + (1 - tau) * (b - a));
    x2 = (a + tau * (b - a));

    X1[node_id] = X[node_id] + x1 * direction;
    X2[node_id] = X[node_id] + x2 * direction;
    cg::sync(block);

    real_t f1 = energy(X1, node_id, dat, constants, reduction_array, N);
    real_t f2 = energy(X2, node_id, dat, constants, reduction_array, N);

    for (uint8_t i = 0; i < 40; i++){
        if (f1 > f2){
            a = x1;
            x1 = x2;
            f1 = f2;
            x2 = a + tau * (b - a);
            X2[node_id] = X[node_id] + x2 * direction;
            __threadfence();
            f2 = energy(X2, node_id, dat, constants, reduction_array, N);
        }else
        {
            b = x2;
            x2 = x1;
            f2 = f1;
            x1 = a + (1 - tau) * (b - a);
            X1[node_id] = X[node_id] + x1 * direction;
            __threadfence();
            f1 = energy(X1, node_id, dat, constants, reduction_array, N);
        }
    }
    //Line search coefficient
    real_t alfa = (a+b)/2;
    X[node_id] = X[node_id] + alfa*direction;
    new_direction = gradient(X,node_id,dat, constants);
    new_direction = -new_direction;
}

__device__ Constants compute_constants(BookkeepingData &dat, node_t node_id){
    real_t* r0 = new real_t[3]; real_t* angle0 = new real_t[3]; real_t* inner_dih0 = new real_t[3]; real_t* outer_angle_m0 = new real_t[3]; real_t* outer_angle_p0 = new real_t[3]; real_t* outer_dih0 = new real_t[3];
    real_t* f_bond = new real_t[3]; real_t* f_inner_angle = new real_t[3]; real_t* f_inner_dihedral = new real_t[3]; real_t* f_outer_angle_m = new real_t[3]; real_t* f_outer_angle_p = new real_t[3]; real_t* f_outer_dihedral = new real_t[3];

    for (uint8_t j = 0; j < 3; j++) {
        uint8_t f_r = dat.face_right[node_id * 3 + j] - 5;
        uint8_t f_l = dat.face_right[node_id * 3 + (2 + j)%3] - 5;

        uint8_t face_sum = dat.face_right[node_id * 3] - 5 + dat.face_right[node_id * 3 + 1] - 5 + dat.face_right[node_id * 3 + 2] - 5;
        uint8_t dihedral_face_sum = dat.face_right[dat.neighbours[node_id * 3 + j] * 3]-5 + dat.face_right[dat.neighbours[node_id * 3 + j] * 3 + 1]-5 +  dat.face_right[dat.neighbours[node_id * 3 + j] * 3 + 2]-5;

        //Load equillibirium distance, angles and dihedral angles from face information.
        r0[j] = optimal_bond_lengths[ f_l + f_r ];
        angle0[j] = optimal_corner_cos_angles[ f_r ];
        inner_dih0[j] = optimal_dih_cos_angles[ face_sum ];
        outer_angle_m0[j] = optimal_corner_cos_angles[ f_l ];
        outer_angle_p0[j] = optimal_corner_cos_angles[ f_r ];
        outer_dih0[j] = optimal_dih_cos_angles[ dihedral_face_sum ];

        //Load force constants from neighbouring face information.
        f_bond[j] = bond_forces[ f_l + f_r ];
        f_inner_angle[j] = angle_forces[ f_l ];
        f_inner_dihedral[j] = dih_forces[ face_sum];
        f_outer_angle_m[j] = angle_forces[ f_r ];
        f_outer_angle_p[j] = angle_forces[ f_l ];
        f_outer_dihedral[j] = dih_forces[ dihedral_face_sum];
    }

    return Constants::Constants(f_bond,f_inner_angle,f_inner_dihedral, f_outer_angle_m, f_outer_angle_p, f_outer_dihedral, r0, angle0, outer_angle_m0, outer_angle_p0, inner_dih0, outer_dih0);
}

__global__ void conjugate_gradient(coord3d* d_X, coord3d* d_X_temp, coord3d* d_X1, coord3d* d_X2, node_t* d_neighbours, node_t* d_next_on_face, node_t* d_prev_on_face, uint8_t* d_face_right, size_t N){
    real_t __shared__ reduction_array[1024];
    coord3d __shared__ sX[400];
    coord3d __shared__ sX_temp[400];
    coord3d __shared__ sX1[400];
    coord3d __shared__ sX2[400];

    coord3d delta_x0, delta_x1, direction;

    size_t iter_count = 0;
    size_t max_iter = N*10;
    size_t gradient_evals = 0;
    size_t energy_evals = 0;
    
    real_t beta, dnorm, r0_norm, direction_norm;
    beta = dnorm = r0_norm = direction_norm = 0.0;

    node_t block_offset = 0;
    node_t node_id = threadIdx.x - block_offset;
    size_t offset = blockIdx.x * blockDim.x + block_offset;
    size_t tid = blockIdx.x * blockDim.x + threadIdx.x;
    cg::thread_block block = cg::this_thread_block();
    
    real_t* local_reduction = &reduction_array[block_offset];
    
    coord3d* X = &d_X[offset];
    coord3d* X_temp = &d_X_temp[offset];
    coord3d* X1 = &d_X1[offset];
    coord3d* X2 = &d_X2[offset];
    
    const node_t* neighbours = &d_neighbours[3*offset];
    const uint8_t* face_right = &d_face_right[3*offset];
    const node_t* next_on_face = &d_next_on_face[3*offset];
    const node_t* prev_on_face = &d_prev_on_face[3*offset];
    //Pre-compute force constants and store in registers.
    ///===============================================///
    
    BookkeepingData local_bookkeeping = BookkeepingData::BookkeepingData(neighbours,face_right,next_on_face,prev_on_face);   
    DevicePointers local_aben = DevicePointers::DevicePointers(X,X_temp,X1,X2);
    Constants constants = compute_constants(local_bookkeeping, node_id);

    cg::sync(block);
    direction = gradient(X, node_id ,local_bookkeeping, constants);
    gradient_evals ++;
    
    
    local_reduction[node_id] = dot(direction,direction);
    reduction(local_reduction, N, node_id, cg::this_thread_block());
    dnorm = sqrtf(local_reduction[0]);
    direction = -direction/dnorm;
    
    sX[node_id] = X[node_id]; sX_temp[node_id] = sX[node_id];
    delta_x0 = direction;
    real_t test = energy(sX_temp, node_id, local_bookkeeping, constants, local_reduction, N);
    if (tid == 0)
    {
        printf("%e \n", test);/* code */
    }
    
    
    for (node_t i = 0; i < (node_t)2.6*N; i++)
    {   
        beta = 0.0; direction_norm = 0.0; dnorm=0.0; r0_norm = 0.0;
        cg::sync(block);
        golden_section_search(sX_temp, direction, delta_x1, sX1, sX2,local_reduction, 0, 1, node_id, N, local_bookkeeping, constants);
        __threadfence();
        gradient_evals++;
        energy_evals += 42;
        //Polak Ribiere method
        local_reduction[node_id] = dot(delta_x0, delta_x0); reduction(local_reduction, N, node_id, cg::this_thread_block()); r0_norm = local_reduction[0];
        local_reduction[node_id] = dot(delta_x1, (delta_x1 - delta_x0)); reduction(local_reduction, N, node_id, cg::this_thread_block()); beta = local_reduction[0] / r0_norm;

    
        if (energy(sX_temp, node_id, local_bookkeeping, constants, local_reduction, N) > energy(sX, node_id, local_bookkeeping, constants, local_reduction, N))
        {   
            sX_temp[node_id] =  sX[node_id];
            delta_x1 =  delta_x0;
            beta = 0.0;
        }
        else
        {   
            sX[node_id] = sX_temp[node_id];
            delta_x0 = delta_x1;
        }
        direction = delta_x1 + beta*direction;
        

        //Calculate gradient and residual gradient norms..

        local_reduction[node_id] = dot(direction,direction); reduction(local_reduction, N, node_id, cg::this_thread_block()); direction_norm = sqrtf(local_reduction[0]);
        local_reduction[node_id] = dot(delta_x1,delta_x1); reduction(local_reduction, N, node_id, cg::this_thread_block()); dnorm = sqrtf(local_reduction[0]);

        //Normalize gradient.
        direction /= direction_norm;
        if (tid == 0) {
            //printf("Norm: %e \n ", dnorm);
        }
        
        if (iter_count > N*10)
        {
            if (tid == 0){
            printf("Conjugated Gradient Terminated Due to Max Iterations : %d \n", N*10);
            printf("Gradient Evaluations: %d \n ", gradient_evals);
            printf("Energy Evaluations: %d \n", energy_evals);
            }
            break;
        }
        iter_count++;
    }
    cg::sync(block);
    real_t end_energy = energy(sX, node_id, local_bookkeeping, constants, local_reduction, N);
    if (tid == 0){
        printf("Conjugated Gradient Terminated Due to Max Iterations : %d \n", iter_count);
        printf("Gradient Evaluations: %d \n ", gradient_evals);
        printf("Energy Evaluations: %d \n", energy_evals);
        printf("Energy At End: %e \n", end_energy);
    }

}

int main(){

    size_t N = 60;

    size_t* d_N;

    coord3d* h_X = new coord3d[N];

    coord3d* d_X;
    coord3d* d_X_temp;
    coord3d* d_X1;
    coord3d* d_X2;
    coord3d* d_delta_x0;
    coord3d* d_delta_x1;
    coord3d* d_direction;


    node_t* d_neighbours;
    uint8_t* d_face_right;
    node_t* d_next_on_face;
    node_t* d_prev_on_face;

    cudaMalloc(&d_X, sizeof(real_t)*3*N);
    cudaMalloc(&d_X_temp, sizeof(real_t)*3*N);
    cudaMalloc(&d_X1, sizeof(real_t)*3*N);
    cudaMalloc(&d_X2, sizeof(real_t)*3*N);
    
    cudaMalloc(&d_neighbours, sizeof(node_t)*3*N);
    cudaMalloc(&d_next_on_face, sizeof(node_t)*3*N);
    cudaMalloc(&d_prev_on_face, sizeof(node_t)*3*N);
    cudaMalloc(&d_face_right, sizeof(uint8_t)*3*N);
    cudaMalloc(&d_N, sizeof(size_t)); cudaMemcpy(d_N, &N, sizeof(size_t), cudaMemcpyHostToDevice);

    cudaMemcpy(d_X, &X_60, sizeof(real_t)*3*N , cudaMemcpyHostToDevice);
    cudaMemcpy(d_neighbours, &cubic_neighbours_60, sizeof(node_t)*3*N, cudaMemcpyHostToDevice);
    cudaMemcpy(d_next_on_face, &next_on_face_60, sizeof(node_t)*3*N, cudaMemcpyHostToDevice);
    cudaMemcpy(d_prev_on_face, &prev_on_face_60, sizeof(node_t)*3*N, cudaMemcpyHostToDevice);
    cudaMemcpy(d_face_right, &face_right_60, sizeof(uint8_t)*3*N, cudaMemcpyHostToDevice);

    DevicePointers dpointers = DevicePointers::DevicePointers(d_X, d_X_temp, d_X1, d_X2);
    BookkeepingData bpointers = BookkeepingData::BookkeepingData(d_neighbours, d_face_right, d_next_on_face, d_prev_on_face);

    void *kernelArgs[] = {
        (void*)&d_X,
        (void*)&d_X_temp,
        (void*)&d_X1,
        (void*)&d_X2,
        (void*)&d_neighbours,
        (void*)&d_next_on_face,
        (void*)&d_prev_on_face,
        (void*)&d_face_right,
        (void*)&N
    };
    dim3 dimBlock(N, 1, 1);
    dim3 dimGrid(1, 1, 1);

    auto start = std::chrono::system_clock::now();
    checkCudaErrors(cudaLaunchCooperativeKernel((void*)conjugate_gradient, dimGrid, dimBlock, kernelArgs, NULL, NULL));
    cudaDeviceSynchronize();
    auto end = std::chrono::system_clock::now();
    std::cout << "Elapsed time: " << (end-start)/ 1ms << "ms\n" ;
    getLastCudaError("Failed to launch kernel: ");
    cudaMemcpy(h_X, d_X, sizeof(real_t)*3*N, cudaMemcpyDeviceToHost);



    for (int i = 0; i<N; i++){
       //print_coord(h_X[i]);
    }

    

}


