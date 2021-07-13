#include "fullerenes/gpu/isomerspace_forcefield.hh"
#include "cuda_runtime.h"
#include <cuda_runtime_api.h>
#include <cuda.h>
#include "device_launch_parameters.h"
#include <stdio.h>
#//include <helper_cuda.h>  // TODO: Get rid of this (not in nvcc std include dirs)
#define getLastCudaError(x) 
#include <iostream>
#include <fstream>
#include <chrono>

namespace IsomerspaceForcefield {

typedef device_real_t real_t;
typedef device_node_t node_t;

#include "C60ih.cu"
#include "coord3d.cu"
#include "helper_functions.cu"

using namespace std::literals;
namespace cg = cooperative_groups;



__device__ struct ArcData{
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
    
    __device__ ArcData(const node_t a, const uint8_t j, const node_t* neighbours, const coord3d* X, const uint8_t* face_right, const node_t* next_on_face, const node_t* prev_on_face){   
        real_t r_rmp;
        coord3d ap, am, ab, ac, ad, mp;
        //printf("Index: %d \n", a*3 + j);

        //Compute the arcs ab, ac, ad, bp, bm, ap, am, mp, bc and cd
        ab = (X[neighbours[a*3 + j]] - X[a]);  r_rab = bond_length(ab); ab_hat = r_rab * ab;
        ac = (X[neighbours[a*3 + (j+1)%3]] - X[a]); r_rac = bond_length(ac); ac_hat = r_rac * ac;
        ad = (X[neighbours[a*3 + (j+2)%3]] - X[a]); r_rad = bond_length(ad); ad_hat = r_rad * ad;
        
        coord3d bp = (X[next_on_face[a*3 + j]] - X[neighbours[a*3 + j]]); bp_hat = unit_vector(bp);
        coord3d bm = (X[prev_on_face[a*3 + j]] - X[neighbours[a*3 + j]]); bm_hat = unit_vector(bm);

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
    __device__ real_t harmonic_energy(const real_t p0, const real_t p) const{
        return (real_t)0.5*(p-p0)*(p-p0);
    }
    __device__ __forceinline__ coord3d  harmonic_energy_gradient(const real_t p0, const real_t p, const coord3d gradp) const{
        return (p-p0)*gradp;     
    }

    __device__ real_t bond() const {return (real_t)1.0/r_rab;}
    __device__ real_t angle() const {return dot(ab_hat,ac_hat);}
    //Returns the inner dihedral angle for the current arc. Used here only for energy calculation, 
    //otherwise embedded in dihedral computation because the planes and angles that make up the dihedral angle computation are required for derivative computation.
    __device__ __forceinline__ real_t dihedral() const 
    { 
        coord3d nabc, nbcd; real_t cos_b, cos_c, r_sin_b, r_sin_c;
        cos_b = dot(ba_hat,bc_hat); r_sin_b = rsqrt((real_t)1.0 - cos_b*cos_b); nabc = cross(ba_hat, bc_hat) * r_sin_b;
        cos_c = dot(-bc_hat,cd_hat); r_sin_c = rsqrt((real_t)1.0 - cos_c*cos_c); nbcd = cross(-bc_hat,cd_hat) * r_sin_c;
        return dot(nabc, nbcd);
    }
    
    // Chain rule terms for angle calculation
    //Computes gradient related to bending term. ~24 FLOPs
    __device__ coord3d inner_angle_gradient() const
    {
        real_t cos_angle = angle(); //Inner angle of arcs ab,ac.
        coord3d grad = cos_angle * (ab_hat * r_rab + ac_hat * r_rac) - ab_hat * r_rac - ac_hat* r_rab; //Derivative of inner angle: Eq. 21. 
        return f_inner_angle * harmonic_energy_gradient(angle0, cos_angle, grad); //Harmonic Energy Gradient: Eq. 21. multiplied by harmonic term.
    }
    //Computes gradient related to bending of outer angles. ~20 FLOPs
    __device__ coord3d outer_angle_gradient_m() const
    {
        real_t cos_angle = -dot(ab_hat, bm_hat); //Compute outer angle. ab,bm
        coord3d grad = (bm_hat + ab_hat * cos_angle) * r_rab; //Derivative of outer angles Eq. 30. Buster Thesis
        return f_outer_angle_m * harmonic_energy_gradient(outer_angle_m0,cos_angle,grad); //Harmonic Energy Gradient: Eq. 30 multiplied by harmonic term.
    }
    __device__ coord3d outer_angle_gradient_p() const
    {
        real_t cos_angle = -dot(ab_hat, bp_hat); //Compute outer angle. ab,bp
        coord3d grad = (bp_hat + ab_hat * cos_angle) * r_rab; //Derivative of outer angles Eq. 28. Buster Thesis
        return f_outer_angle_p * harmonic_energy_gradient(outer_angle_p0,cos_angle,grad); //Harmonic Energy Gradient: Eq. 28 multiplied by harmonic term.
    }
    // Chain rule terms for dihedral calculation
    //Computes gradient related to dihedral/out-of-plane term. ~75 FLOPs
    __device__ coord3d inner_dihedral_gradient() const
    {
        coord3d nabc, nbcd; real_t cos_b, cos_c, r_sin_b, r_sin_c;
        cos_b = dot(ba_hat,bc_hat); r_sin_b = rsqrtf((real_t)1.0 - cos_b*cos_b); nabc = cross(ba_hat, bc_hat) * r_sin_b;
        cos_c = dot(-bc_hat,cd_hat); r_sin_c = rsqrtf((real_t)1.0 - cos_c*cos_c); nbcd = cross(-bc_hat,cd_hat) * r_sin_c;

        real_t cos_beta = dot(nabc, nbcd); //Inner dihedral angle from planes abc,bcd.
        
        real_t cot_b = cos_b * r_sin_b * r_sin_b; //cos(b)/sin(b)^2

        //Derivative w.r.t. inner dihedral angle F and G in Eq. 26
        coord3d grad = cross(bc_hat, nbcd) * r_sin_b * r_rab - ba_hat * cos_beta * r_rab + (cot_b * cos_beta * r_rab) * (bc_hat - ba_hat * cos_b);

        return f_inner_dihedral * harmonic_energy_gradient(inner_dih0, cos_beta, grad); //Eq. 26.
    }

    //Computes gradient from dihedral angles constituted by the planes bam, amp ~162 FLOPs
    __device__ coord3d outer_a_dihedral_gradient() const
    {
        coord3d nbam_hat, namp_hat; real_t cos_a, cos_m, r_sin_a, r_sin_m;

        cos_a = dot(ab_hat,am_hat); r_sin_a = rsqrtf((real_t)1.0 - cos_a*cos_a); nbam_hat = cross(ab_hat,am_hat) * r_sin_a;
        cos_m = dot(-am_hat,mp_hat); r_sin_m = rsqrtf((real_t)1.0 - cos_m*cos_m); namp_hat = cross(-am_hat,mp_hat) * r_sin_m;
        
        real_t cos_beta = dot(nbam_hat, namp_hat); //Outer Dihedral angle bam, amp
        real_t cot_a = cos_a * r_sin_a * r_sin_a;
        real_t cot_m = cos_m * r_sin_m * r_sin_m;

        //Derivative w.r.t. outer dihedral angle, factorized version of Eq. 31.
        coord3d grad = cross(mp_hat,nbam_hat)*r_ram*r_sin_m - (cross(namp_hat,ab_hat)*r_ram + cross(am_hat,namp_hat)*r_rab)*r_sin_a +
                        cos_beta*(ab_hat*r_rab + r_ram * ((real_t)2.0*am_hat + cot_m*(mp_hat+cos_m*am_hat)) - cot_a*(r_ram*(ab_hat - am_hat*cos_a) + r_rab*(am_hat-ab_hat*cos_a)));
        
        //Eq. 31 multiplied by harmonic term.
        return f_outer_dihedral * harmonic_energy_gradient(outer_dih0_a, cos_beta, grad);
    }

    //Computes gradient from dihedral angles constituted by the planes nbmp, nmpa ~92 FLOPs
    __device__ coord3d outer_m_dihedral_gradient() const
    {
        coord3d nbmp_hat, nmpa_hat; real_t cos_m, cos_p, r_sin_m, r_sin_p;
        cos_m = dot(mb_hat,mp_hat); r_sin_m = rsqrtf((real_t)1.0 - cos_m*cos_m); nbmp_hat = cross(mb_hat,mp_hat) * r_sin_m;
        cos_p = dot(-mp_hat,pa_hat); r_sin_p = rsqrtf((real_t)1.0 - cos_p*cos_p); nmpa_hat = cross(-mp_hat,pa_hat) * r_sin_p;
        
        //Cosine to the outer dihedral angle constituted by the planes bmp and mpa
        real_t cos_beta = dot(nbmp_hat, nmpa_hat); //Outer dihedral angle bmp,mpa.
        real_t cot_p = cos_p * r_sin_p * r_sin_p;
        
        //Derivative w.r.t. outer dihedral angle, factorized version of Eq. 32.
        coord3d grad = r_rap * (cot_p*cos_beta * (-mp_hat - pa_hat*cos_p) - cross(nbmp_hat, mp_hat)*r_sin_p - pa_hat*cos_beta );

        //Eq. 32 multiplied by harmonic term.
        return f_outer_dihedral * harmonic_energy_gradient(outer_dih0_m, cos_beta, grad);
    }

    //Computes gradient from dihedral angles constituted by the planes bpa, pam ~162 FLOPs
    __device__ coord3d outer_p_dihedral_gradient() const
    {
        coord3d nbpa_hat, npam_hat; real_t cos_p, cos_a, r_sin_p, r_sin_a;
        cos_a = dot(ap_hat,am_hat); r_sin_a = rsqrtf((real_t)1.0 - cos_a*cos_a); npam_hat = cross(ap_hat,am_hat) * r_sin_a;
        cos_p = dot(pb_hat,-ap_hat); r_sin_p = rsqrtf((real_t)1.0 - cos_p*cos_p); nbpa_hat = cross(pb_hat,-ap_hat) * r_sin_p;

        real_t cos_beta = dot(nbpa_hat, npam_hat); //Outer dihedral angle bpa, pam.
        real_t cot_p = cos_p * r_sin_p * r_sin_p;
        real_t cot_a = cos_a * r_sin_a * r_sin_a;

        //Derivative w.r.t. outer dihedral angle, factorized version of Eq. 33.
        coord3d grad = cross(npam_hat,pb_hat)*r_rap*r_sin_p - (cross(am_hat,nbpa_hat)*r_rap + cross(nbpa_hat,ap_hat)*r_ram)*r_sin_a +
                        cos_beta*(am_hat*r_ram + r_rap * ((real_t)2.0*ap_hat + cot_p*(pb_hat+cos_p*ap_hat)) - cot_a*(r_rap*(am_hat - ap_hat*cos_a) + r_ram*(ap_hat-am_hat*cos_a)));
        
        //Eq. 33 multiplied by harmonic term.
        return f_outer_dihedral * harmonic_energy_gradient(outer_dih0_p, cos_beta, grad);
    }
    // Internal coordinate gradients
    __device__ coord3d bond_length_gradient() const { return - f_bond * harmonic_energy_gradient(r0,bond(),ab_hat);}
    //Sum of angular gradient components.
    __device__ coord3d angle_gradient() const { return inner_angle_gradient() + outer_angle_gradient_p() + outer_angle_gradient_m();}
    //Sum of inner and outer dihedral gradient components.
    __device__ coord3d dihedral_gradient() const { return inner_dihedral_gradient() + outer_a_dihedral_gradient() + outer_m_dihedral_gradient() + outer_p_dihedral_gradient();}
    //coord3d flatness()             const { return ;  }   
    

    //Harmonic energy contribution from bond stretching, angular bending and dihedral angle bending.
    __device__ __forceinline__ real_t energy() const {return ((real_t)0.5)*f_bond *harmonic_energy(bond(),r0)+f_inner_angle* harmonic_energy(angle(),angle0)+f_inner_dihedral* harmonic_energy(dihedral(),inner_dih0);}
    //Sum of bond, angular and dihedral gradient components.
    __device__ coord3d gradient() const{ return bond_length_gradient()+ angle_gradient() + dihedral_gradient();}


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

__device__ coord3d gradient(const coord3d* X, const node_t node_id, const BookkeepingData &dat) {
    coord3d grad = {0.0, 0.0, 0.0};
    for (uint8_t j = 0; j < 3; j++ ){
        ArcData arc = ArcData(node_id, j, dat.neighbours, X, dat.face_right, dat.next_on_face, dat.prev_on_face);
        grad += arc.gradient();
    }
    return grad;
}

__device__ real_t energy(const coord3d* X, const node_t node_id, const BookkeepingData &dat, real_t* reduction_array, node_t N) {
    real_t node_energy = (real_t)0.0;

    for (uint8_t j = 0; j < 3; j++ ){
        ArcData arc = ArcData(node_id, j, dat.neighbours, X, dat.face_right, dat.next_on_face, dat.prev_on_face);
        node_energy += arc.energy();
    }

    cg::sync(cg::this_thread_block());
    reduction_array[threadIdx.x] = node_energy;
    reduction(reduction_array, N);
    return reduction_array[0];
}

__device__ void golden_section_search(coord3d* X, coord3d* direction, coord3d* new_direction,coord3d* X1, coord3d* X2, real_t* reduction_array, real_t a, real_t b, const node_t node_id, const node_t N, const BookkeepingData &dat){
    real_t tau = (sqrtf(5) - 1) / 2;
    cg::thread_block block = cg::this_thread_block();
    //Actual coordinates resulting from each traversal 
    //Line search x - values;
    real_t x1,  x2, dfc;
    x1 = (a + (1 - tau) * (b - a));
    x2 = (a + tau * (b - a));

    X1[node_id] = X[node_id] + x1 * direction[node_id];
    X2[node_id] = X[node_id] + x2 * direction[node_id];
    cg::sync(block);

    real_t f1 = energy(X1, node_id, dat, reduction_array, N);
    real_t f2 = energy(X2, node_id, dat, reduction_array, N);

    for (uint8_t i = 0; i < 20; i++){
        if (f1 > f2){
            a = x1;
            x1 = x2;
            f1 = f2;
            x2 = a + tau * (b - a);
            cg::sync(block);
            X2[node_id] = X[node_id] + x2 * direction[node_id];
            cg::sync(block);
            f2 = energy(X2, node_id, dat, reduction_array, N);
        }else
        {
            b = x2;
            x2 = x1;
            f2 = f1;
            x1 = a + (1 - tau) * (b - a);
            cg::sync(block);
            X1[node_id] = X[node_id] + x1 * direction[node_id];
            cg::sync(block);
            f1 = energy(X1, node_id, dat, reduction_array, N);
        }
    }
    //Line search coefficient
    real_t alfa = (a+b)/2;
    cg::sync(block);
    X1[node_id] = X[node_id] + alfa*direction[node_id];
    cg::sync(block);
    new_direction[node_id] = -gradient(X1,node_id,dat);
}

__global__ void conjugate_gradient(coord3d* d_X, coord3d* d_X_temp, coord3d* d_X1, coord3d* d_X2, coord3d* d_delta_x0, coord3d* d_delta_x1, coord3d* d_direction, node_t* d_neighbours, node_t* d_next_on_face, node_t* d_prev_on_face, uint8_t* d_face_right, size_t N){
    extern __shared__ real_t reduction_array[];

    size_t iter_count = 0;
    size_t max_iter = N*3;
    real_t beta = 0.0;
    real_t dnorm = 0;
    real_t r0_norm;
    real_t direction_norm = 0.0;
    size_t gradient_evals = 0;
    size_t energy_evals = 0;
    size_t node_id = threadIdx.x;

    size_t offset = blockIdx.x * blockDim.x;
    size_t tid = blockIdx.x * blockDim.x + threadIdx.x;
    cg::thread_block block = cg::this_thread_block();
    cg::grid_group grid = cg::this_grid();
    
    coord3d* X = &d_X[offset];
    coord3d* X_temp = &d_X_temp[offset];
    coord3d* X1 = &d_X1[offset];
    coord3d* X2 = &d_X2[offset];
    coord3d* delta_x0 = &d_delta_x0[offset];
    coord3d* delta_x1 = &d_delta_x1[offset];
    coord3d* direction = &d_direction[offset];
    
    const node_t* neighbours = &d_neighbours[3*offset];
    const uint8_t* face_right = &d_face_right[3*offset];
    const node_t* next_on_face = &d_next_on_face[3*offset];
    const node_t* prev_on_face = &d_prev_on_face[3*offset];

    BookkeepingData local_bookkeeping = BookkeepingData(neighbours,face_right,next_on_face,prev_on_face);   
    direction[node_id] = gradient(X, node_id ,local_bookkeeping);
    
    gradient_evals ++;
    
    
    reduction_array[threadIdx.x] = dot(direction[node_id],direction[node_id]);
    reduction(reduction_array, N);
    dnorm = sqrtf(reduction_array[0]);
    direction[node_id] = -direction[node_id]/dnorm;
    
    X_temp[node_id] = X[node_id];
    delta_x0[node_id] = direction[node_id];

    for (node_t i = 0; i < max_iter; i++)
    {   
        beta = 0.0; direction_norm = 0.0; dnorm=0.0; r0_norm = 0.0;
        cg::sync(grid);
        golden_section_search(X, direction, delta_x1, X_temp, X2,reduction_array, 0, 1, node_id, N, local_bookkeeping);

        gradient_evals++;
        energy_evals += 22;
        //Polak Ribiere method
        reduction_array[threadIdx.x] = dot(delta_x0[node_id], delta_x0[node_id]); reduction(reduction_array, N); r0_norm = reduction_array[0];
        cg::sync(block);
        reduction_array[threadIdx.x] = dot(delta_x1[node_id], (delta_x1[node_id] - delta_x0[node_id])); reduction(reduction_array, N); beta = reduction_array[0] / r0_norm;
        cg::sync(block);
        if (energy(X_temp, node_id, local_bookkeeping, reduction_array, N) > energy(X, node_id, local_bookkeeping, reduction_array, N))
        {   
            X_temp[node_id] =  X[node_id];
            delta_x1[node_id] =  delta_x0[node_id];
            beta = 0.0;
        }
        else
        {   
            X[node_id] = X_temp[node_id];
            delta_x0[node_id] = delta_x1[node_id];
        }
        direction[node_id] = delta_x1[node_id] + beta*direction[node_id];
        

        //Calculate gradient and residual gradient norms..
        cg::sync(block);
        reduction_array[threadIdx.x] = dot(direction[node_id],direction[node_id]); reduction(reduction_array, N); direction_norm = sqrtf(reduction_array[0]);
        cg::sync(block);
        reduction_array[threadIdx.x] = dot(delta_x1[node_id],delta_x1[node_id]);  reduction(reduction_array, N); dnorm = sqrtf(reduction_array[0]);
        cg::sync(block);
        //Normalize gradient.
        direction[node_id] /= direction_norm;
        iter_count++;
    }
}

size_t computeBatchSize(size_t N){
    cudaDeviceProp properties;
    cudaGetDeviceProperties(&properties,0);

    /** Compiling with --maxrregcount=64   is necessary to easily (singular blocks / fullerene) parallelize fullerenes of size 20-1024 !**/
    int fullerenes_per_block;
    
    /** Needs 3 storage arrays for coordinates and 1 for reductions **/
    int sharedMemoryPerBlock = sizeof(coord3d)* 3 * (N + 1) + sizeof(real_t)*N;

    /** Calculates maximum number of resident fullerenes on a single Streaming Multiprocessor, multiply with multi processor count to get total batch size**/
    cudaOccupancyMaxActiveBlocksPerMultiprocessor(&fullerenes_per_block, conjugate_gradient, N, sharedMemoryPerBlock);

    return (size_t)(properties.multiProcessorCount*fullerenes_per_block);
}

void OptimizeBatch(real_t* h_X, node_t* h_cubic_neighbours, node_t* h_next_on_face, node_t* h_prev_on_face, uint8_t* h_face_right, const size_t N, const size_t batch_size){
    bool concurrent_kernels = false;
    bool single_block_fullerenes = true;
    dim3 dimBlock = dim3(N, 1, 1);
    dim3 dimGrid = dim3(batch_size, 1, 1);

    size_t* d_N;
    bool* d_single_block_fullerenes;


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
    real_t* d_gdata;

    cudaError_t error;
    error = cudaMalloc(&d_X, sizeof(coord3d)*N*batch_size);
    error = cudaMalloc(&d_X_temp, sizeof(coord3d)*N*batch_size);
    error = cudaMalloc(&d_X1, sizeof(coord3d)*N*batch_size);
    error = cudaMalloc(&d_X2, sizeof(coord3d)*N*batch_size);
    error = cudaMalloc(&d_delta_x0, sizeof(coord3d)*N*batch_size);
    error = cudaMalloc(&d_delta_x1, sizeof(coord3d)*N*batch_size);
    error = cudaMalloc(&d_direction, sizeof(coord3d)*N*batch_size);
    
    error = cudaMalloc(&d_neighbours, sizeof(node_t)*3*N*batch_size);
    error = cudaMalloc(&d_next_on_face, sizeof(node_t)*3*N*batch_size);
    error = cudaMalloc(&d_prev_on_face, sizeof(node_t)*3*N*batch_size);
    error = cudaMalloc(&d_face_right, sizeof(uint8_t)*3*N*batch_size);
    error = cudaMalloc(&d_gdata, sizeof(real_t)*dimGrid.x);
    error = cudaMalloc(&d_N, sizeof(size_t)); cudaMemcpy(d_N, &N, sizeof(size_t), cudaMemcpyHostToDevice);
    error = cudaMalloc(&d_single_block_fullerenes, sizeof(bool)); cudaMemcpy(d_single_block_fullerenes, &single_block_fullerenes, sizeof(bool), cudaMemcpyHostToDevice);

    error = cudaMemcpy(d_X, h_X, sizeof(coord3d)*N*batch_size , cudaMemcpyHostToDevice);
    error = cudaMemcpy(d_neighbours, h_cubic_neighbours, sizeof(node_t)*3*N*batch_size, cudaMemcpyHostToDevice);
    error = cudaMemcpy(d_next_on_face, h_next_on_face, sizeof(node_t)*3*N*batch_size, cudaMemcpyHostToDevice);
    error = cudaMemcpy(d_prev_on_face, h_prev_on_face, sizeof(node_t)*3*N*batch_size, cudaMemcpyHostToDevice);
    error = cudaMemcpy(d_face_right, h_face_right, sizeof(uint8_t)*3*N*batch_size, cudaMemcpyHostToDevice);

    auto start = std::chrono::system_clock::now();

    if (!concurrent_kernels)
    {
        void* kernelArgs[] = {
        (void*)&d_X,
        (void*)&d_X_temp,
        (void*)&d_X1,
        (void*)&d_X2,
        (void*)&d_delta_x0, 
        (void*)&d_delta_x1,
        (void*)&d_direction,
        (void*)&d_neighbours,
        (void*)&d_next_on_face,
        (void*)&d_prev_on_face,
        (void*)&d_face_right,
        (void*)&N,
        };
        cudaLaunchCooperativeKernel((void*)conjugate_gradient, dimGrid, dimBlock, kernelArgs, sizeof(coord3d)*3*(N+1) + sizeof(real_t)*N, NULL);
    } 
    cudaDeviceSynchronize();
    auto end = std::chrono::system_clock::now();
    std::cout << "Elapsed time: " << (end-start)/ 1ms << "ms\n" ;
    std::cout << "Estimated Performance " << ((real_t)(411*N*batch_size*3*N*22  + 2106*N*batch_size*3*N)/(std::chrono::duration_cast<std::chrono::microseconds>(end-start)).count()) * 1.0e6 << "FLOP/s \n";


    cudaMemcpy(h_X, d_X, sizeof(coord3d)*batch_size*N, cudaMemcpyDeviceToHost);
    cudaDeviceSynchronize();

    cudaFree(d_X); cudaFree(d_X2); cudaFree(d_neighbours); cudaFree(d_next_on_face); cudaFree(d_prev_on_face);
    cudaFree(d_X_temp); cudaFree(d_face_right); cudaFree(d_gdata); cudaFree(d_delta_x0); cudaFree(d_delta_x1); cudaFree(d_direction);
}

};

