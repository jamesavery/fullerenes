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



struct ArcData{
    //124 FLOPs;
    __device__ ArcData(const node_t a, const uint8_t j, const coord3d* __restrict__ X, const BookkeepingData& bdat){   
        this->j = j;   

        real_t r_rmp;
        coord3d ap, am, ab, ac, ad, mp;
        coord3d X_a = X[a]; coord3d X_b = X[bdat.neighbours[j]];
        //printf("Index: %d \n", a*3 + j);

        //Compute the arcs ab, ac, ad, bp, bm, ap, am, mp, bc and cd
        ab = (X_b - X_a);  r_rab = bond_length(ab); ab_hat = r_rab * ab;
        ac = (X[bdat.neighbours[(j+1)%3]] - X_a); r_rac = bond_length(ac); ac_hat = r_rac * ac;
        ad = (X[bdat.neighbours[(j+2)%3]] - X_a); r_rad = bond_length(ad); ad_hat = r_rad * ad;
        
        coord3d bp = (X[bdat.next_on_face[j]] - X_b); bp_hat = unit_vector(bp);
        coord3d bm = (X[bdat.prev_on_face[j]] - X_b); bm_hat = unit_vector(bm);

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

    //3 FLOPs
    INLINE real_t harmonic_energy(const real_t p0, const real_t p) const{
        return (real_t)0.5*(p-p0)*(p-p0);
    }
    //4 FLOPs
    INLINE coord3d  harmonic_energy_gradient(const real_t p0, const real_t p, const coord3d gradp) const{
        return (p-p0)*gradp;     
    }

    //1 FLOP
    INLINE real_t bond() const {return (real_t)1.0/r_rab;}

    //5 FLOPs
    INLINE real_t angle() const {return dot(ab_hat,ac_hat);}

    //Returns the inner dihedral angle for the current arc. Used here only for energy calculation, 
    //otherwise embedded in dihedral computation because the planes and angles that make up the dihedral angle computation are required for derivative computation.
    //50 FLOPs
    INLINE real_t dihedral() const 
    { 
        coord3d nabc, nbcd; real_t cos_b, cos_c, r_sin_b, r_sin_c;
        cos_b = dot(ba_hat,bc_hat);  r_sin_b = rsqrt((real_t)1.0 - cos_b*cos_b); nabc = cross(ba_hat, bc_hat) * r_sin_b;
        cos_c = dot(-bc_hat,cd_hat); r_sin_c = rsqrt((real_t)1.0 - cos_c*cos_c); nbcd = cross(-bc_hat,cd_hat) * r_sin_c;
        return dot(nabc, nbcd);
    }
    
    // Chain rule terms for angle calculation
    //Computes gradient related to bending term. ~24 FLOPs
    INLINE coord3d inner_angle_gradient(const Constants& c) const
    {   
        real_t cos_angle = angle(); //Inner angle of arcs ab,ac.
        coord3d grad = cos_angle * (ab_hat * r_rab + ac_hat * r_rac) - ab_hat * r_rac - ac_hat* r_rab; //Derivative of inner angle: Eq. 21. 
        return get(c.f_inner_angle,j) * harmonic_energy_gradient(get(c.angle0,j), cos_angle, grad); //Harmonic Energy Gradient: Eq. 21. multiplied by harmonic term.
    }
    //Computes gradient related to bending of outer angles. ~20 FLOPs
    INLINE coord3d outer_angle_gradient_m(const Constants& c) const
    {
        real_t cos_angle = -dot(ab_hat, bm_hat); //Compute outer angle. ab,bm
        coord3d grad = (bm_hat + ab_hat * cos_angle) * r_rab; //Derivative of outer angles Eq. 30. Buster Thesis
        return get(c.f_outer_angle_m,j) * harmonic_energy_gradient(get(c.outer_angle_m0,j),cos_angle,grad); //Harmonic Energy Gradient: Eq. 30 multiplied by harmonic term.
    }
    INLINE coord3d outer_angle_gradient_p(const Constants& c) const
    {
        real_t cos_angle = -dot(ab_hat, bp_hat); //Compute outer angle. ab,bp
        coord3d grad = (bp_hat + ab_hat * cos_angle) * r_rab; //Derivative of outer angles Eq. 28. Buster Thesis
        return get(c.f_outer_angle_p,j) * harmonic_energy_gradient(get(c.outer_angle_p0,j),cos_angle,grad); //Harmonic Energy Gradient: Eq. 28 multiplied by harmonic term.
    }
    // Chain rule terms for dihedral calculation
    //Computes gradient related to dihedral/out-of-plane term. ~75 FLOPs
    INLINE coord3d inner_dihedral_gradient(const Constants& c) const
    {
        coord3d nabc, nbcd; real_t cos_b, cos_c, r_sin_b, r_sin_c;
        cos_b = dot(ba_hat,bc_hat); r_sin_b = rsqrtf((real_t)1.0 - cos_b*cos_b); nabc = cross(ba_hat, bc_hat) * r_sin_b;
        cos_c = dot(-bc_hat,cd_hat); r_sin_c = rsqrtf((real_t)1.0 - cos_c*cos_c); nbcd = cross(-bc_hat,cd_hat) * r_sin_c;

        real_t cos_beta = dot(nabc, nbcd); //Inner dihedral angle from planes abc,bcd.
        real_t cot_b = cos_b * r_sin_b * r_sin_b; //cos(b)/sin(b)^2

        //Derivative w.r.t. inner dihedral angle F and G in Eq. 26
        coord3d grad = cross(bc_hat, nbcd) * r_sin_b * r_rab - ba_hat * cos_beta * r_rab + (cot_b * cos_beta * r_rab) * (bc_hat - ba_hat * cos_b);

        return get(c.f_inner_dihedral,j) * harmonic_energy_gradient(get(c.inner_dih0,j), cos_beta, grad); //Eq. 26.
    }

    //Computes gradient from dihedral angles constituted by the planes bam, amp ~162 FLOPs
    INLINE coord3d outer_a_dihedral_gradient(const Constants& c) const
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
        return get(c.f_outer_dihedral,j) * harmonic_energy_gradient(get(c.outer_dih0_a,j), cos_beta, grad);
    }

    //Computes gradient from dihedral angles constituted by the planes nbmp, nmpa ~92 FLOPs
    INLINE coord3d outer_m_dihedral_gradient(const Constants& c) const
    {
        coord3d nbmp_hat, nmpa_hat; real_t cos_m, cos_p, r_sin_m, r_sin_p;
        cos_m = dot(mb_hat,mp_hat);  r_sin_m = rsqrtf((real_t)1.0 - cos_m*cos_m); nbmp_hat = cross(mb_hat,mp_hat) * r_sin_m;
        cos_p = dot(-mp_hat,pa_hat); r_sin_p = rsqrtf((real_t)1.0 - cos_p*cos_p); nmpa_hat = cross(-mp_hat,pa_hat) * r_sin_p;
        
        //Cosine to the outer dihedral angle constituted by the planes bmp and mpa
        real_t cos_beta = dot(nbmp_hat, nmpa_hat); //Outer dihedral angle bmp,mpa.
        real_t cot_p = cos_p * r_sin_p * r_sin_p;
        
        //Derivative w.r.t. outer dihedral angle, factorized version of Eq. 32.
        coord3d grad = r_rap * (cot_p*cos_beta * (-mp_hat - pa_hat*cos_p) - cross(nbmp_hat, mp_hat)*r_sin_p - pa_hat*cos_beta );

        //Eq. 32 multiplied by harmonic term.
        return get(c.f_outer_dihedral,j) * harmonic_energy_gradient(get(c.outer_dih0_m,j), cos_beta, grad);
    }

    //Computes gradient from dihedral angles constituted by the planes bpa, pam ~162 FLOPs
    INLINE coord3d outer_p_dihedral_gradient(const Constants& c) const
    {
        coord3d nbpa_hat, npam_hat; real_t cos_p, cos_a, r_sin_p, r_sin_a;
        cos_a = dot(ap_hat,am_hat);  r_sin_a = rsqrtf((real_t)1.0 - cos_a*cos_a); npam_hat = cross(ap_hat,am_hat)  * r_sin_a;
        cos_p = dot(pb_hat,-ap_hat); r_sin_p = rsqrtf((real_t)1.0 - cos_p*cos_p); nbpa_hat = cross(pb_hat,-ap_hat) * r_sin_p;

        real_t cos_beta = dot(nbpa_hat, npam_hat); //Outer dihedral angle bpa, pam.
        real_t cot_p = cos_p * r_sin_p * r_sin_p;
        real_t cot_a = cos_a * r_sin_a * r_sin_a;

        //Derivative w.r.t. outer dihedral angle, factorized version of Eq. 33.
        coord3d grad = cross(npam_hat,pb_hat)*r_rap*r_sin_p - (cross(am_hat,nbpa_hat)*r_rap + cross(nbpa_hat,ap_hat)*r_ram)*r_sin_a +
                        cos_beta*(am_hat*r_ram + r_rap * ((real_t)2.0*ap_hat + cot_p*(pb_hat+cos_p*ap_hat)) - cot_a*(r_rap*(am_hat - ap_hat*cos_a) + r_ram*(ap_hat-am_hat*cos_a)));
        
        //Eq. 33 multiplied by harmonic term.
        return get(c.f_outer_dihedral,j) * harmonic_energy_gradient(get(c.outer_dih0_p,j), cos_beta, grad);
    }
    // Internal coordinate gradients
    INLINE coord3d bond_length_gradient(const Constants& c) const { return - get(c.f_bond,j) * harmonic_energy_gradient(get(c.r0,j),bond(),ab_hat);}
    //Sum of angular gradient components.
    INLINE coord3d angle_gradient(const Constants& c) const { return inner_angle_gradient(c) + outer_angle_gradient_p(c) + outer_angle_gradient_m(c);}
    //Sum of inner and outer dihedral gradient components.
    INLINE coord3d dihedral_gradient(const Constants& c) const { return inner_dihedral_gradient(c) + outer_a_dihedral_gradient(c) + outer_m_dihedral_gradient(c) + outer_p_dihedral_gradient(c);}
    //coord3d flatness()             const { return ;  }   
    
    //Harmonic energy contribution from bond stretching, angular bending and dihedral angle bending.
    //71 FLOPs
    INLINE real_t energy(const Constants& c) const {return (real_t)0.5 *get(c.f_bond,j) *harmonic_energy(bond(),get(c.r0,j))+ get(c.f_inner_angle,j)* harmonic_energy(angle(),get(c.angle0,j)) + get(c.f_inner_dihedral,j)* harmonic_energy(dihedral(),get(c.inner_dih0,j));}
    //Sum of bond, angular and dihedral gradient components.
    INLINE coord3d gradient(const Constants& c) const{return bond_length_gradient(c) + angle_gradient(c) + dihedral_gradient(c);}
    
    uint8_t j;

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

__device__ coord3d gradient(const coord3d* __restrict__ X, const node_t node_id, const BookkeepingData &dat, const Constants &constants) {
    coord3d grad = {0.0, 0.0, 0.0};

    for (uint8_t j = 0; j < 3; j++ ){
        ArcData arc = ArcData(node_id, j, X, dat);
        grad += arc.gradient(constants);
    }
    return grad;
}

__device__ real_t energy(const coord3d* __restrict__ X, const node_t node_id, const BookkeepingData &dat, const Constants &constants, real_t* __restrict__ reduction_array, real_t* __restrict__ gdata, const node_t N, bool single_block_fullerenes) {
    real_t arc_energy = (real_t)0.0;

    //(71 + 124) * 3 * N  = 585*N FLOPs
    for (uint8_t j = 0; j < 3; j++ ){
        ArcData arc = ArcData(node_id, j, X, dat);
        arc_energy += arc.energy(constants);
    }
    cg::sync(cg::this_thread_block());
    reduction_array[threadIdx.x] = arc_energy;
    // (/N // 32) * log2(32) = N//32  * 5 FLOPs 
    reduction(reduction_array,gdata, N, single_block_fullerenes); 
     return reduction_array[0];
}

__device__ void golden_section_search(coord3d* __restrict__ X, coord3d& direction, coord3d& new_direction,coord3d* __restrict__ X1, coord3d* __restrict__ X2, real_t* __restrict__ reduction_array, real_t* __restrict__ gdata, const node_t node_id, const node_t N, const BookkeepingData& dat, const Constants& constants, cg::thread_group sync_group, bool single_block_fullerenes){
    real_t tau = (sqrtf(5) - 1) / 2;
    //Actual coordinates resulting from each traversal 
    //Line search x - values;
    real_t a = 0.0; real_t b = 1.0;
    real_t x1,  x2;/* , dfc; */	/* TODO: dfc is not used: should it be? */


    x1 = (a + (1 - tau) * (b - a));
    x2 = (a + tau * (b - a));

    X1[node_id] = X[node_id] + x1 * direction;
    X2[node_id] = X[node_id] + x2 * direction;


    cg::sync(sync_group);

    real_t f1 = energy(X1, node_id, dat, constants, reduction_array, gdata, N, single_block_fullerenes);
    real_t f2 = energy(X2, node_id, dat, constants, reduction_array, gdata, N, single_block_fullerenes);

    for (uint8_t i = 0; i < 20; i++){
        if (f1 > f2){
            a = x1;
            x1 = x2;
            f1 = f2;
            x2 = a + tau * (b - a);
            cg::sync(sync_group);
            X2[node_id] = X[node_id] + x2 * direction;
            cg::sync(sync_group);
            f2 = energy(X2, node_id, dat, constants, reduction_array, gdata, N, single_block_fullerenes);
        }else
        {
            b = x2;
            x2 = x1;
            f2 = f1;
            x1 = a + (1 - tau) * (b - a);
            cg::sync(sync_group);
            X1[node_id] = X[node_id] + x1 * direction;
            cg::sync(sync_group);
            f1 = energy(X1, node_id, dat, constants, reduction_array, gdata, N, single_block_fullerenes);
        }
    }
    //Line search coefficient
    real_t alfa = (a+b)/2;
    cg::sync(sync_group);
    X1[node_id] = X[node_id] + alfa*direction;
    cg::sync(sync_group);
    new_direction = -gradient(X1,node_id,dat, constants);
    
}

__global__ void conjugate_gradient(coord3d* d_X, coord3d* d_X_temp, coord3d* d_X2, const node_t* d_neighbours, const node_t* d_next_on_face, const node_t* d_prev_on_face, const uint8_t* d_face_right, real_t* gdata, const size_t N, const bool single_block_fullerenes, const size_t MaxIter){
    extern __shared__ real_t smem[];
    

    cg::grid_group grid = cg::this_grid();
    cg::thread_block block = cg::this_thread_block();

    coord3d delta_x0, delta_x1, direction;
    node_t node_id;
    size_t offset;
    size_t gradient_evals = 0;
    size_t energy_evals = 0;
    
    real_t beta, dnorm, r0_norm, direction_norm;
    beta = dnorm = r0_norm = direction_norm = 0.0;

    //If fullerenes are localized to individual blocks then use block threadIdx else use grid ID and use 1 grid per fullerene with concurrent launches.
    //If fullerenes are localized to individual blocks then use block size to determine pointer in array, else assume concurrent kernel launches with pointers to individual fullerenes.
    if (single_block_fullerenes){
        node_id = threadIdx.x;
        offset = blockIdx.x * blockDim.x;
    } else
    {
        node_id = blockDim.x * blockIdx.x + threadIdx.x;
        offset = 0;
    }
    
    coord3d* X = &d_X[offset];
    coord3d* X_temp = &d_X_temp[offset];
    coord3d* X2 = &d_X2[offset];
    
    if (single_block_fullerenes)
    {
        coord3d* sX =&reinterpret_cast<coord3d*>(smem)[(int)ceil(N/3) + 1];
        coord3d* sX_temp =&reinterpret_cast<coord3d*>(smem)[(int)ceil(N/3) + 2 + N];
        coord3d* sX2 =&reinterpret_cast<coord3d*>(smem)[(int)ceil(N/3)+ 3 +2*N];  
        sX[node_id] = X[node_id];
        sX_temp[node_id] = sX[node_id];

        X = &sX[0]; X_temp = &sX_temp[0]; X2 = &sX2[0];
    } else {
        X_temp[node_id] = X[node_id];
    }

    
    //Pre-compute force constants and store in registers.
    BookkeepingData bookit = BookkeepingData(&d_neighbours[3*offset],&d_face_right[3*offset],&d_next_on_face[3*offset],&d_prev_on_face[3*offset]);
    Constants constants = compute_constants(bookit, node_id);

    //Load constant bookkeeping data into registers.
    const node_t neighbours[3] = {d_neighbours[3*(offset+node_id)],d_neighbours[3*(offset+node_id) + 1],d_neighbours[3*(offset+node_id) + 2]};
    const uint8_t face_right[3] = {d_face_right[3*(offset+node_id)],d_face_right[3*(offset+node_id) + 1],d_face_right[3*(offset+node_id) + 2]};;
    const node_t next_on_face[3] = {d_next_on_face[3*(offset+node_id)],d_next_on_face[3*(offset+node_id) + 1],d_next_on_face[3*(offset+node_id) + 2]};
    const node_t prev_on_face[3] = {d_prev_on_face[3*(offset+node_id)],d_prev_on_face[3*(offset+node_id) + 1],d_prev_on_face[3*(offset+node_id) + 2]};
    BookkeepingData bookkeeping = BookkeepingData(&neighbours[0],&face_right[0],&next_on_face[0],&prev_on_face[0]);   

    cg::sync(grid);
    direction = gradient(X, node_id ,bookkeeping, constants);
    gradient_evals ++;
    
    smem[threadIdx.x] = dot(direction,direction);

    reduction(smem,gdata,N,single_block_fullerenes);
    dnorm = sqrtf(smem[0]);
    direction = -direction/dnorm;
    

    delta_x0 = direction;
    
    for (size_t i = 0; i < MaxIter; i++)
    {   
        beta = 0.0; direction_norm = 0.0; dnorm=0.0; r0_norm = 0.0;
        cg::sync(block);
        if (single_block_fullerenes){golden_section_search(X, direction, delta_x1, X_temp, X2, smem, gdata, node_id, N, bookkeeping, constants, block, single_block_fullerenes);} 
        else { golden_section_search(X, direction, delta_x1, X_temp, X2, smem, gdata, node_id, N, bookkeeping, constants, grid, single_block_fullerenes);}
        
        
        cg::sync(block);

        gradient_evals++;
        energy_evals += 22;
        //Polak Ribiere method
        
        smem[threadIdx.x] = dot(delta_x0, delta_x0); reduction(smem,gdata,N,single_block_fullerenes); r0_norm = smem[0];
        cg::sync(block);
        smem[threadIdx.x] = dot(delta_x1, (delta_x1 - delta_x0)); reduction(smem,gdata,N,single_block_fullerenes); beta = smem[0] / r0_norm;
        cg::sync(block);
        real_t E1 = energy(X_temp, node_id, bookkeeping, constants, smem, gdata, N, single_block_fullerenes);
        cg::sync(block);
        real_t E2 = energy(X, node_id, bookkeeping, constants, smem, gdata, N, single_block_fullerenes);
        cg::sync(block);
        if (E1> E2)
        {   
            X_temp[node_id] =  X[node_id];
            delta_x1 =  delta_x0;
            beta = 0.0;
        }
        else
        {   
            X[node_id] = X_temp[node_id];
            delta_x0 = delta_x1;
        }
        direction = delta_x1 + beta*direction;
        //Calculate gradient and residual gradient norms..
        cg::sync(block);
        smem[threadIdx.x] = dot(direction,direction); reduction(smem,gdata,N,single_block_fullerenes); direction_norm = sqrtf(smem[0]);
        cg::sync(block);
        smem[threadIdx.x] = dot(delta_x1,delta_x1); reduction(smem,gdata,N,single_block_fullerenes);dnorm = sqrtf(smem[0]);
        cg::sync(block);
        //Normalize gradient.
        direction /= direction_norm;
    }
    d_X[offset + threadIdx.x] = X[threadIdx.x];

    //Check convergence using bond lengths:
    real_t AvgBondLength = AverageBondLength(smem,X,neighbours,node_id,N);
    if (threadIdx.x == 0)
    {
        //printf("Avg Bond Length: %e \n", AvgBondLength);
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

void AllocateDevicePointers(DevicePointers& p, size_t N, size_t batch_size){
    cudaMalloc(&p.d_X, sizeof(coord3d)*N*batch_size);
    cudaMalloc(&p.d_X_temp, sizeof(coord3d)*N*batch_size);
    cudaMalloc(&p.d_X2, sizeof(coord3d)*N*batch_size);
    cudaMalloc(&p.d_neighbours, sizeof(node_t)*3*N*batch_size);
    cudaMalloc(&p.d_next_on_face, sizeof(node_t)*3*N*batch_size);
    cudaMalloc(&p.d_prev_on_face, sizeof(node_t)*3*N*batch_size);
    cudaMalloc(&p.d_face_right, sizeof(uint8_t)*3*N*batch_size);
    cudaMalloc(&p.d_gdata, sizeof(real_t)*batch_size);
}

void FreePointers(DevicePointers& p){
    cudaFree(p.d_X);
    cudaFree(p.d_X_temp);
    cudaFree(p.d_X2);
    cudaFree(p.d_neighbours);
    cudaFree(p.d_next_on_face);
    cudaFree(p.d_prev_on_face);
    cudaFree(p.d_face_right);
    cudaFree(p.d_gdata);
    cudaDeviceReset();
}

void OptimizeBatch(DevicePointers& p, real_t* h_X, node_t* h_cubic_neighbours, node_t* h_next_on_face, node_t* h_prev_on_face, uint8_t* h_face_right, const size_t N, const size_t batch_size, const size_t MaxIter){
    bool single_block_fullerenes = true;

    getLastCudaError("One or more Mallocs Failed! \n");
    cudaMemcpy(p.d_X, h_X, sizeof(coord3d)*N*batch_size , cudaMemcpyHostToDevice);
    cudaMemcpy(p.d_neighbours, h_cubic_neighbours, sizeof(node_t)*3*N*batch_size, cudaMemcpyHostToDevice);
    cudaMemcpy(p.d_next_on_face, h_next_on_face, sizeof(node_t)*3*N*batch_size, cudaMemcpyHostToDevice);
    cudaMemcpy(p.d_prev_on_face, h_prev_on_face, sizeof(node_t)*3*N*batch_size, cudaMemcpyHostToDevice);
    cudaMemcpy(p.d_face_right, h_face_right, sizeof(uint8_t)*3*N*batch_size, cudaMemcpyHostToDevice);
    getLastCudaError("Memcpy Failed! \n");
    auto start = std::chrono::system_clock::now();
    void* kernelArgs[] = {
    (void*)&(p.d_X),
    (void*)&(p.d_X_temp),
    (void*)&(p.d_X2),
    (void*)&(p.d_neighbours),
    (void*)&(p.d_next_on_face),
    (void*)&(p.d_prev_on_face),
    (void*)&(p.d_face_right),
    (void*)&(p.d_gdata),
    (void*)&N,
    (void*)&single_block_fullerenes,
    (void*)&MaxIter
    };

    cudaLaunchCooperativeKernel((void*)conjugate_gradient, dim3(batch_size, 1, 1), dim3(N, 1, 1), kernelArgs, sizeof(coord3d)*3*(N+1) + sizeof(real_t)*N, NULL);
    cudaDeviceSynchronize();
    auto end = std::chrono::system_clock::now();
    getLastCudaError("Failed to launch kernel: ");
    
    cudaMemcpy(h_X, p.d_X, sizeof(coord3d)*N*batch_size , cudaMemcpyDeviceToHost);
    getLastCudaError("Failed to copy back: ");
    
    std::cout << "Elapsed time: " << (end-start)/ 1ms << "ms\n" ;
    std::cout << "Estimated Performance " << ((real_t)(411*N*batch_size*MaxIter*22  + 2106*N*batch_size*MaxIter)/(std::chrono::duration_cast<std::chrono::microseconds>(end-start)).count()) * 1.0e6 << "FLOP/s \n";
}

};
