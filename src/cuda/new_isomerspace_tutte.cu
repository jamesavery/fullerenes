#include <cuda.h>
#include "fullerenes/gpu/kernels.hh"
namespace gpu_kernels{
namespace isomerspace_tutte{
#include "cuda_runtime.h"
#include <cuda_runtime_api.h>
#include "cooperative_groups.h"
#include "fullerenes/gpu/cuda_definitions.h"
#include "misc_cuda.cu"
#include "reductions.cu"
#include "coord3d.cuh"
#include "coord2d.cuh"
#include "forcefield_structs.cu"

__global__
void __tutte_layout(IsomerBatch B, const size_t iterations){
    DEVICE_TYPEDEFS
    extern __shared__  real_t sharedmem[];
    clear_cache(sharedmem, Block_Size_Pow_2);

    for (size_t isomer_idx = blockIdx.x; isomer_idx < B.isomer_capacity; isomer_idx+= gridDim.x)
    {
    if (B.statuses[blockIdx.x] == NOT_CONVERGED){
    size_t offset = blockIdx.x * blockDim.x;
    DeviceFullereneGraph FG(&B.neighbours[offset*3]); 
    real_t* base_pointer        = sharedmem + Block_Size_Pow_2;
    coord2d* xys        = reinterpret_cast<coord2d*>(base_pointer);
    coord2d* newxys     = reinterpret_cast<coord2d*>(base_pointer) + blockDim.x;


    node3 ns            = (reinterpret_cast<node3*>(B.neighbours) + offset)[threadIdx.x];
    xys[threadIdx.x]    = {real_t(0.0), real_t(0.0)};
    device_node_t outer_face[6];
    device_node_t outer_face_vertex   = 0;
    uint8_t Nface = FG.get_face_oriented(0,FG.neighbours[0], outer_face);    
    reinterpret_cast<bool*>(sharedmem)[threadIdx.x] =  false; BLOCK_SYNC;
    if(threadIdx.x < Nface){
      outer_face_vertex = outer_face[threadIdx.x];
      reinterpret_cast<bool*>(sharedmem)[outer_face_vertex] =  true; 
    }
    BLOCK_SYNC;
    bool fixed = reinterpret_cast<bool*>(sharedmem)[threadIdx.x];

    if(threadIdx.x < Nface) xys[outer_face_vertex] = {sin(threadIdx.x*2*real_t(M_PI)/double(Nface)),cos(threadIdx.x*2*real_t(M_PI)/double(Nface))};
    BLOCK_SYNC
    bool converged          = false;
    real_t max_change       = real_t(0.0);
    if(fixed) newxys[threadIdx.x] = xys[threadIdx.x];

    for (size_t i = 0; i < iterations && !converged; i++)
    {   
        max_change = real_t(0.0);
        BLOCK_SYNC
        coord2d neighbour_sum   = {real_t(0.0),real_t(0.0)};    
        for (uint8_t j = 0; j < 3; j++) neighbour_sum += xys[d_get(ns,j)];

        if(!fixed) newxys[threadIdx.x] = xys[threadIdx.x]*real_t(0.15) + (neighbour_sum/3)*real_t(0.85);
        real_t neighbour_dist = 0.0f;

        for (uint8_t j = 0; j < 3; j++) neighbour_dist += norm(xys[threadIdx.x] - xys[d_get(ns,j)])/3;
        
        BLOCK_SYNC
        real_t relative_change = 0.0f;
        if (neighbour_dist > 0.0f && !fixed){ 
            relative_change = norm(xys[threadIdx.x] - newxys[threadIdx.x])/neighbour_dist;
        }

        real_t iteration_max = reduction_max(sharedmem, relative_change);
        if (iteration_max > max_change) max_change = iteration_max;
        converged = max_change <= 5e-4;

        xys[threadIdx.x] = newxys[threadIdx.x];
    }
    BLOCK_SYNC
    (reinterpret_cast<coord2d*>(B.xys) + offset )[threadIdx.x] = xys[threadIdx.x];
    B.statuses[blockIdx.x] = CONVERGED;
    }
    }
}

cudaError_t tutte_layout(IsomerBatch& B, const size_t max_iterations, const cudaStream_t stream){
    size_t smem = sizeof(device_coord2d)*B.n_atoms*2 + sizeof(device_real_t)*Block_Size_Pow_2;
    static LaunchDims dims((void*)__tutte_layout, B.n_atoms, smem);
    dims.update_dims((void*)__tutte_layout, B.n_atoms, smem);
    void* kargs[]{(void*)&B,(void*)&max_iterations};
    return safeCudaKernelCall((void*)__tutte_layout, dims.get_grid(), dims.get_block(), kargs, smem, stream);  
}

}}