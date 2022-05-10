#include <cuda.h>
#include "fullerenes/gpu/kernels.hh"
#include "fullerenes/gpu/cuda_execution.hh"
namespace gpu_kernels{

namespace isomerspace_X0{
#include "cuda_runtime.h"
#include <cuda_runtime_api.h>
#include "cooperative_groups.h"
#include "fullerenes/gpu/cuda_definitions.h"
#include "misc_cuda.cu"
#include "reductions.cu"
#include "coord3d.cuh"
#include "coord2d.cuh"
#include "forcefield_structs.cu"
#include "print_functions.cu"


__device__
device_node_t multiple_source_shortest_paths(const IsomerBatch& B, device_node_t* distances, const size_t isomer_idx){
    DEVICE_TYPEDEFS
    
    DeviceFullereneGraph FG = DeviceFullereneGraph(&B.neighbours[isomer_idx*blockDim.x*3]);
    node_t outer_face[6];
    uint8_t Nface = FG.get_face_oriented(0, FG.neighbours[0],outer_face);
    distances[threadIdx.x] = node_t(NODE_MAX);    
    BLOCK_SYNC
    if (threadIdx.x < Nface)  distances[outer_face[threadIdx.x]] = 0;
    BLOCK_SYNC
    if (threadIdx.x == 0){
        CuDeque<node_t> queue = CuDeque<node_t>(distances + blockDim.x, blockDim.x);
        for (size_t i = 0; i < Nface; i++) queue.push_back(outer_face[i]);
        while (!queue.empty())
        {   
            node_t v = queue.pop_front();
            for (size_t i = 0; i < 3; i++)
            {   
                node_t w = FG.neighbours[v*3 + i];
                if(distances[w] == NODE_MAX) {
                distances[w] = distances[v]+1;
                queue.push_back(w);
                }
            }
        }
    }
    BLOCK_SYNC
    device_node_t distance = distances[threadIdx.x];
    BLOCK_SYNC
    return distance;
}


__device__
device_coord2d spherical_projection(const IsomerBatch& B, device_node_t* sdata, const size_t isomer_idx){
    DEVICE_TYPEDEFS

    node_t distance =  multiple_source_shortest_paths(B,reinterpret_cast<node_t*>(sdata), isomer_idx);
    BLOCK_SYNC
    clear_cache(reinterpret_cast<real_t*>(sdata), Block_Size_Pow_2); 
    node_t d_max = reduction_max(sdata, distance);

    clear_cache(reinterpret_cast<real_t*>(sdata), Block_Size_Pow_2); 
    ordered_atomic_add(&reinterpret_cast<real_t*>(sdata)[distance],real_t(1.0)); 
    BLOCK_SYNC
    node_t num_of_same_dist = node_t(reinterpret_cast<real_t*>(sdata)[distance]); 
    BLOCK_SYNC
    clear_cache(reinterpret_cast<real_t*>(sdata), Block_Size_Pow_2);
    BLOCK_SYNC
    coord2d xys = reinterpret_cast<coord2d*>(B.xys)[isomer_idx*blockDim.x + threadIdx.x]; BLOCK_SYNC
    ordered_atomic_add(&reinterpret_cast<real_t*>(sdata)[distance*2], xys.x); 
    ordered_atomic_add(&reinterpret_cast<real_t*>(sdata)[distance*2+1], xys.y); BLOCK_SYNC
    coord2d centroid = reinterpret_cast<coord2d*>(sdata)[distance] / num_of_same_dist; BLOCK_SYNC    
    coord2d xy = xys - centroid;
    real_t dtheta = real_t(M_PI)/real_t(d_max+1); 
    real_t phi = dtheta*(distance + real_t(0.5)); 
    real_t theta = atan2(xy.x, xy.y); 
    coord2d spherical_layout = {theta, phi};
    

    return spherical_layout;
}

__global__
void zero_order_geometry_(IsomerBatch B, device_real_t scalerad, int offset){
    DEVICE_TYPEDEFS
    
    extern __shared__  device_real_t sdata[];
    clear_cache(sdata, Block_Size_Pow_2);
    size_t isomer_idx = blockIdx.x + offset;
    if (isomer_idx < B.isomer_capacity && B.statuses[isomer_idx] != EMPTY){
    NodeGraph node_graph = NodeGraph(B, isomer_idx); 
    coord2d angles = spherical_projection(B,reinterpret_cast<device_node_t*>(sdata), isomer_idx);
    real_t theta = angles.x; real_t phi = angles.y;
    real_t x = cos(theta)*sin(phi), y = sin(theta)*sin(phi), z = cos(phi);
    coord3d coordinate = {x, y ,z};

    clear_cache(sdata, Block_Size_Pow_2);
    x = reduction(sdata, coordinate.x); y = reduction(sdata, coordinate.y); z = reduction(sdata,coordinate.z);
    coord3d cm = {x, y, z};
    cm /= blockDim.x;
    coordinate -= cm;
    real_t Ravg = real_t(0.0);
    clear_cache(sdata, Block_Size_Pow_2);
    real_t* base_pointer = sdata + Block_Size_Pow_2; 
    coord3d* X = reinterpret_cast<coord3d*>(base_pointer);
    X[threadIdx.x] = coordinate;
    BLOCK_SYNC
    real_t local_Ravg = real_t(0.0);
    for (uint8_t i = 0; i < 3; i++) {local_Ravg += norm(X[threadIdx.x] - X[d_get(node_graph.neighbours,i)]);}
    Ravg = reduction(sdata, local_Ravg);
    Ravg /= real_t(3*blockDim.x);
    coordinate *= scalerad*1.5/Ravg;
    reinterpret_cast<coord3d*>(B.X)[blockDim.x*isomer_idx + threadIdx.x] = coordinate;
    B.statuses[isomer_idx] = CONVERGED;
    }
}

float kernel_time = 0.0;
float time_spent(){
    return kernel_time;
}

cudaError_t zero_order_geometry(IsomerBatch& B, const device_real_t scalerad, const LaunchCtx& ctx, const LaunchPolicy policy){
    //Need a way of telling whether the kernel has been called previously.
    static bool first_call = true;
    static cudaEvent_t start, stop;
    float single_kernel_time = 0.0;
    //Construct events only once
    if(first_call) {cudaEventCreate(&start); cudaEventCreate(&stop);}

    //Records time from previous kernel call
    cudaEventElapsedTime(&single_kernel_time, start, stop);
    kernel_time += single_kernel_time;

    //If launch ploicy is synchronous then wait.
    if (policy == LaunchPolicy::SYNC) ctx.wait();
    cudaSetDevice(ctx.get_device_id());
    size_t smem =  sizeof(device_coord3d)*B.n_atoms + sizeof(device_real_t)*Block_Size_Pow_2;
    
    //Compute best grid dimensions once.
    static LaunchDims dims((void*)zero_order_geometry_, B.n_atoms, smem, B.isomer_capacity);
    dims.update_dims((void*)zero_order_geometry_, B.n_atoms, smem, B.isomer_capacity);
    cudaError_t error;

    //Note: some memory bug exists when using grid-stride for loops inside the kernel launches
    cudaEventRecord(start, ctx.stream);
    for (int i = 0; i < B.isomer_capacity + (dims.get_grid().x - B.isomer_capacity % dims.get_grid().x ); i += dims.get_grid().x)
    {
        void* kargs[]{(void*)&B, (void*)&scalerad, (void*)&i};
        error = safeCudaKernelCall((void*)zero_order_geometry_, dims.get_grid(), dims.get_block(), kargs, smem, ctx.stream);
    }
    cudaEventRecord(stop, ctx.stream);
    
    if (policy == LaunchPolicy::SYNC) ctx.wait();
    printLastCudaError("Zero order geometry:");
    first_call = false;
    return error;
}

}}