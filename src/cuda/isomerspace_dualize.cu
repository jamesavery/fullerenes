#include <cuda.h>
#include "fullerenes/gpu/kernels.hh"
namespace gpu_kernels{
namespace isomerspace_dual{
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
#include "chrono"

__global__
void cubic_layout_(IsomerBatch B){
    unsigned int power2 = B.n_faces; // compute the next highest power of 2 of 32-bit power2
    power2--;
    power2 |= power2 >> 1;
    power2 |= power2 >> 2;
    power2 |= power2 >> 4;
    power2 |= power2 >> 8;
    power2 |= power2 >> 16;
    power2++;

    DEVICE_TYPEDEFS
    extern __shared__  real_t sharedmem[];
    device_node_t* triangle_numbers = reinterpret_cast<device_node_t*>(sharedmem);
    device_node_t* cached_neighbours = triangle_numbers + B.n_faces*6;
    uint8_t* cached_degrees = reinterpret_cast<uint8_t*>(cached_neighbours+ B.n_faces*6);
    device_node_t* smem = reinterpret_cast<device_node_t*>(cached_neighbours) + B.n_faces*8;
    clear_cache(reinterpret_cast<real_t*>(smem),power2);
    for (int isomer_idx = blockIdx.x; isomer_idx < B.isomer_capacity; isomer_idx += gridDim.x ){
    if (B.statuses[isomer_idx] != IsomerStatus::EMPTY ){

    BLOCK_SYNC
    auto thid = threadIdx.x;
    if (thid < B.n_faces){
        reinterpret_cast<node6*>(cached_neighbours)[thid] = reinterpret_cast<node6*>(B.dual_neighbours)[thid + isomer_idx*B.n_faces];
        cached_degrees[thid] = B.face_degrees[thid + isomer_idx*B.n_faces];
    }
    DeviceFullereneDual FD(cached_neighbours, cached_degrees);
    device_node_t cannon_arcs[6]; memset(cannon_arcs, UINT16_MAX,sizeof(device_node_t)*6);
    int represent_count = 0;
    BLOCK_SYNC
    if(thid < B.n_faces){
        for (int i = 0; i < FD.face_degrees[thid]; ++i){
            device_node2 cannon_arc = FD.get_cannonical_triangle_arc(thid,FD.dual_neighbours[thid*6 + i]);
            if (cannon_arc.x == thid) {
                cannon_arcs[i] = cannon_arc.y;
                represent_count++;
            }
        }
    }

    ex_scan<device_node_t>(smem,represent_count,B.n_faces);
    
    if (thid < B.n_faces){
        int arc_count = 0;
        for (size_t i = 0; i < FD.face_degrees[thid]; ++i){
            if (cannon_arcs[i] != UINT16_MAX) {
                triangle_numbers[thid*6 + i] = smem[thid] + arc_count;
                ++arc_count;
            }
        }
    }

    BLOCK_SYNC
    device_node2* representative_arc_list = reinterpret_cast<device_node2*>(smem);
    if (thid < B.n_faces){
        for (size_t i = 0; i < FD.face_degrees[thid]; i++){
            if(cannon_arcs[i] != UINT16_MAX){
                auto idx = triangle_numbers[thid*6 + i];
                representative_arc_list[idx] = {device_node_t(thid), cannon_arcs[i]}; 
            }
        }
        
    }
    BLOCK_SYNC
    auto [u, v] = representative_arc_list[thid];
    device_node_t w = FD.next(u,v);

    device_node2 edge_b = FD.get_cannonical_triangle_arc(v, u); B.cubic_neighbours[isomer_idx*B.n_atoms*3 + thid*3 + 0] = triangle_numbers[edge_b.x * 6 + FD.dedge_ix(edge_b.x, edge_b.y)];
    device_node2 edge_c = FD.get_cannonical_triangle_arc(w, v); B.cubic_neighbours[isomer_idx*B.n_atoms*3 + thid*3 + 1] = triangle_numbers[edge_c.x * 6 + FD.dedge_ix(edge_c.x, edge_c.y)];
    device_node2 edge_d = FD.get_cannonical_triangle_arc(u, w); B.cubic_neighbours[isomer_idx*B.n_atoms*3 + thid*3 + 2] = triangle_numbers[edge_d.x * 6 + FD.dedge_ix(edge_d.x, edge_d.y)];

}}}

float kernel_time = 0.0;
std::chrono::microseconds time_spent(){
    return std::chrono::microseconds((int) (kernel_time*1000.f));
}

cudaError_t cubic_layout(IsomerBatch& B, const LaunchCtx& ctx, const LaunchPolicy policy){
    cudaSetDevice(ctx.get_device_id());
    static std::vector<bool> first_call(16, true);
    static cudaEvent_t start[16], stop[16];
    float single_kernel_time = 0.0;
    auto dev = ctx.get_device_id();
    if(first_call[dev]) {cudaEventCreate(&start[dev]); cudaEventCreate(&stop[dev]);}

    if(policy == LaunchPolicy::SYNC) {ctx.wait();}
    else if(policy == LaunchPolicy::ASYNC && !first_call[dev]){
        cudaEventElapsedTime(&single_kernel_time, start[dev], stop[dev]);
        kernel_time += single_kernel_time;
    }

    size_t smem = sizeof(device_coord3d)*B.n_atoms*3 + sizeof(device_real_t)*Block_Size_Pow_2;
    static LaunchDims dims((void*)cubic_layout_, B.n_atoms, smem, B.isomer_capacity);
    dims.update_dims((void*)cubic_layout_, B.n_atoms, smem, B.isomer_capacity);
    cudaError_t error;
    void* kargs[]{(void*)&B};
    cudaEventRecord(start[dev], ctx.stream);
    error = safeCudaKernelCall((void*)cubic_layout_, dims.get_grid(), dims.get_block(), kargs, smem, ctx.stream);  
    cudaEventRecord(stop[dev], ctx.stream);
    
    if(policy == LaunchPolicy::SYNC) {
        ctx.wait();
        cudaEventElapsedTime(&single_kernel_time, start[dev], stop[dev]);
        kernel_time += single_kernel_time;
    }
    printLastCudaError("Cubic Layout: ");
    first_call[dev] = false;
    return error;
}

}}