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

__global__
void cubic_layout_(IsomerBatch B){
    DEVICE_TYPEDEFS
    extern __shared__  real_t sharedmem[];
    __shared__ device_node_t* triangle_numbers;

    if(threadIdx.x == 0) triangle_numbers = (device_node_t*)malloc(B.n_faces*6*sizeof(device_node_t));

    BLOCK_SYNC
    for (int isomer_idx = blockIdx.x; isomer_idx < B.isomer_capacity; isomer_idx += gridDim.x ){
    if (B.statuses[isomer_idx] != EMPTY ){
    device_node_t* smem = reinterpret_cast<device_node_t*>(sharedmem);
    DeviceFullereneDual FD(B.dual_neighbours + isomer_idx * B.n_faces * 6, B.face_degrees + isomer_idx * B.n_faces);
    auto thid = threadIdx.x;
    device_node_t cannon_arcs[6]; memset(cannon_arcs,0,sizeof(device_node_t)*6);
    bool is_representative[6]; 
    int represent_count = 0;
    if(thid < B.n_faces){
        for (int i = 0; i < FD.face_degrees[thid]; ++i){
            device_node2 cannon_arc = FD.get_cannonical_triangle_arc(thid,FD.dual_neighbours[thid*6 + i]);
            if (cannon_arc.x == thid) {
                cannon_arcs[i] = cannon_arc.y;
                is_representative[i] = true;
                represent_count++;
            } else {
                is_representative[i] = false;
            }
        }
    }
    exclusive_scan(smem, represent_count);

    if (thid < B.n_faces){
        int arc_count = 0;
        for (size_t i = 0; i < FD.face_degrees[thid]; ++i){
            if (is_representative[i]) {
                triangle_numbers[thid*6 + i] = smem[thid] + arc_count;
                ++arc_count;
            }
        }
    }
    BLOCK_SYNC
    device_node2* representative_arc_list = reinterpret_cast<device_node2*>(sharedmem);
    device_node_t node_id;
    if (thid < B.n_faces){
        int arc_count = 0;
        for (size_t i = 0; i < FD.face_degrees[thid]; i++){
            if(is_representative[i]){
                auto idx = triangle_numbers[thid*6 + i] + arc_count;
                representative_arc_list[idx] = {device_node_t(thid), cannon_arcs[i]}; 
            }
        }
        
    }
    BLOCK_SYNC
    /*
    if(thid < B.n_faces){
        for (int i = 0; i < FD.face_degrees[thid]; ++i){
        device_node2 cannon_arc = FD.get_cannonical_triangle_arc(thid,FD.dual_neighbours[thid*6 + i]);
        if (cannon_arc.x == thid) {
            device_node_t v(FD.dual_neighbours[thid*6 + i]);
            device_node_t w = FD.next(thid,v);
            device_node2 edge_b = FD.get_cannonical_triangle_arc(v, thid); B.cubic_neighbours[isomer_idx*B.n_atoms*3 + triangle_numbers[thid*6 + i]*3 + 0] = triangle_numbers[edge_b.x * 6 + FD.dedge_ix(edge_b.x, edge_b.y)];
            device_node2 edge_c = FD.get_cannonical_triangle_arc(w, v); B.cubic_neighbours[isomer_idx*B.n_atoms*3 + triangle_numbers[thid*6 + i]*3 + 1] = triangle_numbers[edge_c.x * 6 + FD.dedge_ix(edge_c.x, edge_c.y)];
            device_node2 edge_d = FD.get_cannonical_triangle_arc(thid, w); B.cubic_neighbours[isomer_idx*B.n_atoms*3 + triangle_numbers[thid*6 + i]*3 + 2] = triangle_numbers[edge_d.x * 6 + FD.dedge_ix(edge_d.x, edge_d.y)];
        };
        }
    }*/
    BLOCK_SYNC
    auto [u, v] = representative_arc_list[thid];
    device_node_t w = FD.next(u,v);

    device_node2 edge_b = FD.get_cannonical_triangle_arc(v, u); B.cubic_neighbours[isomer_idx*B.n_atoms*3 + thid*3 + 0] = triangle_numbers[edge_b.x * 6 + FD.dedge_ix(edge_b.x, edge_b.y)];
    device_node2 edge_c = FD.get_cannonical_triangle_arc(w, v); B.cubic_neighbours[isomer_idx*B.n_atoms*3 + thid*3 + 1] = triangle_numbers[edge_c.x * 6 + FD.dedge_ix(edge_c.x, edge_c.y)];
    device_node2 edge_d = FD.get_cannonical_triangle_arc(u, w); B.cubic_neighbours[isomer_idx*B.n_atoms*3 + thid*3 + 2] = triangle_numbers[edge_d.x * 6 + FD.dedge_ix(edge_d.x, edge_d.y)];

}}}

float kernel_time = 0.0;
float time_spent(){
    return kernel_time;
}

cudaError_t cubic_layout(IsomerBatch& B, const LaunchCtx& ctx, const LaunchPolicy policy){
    static bool first_call = true;
    static cudaEvent_t start, stop;
    float single_kernel_time = 0.0;
    if(first_call) {cudaEventCreate(&start); cudaEventCreate(&stop);}
    cudaEventElapsedTime(&single_kernel_time, start, stop);
    kernel_time += single_kernel_time;

    if(policy == LaunchPolicy::SYNC) ctx.wait();
    size_t smem = sizeof(device_coord2d)*B.n_atoms*2 + sizeof(device_real_t)*Block_Size_Pow_2;
    static LaunchDims dims((void*)cubic_layout_, B.n_atoms, smem, B.isomer_capacity);
    dims.update_dims((void*)cubic_layout_, B.n_atoms, smem, B.isomer_capacity);
    cudaError_t error;
    cudaEventRecord(start, ctx.stream);

    void* kargs[]{(void*)&B};
    error = safeCudaKernelCall((void*)cubic_layout_, dims.get_grid(), dims.get_block(), kargs, smem, ctx.stream);  
    cudaDeviceSynchronize();
    
    cudaEventRecord(stop, ctx.stream);
    
    if(policy == LaunchPolicy::SYNC) ctx.wait();
    printLastCudaError("Cubic Layout: ");
    first_call = false;
    return error;
}

}}