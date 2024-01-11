#include <cuda.h>
#include <cooperative_groups.h>
#include <cooperative_groups/reduce.h>
#include <cooperative_groups/scan.h>
#include "fullerenes/gpu/kernels.hh"
namespace gpu_kernels{
namespace isomerspace_dual{
#include "device_includes.cu"
template cudaError_t dualize        <GPU, uint16_t>(IsomerBatch<GPU>& B, const LaunchCtx& ctx, const LaunchPolicy policy);
template cudaError_t dualize_2      <GPU, uint16_t>(IsomerBatch<GPU>& B, const LaunchCtx& ctx, const LaunchPolicy policy);
template int optimal_batch_size     <GPU, uint16_t>(const int N, const int device_id);
template int optimal_batch_size_2   <GPU, uint16_t>(const int N, const int device_id);

template <Device U, typename K> __global__
void cubic_layout_(IsomerBatch<U> B){
    unsigned int power2 = B.n_faces; // compute the next highest power of 2 of 32-bit power2
    power2--;
    power2 |= power2 >> 1;
    power2 |= power2 >> 2;
    power2 |= power2 >> 4;
    power2 |= power2 >> 8;
    power2 |= power2 >> 16;
    power2++;

    INT_TYPEDEFS(K);
    extern __shared__  float sharedmem[];
    node_t* triangle_numbers = reinterpret_cast<node_t*>(sharedmem);
    node_t* cached_neighbours = triangle_numbers + B.n_faces*6;
    uint8_t* cached_degrees = reinterpret_cast<uint8_t*>(cached_neighbours+ B.n_faces*6);
    node_t* smem = reinterpret_cast<node_t*>(cached_neighbours) + B.n_faces*8;
    //clear_cache(reinterpret_cast<real_t*>(smem),power2);
    if(threadIdx.x == 0) memset(smem, 0, power2 * sizeof(int));
    for (int isomer_idx = blockIdx.x; isomer_idx < B.isomer_capacity; isomer_idx += gridDim.x ){
    if (B.statuses[isomer_idx] != IsomerStatus::EMPTY ){

    BLOCK_SYNC
    auto thid = threadIdx.x;
    if (thid < B.n_faces){
        reinterpret_cast<node6*>(cached_neighbours)[thid] = reinterpret_cast<node6*>(B.dual_neighbours)[thid + isomer_idx*B.n_faces];
        cached_degrees[thid] = B.face_degrees[thid + isomer_idx*B.n_faces];
    }
    DeviceDualGraph FD(cached_neighbours, cached_degrees);
    node_t cannon_arcs[6]; memset(cannon_arcs, UINT16_MAX,sizeof(node_t)*6);
    int represent_count = 0;
    BLOCK_SYNC
    if(thid < B.n_faces){
        for (int i = 0; i < FD.face_degrees[thid]; ++i){
            node2 cannon_arc = FD.get_cannonical_triangle_arc(thid,FD.dual_neighbours[thid*6 + i]);
            if (cannon_arc[0] == thid) {
                cannon_arcs[i] = cannon_arc[1];
                represent_count++;
            }
        }
    }

    ex_scan<node_t>(smem,represent_count,B.n_faces);
    
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
    node2* representative_arc_list = reinterpret_cast<node2*>(smem);
    if (thid < B.n_faces){
        for (size_t i = 0; i < FD.face_degrees[thid]; i++){
            if(cannon_arcs[i] != UINT16_MAX){
                auto idx = triangle_numbers[thid*6 + i];
                representative_arc_list[idx] = {node_t(thid), cannon_arcs[i]}; 
            }
        }
        
    }
    BLOCK_SYNC
    auto [u, v] = representative_arc_list[thid];
    node_t w = FD.next(u,v);

    node2 edge_b = FD.get_cannonical_triangle_arc(v, u); B.cubic_neighbours[isomer_idx*B.n_atoms*3 + thid*3 + 0] = triangle_numbers[edge_b[0] * 6 + FD.arc_ix(edge_b[0], edge_b[1])];
    node2 edge_c = FD.get_cannonical_triangle_arc(w, v); B.cubic_neighbours[isomer_idx*B.n_atoms*3 + thid*3 + 1] = triangle_numbers[edge_c[0] * 6 + FD.arc_ix(edge_c[0], edge_c[1])];
    node2 edge_d = FD.get_cannonical_triangle_arc(u, w); B.cubic_neighbours[isomer_idx*B.n_atoms*3 + thid*3 + 2] = triangle_numbers[edge_d[0] * 6 + FD.arc_ix(edge_d[0], edge_d[1])];

}}}

template <Device U, typename K> __global__
void dualize_2_(IsomerBatch<U> B){
    /* unsigned int power2 = B.n_faces; // compute the next highest power of 2 of 32-bit power2
    power2--;
    power2 |= power2 >> 1;
    power2 |= power2 >> 2;
    power2 |= power2 >> 4;
    power2 |= power2 >> 8;
    power2 |= power2 >> 16;
    power2++;
 */
    INT_TYPEDEFS(K);
    extern __shared__  float sharedmem[];
    node_t* triangle_numbers = reinterpret_cast<node_t*>(sharedmem);
    node_t* cached_neighbours = triangle_numbers + B.n_faces*6;
    uint8_t* cached_degrees = reinterpret_cast<uint8_t*>(cached_neighbours+ B.n_faces*6);
    node_t* smem = reinterpret_cast<node_t*>(cached_neighbours) + B.n_faces*8;
    if(threadIdx.x == 0) memset(smem, 0, B.n_faces * sizeof(int));
    for (int isomer_idx = blockIdx.x; isomer_idx < B.isomer_capacity; isomer_idx += gridDim.x ){
    if (B.statuses[isomer_idx] != IsomerStatus::EMPTY ){

    BLOCK_SYNC
    auto thid = threadIdx.x;
    if (thid < B.n_faces){
        reinterpret_cast<node6*>(cached_neighbours)[thid] = reinterpret_cast<node6*>(B.dual_neighbours)[thid + isomer_idx*B.n_faces];
        cached_degrees[thid] = B.face_degrees[thid + isomer_idx*B.n_faces];
    }

    DeviceDualGraph FD(cached_neighbours, cached_degrees);
    node_t cannon_arcs[6]; memset(cannon_arcs, UINT16_MAX,sizeof(node_t)*6);
    int represent_count = 0;
    BLOCK_SYNC
    if(thid  < B.n_faces)
    for (int i = 0; i < FD.face_degrees[thid]; ++i){
        node2 cannon_arc = FD.get_cannonical_triangle_arc(thid,FD.dual_neighbours[thid*6 + i]);
        if (cannon_arc[0] == thid) {
            cannon_arcs[i] = cannon_arc[1];
            represent_count++;
        }
    }

    ex_scan<node_t>(smem,represent_count,B.n_faces);

    if(thid  < B.n_faces){
    int arc_count = 0;
    for (size_t i = 0; i < FD.face_degrees[thid]; ++i){
        if (cannon_arcs[i] != UINT16_MAX) {
            triangle_numbers[thid*6 + i] = smem[thid] + arc_count;
            ++arc_count;
        }
    }}

    BLOCK_SYNC
    node2* representative_arc_list = reinterpret_cast<node2*>(smem);
    if(thid < B.n_faces){
    for (size_t i = 0; i < FD.face_degrees[thid]; i++){
        if(cannon_arcs[i] != UINT16_MAX){
            auto idx = triangle_numbers[thid*6 + i];
            representative_arc_list[idx] = {node_t(thid), cannon_arcs[i]}; 
        }
    }}
        
    BLOCK_SYNC
    for (auto tix = thid; tix < B.n_atoms; tix += blockDim.x){
        auto [u, v] = representative_arc_list[tix];
        node_t w = FD.next(u,v);
        node2 edge_b = FD.get_cannonical_triangle_arc(v, u); B.cubic_neighbours[isomer_idx*B.n_atoms*3 + tix*3 + 0] = triangle_numbers[edge_b[0] * 6 + FD.arc_ix(edge_b[0], edge_b[1])];
        node2 edge_c = FD.get_cannonical_triangle_arc(w, v); B.cubic_neighbours[isomer_idx*B.n_atoms*3 + tix*3 + 1] = triangle_numbers[edge_c[0] * 6 + FD.arc_ix(edge_c[0], edge_c[1])];
        node2 edge_d = FD.get_cannonical_triangle_arc(u, w); B.cubic_neighbours[isomer_idx*B.n_atoms*3 + tix*3 + 2] = triangle_numbers[edge_d[0] * 6 + FD.arc_ix(edge_d[0], edge_d[1])];
    }


}}}

float kernel_time = 0.0;
std::chrono::microseconds time_spent(){
    return std::chrono::microseconds((int) (kernel_time*1000.f));
}

void reset_time(){
    kernel_time = 0.0;
}

template <Device U, typename K>
int optimal_batch_size(const int N, const int device_id) {
    cudaSetDevice(device_id);
    static size_t smem = sizeof(device_coord3d)*3*N + sizeof(device_real_t)*Block_Size_Pow_2;
    static LaunchDims dims((void*)cubic_layout_<GPU,K>, N, smem);
    dims.update_dims((void*)cubic_layout_<GPU,K>, N, smem);
    return dims.get_grid().x;
}

template <Device U, typename K>
int optimal_batch_size_2(const int N, const int device_id) {
    cudaSetDevice(device_id);
    static size_t smem = sizeof(node_t)*N*9;
    static LaunchDims dims((void*)dualize_2_<GPU,K>, N, smem);
    dims.update_dims((void*)dualize_2_<GPU,K>, N, smem);
    return dims.get_grid().x;
}
template <Device U, typename K>
cudaError_t dualize(IsomerBatch<U>& B, const LaunchCtx& ctx, const LaunchPolicy policy){
    cudaSetDevice(B.get_device_id());
    static std::vector<bool> first_call(16, true);
    static cudaEvent_t start[16], stop[16];
    float single_kernel_time = 0.0;
    auto dev = B.get_device_id();
    if(first_call[dev]) {cudaEventCreate(&start[dev]); cudaEventCreate(&stop[dev]);}

    if(policy == LaunchPolicy::SYNC) {ctx.wait();}
    else if(policy == LaunchPolicy::ASYNC && !first_call[dev]){
        cudaEventElapsedTime(&single_kernel_time, start[dev], stop[dev]);
        kernel_time += single_kernel_time;
    }

    size_t smem = sizeof(device_coord3d)*B.n_atoms*3 + sizeof(device_real_t)*Block_Size_Pow_2;
    static LaunchDims dims((void*)cubic_layout_<U,K>, B.n_atoms, smem, B.isomer_capacity);
    dims.update_dims((void*)cubic_layout_<U,K>, B.n_atoms, smem, B.isomer_capacity);
    cudaError_t error;
    void* kargs[]{(void*)&B};
    cudaEventRecord(start[dev], ctx.stream);
    error = safeCudaKernelCall((void*)cubic_layout_<U,K>, dims.get_grid(), dims.get_block(), kargs, smem, ctx.stream);  
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

template <Device U, typename K>
cudaError_t dualize_2(IsomerBatch<U>& B, const LaunchCtx& ctx, const LaunchPolicy policy){
    cudaSetDevice(B.get_device_id());
    static std::vector<bool> first_call(16, true);
    static cudaEvent_t start[16], stop[16];
    float single_kernel_time = 0.0;
    auto dev = B.get_device_id();
    if(first_call[dev]) {cudaEventCreate(&start[dev]); cudaEventCreate(&stop[dev]);}

    if(policy == LaunchPolicy::SYNC) {ctx.wait();}
    else if(policy == LaunchPolicy::ASYNC && !first_call[dev]){
        cudaEventElapsedTime(&single_kernel_time, start[dev], stop[dev]);
        kernel_time += single_kernel_time;
    }

    size_t smem = sizeof(device_node_t)*B.n_faces*9;
    int threads_rounded_to_next_multiple_of_32 = (B.n_faces + 31) & ~31;
    static LaunchDims dims((void*)dualize_2_<GPU,K>, threads_rounded_to_next_multiple_of_32, smem, B.isomer_capacity);
    dims.update_dims((void*)dualize_2_<GPU,K>, threads_rounded_to_next_multiple_of_32, smem, B.isomer_capacity);
    cudaError_t error;
    void* kargs[]{(void*)&B};
    cudaEventRecord(start[dev], ctx.stream);
    error = safeCudaKernelCall((void*)dualize_2_<GPU,K>, dims.get_grid(), dims.get_block(), kargs, smem, ctx.stream);  
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