#include <cuda.h>
#include <cooperative_groups.h>
#include <cooperative_groups/reduce.h>
#include <cooperative_groups/scan.h>
#include "fullerenes/gpu/kernels.hh"
namespace gpu_kernels{
namespace isomerspace_tutte{
#include "device_includes.cu"

template cudaError_t tutte_layout<GPU, float, uint16_t>(IsomerBatch<GPU>& B, const size_t max_iterations, const LaunchCtx& ctx, const LaunchPolicy policy);
template cudaError_t tutte_layout<GPU, double, uint16_t>(IsomerBatch<GPU>& B, const size_t max_iterations, const LaunchCtx& ctx, const LaunchPolicy policy);

//WIP: Lets try to find some pentagon-pentagon distances
template <Device U> __device__
std::array<device_node_t, 12> multiple_source_shortest_paths(const IsomerBatch<U>& B, device_node_t* distances, const size_t isomer_idx){
   /*  DEVICE_TYPEDEFS;
    
    DeviceCubicGraph FG = DeviceCubicGraph(&B.cubic_neighbours[isomer_idx*blockDim.x*3]);
    node_t outer_face[6]; memset(outer_face, 0, sizeof(node_t)*6); //Do not rely on uninitialized memory it will only be zero on first touch.
    uint8_t Nface = FG.get_face_oriented(0, FG.cubic_neighbours[0],outer_face);
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
                node_t w = FG.cubic_neighbours[v*3 + i];
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
    return distance; */
}


template <Device U, typename T, typename K> __global__
void tutte_layout_(IsomerBatch<U> B, const size_t iterations){
    TEMPLATE_TYPEDEFS(T,K);
    SMEM(T);

    clear_cache(smem, Block_Size_Pow_2);
    for (int isomer_idx = blockIdx.x; isomer_idx < B.isomer_capacity; isomer_idx+= gridDim.x){
    if (B.statuses[isomer_idx] != IsomerStatus::EMPTY){
    size_t offset = isomer_idx * blockDim.x;

    DeviceCubicGraph FG(&B.cubic_neighbours[offset*3]); 
    real_t* base_pointer        = smem + Block_Size_Pow_2;
    coord2d* xys        = reinterpret_cast<coord2d*>(base_pointer);
    coord2d* newxys     = reinterpret_cast<coord2d*>(base_pointer) + blockDim.x;


    node3 ns            = (reinterpret_cast<node3*>(B.cubic_neighbours) + offset)[threadIdx.x];

    xys[threadIdx.x]    = {real_t(0.0), real_t(0.0)};

    node_t outer_face[6];
    node_t outer_face_vertex   = 0;
    uint8_t Nface = FG.get_face_oriented(0,FG.cubic_neighbours[0], outer_face);    
    /* node_t* int_smem = reinterpret_cast<node_t*>(smem);
    uint8_t local_fsize = FG.face_size(threadIdx.x, FG.cubic_neighbours[threadIdx.x*3]);
    node_t is_pentagon = node_t(local_fsize == uint8_t(5) ? 1 : 0);
    ex_scan(int_smem, is_pentagon, blockDim.x);
    if(threadIdx.x > 0 && (int_smem[threadIdx.x] == node_t(1)) && (int_smem[threadIdx.x - 1] == node_t(0))) int_smem[0] = node_t(threadIdx.x - 1);
    BLOCK_SYNC;
    node_t lowest_index_pentagon = int_smem[0];
    uint8_t Nface = FG.get_face_oriented(lowest_index_pentagon, FG.cubic_neighbours[lowest_index_pentagon*3],outer_face);
    BLOCK_SYNC; */ //This code finds the first pentagon and then uses that as the outer face. 
    reinterpret_cast<bool*>(smem)[threadIdx.x] =  false; BLOCK_SYNC;
    if(threadIdx.x < Nface){
      outer_face_vertex = outer_face[threadIdx.x];
      reinterpret_cast<bool*>(smem)[outer_face_vertex] =  true; 
    }
    BLOCK_SYNC;
    bool fixed = reinterpret_cast<bool*>(smem)[threadIdx.x];

    if(threadIdx.x < Nface) xys[outer_face_vertex] = {SIN((real_t)threadIdx.x*(real_t)2.0*real_t(M_PI)/real_t(Nface)),COS((real_t)threadIdx.x*(real_t)2.0*real_t(M_PI)/real_t(Nface))};
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

        // Calculate the new position of the point
        if(!fixed) newxys[threadIdx.x] = xys[threadIdx.x]*real_t(0.15) + (neighbour_sum/real_t(3.))*real_t(0.85);
        real_t neighbour_dist = 0.0f;

        // Calculate the distance between neighbours
        for (uint8_t j = 0; j < 3; j++) neighbour_dist += norm(xys[threadIdx.x] - xys[d_get(ns,j)])/real_t(3);
        
        BLOCK_SYNC
        real_t relative_change = 0.0f;

        // Calculate the relative change
        if (neighbour_dist > (real_t)0.0f && !fixed){ 
            relative_change = norm(xys[threadIdx.x] - newxys[threadIdx.x])/neighbour_dist;
        }

        // Reduce the relative change to find the maximum change
        real_t iteration_max = reduction_max(smem, relative_change);
        if (iteration_max > max_change) max_change = iteration_max;

        converged = max_change <= 100*numeric_limits<real_t>::epsilon();

        // Update the position of the point
        xys[threadIdx.x] = newxys[threadIdx.x];
    }
    BLOCK_SYNC;
    (reinterpret_cast<std::array<real_t,2>*>(B.xys) + offset )[threadIdx.x]  =  xys[threadIdx.x];
    }
    }
}

float kernel_time = 0.0;
std::chrono::microseconds time_spent(){
    return std::chrono::microseconds((int) (kernel_time*1000.f));
}

void reset_time(){
    kernel_time = 0.0;
}

template <Device U, typename T, typename K>
cudaError_t tutte_layout(IsomerBatch<U>& B, const size_t max_iterations, const LaunchCtx& ctx, const LaunchPolicy policy){
    TEMPLATE_TYPEDEFS(T,K);

    cudaSetDevice(B.get_device_id());
    static std::vector<bool> first_call(16, true);
    static cudaEvent_t start[16], stop[16];
    float single_kernel_time = 0.0;
    auto dev = B.get_device_id();
    if(first_call[dev]) {cudaEventCreate(&start[dev]); cudaEventCreate(&stop[dev]);}

    //If launch ploicy is synchronous then wait.
    if(policy == LaunchPolicy::SYNC){ ctx.wait();}
    else if(policy == LaunchPolicy::ASYNC && !first_call[dev]){
        //Records time from previous kernel call
        cudaEventElapsedTime(&single_kernel_time, start[dev], stop[dev]);
        kernel_time += single_kernel_time;
    }
    size_t smem = sizeof(coord2d)*B.n_atoms*2 + sizeof(real_t)*Block_Size_Pow_2;
    static LaunchDims dims((void*)tutte_layout_<U,T,K>, B.n_atoms, smem, B.isomer_capacity);
    dims.update_dims((void*)tutte_layout_<U,T,K>, B.n_atoms, smem, B.isomer_capacity);
    void* kargs[]{(void*)&B,(void*)&max_iterations};

    cudaEventRecord(start[dev], ctx.stream);
    cudaError_t error = safeCudaKernelCall((void*)tutte_layout_<U,T,K>, dims.get_grid(), dims.get_block(), kargs, smem, ctx.stream);  
    cudaEventRecord(stop[dev], ctx.stream);
    
    if(policy == LaunchPolicy::SYNC) {
        ctx.wait();
        cudaEventElapsedTime(&single_kernel_time, start[dev], stop[dev]);
        kernel_time += single_kernel_time;
    }
    printLastCudaError("Tutte: ");
    first_call[dev] = false;
    return error;
}

}}
