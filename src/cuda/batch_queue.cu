#include "fullerenes/gpu/cuda_definitions.h"
#include <inttypes.h>
#include <queue>
#include "chrono"
#include <fullerenes/polyhedron.hh>
#include "fullerenes/gpu/isomer_batch.hh"
#include "fullerenes/gpu/cuda_execution.hh"

namespace cuda_io{
    #include "misc_cuda.cu"
    #include "reductions.cu"

}

#include "cuda_execution.cu"
#include "cuda_io.cu"
#include "fullerenes/gpu/batch_queue.hh"
namespace cuda_io{
#include "print_functions.cu"

static auto casting_coord3d = [](coord3d in){return device_coord3d{static_cast<device_real_t>(in[0]), static_cast<device_real_t>(in[1]), static_cast<device_real_t>(in[2])};};
static auto casting_node3   = [](std::vector<node_t> in){return device_node3{static_cast<device_node_t>(in[0]), static_cast<device_node_t>(in[1]), static_cast<device_node_t>(in[2])};};
static auto casting_coord2d = [](coord2d in){return device_coord2d{static_cast<device_real_t>(in.first), static_cast<device_real_t>(in.second)};};


__global__ void refill_batch_(IsomerBatch B, IsomerBatch Q_B, BatchQueue::QueueProperties queue){

    __shared__ size_t queue_index;
    auto old_front = *queue.front;
    if ((threadIdx.x + blockIdx.x) == 0) {
        *queue.requests = 0; *queue.capacity = Q_B.isomer_capacity;}
    GRID_SYNC
    for (int isomer_idx = blockIdx.x; isomer_idx < B.isomer_capacity; isomer_idx+= gridDim.x){
    if (B.statuses[isomer_idx] != NOT_CONVERGED){
        if(threadIdx.x == 0){
            queue_index = atomicAdd(queue.front, 1);
            queue_index = queue_index % *queue.capacity;
            atomicAdd(queue.requests, 1);
        } 
        cg::sync(cg::this_thread_block());
        assert(queue_index < *queue.capacity);
        size_t queue_array_idx    = queue_index*blockDim.x+threadIdx.x;
        size_t global_idx         = blockDim.x*isomer_idx + threadIdx.x;
        reinterpret_cast<device_coord3d*>(B.X)[global_idx]           = reinterpret_cast<device_coord3d*>(Q_B.X)[queue_array_idx];  
        reinterpret_cast<device_node3*>(B.neighbours)[global_idx]    = reinterpret_cast<device_node3*>(Q_B.neighbours)[queue_array_idx];  
        if (threadIdx.x == 0){
            B.IDs[isomer_idx] = Q_B.IDs[queue_index];
            B.iterations[isomer_idx] = 0;
            B.statuses[isomer_idx] = Q_B.statuses[queue_index];
            Q_B.statuses[queue_index] = EMPTY;
        }
    }
    }
    GRID_SYNC
    if ((threadIdx.x + blockIdx.x) == 0) {
        //Main thread checks if the number of requests for new isomers is greater than the queue size, then either sets queuesize to 0 by subtracting queue size or subtracts the number of requests.
        bool enough_left_in_queue = *queue.size >= *queue.requests;
        *queue.size -= enough_left_in_queue ? *queue.requests : *queue.size;
        if (*queue.size == 0) {*queue.front = -1; *queue.back =-1;} else{
            *queue.front = enough_left_in_queue ? (old_front + *queue.requests) % *queue.capacity  : *queue.back;
        }
    }
}

__global__ void insert_batch_(IsomerBatch B, IsomerBatch Q_B, BatchQueue::QueueProperties queue){
    for (size_t isomer_idx = blockIdx.x; isomer_idx < B.isomer_capacity; isomer_idx+=gridDim.x){
    size_t queue_idx = (*queue.back + isomer_idx) % *queue.capacity;
    reinterpret_cast<device_coord3d*>(Q_B.X)[queue_idx*blockDim.x + threadIdx.x]        = reinterpret_cast<device_coord3d*>(B.X)[isomer_idx*blockDim.x + threadIdx.x];
    reinterpret_cast<device_node3*>(Q_B.neighbours)[queue_idx*blockDim.x + threadIdx.x] = reinterpret_cast<device_node3*>(B.neighbours)[isomer_idx*blockDim.x + threadIdx.x];
    if (threadIdx.x == 0)
    {
        Q_B.IDs[queue_idx]          = B.IDs[isomer_idx];
        Q_B.iterations[queue_idx]   = 0;
        Q_B.statuses[queue_idx]     = B.statuses[isomer_idx];
    }
    }
    GRID_SYNC
    if ((threadIdx.x + blockIdx.x) == 0) {
        *queue.size += B.isomer_capacity; *queue.back = (*queue.back + B.isomer_capacity) % *queue.capacity;
    }
}

BatchQueue::BatchQueue(const size_t N): N(N){
    cudaMallocManaged(&props.back,    sizeof(size_t));
    cudaMallocManaged(&props.front,   sizeof(size_t));
    cudaMallocManaged(&props.capacity,sizeof(size_t));
    cudaMallocManaged(&props.size,    sizeof(size_t));
    cudaMallocManaged(&props.requests,    sizeof(size_t));
    *props.capacity = 1;
    *props.front = -1;
    *props.back = -1;
    *props.size = 0;
}

BatchQueue::~BatchQueue(){
    cudaFree(props.back);
    cudaFree(props.front);
    cudaFree(props.capacity);
    cudaFree(props.size);
    cudaFree(props.requests);
}

cudaError_t BatchQueue::to_host(const LaunchCtx& ctx){
    if(!is_host_updated){
        copy(host_batch, device_batch, ctx);
        is_host_updated = true;
    }
    return cudaGetLastError();
}

cudaError_t BatchQueue::to_device(const LaunchCtx& ctx){
    if(!is_device_updated){
        copy(device_batch, host_batch, ctx);
        is_device_updated = true;
    }
    return cudaGetLastError();
}

cudaError_t BatchQueue::refill_batch(IsomerBatch& batch, const LaunchCtx& ctx, const LaunchPolicy policy){
    if (policy == LaunchPolicy::SYNC) ctx.wait();
    static LaunchDims dims((void*)refill_batch_, N, 0, batch.isomer_capacity);
    to_device(ctx);
    void* kargs[] = {(void*)&batch, (void*)&device_batch, (void*)&props};
    is_host_updated = false;
    cudaError_t error = safeCudaKernelCall((void*)refill_batch_, dims.get_grid(), dims.get_block(), kargs, 0, ctx.stream);
    if (policy == LaunchPolicy::SYNC) ctx.wait();
    printLastCudaError("Refill failed");
    return error;
}

cudaError_t BatchQueue::resize(const size_t new_capacity,const LaunchCtx& ctx, const LaunchPolicy policy){
    //Lambda function because it is only called in here and avoids duplication of code because these transformations need to be applied to both host and device batches.
    auto queue_resize_batch = [&](IsomerBatch& batch){
        //Construct a tempory batch: allocates the needed amount of memory.
        IsomerBatch temp_batch = IsomerBatch(batch.n_atoms, new_capacity, batch.buffer_type);
        //Copy contents of old batch into newly allocated memory.
        copy(temp_batch,batch,ctx,policy,{0, batch.isomer_capacity- *props.front}, {*props.front,batch.isomer_capacity});
        copy(temp_batch,batch,ctx,policy,{batch.isomer_capacity- *props.front, batch.isomer_capacity},{0,*props.back});
        for (int i = 0; i < batch.pointers.size(); i++)
        {
            void* temp_ptr = *get<1>(batch.pointers[i]);
            printLastCudaError("Free failed");
            //Reassign pointers of the input batch, to the new memory
            *get<1>(batch.pointers[i]) = *get<1>(temp_batch.pointers[i]);
            //Assign old pointers to temporary object, let destructor take care of cleanup.
            *get<1>(temp_batch.pointers[i]) = temp_ptr;
        }
        batch.isomer_capacity = temp_batch.isomer_capacity;
    };
    if (*props.back >= *props.front){
        cuda_io::resize(host_batch, new_capacity, ctx, policy);
        cuda_io::resize(device_batch, new_capacity, ctx, policy);
    }else{
        queue_resize_batch(host_batch);
        queue_resize_batch(device_batch);
        *props.front = 0;
        *props.back  = *props.size - 1;
    }
    return cudaGetLastError();
}

cudaError_t BatchQueue::insert(IsomerBatch& input_batch, const LaunchCtx& ctx, const LaunchPolicy policy, const bool insert_2d){
    if(policy == LaunchPolicy::SYNC) ctx.wait(); 
    static LaunchDims dims((void*)insert_batch_, N, 0, input_batch.isomer_capacity);
    input_batch.buffer_type == DEVICE_BUFFER ? to_device(ctx) : to_host(ctx); 
    if (*props.capacity < (*props.size + input_batch.isomer_capacity))
    {
        *props.capacity += input_batch.isomer_capacity;
        resize(*props.capacity, ctx, policy);
    }
    cudaError_t error;
    if (input_batch.buffer_type == DEVICE_BUFFER){
        void* kargs[]{(void*)&input_batch, (void*)&device_batch, (void*)&props};
        error = safeCudaKernelCall((void*)insert_batch_,dims.get_grid(), dims.get_block(), kargs, 0, ctx.stream);
    } else{
        std::cout << "Failed to insert: Batch insertion is intended for device batches only" << std::endl;
    }   
    if(policy == LaunchPolicy::SYNC) ctx.wait();
    printLastCudaError("Batch insertion failed");
    is_host_updated = false;
    return error;
}

cudaError_t BatchQueue::insert(const PlanarGraph& in, const size_t ID, const LaunchCtx& ctx, const LaunchPolicy policy, const bool insert_2d){
    to_host(ctx);
    if (*props.capacity == *props.size)
    {
        *props.capacity *= 2;
        resize(*props.capacity, ctx, policy);
    }
    if(*props.back == -1){
        *props.front = 0;
        *props.back = 0;
    } else{
    *props.back = *props.back + 1 % *props.capacity;}
    size_t offset = *props.back * N;
    for(int i = 0; i < in.neighbours.size(); i++){
        reinterpret_cast<device_node3*>  (host_batch.neighbours)[offset +i]             = casting_node3(in.neighbours[i]);
        if(insert_2d)   reinterpret_cast<device_coord2d*>(host_batch.xys)[offset + i]   = casting_coord2d(in.layout2d[i]);
    }
    host_batch.statuses[*props.back] = NOT_CONVERGED;
    host_batch.IDs[*props.back] = ID;
    host_batch.iterations[*props.back] = 0;
    *props.size += 1;
    is_device_updated = false;
    return cudaGetLastError();
}

cudaError_t BatchQueue::insert(const Polyhedron& in, const size_t ID, const LaunchCtx& ctx, const LaunchPolicy policy, const bool insert_2d){
    to_host(ctx);
    if (*props.capacity == *props.size)
    {
        *props.capacity *= 2;
        resize(*props.capacity, ctx, policy);
    }
    if(*props.back == -1){
        *props.front = 0;
        *props.back = 0;
    } else{
    *props.back = *props.back + 1 % *props.capacity;}
    size_t offset = *props.back * N;
    for(int i = 0; i < in.neighbours.size(); i++){
        if(!in.neighbours.empty())              reinterpret_cast<device_node3*>  (host_batch.neighbours)[offset +i]  = casting_node3(in.neighbours[i]);
        if(!in.points.empty())                  reinterpret_cast<device_coord3d*>(host_batch.X)[offset + i]          = casting_coord3d(in.points[i]);
        if(!in.layout2d.empty() && insert_2d)   reinterpret_cast<device_coord2d*>(host_batch.xys)[offset + i]        = casting_coord2d(in.layout2d[i]);
    }
    host_batch.statuses[*props.back] = NOT_CONVERGED;
    host_batch.IDs[*props.back] = ID;
    host_batch.iterations[*props.back] = 0;

    *props.size += 1;
    is_device_updated = false;
    return cudaGetLastError();
}
}
