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

//Casting functions to convert Graph objects and its' derivatives to IsomerBatch storage i.e flat memory representation.
static auto casting_coord3d = [](coord3d in){return device_coord3d{static_cast<device_real_t>(in[0]), static_cast<device_real_t>(in[1]), static_cast<device_real_t>(in[2])};};
static auto casting_node3   = [](std::vector<node_t> in){return device_node3{static_cast<device_node_t>(in[0]), static_cast<device_node_t>(in[1]), static_cast<device_node_t>(in[2])};};
static auto casting_coord2d = [](coord2d in){return device_coord2d{static_cast<device_real_t>(in.first), static_cast<device_real_t>(in.second)};};


__global__ void refill_batch_(IsomerBatch B, IsomerBatch Q_B, IsomerQueue::QueueProperties queue){

    //Create a shared queue index among all threads in a block. 
    __shared__ size_t queue_index;

    //Store a copy of the front before incrementing the counters, 
    //This is to allow for unordered incrementing of the queue counters avoiding the use of mutexes or similar techniques. 
    auto old_front = *queue.front;

    //Reset the queue requests at each call to this function.
    if ((threadIdx.x + blockIdx.x) == 0) {
        *queue.requests = 0;}
    GRID_SYNC
    //Grid stride for loop, allows for handling of any batch size.
    for (int isomer_idx = blockIdx.x; isomer_idx < B.isomer_capacity + gridDim.x - B.isomer_capacity%gridDim.x; isomer_idx+= gridDim.x){
    GRID_SYNC
    //Checks if the isomer owned by a given block is finished or empty, if it is we want to replace it with a new one.
    if (isomer_idx < B.isomer_capacity)
    if (B.statuses[isomer_idx] != NOT_CONVERGED){
        if(threadIdx.x == 0){
            //Increment front counter atomically and calculate queue index. 
            queue_index = atomicAdd(queue.front, 1);
            queue_index = queue_index % *queue.capacity;
            //Increment number of queue requests.
            atomicAdd(queue.requests, 1);
        } 
        cg::sync(cg::this_thread_block());
        assert(queue_index < *queue.capacity);

        //Given the queue index, copy data from the queue (container Q_B) to the target batch B.
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
            //Here we correct the front index by considering, what the front was at the start, how many requests were made and, what the capacity of the queue is.
            *queue.front = enough_left_in_queue ? (old_front + *queue.requests) % *queue.capacity  : *queue.back;
        }
    }
}

//This function simply copys data to the end of the batch, could be achieved with the cuda_io::copy() function using ranges.
//TODO: delete this function and replace with a copy call.
__global__ void insert_batch_(IsomerBatch B, IsomerBatch Q_B, IsomerQueue::QueueProperties queue){
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

//Allocate managed (unified) memory for queue counters.
IsomerQueue::IsomerQueue(const size_t N): N(N){
    cudaMallocManaged(&props.back,    sizeof(size_t));
    cudaMallocManaged(&props.front,   sizeof(size_t));
    cudaMallocManaged(&props.capacity,sizeof(size_t));
    cudaMallocManaged(&props.size,    sizeof(size_t));
    cudaMallocManaged(&props.requests,sizeof(size_t));
    *props.capacity = 1;
    *props.front = -1;
    *props.back = -1;
    *props.size = 0;
}

//Free counters
IsomerQueue::~IsomerQueue(){
    cudaFree(props.back);
    cudaFree(props.front);
    cudaFree(props.capacity);
    cudaFree(props.size);
    cudaFree(props.requests);
}

//Copies data to host if data is not up to date on host
cudaError_t IsomerQueue::to_host(const LaunchCtx& ctx){
    if(!is_host_updated){
        copy(host_batch, device_batch, ctx);
        is_host_updated = true;
    }
    return cudaGetLastError();
}

//Copies data to device if data is not up to date on device
cudaError_t IsomerQueue::to_device(const LaunchCtx& ctx){
    if(!is_device_updated){
        copy(device_batch, host_batch, ctx);
        is_device_updated = true;
    }
    return cudaGetLastError();
}

//Kernel wrapper for refill_batch_ function
cudaError_t IsomerQueue::refill_batch(IsomerBatch& batch, const LaunchCtx& ctx, const LaunchPolicy policy){
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

//Resizes the underlying containers (host and device batches) and updates the queue counters accordingly
cudaError_t IsomerQueue::resize(const size_t new_capacity,const LaunchCtx& ctx, const LaunchPolicy policy){
    //Lambda function because it is only called in here and avoids duplication of code because these transformations need to be applied to both host and device batches.
    auto queue_resize_batch = [&](IsomerBatch& batch){
        //Construct a tempory batch: allocates the needed amount of memory.
        IsomerBatch temp_batch = IsomerBatch(batch.n_atoms, new_capacity, batch.buffer_type);
        //Copy contents of old batch into newly allocated memory.
        copy(temp_batch,batch,ctx,policy,{0, batch.isomer_capacity- *props.front}, {*props.front,batch.isomer_capacity});
        copy(temp_batch,batch,ctx,policy,{batch.isomer_capacity- *props.front, batch.isomer_capacity - *props.front + *props.back},{0,*props.back});
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

//Inserts entire IsomerBatch in the queue.
cudaError_t IsomerQueue::insert(IsomerBatch& input_batch, const LaunchCtx& ctx, const LaunchPolicy policy, const bool insert_2d){
    //Wait for all processes on the stream to finish if policy is synchronous otherwise just proceed.
    if(policy == LaunchPolicy::SYNC) ctx.wait(); 

    //Static storage specifier of LaunchDims object allows us to perform the calculations on first function call,
    // and store the result for future calls. Thus avoiding having to perform this calculation outside and pass the state to the function.
    static LaunchDims dims((void*)insert_batch_, N, 0, input_batch.isomer_capacity);

    input_batch.buffer_type == DEVICE_BUFFER ? to_device(ctx) : to_host(ctx); 
    if (*props.capacity < (*props.size + input_batch.isomer_capacity))
    {
        *props.capacity += input_batch.isomer_capacity;
        resize(*props.capacity, ctx, policy);
    }

    //TODO: Maybe we want to implement some inefficient way of copying a host batch to the device.
    if (input_batch.buffer_type == DEVICE_BUFFER){
        void* kargs[]{(void*)&input_batch, (void*)&device_batch, (void*)&props};
        safeCudaKernelCall((void*)insert_batch_,dims.get_grid(), dims.get_block(), kargs, 0, ctx.stream);
    } else{
        std::cout << "Failed to insert: Batch insertion is intended for device batches only" << std::endl;
    }   
    if(policy == LaunchPolicy::SYNC) ctx.wait();
    printLastCudaError("Batch insertion failed");
    is_host_updated = false;
    return cudaGetLastError();
}

cudaError_t IsomerQueue::insert(const PlanarGraph& in, const size_t ID, const LaunchCtx& ctx, const LaunchPolicy policy, const bool insert_2d){
    //Before inserting a new isomer, make sure that the host batch is up to date with the device version.
    to_host(ctx);

    //If the queue is full, double the size, same beahvior as dynamically allocated containers in the std library. 
    if (*props.capacity == *props.size){
        *props.capacity *= 2;
        resize(*props.capacity, ctx, policy);
    }

    //Edge case: queue has no members 
    if(*props.back == -1){
        *props.front = 0;
        *props.back = 0;
    } else{
    *props.back = *props.back + 1 % *props.capacity;}

    //Extract the graph information (neighbours) from the PlanarGraph object and insert it at the appropriate location in the queue.
    size_t offset = *props.back * N;

    for(node_t u=0;u<N;u++){
        for(int j=0;j<3;j++)
            host_batch.neighbours[3*(offset+u)+j] = in.neighbours[u][j];

        if(insert_2d){
            host_batch.xys[2*(offset+u) + 0] = in.layout2d[u].first;
            host_batch.xys[2*(offset+u) + 1] = in.layout2d[u].second; 
        }
    }


    //Assign metadata.
    host_batch.statuses[*props.back] = NOT_CONVERGED;
    host_batch.IDs[*props.back] = ID;
    host_batch.iterations[*props.back] = 0;
    *props.size += 1;

    //Since the host batch has been updated the device is no longer up to date.
    is_device_updated = false;
    return cudaGetLastError();
}

cudaError_t IsomerQueue::insert(const Polyhedron& in, const size_t ID, const LaunchCtx& ctx, const LaunchPolicy policy, const bool insert_2d){
    //Before inserting a new isomer, make sure that the host batch is up to date with the device version.
    to_host(ctx);

    //If the queue is full, double the size, same beahvior as dynamically allocated containers in the std library. 
    if (*props.capacity == *props.size){
        *props.capacity *= 2;
        resize(*props.capacity, ctx, policy);
    }
    //Edge case: queue has no members 
    if(*props.back == -1){
        *props.front = 0;
        *props.back = 0;
    } else{
    *props.back = *props.back + 1 % *props.capacity;}

    //Extract the graph information (neighbours) from the PlanarGraph object and insert it at the appropriate location in the queue.
    size_t offset = *props.back * N;
    for(int i = 0; i < in.neighbours.size(); i++){
        if(!in.neighbours.empty())              reinterpret_cast<device_node3*>  (host_batch.neighbours)[offset +i]  = casting_node3(in.neighbours[i]);
        if(!in.points.empty())                  reinterpret_cast<device_coord3d*>(host_batch.X)[offset + i]          = casting_coord3d(in.points[i]);
        if(!in.layout2d.empty() && insert_2d)   reinterpret_cast<device_coord2d*>(host_batch.xys)[offset + i]        = casting_coord2d(in.layout2d[i]);
    }
    //Assign metadata.
    host_batch.statuses[*props.back] = NOT_CONVERGED;
    host_batch.IDs[*props.back] = ID;
    host_batch.iterations[*props.back] = 0;

    *props.size += 1;
    //Since the host batch has been updated the device is no longer up to date.
    is_device_updated = false;
    return cudaGetLastError();
}
}
