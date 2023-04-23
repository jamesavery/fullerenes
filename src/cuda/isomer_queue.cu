#include "fullerenes/gpu/cuda_definitions.h"
#include <cooperative_groups.h>
#include <cooperative_groups/reduce.h>
#include <cooperative_groups/scan.h>
#include <inttypes.h>
#include <queue>
#include "chrono"
#include <fullerenes/polyhedron.hh>
#include "fullerenes/gpu/misc_cuda.cuh"
#include "fullerenes/gpu/launch_ctx.hh"
#include "launch_dims.cu"
#include "fullerenes/gpu/isomer_batch.hh"
#include "fullerenes/gpu/cuda_io.hh"
#include "isomer_batch.cu"
#include "cub/device/device_scan.cuh"

namespace cuda_io{
    
    #include "reductions.cu"

}

#include "launch_ctx.cu"
#include "cuda_io.cu"
#include "fullerenes/gpu/isomer_queue.hh"
namespace cuda_io{
#include "print_functions.cu"

//Casting functions to convert Graph objects and its' derivatives to IsomerBatch storage i.e flat memory representation.
static auto casting_coord3d = [](coord3d in){return device_coord3d{static_cast<float>(in[0]), static_cast<float>(in[1]), static_cast<float>(in[2])};};
static auto casting_node3   = [](std::vector<node_t> in){return device_node3{static_cast<device_node_t>(in[0]), static_cast<device_node_t>(in[1]), static_cast<device_node_t>(in[2])};};
static auto casting_coord2d = [](coord2d in){return device_coord2d{static_cast<float>(in.first), static_cast<float>(in.second)};};


__global__ void refill_batch_(IsomerBatch B, IsomerBatch Q_B, IsomerQueue::QueueProperties queue,  int* scan_array){
    auto num_inside_circular_range = [](const int begin, const int end, const int test_val){
        if(begin < 0 || end < 0) {return false;}
        if (begin <= end) {return (test_val >= begin && test_val <= end);}
        else {return (test_val >= begin || test_val <= end);}
    };

    //Must ensure that all writes to queue counters from the host are visible to the device threads before reading them.
    __threadfence_system();
    DEVICE_TYPEDEFS
    extern __shared__ int smem[];
    auto Nf = B.n_faces; // Number of faces
    auto queue_requests = 0;
    //Grid stride for loop, allows for handling of any batch size.
    auto limit = ((B.isomer_capacity + gridDim.x - 1) / gridDim.x ) * gridDim.x;  //Fast ceiling integer division.
    for (int isomer_idx = blockIdx.x; isomer_idx < limit; isomer_idx+= gridDim.x){
    GRID_SYNC
    bool access_queue = false;
    if (isomer_idx < B.isomer_capacity) access_queue = B.statuses[isomer_idx] != IsomerStatus::NOT_CONVERGED;
    
    //Perform a grid scan over the array of 0s and 1s (1 if a block needs to access the queue and 0 otherwise).
    //If the grid is larger than the batch capacity, we only want to scan over the interval [0,B.isomer_capacity-1]
    grid_ex_scan<int>(scan_array, smem, (int)access_queue, min((int)B.isomer_capacity, (int)gridDim.x) + 1);
    
    //The block gets its uniqueue index from the back of the queue plust its position in the scan array.
    auto queue_index = (*queue.front + scan_array[blockIdx.x] + queue_requests) % *queue.capacity; 
    
    //Check that the index is inside the queue 
    access_queue &= num_inside_circular_range(*queue.front, *queue.back, queue_index);
    
    if (access_queue){
        //Given the queue index, copy data from the queue (container Q_B) to the target batch B.
        size_t queue_array_idx    = queue_index*blockDim.x+threadIdx.x;
        size_t global_idx         = blockDim.x*isomer_idx + threadIdx.x;
        reinterpret_cast<coord3d*>(B.X)[global_idx]                 = reinterpret_cast<coord3d*>(Q_B.X)[queue_array_idx];  
        reinterpret_cast<node3*>(B.cubic_neighbours)[global_idx]    = reinterpret_cast<node3*>(Q_B.cubic_neighbours)[queue_array_idx];
        reinterpret_cast<coord2d*>(B.xys)[global_idx]               = reinterpret_cast<coord2d*>(Q_B.xys)[queue_array_idx];
        
        //Face parallel copying 
        if (threadIdx.x < Nf){
            size_t queue_face_idx     = queue_index* Nf + threadIdx.x; 
            size_t output_face_idx    = isomer_idx * Nf + threadIdx.x;
            reinterpret_cast<node6*>(B.dual_neighbours)[output_face_idx]               = reinterpret_cast<node6*>(Q_B.dual_neighbours)[queue_face_idx];  
            reinterpret_cast<uint8_t*>(B.face_degrees)[output_face_idx]                = reinterpret_cast<uint8_t*>(Q_B.face_degrees)[queue_face_idx];  
        }
        //Per isomer meta data
        if (threadIdx.x == 0){
            B.IDs[isomer_idx] = Q_B.IDs[queue_index];
            B.iterations[isomer_idx] = 0;
            B.statuses[isomer_idx] = Q_B.statuses[queue_index];
            Q_B.statuses[queue_index] = IsomerStatus::EMPTY;
        }
    }
    queue_requests += scan_array[min((int)B.isomer_capacity, (int)gridDim.x)]; //The last element of the scanned array will tell us how many blocks, accessed the queue.
    }
    GRID_SYNC
    if ((threadIdx.x + blockIdx.x) == 0) {
        //Main thread checks if the number of requests for new isomers is greater than the queue size, then either sets queuesize to 0 by subtracting queue size or subtracts the number of requests.
        bool enough_left_in_queue = *queue.size >= queue_requests;
        atomicSub_system(queue.size, enough_left_in_queue ? queue_requests : *queue.size);
        if (*queue.size == 0) {
            atomicExch_system(queue.front,-1);
            atomicExch_system(queue.back, -1);
            } else{
            //Here we correct the front index by considering, what the front was at the start, how many requests were made and, what the capacity of the queue is.
            atomicExch_system(queue.front, enough_left_in_queue ? (*queue.front + queue_requests) % *queue.capacity  : *queue.back);
        }
    }
}

__global__ void push_(IsomerBatch B, IsomerBatch Q_B, IsomerQueue::QueueProperties queue,  int* scan_array){
    auto num_inside_capacity = [](const int size, const int capacity, const int test_val){
        if(size < 0 || capacity < 0 || test_val < 0) {return false;}
        return size + test_val <= capacity;
    };

    //Must ensure that all writes to queue counters from the host are visible to the device threads before reading them.
    __threadfence_system();
    DEVICE_TYPEDEFS
    extern __shared__ int smem[];
 
    auto Nf = B.n_faces; // Number of faces
    auto queue_requests = 0;
    //Grid stride for loop, allows for handling of any batch size.
    auto limit = ((B.isomer_capacity + gridDim.x - 1) / gridDim.x ) * gridDim.x;  //Fast ceiling integer division.
    for (int isomer_idx = blockIdx.x; isomer_idx < limit; isomer_idx+= gridDim.x){
    GRID_SYNC
    bool access_queue = false;
    if (isomer_idx < B.isomer_capacity) access_queue = (B.statuses[isomer_idx] == IsomerStatus::CONVERGED) || (B.statuses[isomer_idx] == IsomerStatus::FAILED);
    
    grid_ex_scan<int>(scan_array, smem, (int)access_queue, min((int)B.isomer_capacity, (int)gridDim.x) + 1);
    //Checks if the isomer owned by a given block is finished or empty, if it is we want to replace it with a new one.


    int queue_index = *queue.back < 0 ? (scan_array[blockIdx.x] + queue_requests) : (*queue.back + 1 + queue_requests + scan_array[blockIdx.x]) % *queue.capacity;
    access_queue &= num_inside_capacity(*queue.size, *queue.capacity, queue_requests + scan_array[blockIdx.x]);
    if (access_queue){
        
        //Given the queue index, copy data from the queue (container Q_B) to the target batch B.
        size_t queue_array_idx    = queue_index*blockDim.x+threadIdx.x;
        size_t global_idx         = blockDim.x*isomer_idx + threadIdx.x;
        reinterpret_cast<coord3d*>(Q_B.X)[queue_array_idx]              = reinterpret_cast<coord3d*>(B.X)[global_idx];
        reinterpret_cast<node3*>(Q_B.cubic_neighbours)[queue_array_idx] = reinterpret_cast<node3*>(B.cubic_neighbours)[global_idx];
        reinterpret_cast<coord2d*>(Q_B.xys)[queue_array_idx]            = reinterpret_cast<coord2d*>(B.xys)[global_idx];
        
        //Face parallel copying 
        if (threadIdx.x < Nf){
            size_t queue_face_idx     = queue_index* Nf + threadIdx.x; 
            size_t output_face_idx    = isomer_idx * Nf + threadIdx.x;
            reinterpret_cast<node6*>(Q_B.dual_neighbours)[queue_face_idx] = reinterpret_cast<node6*>(B.dual_neighbours)[output_face_idx];
            reinterpret_cast<uint8_t*>(Q_B.face_degrees)[queue_face_idx]  = reinterpret_cast<uint8_t*>(B.face_degrees)[output_face_idx];  
        }
        //Per isomer meta data
        if (threadIdx.x == 0){
            Q_B.IDs[queue_index] = B.IDs[isomer_idx];
            Q_B.iterations[queue_index] = B.iterations[isomer_idx];
            Q_B.statuses[queue_index] = B.statuses[isomer_idx];
            B.statuses[isomer_idx] = IsomerStatus::EMPTY;
        }
    }
    queue_requests += scan_array[min((int)B.isomer_capacity, (int)gridDim.x)]; //The last element of the scanned array will tell us how many blocks, accessed the queue.
    }
    GRID_SYNC
    if ((threadIdx.x + blockIdx.x) == 0) {
        if (*queue.size == 0 && queue_requests > 0) {
            atomicExch_system(queue.front, 0); atomicExch_system(queue.back, queue_requests-1);} 
        else{
            atomicExch_system(queue.back, (*queue.back + queue_requests) % *queue.capacity);
        }
        //Pushing to queue simply increases the size by the number of push requests.
        atomicAdd_system(queue.size, queue_requests);
    }
}



//This function simply copys data to the end of the batch, could be achieved with the cuda_io::copy() function using ranges.
//TODO: delete this function and replace with a copy call.
__global__ void insert_batch_(IsomerBatch B, IsomerBatch Q_B, IsomerQueue::QueueProperties queue){
    __threadfence_system();
    //Grid stride for loop, allows for handling of any batch size.
    auto limit = ((B.isomer_capacity + gridDim.x - 1) / gridDim.x ) * gridDim.x;  //Fast ceiling integer division.
    for (int isomer_idx = blockIdx.x; isomer_idx < limit; isomer_idx+= gridDim.x){
    GRID_SYNC
        if (isomer_idx < B.isomer_capacity){
            int queue_idx = *queue.back < 0 ? isomer_idx% *queue.capacity : (*queue.back + 1 + isomer_idx) % *queue.capacity;

            reinterpret_cast<device_coord3d*>(Q_B.X)[queue_idx*blockDim.x + threadIdx.x]              = reinterpret_cast<device_coord3d*>(B.X)[isomer_idx*blockDim.x + threadIdx.x];
            reinterpret_cast<device_node3*>(Q_B.cubic_neighbours)[queue_idx*blockDim.x + threadIdx.x] = reinterpret_cast<device_node3*>(B.cubic_neighbours)[isomer_idx*blockDim.x + threadIdx.x];
            auto Nf = B.n_faces;
            if (threadIdx.x < Nf){
                reinterpret_cast<device_node6*>(Q_B.dual_neighbours)[queue_idx*Nf + threadIdx.x] = reinterpret_cast<device_node6*>(B.dual_neighbours)[isomer_idx*Nf + threadIdx.x];
                reinterpret_cast<uint8_t*>(Q_B.face_degrees)     [queue_idx*Nf + threadIdx.x]    = reinterpret_cast<uint8_t*>(B.face_degrees)     [isomer_idx*Nf + threadIdx.x];
            }
            if (threadIdx.x == 0)
            {
                Q_B.IDs[queue_idx]          = B.IDs[isomer_idx];
                Q_B.iterations[queue_idx]   = 0;
                Q_B.statuses[queue_idx]     = B.statuses[isomer_idx];
                B.statuses[isomer_idx]      = IsomerStatus::EMPTY;
            }
        }
    }
    GRID_SYNC
    if ((threadIdx.x + blockIdx.x) == 0) {
        //Edge case: queue has no members 
        if(*queue.back == -1){

            atomicExch_system(queue.front,0);
            atomicExch_system(queue.back,B.isomer_capacity - 1);
            atomicExch_system(queue.size,B.isomer_capacity);
        } else{
            atomicAdd_system(queue.size, B.isomer_capacity);
            atomicExch_system(queue.back, (*queue.back + B.isomer_capacity) % *queue.capacity);
        }
    }
}

//Allocate managed (unified) memory for queue counters.
IsomerQueue::IsomerQueue(const size_t N, int device): N(N), m_device(device){
    cudaSetDevice(device);
    cudaDeviceProp d_props;
    cudaGetDeviceProperties(&d_props, device);
    cudaMallocManaged(&props.back,    sizeof(size_t));
    cudaMallocManaged(&props.front,   sizeof(size_t));
    cudaMallocManaged(&props.capacity,sizeof(size_t));
    cudaMallocManaged(&props.size,    sizeof(size_t));
    cudaMallocManaged(&props.requests,sizeof(size_t));
    *props.capacity = 1;
    *props.front = -1;
    *props.back = -1;
    *props.size = 0;
    device_batch = IsomerBatch(N,1,DEVICE_BUFFER,m_device);
    cudaMalloc(&g_scan_array, sizeof(int)*d_props.maxBlocksPerMultiProcessor*d_props.multiProcessorCount * 2);
}

IsomerQueue::IsomerQueue(): N(20), m_device(0){
    std::cout << "Warning: Default constructor for IsomerQueue is not recommended. Use IsomerQueue(const size_t N, int device) instead." << std::endl;
    cudaSetDevice(m_device);
    cudaDeviceProp d_props;
    cudaGetDeviceProperties(&d_props, m_device);
    cudaMallocManaged(&props.back,    sizeof(size_t));
    cudaMallocManaged(&props.front,   sizeof(size_t));
    cudaMallocManaged(&props.capacity,sizeof(size_t));
    cudaMallocManaged(&props.size,    sizeof(size_t));
    cudaMallocManaged(&props.requests,sizeof(size_t));
    *props.capacity = 1;
    *props.front = -1;
    *props.back = -1;
    *props.size = 0;
    device_batch = IsomerBatch(N,1,DEVICE_BUFFER,m_device);
    cudaMalloc(&g_scan_array, sizeof(int)*d_props.maxBlocksPerMultiProcessor*d_props.multiProcessorCount * 2);
}

//Free counters
IsomerQueue::~IsomerQueue(){
    cudaSetDevice(m_device);
    
    cudaFree(props.back);
    cudaFree(props.front);
    cudaFree(props.capacity);
    cudaFree(props.size);
    cudaFree(props.requests);
    cudaFree(g_scan_array);
}

int IsomerQueue::get_size() const {cudaSetDevice(m_device);  return *props.size;};
int IsomerQueue::get_front() const {cudaSetDevice(m_device);  return *props.front;};
int IsomerQueue::get_back() const {cudaSetDevice(m_device);  return *props.back;};
int IsomerQueue::get_capacity() const {cudaSetDevice(m_device);  return *props.capacity;};
int IsomerQueue::get_requests() const {cudaSetDevice(m_device);  return *props.requests;};


//Copies data to host if data is not up to date on host
cudaError_t IsomerQueue::to_host(const LaunchCtx& ctx){
    if(!is_host_updated){
        cudaSetDevice(m_device);
        cuda_io::copy(host_batch, device_batch, ctx);
        is_host_updated = true;
    }
    return cudaGetLastError();
}

//Copies data to device if data is not up to date on device
cudaError_t IsomerQueue::to_device(const LaunchCtx& ctx){
    if(!is_device_updated){
        cudaSetDevice(m_device);
        cuda_io::copy(device_batch, host_batch, ctx);
        is_device_updated = true;
    }
    return cudaGetLastError();
}

//Kernel wrapper for refill_batch_ function
cudaError_t IsomerQueue::refill_batch(IsomerBatch& batch, const LaunchCtx& ctx, const LaunchPolicy policy){
    cudaSetDevice(m_device);
    int smem_bytes = sizeof(int)*N*2;
    if (policy == LaunchPolicy::SYNC) ctx.wait();
    static LaunchDims dims((void*)refill_batch_, N, smem_bytes, batch.isomer_capacity);
    dims.update_dims((void*)refill_batch_, N, smem_bytes, batch.isomer_capacity);
    to_device(ctx);
    void* kargs[] = {(void*)&batch, (void*)&device_batch, (void*)&props, (void*)&g_scan_array};
    is_host_updated = false;
    cudaError_t error = safeCudaKernelCall((void*)refill_batch_, dims.get_grid(), dims.get_block(), kargs, smem_bytes, ctx.stream);
    if (policy == LaunchPolicy::SYNC) ctx.wait();
    cudaMemPrefetchAsync(props.size, sizeof(size_t), cudaCpuDeviceId, ctx.stream);
    printLastCudaError("Refill failed");
    return error;
}

//Kernel wrapper for push_ function
cudaError_t IsomerQueue::push(IsomerBatch& batch, const LaunchCtx& ctx, const LaunchPolicy policy){
    cudaSetDevice(m_device);
    int smem_bytes = sizeof(int)*N*2;
    if (policy == LaunchPolicy::SYNC) ctx.wait();
    static LaunchDims dims((void*)push_, N, smem_bytes, batch.isomer_capacity);
    dims.update_dims((void*)push_, N, smem_bytes, batch.isomer_capacity);
    to_device(ctx);
    if(*props.capacity < *props.size + batch.isomer_capacity)  resize(*props.size+batch.isomer_capacity, ctx, policy);
    void* kargs[] = {(void*)&batch, (void*)&device_batch, (void*)&props, (void*)&g_scan_array};
    cudaError_t error = safeCudaKernelCall((void*)push_, dims.get_grid(), dims.get_block(), kargs, smem_bytes, ctx.stream);
    is_host_updated = false;
    if (policy == LaunchPolicy::SYNC) ctx.wait();
    //cudaMemPrefetchAsync(props.size, sizeof(size_t), cudaCpuDeviceId, ctx.stream);
    //cudaMemPrefetchAsync(props.back, sizeof(size_t), cudaCpuDeviceId, ctx.stream);
    //cudaMemPrefetchAsync(props.front, sizeof(size_t), cudaCpuDeviceId, ctx.stream);
    printLastCudaError("Push failed");
    return error;
}

Polyhedron IsomerQueue::pop(const LaunchCtx& ctx, const LaunchPolicy policy){
    cudaSetDevice(m_device);
    to_host(ctx);
    if(*props.size == 0) throw std::runtime_error("Queue is empty");
    neighbours_t out_neighbours(N, std::vector<int>(3));
    for(int i = 0; i < N; i++){
        out_neighbours[i][0] = host_batch.cubic_neighbours[(*props.front)*N*3 + i*3 + 0];
        out_neighbours[i][1] = host_batch.cubic_neighbours[(*props.front)*N*3 + i*3 + 1];
        out_neighbours[i][2] = host_batch.cubic_neighbours[(*props.front)*N*3 + i*3 + 2];
    }
    std::vector<coord3d> out_coords(N);
    for(int i = 0; i < N; i++){
        out_coords[i][0] = host_batch.X[(*props.front)*N*3 + i*3 + 0];
        out_coords[i][1] = host_batch.X[(*props.front)*N*3 + i*3 + 1];
        out_coords[i][2] = host_batch.X[(*props.front)*N*3 + i*3 + 2];
    }
    *props.front = (*props.front + 1) % *props.capacity;
    *props.size -= 1;
    return Polyhedron(Graph(out_neighbours,true),out_coords);
}

//Resizes the underlying containers (host and device batches) and updates the queue counters accordingly
cudaError_t IsomerQueue::resize(const size_t new_capacity,const LaunchCtx& ctx, const LaunchPolicy policy){
    cudaSetDevice(m_device);
    //Lambda function because it is only called in here and avoids duplication of code because these transformations need to be applied to both host and device batches.
    auto queue_resize_batch = [&](IsomerBatch& batch){
        //Construct a tempory batch: allocates the needed amount of memory.
        IsomerBatch temp_batch = IsomerBatch(batch.n_atoms, new_capacity, batch.buffer_type, batch.get_device_id());
        //Copy contents of old batch into newly allocated memory.
        cuda_io::copy(temp_batch,batch,ctx,policy,{0, batch.isomer_capacity - *props.front}, {*props.front,batch.isomer_capacity});
        cuda_io::copy(temp_batch,batch,ctx,policy,{batch.isomer_capacity- *props.front, batch.isomer_capacity - *props.front + *props.back},{0,*props.back});
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
    *props.capacity = new_capacity;
    is_host_updated = true;
    is_device_updated = true;
    return cudaGetLastError();
}

//Inserts entire IsomerBatch in the queue.
cudaError_t IsomerQueue::insert(IsomerBatch& input_batch, const LaunchCtx& ctx, const LaunchPolicy policy, const bool insert_2d){
    cudaSetDevice(m_device);
    //Wait for all processes on the stream to finish if policy is synchronous otherwise just proceed.
    if(policy == LaunchPolicy::SYNC) ctx.wait(); 

    //Static storage specifier of LaunchDims object allows us to perform the calculations on first function call,
    // and store the result for future calls. Thus avoiding having to perform this calculation outside and pass the state to the function.
    static LaunchDims dims((void*)insert_batch_, N, 0, input_batch.isomer_capacity);
    dims.update_dims((void*)insert_batch_, N, 0, input_batch.isomer_capacity);

    input_batch.buffer_type == DEVICE_BUFFER ? to_device(ctx) : to_host(ctx); 
    if (*props.capacity < (*props.size + input_batch.isomer_capacity))
    {
        auto new_capacity = *props.size + input_batch.isomer_capacity;
        resize(new_capacity, ctx, policy);
    }

    //TODO: Maybe we want to implement some inefficient way of copying a host batch to the device.
    if (input_batch.buffer_type == DEVICE_BUFFER){
        void* kargs[]{(void*)&input_batch, (void*)&device_batch, (void*)&props};
        safeCudaKernelCall((void*)insert_batch_,dims.get_grid(), dims.get_block(), kargs, 0, ctx.stream);
        //cudaMemPrefetchAsync(props.size, sizeof(size_t), cudaCpuDeviceId, ctx.stream);
    } else{
        std::cout << "Failed to insert: Batch insertion is intended for device batches only" << std::endl;
    }   
    if(policy == LaunchPolicy::SYNC) ctx.wait();
    printLastCudaError("Batch insertion failed");
    is_host_updated = false;
    is_device_updated = true;
    return cudaGetLastError();
}

cudaError_t IsomerQueue::insert(const Graph& in, const size_t ID, const LaunchCtx& ctx, const LaunchPolicy policy){
    cudaSetDevice(m_device);
    //Before inserting a new isomer, make sure that the host batch is up to date with the device version.
    to_host(ctx);

    //If the queue is full, double the size, same beahvior as dynamically allocated containers in the std library. 
    if (*props.capacity == *props.size){
        resize(*props.capacity * 2, ctx, policy);
    }

    //Edge case: queue has no members 
    if(*props.back == -1){
        *props.front = 0;
        *props.back = 0;
    } else{
    *props.back = *props.back + 1 % *props.capacity;}

    //Extract the graph information (neighbours) from the PlanarGraph object and insert it at the appropriate location in the queue.
    auto Nf = N / 2 + 2;
    size_t face_offset = *props.back * Nf;

    for(node_t u=0;u<in.neighbours.size();u++){
        host_batch.face_degrees[face_offset + u] = in.neighbours[u].size();
        for(int j=0;j<in.neighbours[u].size();j++){
            host_batch.dual_neighbours[6*(face_offset+u)+j] = in.neighbours[u][j];
        }
    }


    //Assign metadata.
    host_batch.statuses[*props.back] = IsomerStatus::NOT_CONVERGED;
    host_batch.IDs[*props.back] = ID;
    host_batch.iterations[*props.back] = 0;
    *props.size += 1;

    //Since the host batch has been updated the device is no longer up to date.
    is_device_updated = false;
    return cudaGetLastError();
}

cudaError_t IsomerQueue::insert(const PlanarGraph& in, const size_t ID, const LaunchCtx& ctx, const LaunchPolicy policy, const bool insert_2d){
    cudaSetDevice(m_device);
    //Before inserting a new isomer, make sure that the host batch is up to date with the device version.
    to_host(ctx);

    //If the queue is full, double the size, same beahvior as dynamically allocated containers in the std library. 
    if (*props.capacity == *props.size){
        resize(*props.capacity * 2, ctx, policy);
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
        for(int j=0;j<3;j++){
            host_batch.cubic_neighbours[3*(offset+u)+j] = in.neighbours[u][j];
        }
        if(insert_2d){
            host_batch.xys[2*(offset+u) + 0] = in.layout2d[u].first;
            host_batch.xys[2*(offset+u) + 1] = in.layout2d[u].second; 
        }
    }


    //Assign metadata.
    host_batch.statuses[*props.back] = IsomerStatus::NOT_CONVERGED;
    host_batch.IDs[*props.back] = ID;
    host_batch.iterations[*props.back] = 0;
    *props.size += 1;

    //Since the host batch has been updated the device is no longer up to date.
    is_device_updated = false;
    return cudaGetLastError();
}

cudaError_t IsomerQueue::insert(const Polyhedron& in, const size_t ID, const LaunchCtx& ctx, const LaunchPolicy policy, const bool insert_2d){
    cudaSetDevice(m_device);
    //Before inserting a new isomer, make sure that the host batch is up to date with the device version.
    to_host(ctx);

    //If the queue is full, double the size, same beahvior as dynamically allocated containers in the std library. 
    if (*props.capacity == *props.size){
        resize(*props.capacity * 2, ctx, policy);
    }
    //Edge case: queue has no members 
    if(*props.back == -1){
        *props.front =  0;
        *props.back = 0;
    } else{
    *props.back =  *props.back + 1 % *props.capacity;}

    //Extract the graph information (neighbours) from the PlanarGraph object and insert it at the appropriate location in the queue.
    size_t offset = *props.back * N;
    for(int i = 0; i < in.neighbours.size(); i++){
        if(!in.neighbours.empty())              reinterpret_cast<device_node3*>  (host_batch.cubic_neighbours)[offset +i]  = casting_node3(in.neighbours[i]);
        if(!in.points.empty())                  reinterpret_cast<device_coord3d*>(host_batch.X)[offset + i]          = casting_coord3d(in.points[i]);
        if(!in.layout2d.empty() && insert_2d)   reinterpret_cast<device_coord2d*>(host_batch.xys)[offset + i]        = casting_coord2d(in.layout2d[i]);
    }
    //Assign metadata.
    host_batch.statuses[*props.back] = IsomerStatus::NOT_CONVERGED;
    host_batch.IDs[*props.back] = ID;
    host_batch.iterations[*props.back] = 0;

    *props.size += 1;
    //Since the host batch has been updated the device is no longer up to date.
    is_device_updated = false;
    return cudaGetLastError();
}

cudaError_t IsomerQueue::clear(const LaunchCtx& ctx, const LaunchPolicy policy){
    if (policy == LaunchPolicy::SYNC) { cudaStreamSynchronize(ctx.stream);}
    cudaSetDevice(m_device);
    //Before clearing the queue, make sure that the host batch is up to date with the device version.
    //Clear the queue.
    *props.size = 0;
    *props.front = -1;
    *props.back = -1;
    //cudaMemset(host_batch.statuses, int(IsomerStatus::EMPTY), *props.capacity * sizeof(IsomerStatus));
    //cudaMemset(host_batch.iterations, 0, *props.capacity * sizeof(size_t));
    //cudaMemset(device_batch.statuses, int(IsomerStatus::EMPTY), *props.capacity * sizeof(IsomerStatus));
    //cudaMemset(device_batch.iterations, 0, *props.capacity * sizeof(size_t));
    if(policy == LaunchPolicy::SYNC){ cudaStreamSynchronize(ctx.stream); }
    //cudaMemPrefetchAsync(props.size, sizeof(size_t), m_device, ctx.stream);
    //cudaMemPrefetchAsync(props.front, sizeof(size_t), m_device, ctx.stream);
    //cudaMemPrefetchAsync(props.back, sizeof(size_t), m_device, ctx.stream);
    return cudaGetLastError();
}

std::ostream& operator << (std::ostream& os, const IsomerQueue& input){
    cudaSetDevice(input.m_device);
    cudaDeviceSynchronize();
    os << "Queue [Size, Capacity]: (" << input.get_size() << " \t, " << input.get_capacity() << ")\n";
    os << "Queue [Front, Back]   : (" << input.get_front() << " \t, " << input.get_back() << ")\n";
    return os;
}



void IsomerQueue::operator=(const IsomerQueue& other){
    std::cout << "Copying queue from device " << other.m_device << " to device " << m_device << std::endl;
    cudaSetDevice(m_device);
    cudaFree(props.size);
    cudaFree(props.capacity);
    cudaFree(props.front);
    cudaFree(props.back);
    cudaFree(props.requests);
    cudaFree(g_scan_array); 

    cudaSetDevice(other.m_device);
    cudaDeviceProp d_props;
    cudaGetDeviceProperties(&d_props, other.m_device);
    device_batch = other.device_batch;
    host_batch = other.host_batch;
    cudaMallocManaged(&props.size, sizeof(size_t));
    cudaMallocManaged(&props.capacity, sizeof(size_t));
    cudaMallocManaged(&props.front, sizeof(size_t));
    cudaMallocManaged(&props.back, sizeof(size_t));
    cudaMallocManaged(&props.requests, sizeof(size_t));
    cudaMalloc(&g_scan_array, sizeof(int)*d_props.maxBlocksPerMultiProcessor*d_props.multiProcessorCount * 2);
    cudaDeviceSynchronize();
    //Copy the properties of the other queue.
    *props.size = *other.props.size;
    *props.capacity = *other.props.capacity;
    *props.front = *other.props.front;
    *props.back = *other.props.back;

    is_device_updated = true;
    is_host_updated = true;
    m_device = other.m_device;
}


}
