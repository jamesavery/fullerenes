#include "fullerenes/gpu/cuda_definitions.h"
#include <inttypes.h>
#include <queue>
#include <fullerenes/polyhedron.hh>
#include "fullerenes/gpu/isomer_batch.hh"

namespace cuda_io{
    #include "misc_cuda.cu"
}

#include "cuda_io.cu"
#include "fullerenes/gpu/batch_queue.hh"
namespace cuda_io{
__global__ void __refill_batch(IsomerBatch b, IsomerBatch q_b, BatchQueue::QueueProperties queue){

    __shared__ size_t queue_index;
    if ((threadIdx.x + blockIdx.x) == 0) {*queue.requests = 0; *queue.capacity = q_b.isomer_capacity;}
    GRID_SYNC
    for (int isomer_idx = blockIdx.x; isomer_idx < b.isomer_capacity; isomer_idx+= gridDim.x){
    if (b.statuses[isomer_idx] != NOT_CONVERGED){
        if(threadIdx.x == 0){
            queue_index = atomicAdd(queue.front, 1);
            queue_index = queue_index % *queue.capacity;
            atomicAdd(queue.requests, 1);
        } 
        cg::sync(cg::this_thread_block());
        assert(queue_index < *queue.capacity);
        size_t queue_array_idx    = queue_index*blockDim.x+threadIdx.x;
        size_t global_idx         = blockDim.x*isomer_idx + threadIdx.x;
        reinterpret_cast<device_coord3d*>(b.X)[global_idx]           = reinterpret_cast<device_coord3d*>(q_b.X)[queue_array_idx];  
        reinterpret_cast<device_node3*>(b.neighbours)[global_idx]    = reinterpret_cast<device_node3*>(q_b.neighbours)[queue_array_idx];  
        if (threadIdx.x == 0){
            b.IDs[isomer_idx] = q_b.IDs[queue_index];
            b.iterations[isomer_idx] = 0;
            b.statuses[isomer_idx] = q_b.statuses[queue_index];
            q_b.statuses[queue_index] = EMPTY;
        }
    }
    }
    GRID_SYNC
    if ((threadIdx.x + blockIdx.x) == 0) {
        //Main thread checks if the number of requests for new isomers is greater than the queue size, then either sets queuesize to 0 by subtracting queue size or subtracts the number of requests.
        bool enough_left_in_queue = *queue.size >= *queue.requests;
        *queue.size -= enough_left_in_queue ? *queue.requests : *queue.size;
        *queue.front = enough_left_in_queue ? (*queue.back - *queue.size + *queue.capacity) % *queue.capacity : *queue.back;
    }
}

BatchQueue::BatchQueue(const size_t N): N(N){
    cudaMallocManaged(&props.back,    sizeof(size_t));
    cudaMallocManaged(&props.front,   sizeof(size_t));
    cudaMallocManaged(&props.capacity,sizeof(size_t));
    cudaMallocManaged(&props.size,    sizeof(size_t));
    cudaMalloc(&props.requests,    sizeof(size_t));
    *props.capacity = 1;
    *props.front = 0;
    *props.back = 0;
    *props.size = 0;
}

BatchQueue::~BatchQueue(){
    cudaFree(props.back);
    cudaFree(props.front);
    cudaFree(props.capacity);
    cudaFree(props.size);
    cudaFree(props.requests);
}

cudaError_t BatchQueue::to_host(const cudaStream_t stream){
    if(!is_host_updated){
        copy(host_batch, device_batch, stream);
        is_host_updated = true;
    }
    return cudaGetLastError();
}

cudaError_t BatchQueue::to_device(const cudaStream_t stream){
    if(!is_device_updated){
        copy(device_batch, host_batch, stream);
        is_device_updated = true;
    }
    return cudaGetLastError();
}

cudaError_t BatchQueue::refill_batch(IsomerBatch& batch, const cudaStream_t stream){
    static LaunchDims dims((void*)__refill_batch, N, 0);
    to_device(stream);
    void* kargs[] = {(void*)&batch, (void*)&device_batch, (void*)&props};
    is_host_updated = false;
    cudaError_t error = safeCudaKernelCall((void*)__refill_batch, dims.get_grid(), dims.get_block(), kargs, 0, stream);
    printLastCudaError("Refill failed");
    return error;
}

cudaError_t BatchQueue::insert(IsomerBatch& input_batch, const cudaStream_t stream, const bool insert_2d){
    reset_convergence_statuses(input_batch, stream);   
    input_batch.buffer_type == DEVICE_BUFFER ? to_device(stream) : to_host(stream); 
    if (*props.capacity < (*props.size + input_batch.isomer_capacity))
    {
        *props.capacity += input_batch.isomer_capacity;
        resize(host_batch, *props.capacity, stream);
        resize(device_batch, *props.capacity, stream);
    }
    if (input_batch.buffer_type == DEVICE_BUFFER){
        for (size_t i = 0; i < input_batch.pointers.size(); i++)
        {
            auto& [name, insert_ptr, bytes, single_isomer_per_batch] = input_batch.pointers[i];
            size_t num_elements = get<3>(input_batch.pointers[i]) ?  input_batch.n_atoms * input_batch.n_isomers : input_batch.n_isomers;
            cudaMemcpyAsync(((char*)(*get<1>(input_batch.pointers[i]))) + bytes * N * (*props.back), *insert_ptr, bytes * num_elements, cudaMemcpyDeviceToDevice, stream);
        }
    }
    is_host_updated = false;
    return cudaGetLastError();
}

static auto casting_coord3d = [](coord3d in){return device_coord3d{static_cast<device_real_t>(in[0]), static_cast<device_real_t>(in[1]), static_cast<device_real_t>(in[2])};};
static auto casting_node3   = [](std::vector<node_t> in){return device_node3{static_cast<device_node_t>(in[0]), static_cast<device_node_t>(in[1]), static_cast<device_node_t>(in[2])};};
static auto casting_coord2d = [](coord2d in){return device_coord2d{static_cast<device_real_t>(in.first), static_cast<device_real_t>(in.second)};};

cudaError_t BatchQueue::insert(const Polyhedron& in, const size_t ID, const cudaStream_t stream, const bool insert_2d){
    to_host(stream);
    if (*props.capacity == *props.size)
    {
        *props.capacity *= 2;
        resize(host_batch, *props.capacity, stream);
        resize(device_batch, *props.capacity, stream);
    }
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
    *props.back = *props.back == *props.capacity ? 0 : (*props.back + 1);

    is_device_updated = false;
    return cudaGetLastError();
}
}
