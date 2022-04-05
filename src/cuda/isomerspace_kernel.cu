#pragma once
#include "fullerenes/gpu/isomerspace_kernel.hh"
#include "io.cu"
#include "cuda_runtime.h"
#include "cuda.h"
#include "cuda_runtime_api.h"

size_t IsomerspaceKernel::get_batch_capacity(){
    cudaDeviceProp properties;
    int device_count = 0;
    int total_capacity = 0;
    int fullerenes_per_SM;
    cudaGetDeviceCount(&device_count);
    for (size_t i = 0; i < device_count; i++)
    {
        cudaGetDeviceProperties(&properties,i);
        /** Compiling with --maxrregcount=64   is necessary to easily (singular blocks / fullerene) parallelize fullerenes of size 20-1024 !**/
        /** Needs 3 storage arrays for coordinates and 1 for reductions **/
        /** Calculates maximum number of resident fullerenes on a single Streaming Multiprocessor, multiply with multi processor count to d_get total batch size**/
        cudaOccupancyMaxActiveBlocksPerMultiprocessor(&fullerenes_per_SM, kernel_pointer, N, shared_memory_bytes); // How many threads per block
        this->device_capacities[i] = properties.multiProcessorCount*fullerenes_per_SM;
        total_capacity += properties.multiProcessorCount*fullerenes_per_SM;
    }
    return (size_t)total_capacity;
}

IsomerspaceKernel::IsomerspaceKernel(const size_t N, void* kernel){
    cudaGetDeviceCount(&this->device_count);
    kernel_pointer          = kernel;
    global_reduction_arrays = std::vector<device_real_t*>(device_count);
    batch_sizes             = std::vector<int>(device_count);
    device_capacities       = std::vector<int>(device_count);
    this->N                 = N;
    batch_capacity          = get_batch_capacity();    
    index_queue             = std::vector<std::queue<int>>(device_count);
    
    //Create 2 streams that have no implicit synchronization with the default stream.
    cudaStreamCreateWithFlags(&main_stream, cudaStreamNonBlocking);
    cudaStreamCreateWithFlags(&copy_to_host_stream, cudaStreamNonBlocking);

    d_batch = std::vector<IsomerBatch>(device_count);
    h_batch = std::vector<IsomerBatch>(device_count);

    for (size_t i = 0; i < device_count; i++)
    {
        cudaSetDevice(i);
        GPUDataStruct::allocate(d_batch[i]         , N, device_capacities[i], DEVICE_BUFFER);
        GPUDataStruct::allocate(h_batch[i]         , N, device_capacities[i], HOST_BUFFER);

        cudaMalloc(&global_reduction_arrays[i], sizeof(device_real_t)*N*device_capacities[i]);
        batch_sizes[i] = 0;
        for (size_t j = 0; j < device_capacities[i]; j++) index_queue[i].push(j);
        for (size_t j = 0; j < device_capacities[i]; j++) h_batch[i].statuses[j] = EMPTY;

        d_batch[i] <<= h_batch[i];
    }
}

IsomerspaceKernel::~IsomerspaceKernel(){
    for (size_t i = 0; i < device_count; i++)
    {
        cudaSetDevice(i);
        GPUDataStruct::free(h_batch[i]);
        cudaDeviceReset();
    }
}