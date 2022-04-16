#include "fullerenes/gpu/misc_cuda.cuh"
#include <cooperative_groups.h>
#include "iostream"
#include "assert.h"

namespace cg = cooperative_groups;

__device__ void clear_cache(device_real_t* sdata, size_t N){
    BLOCK_SYNC
    for (size_t index = threadIdx.x; index < N; index+=blockDim.x)
    {
        sdata[index] = (device_real_t)0.0;
    }
    BLOCK_SYNC
}

template <typename T>
__device__ void swap_reals(T& a, T& b){
    T temp = a;
    a = b;
    b = temp;
}

__device__ void ordered_atomic_add(device_real_t* data, const device_real_t element){
    for (size_t i = 0; i < blockDim.x; i++)
    {
        BLOCK_SYNC
        if(threadIdx.x == i){
            *data += element;
        }
    }
}

void printLastCudaError(std::string message){
    cudaError_t error = cudaGetLastError();
    if(error != cudaSuccess){
        std::cout << "\n" << message << " :\t";
        std::cout << cudaGetErrorString(error);
        printf("\n");
        assert(false);
    }
}

cudaError_t safeCudaKernelCall(const void* func, const dim3 gridDim, const dim3 blockDim, void** args, const size_t sharedMem, const cudaStream_t stream){
    if (gridDim.x > 0 && gridDim.y == 1 && gridDim.z == 1 && blockDim.x > 0 && blockDim.y == 1 && blockDim.z == 1)
    {
        return cudaLaunchCooperativeKernel(func,gridDim,blockDim,args,sharedMem,stream);
    }
    else
    {
        std::cout << "WARNING: Attempted to launch kernel with 1 or more dimensions <= 0 \n";
        return cudaErrorInvalidValue;
    }
    
}

template <typename T>
__global__ void kernel_fill_array(T* cu_array, size_t size, T fillvalue){
    auto tid = blockIdx.x * blockDim.x + threadIdx.x;
    if(tid < size) cu_array[tid] = fillvalue;

}

template <typename T>
cudaError_t fill_cu_array(T* cu_array, size_t size, T fill_value) {
    size_t Nblocks = size / 64 + 1;
    kernel_fill_array<<<dim3(Nblocks, 1, 1), dim3(64, 1, 1)>>>(cu_array, size, fill_value);
    return cudaDeviceSynchronize();
}
