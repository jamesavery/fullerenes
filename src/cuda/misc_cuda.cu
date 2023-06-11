#include "fullerenes/gpu/misc_cuda.cuh"
#include <cooperative_groups.h>
#include "iostream"
#include "assert.h"
#include "fullerenes/gpu/launch_ctx.hh"

namespace cg = cooperative_groups;

template <typename T>
__device__ __forceinline__ T d_max(const T& a, const T& b){
    return a > b ? a : b; 
}
/**
 * @brief This function clears the cache in preparation for a new operation
 *
 * @param sdata: pointer to the cache to be cleared
 * @param N: number of elements to clear
 */
__device__ void clear_cache(device_real_t* sdata, size_t N){
    BLOCK_SYNC
    for (size_t index = threadIdx.x; index < N; index+=blockDim.x)
    {
        sdata[index] = (device_real_t)0.0;
    }
    BLOCK_SYNC
}

/** 
* @brief This code swaps the values of two variables of type T.
* The names of the variables are a and b.
* @param a: the first variable
* @param b: the second variable
*/
template <typename T>
__device__ void swap_reals(T& a, T& b){
    T temp = a;
    a = b;
    b = temp;
}

/**
* @brief This function takes a data pointer and an element, then atomically adds the
* element to the data in order of thread index. This is useful for summing up floating point
* values in a block, as floating point addition is not associative.
* @param data: pointer to the data to be added to
* @param element: the element to be added    
*/
__device__ void ordered_atomic_add(device_real_t* data, const device_real_t element){
    for (size_t i = 0; i < blockDim.x; i++)
    {
        BLOCK_SYNC
        if(threadIdx.x == i){
            *data += element;
        }
    }
}


/**
 * This function print the last CUDA error along with a message to the console.
 * @param message The error message to print
 */
void printLastCudaError(std::string message){
    cudaError_t error = cudaGetLastError();
    if(error != cudaSuccess){
        std::cout << "\n" << message << " :\t";
        std::cout << cudaGetErrorString(error);
        printf("\n");
        assert(false);
    }
}



/** 
* @brief This function launches a kernel provided that the grid and block dimensions are valid.
* It also checks for errors in the kernel launch
* @param func is a pointer to the kernel function
* @param gridDim is the size of the grid (number of blocks) in each dimension
* @param blockDim is the size of the blocks in each dimension
* @param args is a pointer to the arguments for the kernel
* @param sharedMem is the amount of shared memory to allocate for the kernel
* @param stream is the stream that the kernel should be launched on
* @return cudaError_t is the error code returned by the kernel launch
*/
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

/**
 * @brief Fill an array with a given value
 * @param cu_array pointer to memory
 * @param size number of elements
 * @param fillvalue value to fill array with
 */
template <typename T>
__global__ void kernel_fill_array(T* cu_array, size_t size, T fillvalue){
    auto tid = blockIdx.x * blockDim.x + threadIdx.x;
    if(tid < size) cu_array[tid] = fillvalue;
}

/** @brief Fills an array with a given value.
 *  @param cu_array device array to fill
 *  @param size number of elements to fill
 *  @param fill_value value to fill with
 *  @return cudaError_t error code
 */
template <typename T>
cudaError_t fill_cu_array(T* cu_array, size_t size, T fill_value) {
    size_t Nblocks = size / 64 + 1;
    kernel_fill_array<<<dim3(Nblocks, 1, 1), dim3(64, 1, 1)>>>(cu_array, size, fill_value);
    return cudaDeviceSynchronize();
}

