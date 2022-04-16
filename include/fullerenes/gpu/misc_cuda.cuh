#ifndef MISCEllANEOUS_CUDA_FUN
#define MISCEllANEOUS_CUDA_FUN

#include "cuda_definitions.h"
#include "cuda_runtime_api.h"
#include "string"

//Sets the memory to 0. This exists because L1 cache is reused for different purposes, thus pointers alias.
//Algorithms that depend on cache being 0 initialized to avoid bound checking must be prepared with a call to this function.
__device__ void clear_cache(device_real_t* sdata, size_t N);

//Swaps anything that implements the copy assignment operator.
//TODO: Rename: name should be swap but clashes with std::swap because of using namespace std.
template <typename T>
__device__ void swap_reals(T& real1, T& real2);

//Atomically updates the data at the location *data, in order of threadIdx.x (floating point math is not associative)
__device__ void ordered_atomic_add(device_real_t* data, const device_real_t element);

//Initializes device memory with the provided value.

//Prints the last cuda error prepended by a custom string.
void printLastCudaError(std::string message = "");

//Checks if the dimensions passed are non-zero.
cudaError_t safeCudaKernelCall(const void* func, const dim3 gridDim, const dim3 blockDim, void** args, const size_t sharedMem, const cudaStream_t stream = NULL);

//Initialization function for device arrays.
template <typename T>
cudaError_t fill_cu_array(T* cu_array, size_t size, T fill_value);

#endif