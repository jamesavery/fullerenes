#include "cuda_definitions.h"
#include "cuda_runtime_api.h"

//Sets the memory to 0. This exists because L1 cache is reused for different purposes, thus pointers alias.
//Algorithms that depend on cache being 0 initialized to avoid bound checking must be prepared with a call to this function.
__device__ void clear_cache(GPU_REAL* sdata, size_t N);

//Swaps anything that implements the copy assignment operator.
//TODO: Rename: name should be swap but clashes with std::swap because of using namespace std.
template <typename T>
__device__ void swap_reals(T& real1, T& real2);

//Atomically updates the data at the location *data, in order of threadIdx.x (floating point math is not associative)
__device__ void ordered_atomic_add(GPU_REAL* data, const GPU_REAL element);
