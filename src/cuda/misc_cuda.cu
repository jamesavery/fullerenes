#include "fullerenes/gpu/misc_cuda.cuh"
#include <cooperative_groups.h>

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

