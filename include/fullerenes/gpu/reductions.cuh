#ifndef CUDA_REDUCTIONS
#define CUDA_REDUCTIONS

namespace cg = cooperative_groups;

template <typename T>
__device__ T reduction(T* sdata, const T data);

template <typename T>
__device__ T reduction_max(T* sdata, const T data);

__device__ device_node_t reduction_max(device_node_t* sdata, const device_node_t data);

template <typename T>
__device__ T global_reduction(T* sdata, T* gdata, T data, bool mask = true);

#endif
