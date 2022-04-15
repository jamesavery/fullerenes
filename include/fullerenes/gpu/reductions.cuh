#ifndef CUDA_REDUCTIONS
#define CUDA_REDUCTIONS
#include <cooperative_groups.h>
#include <cooperative_groups/reduce.h>
#include "cuda_definitions.h"
#include "cuda_runtime.h"
namespace cg = cooperative_groups;

__device__ GPU_REAL reduction(GPU_REAL* sdata, const GPU_REAL data);

__device__ GPU_REAL reduction_max(GPU_REAL* sdata, const GPU_REAL data);

__device__ GPU_NODE reduction_max(GPU_NODE* sdata, const GPU_NODE data);

__device__ GPU_REAL global_reduction(GPU_REAL* sdata, GPU_REAL* gdata, GPU_REAL data, bool mask = true);

#endif
