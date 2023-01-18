#ifndef FULLERENE_CUDA_DEFINES
#define FULLERENE_CUDA_DEFINES

#include "inttypes.h"
#include "cuda_runtime.h"


#define BLOCK_SYNC __syncthreads();
#define GRID_SYNC cg::sync(cg::this_grid());
#define INLINE __device__ __forceinline__
#define sin sin		// change to sin for double precision, sinl for long double
#define cos cos 		// change to cos for double precision, cosl for long double
#define NODE_MAX UINT16_MAX
#define MAKE_NODE3 make_ushort3
#define Block_Size_Pow_2 256
#define SEMINARIO_FORCE_CONSTANTS 0
#define USE_MAX_NORM 0
#define REDUCTION_METHOD 3
#define LINESEARCH_METHOD GSS
#define FORCEFIELD_VERSION FLATNESS_ENABLED
#define USE_CONSTANT_INDICES 0
typedef float device_real_t;
typedef uint16_t device_node_t;
typedef float3 device_coord3d;
typedef ushort3 device_node3;
typedef ushort2 device_node2;
typedef float2 device_coord2d;
typedef struct device_node6
{
    device_node_t b;
    device_node_t c;
    device_node_t d;
    device_node_t e;
    device_node_t f;
    device_node_t g;
} device_node6;
#define DEVICE_TYPEDEFS typedef device_coord3d coord3d; typedef device_coord2d coord2d; typedef device_real_t real_t; typedef device_node3 node3; typedef device_node_t node_t; typedef device_node6 node6;

#endif
