#ifndef FULLERENE_CUDA_DEFINES
#define FULLERENE_CUDA_DEFINES

#define GPU_REAL float       
#define GPU_REAL3 float3
#define GPU_REAL2 float2
#define GPU_NODE2 ushort2
#define GPU_NODE3 ushort3
#define sin sinf		// change to sin for double precision, sinl for long double
#define cos cosf 		// change to cos for double precision, cosl for long double
#define NODE_MAX UINT16_MAX
#define MAKE_NODE3 make_ushort3
#define Block_Size_Pow_2 256
#define SEMINARIO_FORCE_CONSTANTS 1
#define USE_MAX_NORM 0
#define REDUCTION_METHOD 0
#define LINESEARCH_METHOD GSS


#endif