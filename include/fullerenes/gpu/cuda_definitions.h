#ifndef FULLERENE_CUDA_DEFINES
#define FULLERENE_CUDA_DEFINES

#include "inttypes.h"
#include "cuda_runtime.h"
#include "array"
#include "type_traits"


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
#define REDUCTION_METHOD 3 //0 is classical shared memory reduction 3 is warp reduction
#define LINESEARCH_METHOD GSS
#define FORCEFIELD_VERSION FLATNESS_ENABLED
#define USE_CONSTANT_INDICES 0
//Half-precision and particularly float16 currently do not seem to provide accurate enough results.
#ifndef FLOAT_TYPE
#define FLOAT_TYPE 2
#endif
#if FLOAT_TYPE == 0
    #include "cuda_fp16.h"
    typedef __half device_real_t;
    typedef std::array<__half,2> device_coord2d;
    typedef std::array<__half,3> device_coord3d;
    #define SQRT hsqrt 
    #define RSQRT hrsqrt
    #define ABS __habs
    #define COS hcos
    #define ACOS(arg) acosf((float)arg)
    #define SIN hsin
    #define ISNAN __hisnan
#elif FLOAT_TYPE == 1
    #include "cuda_bf16.h"
    typedef __nv_bfloat16 device_real_t;
    typedef std::array<__nv_bfloat16,2> device_coord2d;
    typedef std::array<__nv_bfloat16,3> device_coord3d;
    #define SQRT hsqrt 
    #define RSQRT hrsqrt
    #define ABS __habs
    #define COS hcos
    #define ACOS(arg) acosf((float)arg)
    #define SIN hsin
    #define ISNAN __hisnan
#elif FLOAT_TYPE == 2
    typedef float device_real_t;
    typedef std::array<float,2> device_coord2d;
    typedef std::array<float,3> device_coord3d;
    #define SQRT sqrtf
    #define RSQRT rsqrtf
    #define ABS abs
    #define COS cosf
    #define ACOS(arg) acos((double)arg)
    #define SIN sinf
    #define ISNAN isnan
#elif FLOAT_TYPE == 3
    typedef double device_real_t;
    typedef std::array<double,2> device_coord2d;
    typedef std::array<double,3> device_coord3d;
    #define SQRT sqrt
    #define RSQRT rsqrt
    #define ABS abs
    #define COS cos
    #define ACOS acos
    #define SIN sin
    #define ISNAN isnan
#endif
//typedef double device_hpreal_t;	/* High precision real type, for when we need it regardless of speed. */
//#define HPREAL_DOUBLE 1
typedef device_real_t device_hpreal_t;	/* High precision real type, for when we need it regardless of speed. */

typedef uint16_t device_node_t;
typedef std::array<uint16_t,3> device_node3;
typedef std::array<uint16_t,2> device_node2;
typedef std::array<uint16_t,6> device_node6;
#define DEVICE_TYPEDEFS typedef device_coord3d coord3d; typedef device_coord2d coord2d; typedef device_real_t real_t; typedef device_node3 node3; typedef device_node_t node_t; typedef device_node6 node6; typedef device_hpreal_t hpreal_t;
#define FLOAT_TYPEDEFS(T) static_assert(std::is_floating_point<T>::value, "T must be float"); typedef std::array<T,3> coord3d; typedef std::array<T,2> coord2d; typedef T real_t;
#define INT_TYPEDEFS(K) static_assert(std::is_integral<K>::value, "K must be integral type"); typedef std::array<K,3> node3; typedef std::array<K,2> node2; typedef K node_t; typedef std::array<K,6> node6;
#define TEMPLATE_TYPEDEFS(T,K) FLOAT_TYPEDEFS(T) INT_TYPEDEFS(K)
#define SMEM(T) extern __shared__ unsigned char my_smem[]; T* smem = reinterpret_cast<T*>(my_smem);


#endif
