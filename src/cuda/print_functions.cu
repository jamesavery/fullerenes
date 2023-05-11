#ifndef CUDA_PRINT_FUNCTIONS
#define CUDA_PRINT_FUNCTIONS

#include "fullerenes/gpu/cuda_definitions.h"
#include "cuda_runtime_api.h"

/* __device__ void print(const device_coord3d& ab){
    printf("[%.6f,%.6f,%.6f]",ab.x,ab.y,ab.z);
}
 */
__device__ void print(device_real_t a){
    printf("%.6f", a);
}

__device__ void print(bool b){
    printf("%d",int(b));
}

__device__ void print(int a){
    printf("%d",a);
}

__device__ void print(size_t a){
    printf("%u",(unsigned int)a);
}

__device__ void print(const char* a, int thread_id = 0){
    if (threadIdx.x != thread_id) return;
    printf(a);
}

__device__ void print(const device_coord3d& a, int thread_id = 0){
    if (threadIdx.x != thread_id) return;
    printf("[%.6f,%.6f,%.6f]\n",a[0],a[1],a[2]);
}
/* __device__ void print(const device_node2& a){
    printf("[%d,%d]",a.x,a.y);
}

__device__ void print(const device_node3& a){
    printf("[%d,%d,%d]",a.x,a.y,a.z);
}

__device__ void print(const device_node6& a){
    printf("[%d,%d,%d,%d,%d,%d]",a.b, a.c, a.d, a.e, a.f, a.g);
}
__device__ void print(const device_coord2d& a){
    printf("[%.6f,%.6f]",a.x,a.y);
} */
template <typename T>
__device__ void print_single(T data){
    if (threadIdx.x + blockIdx.x == 0) {
        print(data);
    }
}

template <typename T>
__device__ void sequential_print(T data, size_t fullerene_id){
    if (blockIdx.x == fullerene_id)
    {
    if (threadIdx.x == 0) printf("[");
    cg::sync(cg::this_thread_block());
    for (size_t i = 0; i < blockDim.x; i++)
    {
        if (threadIdx.x == i)
        {   
            if (i != blockDim.x-1)
            {
                print(data); printf(",");
            } else{
                print(data);
            }
        }
        cg::sync(cg::this_thread_block());
    }
    if (threadIdx.x == 0) printf("]\n");
    cg::sync(cg::this_thread_block());
    }
}

template <typename T>
__device__ void grid_print(T data){

    if (threadIdx.x + blockIdx.x == 0) printf("[");
    cg::sync(cg::this_grid());
    for (size_t i = 0; i < gridDim.x; i++)
    {   
            if(threadIdx.x == 0){
            if (blockIdx.x == i)
            {
            if (i != gridDim.x-1)
            {
                print(data); printf(",");
            } else{
                print(data);
            }}

        }
        cg::sync(cg::this_grid());
    }
    if(threadIdx.x + blockIdx.x == 0) printf("]\n");
    cg::sync(cg::this_grid());
}

template <typename T>
__device__ void grid_print(T data, bool mask){

    if (threadIdx.x + blockIdx.x == 0) printf("[");
    cg::sync(cg::this_grid());
    for (size_t i = 0; i < gridDim.x; i++)
    {   
            if(threadIdx.x == 0){
            if (blockIdx.x == i)
            {
            if (i != gridDim.x-1 && mask)
            {
                print(data); printf(",");
            } else if(mask){
                print(data);
            }}

        }
        cg::sync(cg::this_grid());
    }
    if(threadIdx.x + blockIdx.x == 0) printf("]\n");
    cg::sync(cg::this_grid());
}

#endif
