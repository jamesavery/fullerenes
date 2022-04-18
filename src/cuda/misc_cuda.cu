#include "fullerenes/gpu/misc_cuda.cuh"
#include <cooperative_groups.h>
#include "iostream"
#include "assert.h"

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

void printLastCudaError(std::string message){
    cudaError_t error = cudaGetLastError();
    if(error != cudaSuccess){
        std::cout << "\n" << message << " :\t";
        std::cout << cudaGetErrorString(error);
        printf("\n");
        assert(false);
    }
}

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

template <typename T>
__global__ void kernel_fill_array(T* cu_array, size_t size, T fillvalue){
    auto tid = blockIdx.x * blockDim.x + threadIdx.x;
    if(tid < size) cu_array[tid] = fillvalue;
}

template <typename T>
cudaError_t fill_cu_array(T* cu_array, size_t size, T fill_value) {
    size_t Nblocks = size / 64 + 1;
    kernel_fill_array<<<dim3(Nblocks, 1, 1), dim3(64, 1, 1)>>>(cu_array, size, fill_value);
    return cudaDeviceSynchronize();
}

//Container for launch dimensions on all devices, can be stored statically in a function.
struct LaunchDims
{   
    //Figures out which device is being used and gets the appropriate grid dimensions.
    dim3 get_grid(){
        int device_id; cudaGetDevice(&device_id);
        return grid_dims[device_id];
    }

    //Figures out which device is being used and gets the appropriate block dimensions.
    dim3 get_block(){
        int device_id; cudaGetDevice(&device_id);
        return block_dims[device_id];
    }

    //Allocates memory according to device count and computes all launch dimensions for all devices.
    LaunchDims(const void* kernel, const int __block_size, const size_t __smem = 0) {
        cudaGetDeviceCount(&device_count);
        block_dims.resize(device_count);
        grid_dims.resize(device_count);
        SMs.resize(device_count);
        //Fetch the multiprocessor count for each device.
        for (size_t i = 0; i < device_count; i++) {
            cudaDeviceProp prop; cudaGetDeviceProperties(&prop, i);
            SMs[i] = prop.multiProcessorCount;
        }
        
        update_dims(kernel, __block_size, __smem);
    }

    //Only recomputes dimensions if a new block size or new amount of dynamic shared memory is specified.
    cudaError_t update_dims(const void* kernel, const int __block_size, const size_t __smem = 0){
        if(__block_size == block_size && __smem == smem) return cudaSuccess;
        block_size = __block_size; smem = __smem;
        int current_device; cudaGetDevice(&current_device);
        int temp_int; //CUDA api is pretty particular about integer types
        for (size_t i = 0; i < device_count; i++)
        {   
            cudaSetDevice(i);
            cudaOccupancyMaxActiveBlocksPerMultiprocessor(&temp_int,kernel,block_size,smem);
            grid_dims[i].x = temp_int * SMs[i];
            block_dims[i].x = block_size;
        }
        cudaSetDevice(current_device);
        cudaDeviceSynchronize();
        return cudaGetLastError();
    }

    private:
        int block_size{};
        size_t smem{};
        int device_count{};
        std::vector<int> SMs;
        std::vector<dim3> grid_dims;
        std::vector<dim3> block_dims;
};

__global__
void __reset_convergence_status(IsomerBatch B){
    if (threadIdx.x == 0) B.statuses[blockIdx.x] = B.statuses[blockIdx.x] != EMPTY ? NOT_CONVERGED : EMPTY;
    if (threadIdx.x == 0) B.iterations[blockIdx.x] = 0;
}

cudaError_t reset_convergence_statuses(IsomerBatch& B, const cudaStream_t stream){
    static LaunchDims dims((void*)__reset_convergence_status, B.n_atoms);
    void* kargs[]{(void*)&B};
    return cudaLaunchCooperativeKernel((void*)__reset_convergence_status, dims.get_grid(), dims.get_block(), kargs, 0, stream);
}

