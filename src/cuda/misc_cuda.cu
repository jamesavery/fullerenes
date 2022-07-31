#include "fullerenes/gpu/misc_cuda.cuh"
#include <cooperative_groups.h>
#include "iostream"
#include "assert.h"
#include "fullerenes/gpu/cuda_execution.hh"

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
        return m_grid_dims[device_id];
    }

    //Figures out which device is being used and gets the appropriate block dimensions.
    dim3 get_block(){
        int device_id; cudaGetDevice(&device_id);
        return m_block_dims[device_id];
    }

    //Allocates memory according to device count and computes all launch dimensions for all devices.
    LaunchDims(const void* kernel, const int block_size_, const size_t smem_ = 0, const int num_tasks = -1) {
        cudaGetDeviceCount(&m_device_count);
        m_block_dims.resize(m_device_count);
        m_grid_dims.resize(m_device_count);
        m_properties.resize(m_device_count);
        //Fetch the multiprocessor count for each device.
        for (size_t i = 0; i < m_device_count; i++) {
            cudaGetDeviceProperties(&m_properties[i], i);
        }
        
        update_dims(kernel, block_size_, smem_, num_tasks);
    }

    //Only recomputes dimensions if a new block size or new amount of dynamic shared memory is specified.
    cudaError_t update_dims(const void* kernel, const int block_size_, const size_t smem_ = 0, const int num_tasks = -1){
        if(block_size_ == m_block_size && smem_ == m_smem) return cudaSuccess;
        m_block_size = block_size_; m_smem = smem_;
        int current_device; cudaGetDevice(&current_device);
        int temp_int; //CUDA api is pretty particular about integer types
        for (size_t i = 0; i < m_device_count; i++)
        {   
            cudaSetDevice(i);
            cudaOccupancyMaxActiveBlocksPerMultiprocessor(&temp_int,kernel,m_block_size,m_smem);
            int SMs = m_properties[i].multiProcessorCount;
            if( num_tasks == -1) {m_grid_dims[i].x = temp_int * SMs;}
            else{ 
                int num_blocks = 0;
                int best_blocks = temp_int*SMs;
                int best_remainder = temp_int*SMs;
                for (num_blocks = temp_int*SMs; num_blocks >= SMs ; num_blocks -= 1)
                {    
                    int remainder = (num_tasks / num_blocks + (num_tasks % num_blocks != 0)) * num_blocks - num_tasks;
                    if (remainder < best_remainder) {
                        best_blocks = num_blocks; 
                        best_remainder = remainder;}
                    if (best_remainder == 0) break;
                }
                m_grid_dims[i].x = best_blocks;
            }
            m_block_dims[i].x = m_block_size;
        }
        cudaSetDevice(current_device);
        cudaDeviceSynchronize();
        return cudaGetLastError();
    }

private:

    int m_block_size{};
    int m_smem{};
    int m_device_count{};
    int m_num_tasks{};
    std::vector<cudaDeviceProp> m_properties;
    std::vector<dim3> m_grid_dims;
    std::vector<dim3> m_block_dims;
};

__global__
void reset_convergence_status_(IsomerBatch B){
    for (int isomer_idx = blockIdx.x; isomer_idx < B.isomer_capacity; isomer_idx+=gridDim.x)
    {
        if (threadIdx.x == 0) B.statuses[isomer_idx] = B.statuses[isomer_idx] != EMPTY ? NOT_CONVERGED : EMPTY;
        if (threadIdx.x == 0) B.iterations[isomer_idx] = 0;
    }
    
}

cudaError_t reset_convergence_statuses(IsomerBatch& B, const LaunchCtx& ctx, const LaunchPolicy policy){
    cudaSetDevice(ctx.get_device_id());
    static LaunchDims dims((void*)reset_convergence_status_, B.n_atoms);
    void* kargs[]{(void*)&B};
    return cudaLaunchCooperativeKernel((void*)reset_convergence_status_, dims.get_grid(), dims.get_block(), kargs, 0, ctx.stream);
}

