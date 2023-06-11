
#include <cooperative_groups.h>
#include <cooperative_groups/reduce.h>
#include <cooperative_groups/scan.h>
#include "fullerenes/gpu/benchmark_functions.hh"
#include "fullerenes/progress_bar.hh"
#include "chrono"
#include <vector>
#include <iostream>
#include <iomanip>
#include <filesystem>
#include <fstream>
#include <numeric>
#include <random>
#include "fullerenes/polyhedron.hh"
#include "fullerenes/gpu/misc_cuda.cuh"
namespace cuda_benchmark {
    #include "fullerenes/gpu/cuda_definitions.h"
    #include "coord3d.cuh"
    #include "reductions.cu"
    __global__ void test_scan_V0_(const int n_times){
        extern __shared__  device_node_t smem[];
        auto val = 1;
        for (int i = 0; i < n_times; i++){
            exclusive_scan(smem, val, blockDim.x);
        }
    }
    __global__ void test_scan_V1_(const int n_times){
        unsigned int power2 = blockDim.x; // compute the next highest power of 2 of 32-bit power2
        power2--;
        power2 |= power2 >> 1;
        power2 |= power2 >> 2;
        power2 |= power2 >> 4;
        power2 |= power2 >> 8;
        power2 |= power2 >> 16;
        power2++;
        extern __shared__  device_node_t smem[];
        auto val = 1;
        for (int i = 0; i < n_times; i++){
            prescan<device_node_t>(smem, val, power2);
        }
    }
    __global__ void test_scan_V2_(const int n_times){
        extern __shared__  device_node_t smem[];
        auto val = 1;
        for (int i = 0; i < n_times; i++){
            ex_scan<device_node_t>(smem, val, blockDim.x);
        }
    }

    __global__ void test_global_scan_V0_(device_node_t* g_out, const int n_times, bool* pass){
        extern __shared__  device_node_t smem[];
        auto val = 0;
        for (int j = 0; j < gridDim.x; j++)
        {
            if(blockIdx.x <= j) val = 1;    
            grid_ex_scan<device_node_t>(g_out, smem, val, gridDim.x + 1);
            if(threadIdx.x + blockIdx.x == 0){
                for (int k = 0; k < gridDim.x + 1; k++){
                    *pass = *pass && g_out[k] == k <= j + 1 ? k - 1 : j;
                }
            }
        }
    }
    __global__ void test_reduction_V0_(const int n_times){
        extern __shared__  device_real_t shmem[];
        for (int i = 0; i < n_times; i++){
            reduction_V0<device_real_t>(shmem, (device_real_t)threadIdx.x);
        }
    }

    __global__ void test_reduction_V1_(const int n_times){
        extern __shared__  device_real_t shmem[];
        for (int i = 0; i < n_times; i++){
            reduction_V1<device_real_t>(shmem, (device_real_t)threadIdx.x);
        }
    }

    __global__ void test_reduction_V2_(const int n_times){
        extern __shared__  device_real_t shmem[];
        for (int i = 0; i < n_times; i++){
            reduction_V2<device_real_t>(shmem, (device_real_t)threadIdx.x);
        }
    }

    __global__ void test_reduction_V3_(const int n_times){
        extern __shared__  device_real_t shmem[];
        for (int i = 0; i < n_times; i++){
            reduction_V3<device_real_t>(shmem, (device_real_t)threadIdx.x);
        }
    }

    __global__ void test_reduction_V4_(const int n_times){
        auto num_warps = (blockDim.x >> 5) + 1;
        int pow2 = 0;
        if (num_warps >= 2) pow2 = 1;
        if (num_warps >= 4) pow2 = 2;
        if (num_warps >= 8) pow2 = 3;
        if (num_warps >= 16) pow2 = 4;
        if (num_warps >= 32) pow2 = 5;

        extern __shared__  device_real_t shmem[];
        for (int i = 0; i < n_times; i++){
            reduction_V4<device_real_t>(shmem, (device_real_t)threadIdx.x, num_warps, pow2);
        }
    }

    __global__ void test_reduction_V5_(const int n_times){
        extern __shared__  device_real_t shmem[];
        __shared__ device_real_t result;
        for (int i = 0; i < n_times; i++){
            reduction_V5<device_real_t>(shmem, (device_real_t)threadIdx.x,&result);
        }
    }

    __global__ void test_reduction_V6_(const int n_times){
        extern __shared__  device_real_t shmem[];
        __shared__ device_real_t result;
        for (int i = 0; i < n_times; i++){
            reduction_V6<device_real_t>(shmem, (device_real_t)threadIdx.x,&result);
        }
    }

    __global__ void test_reduction_V7_(const int n_times){
        auto num_warps = (blockDim.x >> 5) + 1;
        int pow2 = 0;
        if (num_warps >= 2) pow2 = 1;
        if (num_warps >= 4) pow2 = 2;
        if (num_warps >= 8) pow2 = 3;
        if (num_warps >= 16) pow2 = 4;
        if (num_warps >= 32) pow2 = 5;

        extern __shared__  device_real_t shmem[];
        __shared__ device_real_t result;
        for (int i = 0; i < n_times; i++){
            reduction_V7<device_real_t>(shmem, (device_real_t)threadIdx.x, &result, num_warps, pow2);
        }
    }

    __global__ void nothing_kernel_(float* A, float* B, float* C){
        auto tid    = blockDim.x * blockIdx.x + threadIdx.x;
        for(int i = 0; i < 5000; i++){
            C[tid] = (B[tid] + 3.5f*A[tid])/5.6f;
        }
    }

    bool test_global_scan(const int n_elements, const int n_times){
        int n_blocks;
        cudaDeviceProp prop;
        cudaGetDeviceProperties(&prop,0);

        size_t smem = sizeof(device_node_t)*n_elements*2;
        cudaOccupancyMaxActiveBlocksPerMultiprocessor(&n_blocks, test_global_scan_V0_, n_elements, smem);
        n_blocks *= prop.multiProcessorCount;
        device_node_t* g_data;
        bool* pass;
        cudaMalloc(&g_data, sizeof(device_node_t)*n_blocks*2);
        cudaMalloc(&pass, sizeof(bool));

        void* kargs[]{(void*)&g_data, (void*)&n_times, (void*)&pass};
        cudaLaunchCooperativeKernel((void*)test_global_scan_V0_, n_blocks, n_elements, kargs, smem);
        bool h_bool = false;
        device_node_t* h_data;
        cudaMallocHost(&h_data, sizeof(device_node_t)*n_blocks*2);
        cudaMemcpy(&h_bool, pass, sizeof(bool), cudaMemcpyDeviceToHost);
        cudaMemcpy(h_data, g_data, sizeof(device_node_t)*n_blocks*2,cudaMemcpyDeviceToHost);
        
        std::cout << std::vector<device_node_t>(h_data, h_data + n_blocks*2);
        std::cout << n_blocks << endl;
        return h_bool;
    }

    std::chrono::microseconds benchmark_scan(const int n_elements, const int n_times, const int scan_version){

        cudaEvent_t start, stop;
        float single_kernel_time = 0.0;
        cudaEventCreate(&start); cudaEventCreate(&stop);
        cudaEventRecord(start, 0);
        int n_blocks;
        cudaDeviceProp prop;
        cudaGetDeviceProperties(&prop,0);


        size_t smem = sizeof(device_node_t)*n_elements*2;
        cudaOccupancyMaxActiveBlocksPerMultiprocessor(&n_blocks, test_scan_V1_, n_elements, smem);
        n_blocks *= prop.multiProcessorCount;
        void* kargs[]{(void*)&n_times};
        switch (scan_version)
        {
        case 0:
            cudaLaunchCooperativeKernel((void*)test_scan_V0_, n_blocks, n_elements, kargs, smem);  
            break;
        case 1:
            cudaLaunchCooperativeKernel((void*)test_scan_V1_, n_blocks, n_elements, kargs, smem);  
            break;
        case 2:
            cudaLaunchCooperativeKernel((void*)test_scan_V2_, n_blocks, n_elements, kargs, smem);  
        default:
            break;
        }

        cudaEventRecord(stop);
        cudaDeviceSynchronize();
        cudaEventElapsedTime(&single_kernel_time, start, stop);
        

        return 1us* ((int) (single_kernel_time * 1000.0));
    }

    size_t n_blocks(const int n_elements){
        int n_blocks;   
        size_t smem = sizeof(device_real_t)*(n_elements*2 + 32);
        cudaOccupancyMaxActiveBlocksPerMultiprocessor(&n_blocks, test_reduction_V0_, n_elements, smem);
        cudaDeviceProp prop;
        cudaGetDeviceProperties(&prop,0);
        n_blocks *= prop.multiProcessorCount;
        return n_blocks;
    }

    std::chrono::nanoseconds benchmark_reduction(const int n_elements, const int n_times, const int scan_version){

        cudaEvent_t start, stop;
        float single_kernel_time = 0.0;
        cudaEventCreate(&start); cudaEventCreate(&stop);
        cudaEventRecord(start, 0);
       

        size_t smem = sizeof(device_real_t)*(n_elements*4 + 32);
        int M = n_blocks(n_elements);
        void* kargs[]{(void*)&n_times};
        switch (scan_version)
        {
        case 0:
            cudaLaunchCooperativeKernel((void*)test_reduction_V0_, M, n_elements, kargs, smem);  
            break;
        case 1:
            cudaLaunchCooperativeKernel((void*)test_reduction_V1_, M, n_elements, kargs, smem);  
            break;
        case 2:
            cudaLaunchCooperativeKernel((void*)test_reduction_V2_, M, n_elements, kargs, smem);
            break;  
        case 3:
            cudaLaunchCooperativeKernel((void*)test_reduction_V3_, M, n_elements, kargs, smem);  
            break;
        case 4:
            cudaLaunchCooperativeKernel((void*)test_reduction_V4_, M, n_elements, kargs, smem);  
            break;
        case 5:
            cudaLaunchCooperativeKernel((void*)test_reduction_V5_, M, n_elements, kargs, smem);  
            break;
        case 6: 
            cudaLaunchCooperativeKernel((void*)test_reduction_V6_, M, n_elements, kargs, smem);  
            break;
        case 7: 
            cudaLaunchCooperativeKernel((void*)test_reduction_V7_, M, n_elements, kargs, smem);  
            break;
        default:
            break;
        }

        cudaEventRecord(stop);
        cudaDeviceSynchronize();
        cudaEventElapsedTime(&single_kernel_time, start, stop);
        

        return 1ns* ((int) (single_kernel_time * 1e6));
    }

    void warmup_kernel(std::chrono::seconds warmup_time){
        std::cout << "Warming Up... \n\n";
        auto progress_bar = ProgressBar('#',30);
        auto start = std::chrono::high_resolution_clock::now();
        int nBlocks = n_blocks(128);
        
        auto Ndevices = 0;
        cudaGetDeviceCount(&Ndevices);
        cudaStream_t stream[Ndevices];
        float* A[Ndevices]; float* B[Ndevices]; float* C[Ndevices];
        for (size_t i = 0; i < Ndevices; i++)
        {
            cudaSetDevice(i);
            cudaStreamCreate(&stream[i]);
            cudaMalloc(&A[i],nBlocks*128*sizeof(float));
            cudaMalloc(&B[i],nBlocks*128*sizeof(float));
            cudaMalloc(&C[i],nBlocks*128*sizeof(float));
        }
        while (chrono::high_resolution_clock::now() - start < warmup_time){
            for(int i = 0; i < Ndevices; i++){
                void* kargs[]{(void*)&A[i], (void*)&B[i], (void*)&C[i]};
                cudaSetDevice(i);
                cudaLaunchCooperativeKernel((void*)nothing_kernel_,nBlocks,128,kargs,0, stream[i]);
            }
            for(int i = 0; i < Ndevices; i++){
                cudaSetDevice(i);
                cudaStreamSynchronize(stream[i]);
            }
            progress_bar.update_progress((float) ((chrono::high_resolution_clock::now() - start)/1ms) / (float)(warmup_time/1ms));
        }
        for (size_t i = 0; i < Ndevices; i++)
        {
            cudaFree(A[i]); cudaFree(B[i]); cudaFree(C[i]);
        }
        
    }

    /* @brief This function reads a file containing a list of isomers and returns a random one
    * @param path Path to the file containing the isomers
    * @param G Graph to be filled with the isomer
    * @return id of the isomer within the file
    */
    int random_isomer(const std::string &path, Graph& G){ 
        static int id_counter = 0;
        static std::string m_path = path;
        static ifstream isomer_sample(path,std::ios::binary); 
        static auto available_samples = std::filesystem::file_size(path) / (sizeof(device_node_t) * G.N * 6);
        static auto m_fsize = available_samples * sizeof(device_node_t) * G.N * 6;
        static int Nf = G.N;
        static std::vector<device_node_t> input_buffer(m_fsize/sizeof(device_node_t));
        static bool first_time = true;
        static std::vector<int> random_IDs(available_samples);

        if (first_time || G.N != Nf || path != m_path){
            first_time = false;
            m_path = path;
            isomer_sample.close();
            isomer_sample.open(path,std::ios::binary);
            available_samples = std::filesystem::file_size(path) / (sizeof(device_node_t) * G.N * 6);
            Nf = G.N;
            m_fsize = available_samples * sizeof(device_node_t) * G.N * 6;
            input_buffer.resize(m_fsize/sizeof(device_node_t));
            random_IDs.resize(available_samples);

            isomer_sample.read(reinterpret_cast<char*>(input_buffer.data()), m_fsize);
            std::iota(random_IDs.begin(), random_IDs.end(), 0);
            std::shuffle(random_IDs.begin(), random_IDs.end(), std::mt19937{42});
            id_counter = 0;
        }


        

        for (size_t i = 0; i < Nf; i++){
            G.neighbours.at(i).clear();
            for (size_t j = 0; j < 6; j++) {
                auto u = input_buffer.at(random_IDs.at(id_counter)*Nf*6 + i*6 +j);
                if(u != UINT16_MAX) G.neighbours.at(i).push_back(u);
            }
        }

        id_counter = (id_counter + 1) % available_samples;
        return random_IDs.at(id_counter);
    }




}