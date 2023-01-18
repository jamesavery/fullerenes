#pragma once
#include "cuda_runtime.h"
#include <cuda_runtime_api.h>
#include <cuda.h>
#include "device_launch_parameters.h"
#include <stdio.h>
#include <iostream>
#include <fstream>
#include <chrono>
#include "fullerenes/gpu/isomerspace_tutte.hh"

#include "coord2d.cuh"

__global__
void kernel_tutte_layout(IsomerBatch G, const size_t iterations){
    DEVICE_TYPEDEFS
    extern __shared__  real_t sharedmem[];
    clear_cache(sharedmem, Block_Size_Pow_2);

    if (G.statuses[blockIdx.x] == IsomerStatus::NOT_CONVERGED)
    {
    size_t offset = blockIdx.x * blockDim.x;
    CubicGraph FG(&G.cubic_neighbours[offset*3]); 
    real_t* base_pointer        = sharedmem + Block_Size_Pow_2;
    coord2d* xys        = reinterpret_cast<coord2d*>(base_pointer);
    coord2d* newxys     = reinterpret_cast<coord2d*>(base_pointer) + blockDim.x;


    node3 ns            = (reinterpret_cast<node3*>(G.cubic_neighbours) + offset)[threadIdx.x];
    xys[threadIdx.x]    = {real_t(0.0), real_t(0.0)};
    device_node_t outer_face[6];
    device_node_t outer_face_vertex   = 0;
    uint8_t Nface = FG.get_face_oriented(0,FG.cubic_neighbours[0], outer_face);    
    reinterpret_cast<bool*>(sharedmem)[threadIdx.x] =  false; BLOCK_SYNC;
    if(threadIdx.x < Nface){
      outer_face_vertex = outer_face[threadIdx.x];
      reinterpret_cast<bool*>(sharedmem)[outer_face_vertex] =  true; 
    }
    BLOCK_SYNC;
    bool fixed = reinterpret_cast<bool*>(sharedmem)[threadIdx.x];

    if(threadIdx.x < Nface) xys[outer_face_vertex] = {sin(threadIdx.x*2*real_t(M_PI)/double(Nface)),cos(threadIdx.x*2*real_t(M_PI)/double(Nface))};
    BLOCK_SYNC
    bool converged          = false;
    real_t max_change       = real_t(0.0);
    if(fixed) newxys[threadIdx.x] = xys[threadIdx.x];

    for (size_t i = 0; i < iterations && !converged; i++)
    {   
        max_change = real_t(0.0);
        BLOCK_SYNC
        coord2d neighbour_sum   = {real_t(0.0),real_t(0.0)};    
        for (uint8_t j = 0; j < 3; j++) neighbour_sum += xys[d_get(ns,j)];

        if(!fixed) newxys[threadIdx.x] = xys[threadIdx.x]*real_t(0.15) + (neighbour_sum/3)*real_t(0.85);
        real_t neighbour_dist = 0.0f;

        for (uint8_t j = 0; j < 3; j++) neighbour_dist += norm(xys[threadIdx.x] - xys[d_get(ns,j)])/3;
        
        BLOCK_SYNC
        real_t relative_change = 0.0f;
        if (neighbour_dist > 0.0f && !fixed){ 
            relative_change = norm(xys[threadIdx.x] - newxys[threadIdx.x])/neighbour_dist;
        }

        real_t iteration_max = reduction_max(sharedmem, relative_change);
        if (iteration_max > max_change) max_change = iteration_max;
        converged = max_change <= 5e-4;

        xys[threadIdx.x] = newxys[threadIdx.x];
    }
    BLOCK_SYNC
    (reinterpret_cast<coord2d*>(G.xys) + offset )[threadIdx.x] = xys[threadIdx.x];
    G.statuses[blockIdx.x] = IsomerStatus::CONVERGED;
    }
}

void IsomerspaceTutte::tutte_layout(){
    printLastCudaError("Memcpy Failed! \n");
    auto start = std::chrono::system_clock::now();
    for (size_t i = 0; i < device_count; i++) {cudaSetDevice(i); d_batch[i] <<= h_batch[i];}
    for (size_t i = 0; i < device_count; i++)
    {
        cudaSetDevice(i);
        void* kernelArgs[] = {(void*)&d_batch[i],(void*)&TUTTE_MAX_ITERATION};
        safeCudaKernelCall((void*)kernel_tutte_layout, dim3(device_capacities[i], 1, 1), dim3(N, 1, 1), kernelArgs, shared_memory_bytes);
    }
    for (size_t i = 0; i < device_count; i++) {cudaSetDevice(i); h_batch[i] <<= d_batch[i];}
        
    cudaDeviceSynchronize();
    
    auto end = std::chrono::system_clock::now();
    printLastCudaError("Tutte kernel launch failed: ");
    
}
void IsomerspaceTutte::check_batch(){
    batch_size = 0;
    for (size_t i = 0; i < device_count; i++)
    {
        cudaSetDevice(i);
        int num_of_not_converged_isomers = 0;
        h_batch[i] <<= d_batch[i];
        
        IsomerStatus* statuses = h_batch[i].statuses;
        for (int j = 0; j < device_capacities[i]; j++)
        {   
            num_of_not_converged_isomers += (int)(statuses[j] == IsomerStatus::NOT_CONVERGED);
            if (statuses[j] != IsomerStatus::NOT_CONVERGED)
            {   
                index_queue[i].push(j);
            }
            
        }
        batch_sizes[i] = num_of_not_converged_isomers;
        batch_size += num_of_not_converged_isomers;
    }
}


IsomerspaceTutte::IsomerspaceTutte(const size_t N) : IsomerspaceKernel::IsomerspaceKernel(N, (void*)kernel_tutte_layout){
    this->shared_memory_bytes = sizeof(device_coord2d)*N*2 + sizeof(device_real_t)*Block_Size_Pow_2;
    std::cout << "Tutte Capacity: " << this->batch_capacity << "\n";
    printLastCudaError("Tutte kernel class instansiation failed!");
}

IsomerspaceTutte::~IsomerspaceTutte(){}
