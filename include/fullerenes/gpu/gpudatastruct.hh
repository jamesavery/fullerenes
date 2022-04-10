#pragma once
#include <inttypes.h>
#include <stddef.h>
#include <vector>
#include <string>
#include <cuda_runtime.h>

enum BufferType   {HOST_BUFFER, DEVICE_BUFFER};

struct GPUDataStruct{
    bool allocated = false;
    size_t n_isomers        = 0;
    size_t n_atoms          = 0;
    size_t isomer_capacity  = 0;

    BufferType buffer_type;
    GPUDataStruct(){}
    std::vector<std::tuple<std::string,void**,size_t,bool>> pointers;
    static void allocate(GPUDataStruct& G,const size_t N,const size_t batch_size, const BufferType buffer_type);
    static void free(GPUDataStruct& G);
    static void copy(GPUDataStruct& destination, const GPUDataStruct& source, cudaStream_t stream = NULL);
};  