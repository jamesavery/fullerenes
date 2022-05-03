#include "iostream"
#include "fullerenes/gpu/cuda_execution.hh"
#include "cuda_runtime_api.h"

int LaunchCtx::get_device_id() const {
    return m_device_id;
}

bool LaunchCtx::is_finished() const{
    return cudaStreamQuery(stream) == cudaSuccess;
}

bool LaunchCtx::is_default_stream() const{
    return m_unique_stream_idx == -1;
}

void LaunchCtx::wait() const {
    cudaStreamSynchronize(stream);
}

void LaunchCtx::wait_all(){
    for (auto& it: m_all_streams) cudaStreamSynchronize(**it.second);
}   

int LaunchCtx::get_device_count(){
    int count;
    cudaGetDeviceCount(&count);
    return count;
}

int LaunchCtx::get_stream_count(){
    return m_all_streams.size();
}

LaunchCtx::LaunchCtx(){
    cudaGetDeviceCount(&m_device_count);
    if (m_device_count < 1) {
        std::cout << "Error: no CUDA enabled devices found" << std::endl; 
        return;
    }
    stream = cudaStream_t(NULL);
    cudaStream_t* stream_ptr = &stream;
    m_unique_stream_idx = int(-1);
    m_device_id = 0;
    if(default_ctx_created) m_all_streams.insert({m_unique_stream_idx,&stream_ptr});
    default_ctx_created = true;
}

LaunchCtx::LaunchCtx(int device){
    cudaGetDeviceCount(&m_device_count);
    if (m_device_count < device) {std::cout << "Error: requested device was not found" << std::endl; return;}
    int temp_device; cudaGetDevice(&temp_device);
    cudaSetDevice(device);
    cudaStreamCreateWithFlags(&stream, cudaStreamNonBlocking);
    cudaStream_t* stream_ptr = &stream;
    m_unique_stream_idx = m_object_counter++;
    m_all_streams.insert({m_unique_stream_idx,&stream_ptr});
    cudaSetDevice(temp_device);
}

LaunchCtx::~LaunchCtx(){
    //Never destroy the default stream everything will break if you do.
    if (!is_default_stream())
    {
        cudaStreamDestroy(stream);
        m_all_streams.erase(m_unique_stream_idx);
    }
}

