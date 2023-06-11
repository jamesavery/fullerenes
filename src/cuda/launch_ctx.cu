#include "iostream"
#include "fullerenes/gpu/launch_ctx.hh"
#include "cuda_runtime_api.h"

LaunchCtx& LaunchCtx::operator=(const LaunchCtx& other){
    if (this != &other){
        m_device_id = other.m_device_id;
        cudaSetDevice(m_device_id);
        cudaStreamCreateWithFlags(&stream, cudaStreamNonBlocking);
        cudaEventCreateWithFlags(&m_start, cudaEventBlockingSync);
        cudaEventCreateWithFlags(&m_stop, cudaEventBlockingSync);
        
    }
    return *this;
}

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

void LaunchCtx::clear_allocations(){
    for (int i = 0; i < get_device_count(); ++i){
        cudaSetDevice(i);
        cudaDeviceSynchronize();
        cudaDeviceReset();
    }
}

void LaunchCtx::start_timer(){
    cudaEventRecord(m_start, stream);
}

std::chrono::nanoseconds LaunchCtx::stop_timer(){
    float elapsed_time = 0.0f;
    cudaEventRecord(m_stop, stream);
    cudaEventSynchronize(m_stop);
    cudaEventElapsedTime(&elapsed_time, m_start, m_stop); //elapsed_time is in ms
    return std::chrono::nanoseconds((int) (elapsed_time*1e6));
}

LaunchCtx::LaunchCtx(){
    cudaGetDeviceCount(&m_device_count);
    if (m_device_count < 1) {
        std::cout << "Error: no CUDA enabled devices found" << std::endl; 
        return;
    }
    stream = cudaStream_t(NULL);
    cudaStream_t* stream_ptr = &stream;
    cudaEventCreateWithFlags(&m_start, cudaEventBlockingSync);
    cudaEventCreateWithFlags(&m_stop, cudaEventBlockingSync);
    m_unique_stream_idx = int(-1);
    m_device_id = 0;
    if(!default_ctx_created) m_all_streams.insert({m_unique_stream_idx,&stream_ptr});
    default_ctx_created = true;
}

LaunchCtx::LaunchCtx(int device){
    static int m_device_count = get_device_count();
    static std::vector<bool> first_call(16, true);
    if (m_device_count < device) {std::cout << "Error: requested device was not found" << std::endl; return;}
    int temp_device; cudaGetDevice(&temp_device);
    m_device_id = device;
    cudaSetDevice(device);
    if(first_call[device]) cudaSetDeviceFlags(cudaDeviceScheduleBlockingSync);
    cudaStreamCreateWithFlags(&stream, cudaStreamNonBlocking);
    cudaStream_t* stream_ptr = &stream;
    cudaEventCreateWithFlags(&m_start, cudaEventBlockingSync);
    cudaEventCreateWithFlags(&m_stop, cudaEventBlockingSync);
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
    cudaEventDestroy(m_start); cudaEventDestroy(m_stop);
}

