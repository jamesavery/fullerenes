#include "fullerenes/gpu/cu_array.hh"

template <typename T>
CuArray<T>::CuArray(const size_t N, const size_t Nisomers, const BufferType buffer_type) : N(N), Nisomers(Nisomers), buffer_type(buffer_type){
    size_t bytes = sizeof(T) * N * Nisomers;
    buffer_type == DEVICE_BUFFER ? cudaMalloc(&this->data, bytes) : cudaMallocHost(&this->data, bytes);
}

template <typename T>
CuArray<T>::~CuArray(){
    buffer_type == DEVICE_BUFFER ? cudaFree(&this->data) : cudaFreeHost(&this->data);
}