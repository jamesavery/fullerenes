#pragma once
#include "cuda_runtime.h"
#include "cuda.h"
#include "cuda_runtime_api.h"
#include <ostream>

template <typename T>
struct CuArray
{
public:
    CuArray();
    CuArray(const size_t size);
    CuArray(const CuArray& other);
    CuArray(const size_t size, const T& value);
    ~CuArray();
    
    void resize(const size_t capacity);

    //Unsafe indexing
    __device__ __host__ T& operator[] (const size_t i);

    //Copy assignment
    CuArray& operator=(const CuArray& other);

    size_t size();
    void to_device(const int device); //Forces a copy to device
    void to_host(const int device); //Forces a copy to host
    void fill(const T& value);

    T* data;
    size_t size_ = 0;
private:
    int* flag;
    bool initialized_ = false;
    size_t capacity_ = 0;
};

template <typename T>
std::ostream& operator<<(std::ostream& os, const CuArray<T>& input);

