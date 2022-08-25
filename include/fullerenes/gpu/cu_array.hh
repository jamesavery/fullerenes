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

    ~CuArray();
    
    void resize(const size_t capacity);

    //Unsafe indexing
    __device__ __host__ T& operator[] (const size_t i);

    size_t size();

    //Safe indexing with bounds check.
    T& at(const size_t i);

    T* data;
    size_t size_ = 0;
private:
    int* flag;
    size_t capacity_ = 0;
};

template <typename T>
std::ostream& operator<<(std::ostream& os, const CuArray<T>& input);

