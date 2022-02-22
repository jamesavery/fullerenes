#include "gpudatastruct.hh"

template <typename T>
struct CuArray
{
public:
    size_t N = 0;
    size_t Nisomers = 0;
    BufferType buffer_type;

    T* data;
    
    CuArray(){}
    CuArray(T* pointer, size_t stride){
        this->data = pointer;
        this->stride = stride;
    }
    CuArray(const size_t N, const size_t Nisomers, const BufferType buffer_type);
    ~CuArray();
    T* operator[] (const size_t i) {
        return data + i*N;
    }
};
