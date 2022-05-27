#include "fullerenes/gpu/cu_array.hh"
#include "cuda_runtime.h"
#include "fullerenes/gpu/cuda_definitions.h"
#include <ostream>
#include <iostream>

template <typename T>
CuArray<T>::CuArray(const size_t size) : size_(size), capacity_(size){
    cudaMallocManaged(&data, size*sizeof(T));
    //cudaMemSet(data,0, size*sizeof(T));
}

template <typename T>
void CuArray<T>::resize(const size_t size){
    T* new_ptr;
    cudaMallocManaged(&new_ptr, size*sizeof(T));
    cudaMemcpy(new_ptr, data, size_*sizeof(T), cudaMemcpyDeviceToDevice );
    
    cudaFree(data);
    data = new_ptr;
}

template <typename T>
size_t CuArray<T>::size(){
    return size_;
}

template <typename T>
T& CuArray<T>::operator[] (const size_t i){
    return data[i];
}

template <typename T>
std::ostream& operator<<(std::ostream& os, const CuArray<T>& input){
    os << "[" << input.data[0];
    for (size_t i = 1; i < input.size_; ++i){
        os << "," << input.data[i];
    }
    os << "]\n";
    return os;
}

template <typename T>
CuArray<T>::~CuArray(){
    cudaFree(data);
}

//Primitive way of handling the fact that templated code in this translation unit wont be generated unless explicitly instantiated somewhere.
int declare_cu_arrays(){
    CuArray<float> f1(1);           f1[0] = {};   std::cout << f1;
    CuArray<double> f2(1);          f2[0] = {};   std::cout << f2;
    //CuArray<device_coord3d> f3(1);  f3[0] = {};   std::cout << f3;
    //CuArray<device_node3> f4(1);    f4[0] = {};   std::cout << f4;
    //CuArray<device_coord2d> f5(1);  f5[0] = {};   std::cout << f5;
    //CuArray<device_node_t> f6(1);   f6[0] = {};   std::cout << f6;
    CuArray<int> f7(1);             f7[0] = {};   std::cout << f7;
    CuArray<size_t> f8(1);          f8[0] = {};   std::cout << f8;
    CuArray<uint8_t> f9(1);         f9[0] = {};   std::cout << f9;
    CuArray<uint16_t> f10(1);       f10[0] = {};  std::cout << f10;
    CuArray<uint32_t> f11(1);       f11[0] = {};  std::cout << f11;
    CuArray<uint64_t> f12(1);       f12[0] = {};  std::cout << f12;
    CuArray<char> f13(1);           f13[0] = {};  std::cout << f13;
    CuArray<bool> f14(1);           f14[0] = {};  std::cout << f14;
    return 1;
}
