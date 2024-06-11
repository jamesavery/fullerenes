#pragma once
#include <fullerenes/sycl-isomer-batch.hh>


template <typename T = float, typename K = uint16_t>
void dualize(sycl::queue&Q, IsomerBatch<T,K>& batch, const LaunchPolicy policy);

template <typename T = float, typename K = uint16_t>
void dualize_V1(sycl::queue&Q, IsomerBatch<T,K>& batch, const LaunchPolicy policy);

template <int MaxDegIn, int MaxDegOut, typename K = uint16_t>
void dualize_general_for_visualization(sycl::queue&Q, sycl::buffer<K>& G_in, sycl::buffer<K>& Deg_in, sycl::buffer<K>& G_out, sycl::buffer<K>& Deg_out, int Nin, int Nout, LaunchPolicy policy, bool output_intermediate = false);

template <int MaxDegIn, int MaxDegOut, typename K = uint16_t>
void dualize_general(sycl::queue&Q, sycl::buffer<K>& G_in, sycl::buffer<K>& Deg_in, sycl::buffer<K>& G_out, sycl::buffer<K>& Deg_out, int Nin, int Nout, LaunchPolicy policy);

template <typename T = float, typename K = uint16_t>
void tutte_layout(sycl::queue&Q, IsomerBatch<T,K>& batch, const LaunchPolicy policy);

template <typename T = float, typename K = uint16_t>
void spherical_projection(sycl::queue&Q, IsomerBatch<T,K>& batch, const LaunchPolicy policy);

template <ForcefieldType FFT = PEDERSEN, typename T = float, typename K = uint16_t>
void forcefield_optimize(sycl::queue&Q, IsomerBatch<T,K>& batch, const int iterations, const int max_iterations, const LaunchPolicy policy);

template <ForcefieldType FFT = PEDERSEN, typename T = float, typename K = uint16_t>
void compute_hessians(sycl::queue&Q, IsomerBatch<T,K>& batch, sycl::buffer<T,1>& hessians, sycl::buffer<K,1>& cols, const LaunchPolicy policy);

template <typename T = float, typename K = uint16_t>
void nop_kernel(sycl::queue&Q, IsomerBatch<T,K>& batch, const LaunchPolicy policy);

template<class T> sycl::buffer<T,1>& default_object() {
    static sycl::buffer<T,1> x(0);
    return x;
} //For providing default lvalue reference arguments, super toxic hack 
template <EigensolveMode mode, typename T, typename K>
void eigensolve(sycl::queue& ctx, IsomerBatch<T,K> B, sycl::buffer<T,1>& hessians, sycl::buffer<K,1>& cols, sycl::buffer<T,1>& eigenvalues, const LaunchPolicy policy = LaunchPolicy::SYNC,  size_t Nlanczos = 50, sycl::buffer<T,1>& eigenvectors = default_object<T>());