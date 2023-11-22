#pragma once

#include <fullerenes/sycl-isomer-batch.hh>

template <typename T = float, typename K = uint16_t>
void dualise(sycl::queue&Q, IsomerBatch<T,K>& batch, const LaunchPolicy policy);

template <typename T = float, typename K = uint16_t>
void tutte_layout(sycl::queue&Q, IsomerBatch<T,K>& batch, const LaunchPolicy policy);

template <typename T = float, typename K = uint16_t>
void spherical_projection(sycl::queue&Q, IsomerBatch<T,K>& batch, const LaunchPolicy policy);

template <ForcefieldType FFT = PEDERSEN, typename T = float, typename K = uint16_t>
void forcefield_optimise(sycl::queue&Q, IsomerBatch<T,K>& batch, const int iterations, const int max_iterations, const LaunchPolicy policy);

template <ForcefieldType FFT = PEDERSEN, typename T = float, typename K = uint16_t>
void compute_hessians(sycl::queue&Q, IsomerBatch<T,K>& batch, sycl::buffer<T,1>& hessians, sycl::buffer<K,1>& cols, const LaunchPolicy policy);
