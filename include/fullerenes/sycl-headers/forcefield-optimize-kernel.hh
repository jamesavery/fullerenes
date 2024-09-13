#pragma once
#include <fullerenes/sycl-headers/base-kernel.hh>

template <ForcefieldType FFT, typename T, typename K>
struct ForcefieldOptimizeFunctor: public KernelFunctor<ForcefieldOptimizeFunctor<FFT,T,K>> {
    //void operator()(SyclQueue& Q, Fullerene<T,K> fullerene, const int iterations, const int max_iterations, LaunchPolicy policy);
    //void operator()(SyclQueue& Q, FullereneBatchView<T,K> batch, const int iterations, const int max_iterations, LaunchPolicy policy);
    template <typename... Data>
    SyclEvent compute(SyclQueue& Q, Fullerene<T,K> fullerene, size_t iterations, size_t max_iterations, Data&&... data);

    template <typename... Data>
    SyclEvent compute(SyclQueue& Q, FullereneBatchView<T,K> batch, size_t iterations, size_t max_iterations, Data&&... data);


    mutable FunctorArrays<K> indices_;
    mutable FunctorArrays<std::array<T,3>> X1_;
    mutable FunctorArrays<std::array<T,3>> X2_;
    mutable FunctorArrays<std::array<T,3>> g0_;
    mutable FunctorArrays<std::array<T,3>> g1_;
    mutable FunctorArrays<std::array<T,3>> s_;

    inline constexpr auto to_tuple(size_t N, size_t iterations, size_t max_iterations) const {
        return  std::make_tuple(
                std::make_pair(std::ref(indices_), N),
                std::make_pair(std::ref(X1_), N),
                std::make_pair(std::ref(X2_), N),
                std::make_pair(std::ref(g0_), N),
                std::make_pair(std::ref(g1_), N),
                std::make_pair(std::ref(s_), N)
                );
    }
};
