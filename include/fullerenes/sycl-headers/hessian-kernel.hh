#pragma once
#include <fullerenes/sycl-headers/base-kernel.hh>

template<ForcefieldType FFT, typename T, typename K>
struct HessianFunctor : public KernelFunctor<HessianFunctor<FFT, T, K>> {
    template <typename... Data>
    SyclEvent compute(SyclQueue& Q, FullereneBatchView<T, K> batch, Span<T> out_hessian, Span<K> out_cols, Data&&... data);
    template <typename... Data>
    SyclEvent compute(SyclQueue& Q, Fullerene<T, K> batch, Span<T> out_hessian, Span<K> out_cols, Data&&... data);

    mutable FunctorArrays<K> indices_;

    inline constexpr auto to_tuple(size_t N, Span<T> out_hessian, Span<K> out_cols) const {
        return  std::make_tuple(
                std::make_pair(std::ref(indices_), N)
                );
    }
};