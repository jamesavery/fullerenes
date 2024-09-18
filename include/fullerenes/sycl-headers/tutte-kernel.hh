#pragma once
#include <fullerenes/sycl-headers/base-kernel.hh>

template <typename T, typename K>
struct TutteFunctor : public KernelFunctor<TutteFunctor<T,K>> {
    SyclEvent compute(SyclQueue& Q, Fullerene<T,K> fullerene, Span<std::array<T,2>> newxys, Span<bool> fixed, Span<T> max_change);
    SyclEvent compute(SyclQueue& Q, FullereneBatchView<T,K> batch);

    mutable FunctorArrays<std::array<T,2>> newxys_;
    mutable FunctorArrays<bool> fixed_;
    mutable FunctorArrays<T> max_change_;

    inline constexpr auto to_tuple(size_t N) const {
        return  std::make_tuple(
                std::make_pair(std::ref(newxys_),    N), 
                std::make_pair(std::ref(fixed_),     N), 
                std::make_pair(std::ref(max_change_),N));
    }

    inline constexpr auto to_tuple_batch(size_t N) const {
        return std::make_tuple();
    }
};