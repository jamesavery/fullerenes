#pragma once
#include <fullerenes/sycl-headers/base-kernel.hh>

template <typename T, typename K>
struct TutteFunctor : public KernelFunctor<TutteFunctor<T,K>> {
    template <typename... Data>
    SyclEvent compute(SyclQueue& Q, Fullerene<T,K> fullerene, Data&&... data);
    template <typename... Data>
    SyclEvent compute(SyclQueue& Q, FullereneBatchView<T,K> batch, Data&&... data);

    mutable FunctorArrays<std::array<T,2>> newxys_;
    mutable FunctorArrays<bool> fixed_;
    mutable FunctorArrays<T> max_change_;

    inline constexpr auto to_tuple(size_t N) const {
        return  std::make_tuple(
                std::make_pair(std::ref(newxys_),    N), 
                std::make_pair(std::ref(fixed_),     N), 
                std::make_pair(std::ref(max_change_),N));
    }
};