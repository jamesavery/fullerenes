#pragma once
#include <fullerenes/sycl-headers/base-kernel.hh>

template <typename T, typename K>
struct SphericalProjectionFunctor : public KernelFunctor<SphericalProjectionFunctor<T,K>> {
    template <typename... Data>
    SyclEvent compute(SyclQueue& Q, Fullerene<T,K> fullerene, Data&&... data);
    template <typename... Data>
    SyclEvent compute(SyclQueue& Q, FullereneBatchView<T,K> batch, Data&&... data);

    mutable FunctorArrays<K> topological_distances_;
    mutable FunctorArrays<K> reduce_in_;
    mutable FunctorArrays<K> reduce_out_;
    mutable FunctorArrays<K> output_keys_;
    mutable FunctorArrays<std::array<T,2>> sorted_xys_;

    inline constexpr auto to_tuple(size_t N) const {
        return  std::make_tuple(
                std::make_pair(std::ref(topological_distances_), N),
                std::make_pair(std::ref(reduce_in_), N * 4),
                std::make_pair(std::ref(reduce_out_),N * 4),
                std::make_pair(std::ref(output_keys_), N),
                std::make_pair(std::ref(sorted_xys_), N));
    }
};
