#pragma once
#include <fullerenes/sycl-headers/base-kernel.hh>

template <typename T, typename K>
struct DualizeFunctor : public KernelFunctor<DualizeFunctor<T,K>> {
    template <typename... Data>
    SyclEvent compute(SyclQueue& Q, Fullerene<T,K> fullerene, Data&&... data);
    template <typename... Data>
    SyclEvent compute(SyclQueue& Q, FullereneBatchView<T,K> batch, Data&&... data);

    mutable FunctorArrays<K> cannon_ixs_;
    mutable FunctorArrays<K> rep_count_;
    mutable FunctorArrays<K> scan_array_;
    mutable FunctorArrays<K> triangle_numbers_;
    mutable FunctorArrays<K> arc_list_;
    
    inline constexpr auto to_tuple(size_t N) const {
        size_t Nin = (N/2) + 2;
        size_t Nout = N;
        size_t MaxDegree = 6;
        return  std::make_tuple(
                std::make_pair(std::ref(cannon_ixs_),       Nin * MaxDegree), 
                std::make_pair(std::ref(rep_count_),        Nin), 
                std::make_pair(std::ref(scan_array_),       Nin), 
                std::make_pair(std::ref(triangle_numbers_), Nin*MaxDegree),
                std::make_pair(std::ref(arc_list_),         Nout*2 ));
    }
};