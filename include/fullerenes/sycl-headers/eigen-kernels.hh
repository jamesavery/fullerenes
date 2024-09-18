#pragma once
#include <fullerenes/sycl-headers/base-kernel.hh>

template<EigensolveMode mode, typename T, typename K>
struct EigenFunctor : public KernelFunctor<EigenFunctor<mode, T, K>> {
    SyclEvent compute(SyclQueue& Q, FullereneBatchView<T, K> batch, Span<T> hessian, Span<K> cols, size_t n_lanczos, Span<T> eigenvalues, Span<T> eigenvectors, 
                        Span<T> off_diagonal, Span<T> qmat, Span<T> lanczos, Span<T> diag, Span<K> ends_idx);
    SyclEvent compute(SyclQueue& Q, Fullerene<T, K> fullerene, Span<T> hessian, Span<K> cols, size_t n_lanczos, Span<T> eigenvalues, Span<T> eigenvectors,
                        Span<K> indices, Span<T> off_diagonal, Span<T> qmat, Span<T> lanczos, Span<T> diag, Span<K> ends_idx);


    mutable FunctorArrays<K> indices_;
    mutable FunctorArrays<T> off_diagonal_;
    mutable FunctorArrays<T> qmat_;
    mutable FunctorArrays<T> lanczos_;
    mutable FunctorArrays<T> diag_;
    mutable FunctorArrays<K> ends_idx_;

    inline constexpr auto to_tuple(size_t N, Span<T> hessian, Span<K> cols, size_t n_lanczos, Span<T> eigenvalues, Span<T> eigenvectors) const {
        return  std::make_tuple(
                std::make_pair(std::ref(indices_), N),
                std::make_pair(std::ref(off_diagonal_), n_lanczos),
                std::make_pair(std::ref(qmat_), n_lanczos*n_lanczos),
                std::make_pair(std::ref(lanczos_), n_lanczos*N*3),
                std::make_pair(std::ref(diag_), n_lanczos),
                std::make_pair(std::ref(ends_idx_), 2)
                );
    }

    inline constexpr auto to_tuple_batch(size_t N, Span<T> hessian, Span<K> cols, size_t n_lanczos, Span<T> eigenvalues, Span<T> eigenvectors) const {
        return std::make_tuple(
                std::make_pair(std::ref(off_diagonal_), n_lanczos),
                std::make_pair(std::ref(qmat_), n_lanczos*n_lanczos),
                std::make_pair(std::ref(lanczos_), n_lanczos*N*3),
                std::make_pair(std::ref(diag_), n_lanczos),
                std::make_pair(std::ref(ends_idx_), 2)
        );
    }
};  