#pragma once
#include <fullerenes/kernel-headers/base-functor.hh>
#include <fullerenes/batch/batch.hh>
#include <fullerenes/batch/batch_state.hh>
#include <fullerenes/dense_graph.hh>

template<EigensolveMode mode, typename T, typename K>
struct EigenFunctor : public KernelFunctor<EigenFunctor<mode, T, K>> {
    // View-based batch overload (Phase 7).
    // xyz:         capacity*N 3D coordinates (read-only, for deflation against rigid-body modes).
    SyclEvent compute(SyclQueue& Q,
                      std::span<std::array<T,3>> xyz,
                      int N, int capacity,
                      std::span<T> hessian, std::span<K> cols, size_t n_lanczos,
                      std::span<T> eigenvalues, std::span<T> eigenvectors,
                      std::span<T> off_diagonal, std::span<T> qmat,
                      std::span<T> lanczos, std::span<T> diag, std::span<K> ends_idx);


    mutable FunctorArrays<K> indices_;
    mutable FunctorArrays<T> off_diagonal_;
    mutable FunctorArrays<T> qmat_;
    mutable FunctorArrays<T> lanczos_;
    mutable FunctorArrays<T> diag_;
    mutable FunctorArrays<K> ends_idx_;

    template <typename... Args>
    inline constexpr auto to_tuple(size_t N, std::span<T> hessian, std::span<K> cols, size_t n_lanczos, Args&&... args) const {
        return  std::make_tuple(
                std::make_pair(std::ref(indices_), N),
                std::make_pair(std::ref(off_diagonal_), n_lanczos),
                std::make_pair(std::ref(qmat_), n_lanczos*n_lanczos),
                std::make_pair(std::ref(lanczos_), n_lanczos*N*3),
                std::make_pair(std::ref(diag_), n_lanczos),
                std::make_pair(std::ref(ends_idx_), 2)
                );
    }

    template <typename... Args>
    inline constexpr auto to_tuple_batch(size_t N, std::span<T> hessian, std::span<K> cols, size_t n_lanczos, Args&&... args) const {
        return std::make_tuple(
                std::make_pair(std::ref(off_diagonal_), n_lanczos),
                std::make_pair(std::ref(qmat_), n_lanczos*n_lanczos),
                std::make_pair(std::ref(lanczos_), n_lanczos*N*3),
                std::make_pair(std::ref(diag_), n_lanczos),
                std::make_pair(std::ref(ends_idx_), 2)
        );
    }
};  