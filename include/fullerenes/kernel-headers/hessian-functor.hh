#pragma once
#include <fullerenes/kernel-headers/base-functor.hh>
#include <fullerenes/batch/batch.hh>
#include <fullerenes/batch/batch_state.hh>
#include <fullerenes/dense_graph.hh>

template<ForcefieldType FFT, typename T, typename K>
struct HessianFunctor : public KernelFunctor<HessianFunctor<FFT, T, K>> {
    SyclEvent compute(SyclQueue& Q, FullereneBatchView<T, K> batch, Span<T> out_hessian, Span<K> out_cols);
    SyclEvent compute(SyclQueue& Q, Fullerene<T, K> batch, Span<T> out_hessian, Span<K> out_cols, Span<K> indices);

    // View-based batch overload (Phase 7).
    // graph:       cubic-graph batch (dmax==3, N nodes per isomer). Reads adjacency.
    // xyz:         capacity*N 3D coordinates (read-only).
    // state:       per-entry status; honours FULLERENEGRAPH_PREPARED.
    // out_hessian: flattened hessian blocks, must be >= 90*N*capacity.
    // out_cols:    column index arrays,     must be >= 90*N*capacity.
    SyclEvent compute(SyclQueue& Q,
                      batch::BatchView<Spanify::RSRAdjacencyView<K>> graph,
                      Span<std::array<T,3>>                          xyz,
                      batch::BatchStateView                          state,
                      Span<T> out_hessian, Span<K> out_cols);

    mutable FunctorArrays<K> indices_;

    inline constexpr auto to_tuple(size_t N, Span<T> out_hessian, Span<K> out_cols) const {
        return  std::make_tuple(
                std::make_pair(std::ref(indices_), N)
                );
    }

    inline constexpr auto to_tuple_batch(size_t N, Span<T> out_hessian, Span<K> out_cols) const {
        return std::make_tuple();
    }
};