#pragma once
#include <fullerenes/kernel-headers/base-functor.hh>
#include <fullerenes/batch/batch.hh>
#include <fullerenes/batch/batch_state.hh>
#include <fullerenes/dense_graph.hh>

template<ForcefieldType FFT, typename T, typename K>
struct HessianFunctor : public KernelFunctor<HessianFunctor<FFT, T, K>> {
    // View-based batch overload (Phase 7).
    // graph:       cubic-graph batch (dmax==3, N nodes per isomer). Reads adjacency.
    // xyz:         capacity*N 3D coordinates (read-only).
    // state:       per-entry status; honours FULLERENEGRAPH_PREPARED.
    // out_hessian: flattened hessian blocks, must be >= 90*N*capacity.
    // out_cols:    column index arrays,     must be >= 90*N*capacity.
    SyclEvent compute(SyclQueue& Q,
                      batch::BatchView<Spanify::RSRAdjacencyView<K>> graph,
                      std::span<std::array<T,3>>                          xyz,
                      batch::BatchStateView                          state,
                      std::span<T> out_hessian, std::span<K> out_cols);

    mutable FunctorArrays<K> indices_;

    inline constexpr auto to_tuple(size_t N, std::span<T> out_hessian, std::span<K> out_cols) const {
        return  std::make_tuple(
                std::make_pair(std::ref(indices_), N)
                );
    }

    inline constexpr auto to_tuple_batch(size_t N, std::span<T> out_hessian, std::span<K> out_cols) const {
        return std::make_tuple();
    }
};