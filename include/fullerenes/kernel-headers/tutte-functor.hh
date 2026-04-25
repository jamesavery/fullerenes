#pragma once
#include <fullerenes/kernel-headers/base-functor.hh>
#include <fullerenes/batch/batch.hh>
#include <fullerenes/batch/batch_state.hh>
#include <fullerenes/dense_graph.hh>

template <typename T, typename K>
struct TutteFunctor : public KernelFunctor<TutteFunctor<T,K>> {
    // View-based batch overload (Phase 7).
    // graph: cubic-graph batch (dmax==3, N nodes per isomer). Reads adjacency.
    // layout: capacity*N 2D coords, updated in-place with Tutte layout.
    // state: per-entry status; honours FULLERENEGRAPH_PREPARED, sets CONVERGED_2D.
    SyclEvent compute(SyclQueue& Q,
                      batch::BatchView<Spanify::RSRAdjacencyView<K>> graph,
                      std::span<std::array<T,2>>                          layout,
                      batch::BatchStateView                          state);

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