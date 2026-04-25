#pragma once
#include <fullerenes/kernel-headers/base-functor.hh>
#include <fullerenes/batch/batch.hh>
#include <fullerenes/batch/batch_state.hh>
#include <fullerenes/dense_graph.hh>

template <ForcefieldType FFT, typename T, typename K>
struct ForcefieldOptimizeFunctor: public KernelFunctor<ForcefieldOptimizeFunctor<FFT,T,K>> {
    // View-based batch overload (Phase 7).
    // graph:       cubic-graph batch (dmax==3, N nodes per isomer). Reads adjacency.
    // xyz:         capacity*N 3D coordinates, updated in-place.
    // state:       per-entry status; honours NOT_CONVERGED, updates CONVERGED_3D/FAILED_3D.
    //              state.iteration[i] is incremented by batch_iters each call.
    SyclEvent compute(SyclQueue& Q,
                      batch::BatchView<Spanify::RSRAdjacencyView<K>> graph,
                      std::span<std::array<T,3>>                          xyz,
                      batch::BatchStateView                          state,
                      size_t batch_iters, size_t max_iters);


    mutable FunctorArrays<K> indices_;
    mutable FunctorArrays<std::array<T,3>> X1_;
    mutable FunctorArrays<std::array<T,3>> X2_;
    mutable FunctorArrays<std::array<T,3>> g0_;
    mutable FunctorArrays<std::array<T,3>> g1_;
    mutable FunctorArrays<std::array<T,3>> s_;

    inline constexpr auto to_tuple(size_t N, size_t iterations, size_t max_iterations) const {
        return  std::make_tuple(
                std::make_pair(std::ref(indices_), N),
                std::make_pair(std::ref(X1_), N),
                std::make_pair(std::ref(X2_), N),
                std::make_pair(std::ref(g0_), N),
                std::make_pair(std::ref(g1_), N),
                std::make_pair(std::ref(s_), N)
                );
    }

    inline constexpr auto to_tuple_batch(size_t N, size_t iterations, size_t max_iterations) const {
        return std::make_tuple();
    }
};
