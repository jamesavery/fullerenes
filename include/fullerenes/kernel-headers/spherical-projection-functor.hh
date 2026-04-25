#pragma once
#include <fullerenes/kernel-headers/base-functor.hh>
#include <fullerenes/batch/batch.hh>
#include <fullerenes/batch/batch_state.hh>
#include <fullerenes/dense_graph.hh>

template <typename T, typename K>
struct SphericalProjectionFunctor : public KernelFunctor<SphericalProjectionFunctor<T,K>> {
    
    // View-based batch overload (Phase 7).
    // graph:      cubic-graph batch (dmax==3, N nodes per isomer). Reads adjacency.
    // layout_2d:  capacity*N 2D input coordinates (Tutte layout).
    // xyz_3d:     capacity*N 3D output coordinates (spherical projection result).
    // state:      per-entry status; honours FULLERENEGRAPH_PREPARED, sets NOT_CONVERGED.
    SyclEvent compute(SyclQueue& Q,
                      batch::BatchView<Spanify::RSRAdjacencyView<K>> graph,
                      std::span<std::array<T,2>>                          layout_2d,
                      std::span<std::array<T,3>>                          xyz_3d,
                      batch::BatchStateView                          state);

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
                std::make_pair(std::ref(sorted_xys_), N * 2));
    }

    inline constexpr auto to_tuple_batch(size_t N) const {
        return std::make_tuple();
    }
};
