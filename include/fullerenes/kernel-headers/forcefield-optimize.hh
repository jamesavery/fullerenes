#pragma once
#include <fullerenes/sycl-headers/misc-enums.hh>
#include <fullerenes/sycl-headers/workspace.hh>
#include <fullerenes/batch/batch.hh>
#include <fullerenes/batch/batch_state.hh>
#include <fullerenes/dense_graph.hh>

// See dualize_buffer_size — currently returns 0.
template <ForcefieldType FFT, typename T, typename K>
size_t forcefield_optimize_buffer_size(const Device& device, int N, int capacity);

// graph:       cubic-graph batch (dmax==3, N nodes per isomer). Reads adjacency.
// xyz:         capacity*N 3D coordinates, updated in-place.
// state:       per-entry status; honours NOT_CONVERGED, updates CONVERGED_3D/FAILED_3D.
//              state.iteration[i] is incremented by batch_iters each call.
template <ForcefieldType FFT, typename T, typename K>
SyclEvent forcefield_optimize(SyclQueue& Q,
                              batch::BatchView<Spanify::RSRAdjacencyView<K>> graph,
                              std::span<std::array<T,3>>                          xyz,
                              batch::BatchStateView                          state,
                              size_t batch_iters, size_t max_iters,
                              Workspace                                      ws = {});
