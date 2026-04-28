#pragma once
#include <fullerenes/sycl-headers/misc-enums.hh>
#include <fullerenes/sycl-headers/workspace.hh>
#include <fullerenes/batch/batch.hh>
#include <fullerenes/batch/batch_state.hh>
#include <fullerenes/dense_graph.hh>

// See dualize_buffer_size — currently returns 0.
template <ForcefieldType FFT, typename T, typename K>
size_t compute_hessian_buffer_size(const Device& device, int N, int capacity);

// graph:       cubic-graph batch (dmax==3, N nodes per isomer). Reads adjacency.
// xyz:         capacity*N 3D coordinates (read-only).
// state:       per-entry status; honours FULLERENEGRAPH_PREPARED.
// out_hessian: flattened hessian blocks, must be >= 90*N*capacity.
// out_cols:    column index arrays,     must be >= 90*N*capacity.
template <ForcefieldType FFT, typename T, typename K>
SyclEvent compute_hessian(SyclQueue& Q,
                          batch::BatchView<Spanify::RSRAdjacencyView<K>> graph,
                          std::span<std::array<T,3>>                          xyz,
                          batch::BatchStateView                          state,
                          std::span<T> out_hessian, std::span<K> out_cols,
                          Workspace                                      ws = {});