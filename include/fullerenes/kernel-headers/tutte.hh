#pragma once
#include <fullerenes/sycl-headers/workspace.hh>
#include <fullerenes/batch/batch.hh>
#include <fullerenes/batch/batch_state.hh>
#include <fullerenes/dense_graph.hh>

// See dualize_buffer_size — currently returns 0.
template <typename T, typename K>
size_t tutte_layout_buffer_size(const Device& device, int N, int capacity);

// graph: cubic-graph batch (dmax==3, N nodes per isomer). Reads adjacency.
// layout: capacity*N 2D coords, updated in-place with Tutte layout.
// state: per-entry status; honours FULLERENEGRAPH_PREPARED, sets CONVERGED_2D.
template <typename T, typename K>
SyclEvent tutte_layout(SyclQueue& Q,
                       batch::BatchView<Spanify::RSRAdjacencyView<K>> graph,
                       std::span<std::array<T,2>>                          layout,
                       batch::BatchStateView                          state,
                       Workspace                                      ws = {});