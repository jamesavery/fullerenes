#pragma once
#include <fullerenes/sycl-headers/workspace.hh>
#include <fullerenes/batch/batch.hh>
#include <fullerenes/batch/batch_state.hh>
#include <fullerenes/dense_graph.hh>

// See dualize_buffer_size — currently returns 0.
template <typename T, typename K>
size_t spherical_projection_buffer_size(const Device& device, int N, int capacity);

// graph:      cubic-graph batch (dmax==3, N nodes per isomer). Reads adjacency.
// layout_2d:  capacity*N 2D input coordinates (Tutte layout).
// xyz_3d:     capacity*N 3D output coordinates (spherical projection result).
// state:      per-entry status; honours FULLERENEGRAPH_PREPARED, sets NOT_CONVERGED.
template <typename T, typename K>
SyclEvent spherical_projection(SyclQueue& Q,
                               batch::BatchView<Spanify::RSRAdjacencyView<K>> graph,
                               std::span<std::array<T,2>>                          layout_2d,
                               std::span<std::array<T,3>>                          xyz_3d,
                               batch::BatchStateView                          state,
                               Workspace                                      ws = {});
