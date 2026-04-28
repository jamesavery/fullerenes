#pragma once
#include <fullerenes/sycl-headers/workspace.hh>
#include <fullerenes/batch/batch.hh>
#include <fullerenes/batch/batch_state.hh>
#include <fullerenes/dense_graph.hh>

// Bytes the dualize kernel needs from a Workspace for a (N, capacity)
// batch. The current view-based path uses sycl::local_accessor for all
// per-call scratch, so this is 0; the API exists so callers can size a
// shared workspace via std::max(...) over every kernel and not need to
// change when a future split-kernel impl introduces persistent USM
// scratch.
template <typename T, typename K>
size_t dualize_buffer_size(const Device& device, int N, int capacity);

//   src: dual-graph batch (Nv=Nf, dmax<=6). Reads neighbours/deg.
//   dst: cubic-graph batch (Nv=N,  dmax==3). Writes neighbours.
//   state: per-entry status; honours DUAL_INITIALIZED and sets
//          FULLERENEGRAPH_PREPARED on success.
//   faces_cubic: capacity*Nf output of triangle indices per dual vertex.
//   faces_dual : capacity*N  output of dual-vertex triples per triangle.
template <typename T, typename K>
SyclEvent dualize(SyclQueue& Q,
                  batch::BatchView<Spanify::RSRAdjacencyView<K>> src,
                  batch::BatchView<Spanify::RSRAdjacencyView<K>> dst,
                  batch::BatchStateView                          state,
                  std::span<std::array<K,6>>                          faces_cubic,
                  std::span<std::array<K,3>>                          faces_dual,
                  Workspace                                      ws = {});