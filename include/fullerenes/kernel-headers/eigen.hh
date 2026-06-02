#pragma once
#include <fullerenes/sycl-headers/misc-enums.hh>
#include <fullerenes/sycl-headers/workspace.hh>
#include <fullerenes/batch/batch.hh>
#include <fullerenes/batch/batch_state.hh>
#include <fullerenes/dense_graph.hh>

// Bytes the eigen kernel needs from a Workspace for (N, capacity, n_lanczos).
// In a future iteration the eigen entry-point will allocate its scratch
// (off_diagonal, qmat, lanczos, diag, ends_idx, indices) from this
// workspace instead of asking callers to pre-size each Span individually.
// Until then the existing API still takes those Spans externally and this
// returns 0.
template <EigensolveMode mode, typename T, typename K>
size_t eigensolve_buffer_size(const Device& device, int N, int capacity, size_t n_lanczos);

// xyz: capacity*N 3D coordinates (read-only, for deflation against rigid-body modes).
template <EigensolveMode mode, typename T, typename K>
SyclEvent eigensolve(SyclQueue& Q,
                     std::span<std::array<T,3>> xyz,
                     int N, int capacity,
                     std::span<T> hessian, std::span<K> cols, size_t n_lanczos,
                     std::span<T> eigenvalues, std::span<T> eigenvectors,
                     std::span<T> off_diagonal, std::span<T> qmat,
                     std::span<T> lanczos, std::span<T> diag, std::span<K> ends_idx,
                     Workspace ws = {});
