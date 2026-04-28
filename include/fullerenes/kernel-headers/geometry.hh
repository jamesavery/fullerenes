#pragma once
#include <fullerenes/sycl-headers/workspace.hh>
#include <fullerenes/batch/batch.hh>
#include <fullerenes/batch/batch_state.hh>

// See dualize_buffer_size — currently all return 0.
template <typename T, typename K>
size_t eccentricity_buffer_size(const Device& device, int N, int capacity);
template <typename T, typename K>
size_t inertia_buffer_size(const Device& device, int N, int capacity);
template <typename T, typename K>
size_t transform_to_principal_axes_buffer_size(const Device& device, int N, int capacity);
template <typename T, typename K>
size_t surface_area_buffer_size(const Device& device, int N, int capacity);
template <typename T, typename K>
size_t volume_buffer_size(const Device& device, int N, int capacity);

// Free-function entry points.
template <typename T, typename K>
SyclEvent eccentricity(SyclQueue& Q, std::span<std::array<T,3>> xyz, int N, int capacity,
                       batch::BatchStateView state, std::span<T> out_ellipticity,
                       Workspace ws = {});
template <typename T, typename K>
SyclEvent inertia(SyclQueue& Q, std::span<std::array<T,3>> xyz, int N, int capacity,
                  batch::BatchStateView state, std::span<std::array<T,3>> out_inertia,
                  Workspace ws = {});
template <typename T, typename K>
SyclEvent transform_to_principal_axes(SyclQueue& Q, std::span<std::array<T,3>> xyz,
                                      int N, int capacity, batch::BatchStateView state,
                                      Workspace ws = {});
template <typename T, typename K>
SyclEvent surface_area(SyclQueue& Q, std::span<std::array<T,3>> xyz, int N, int capacity,
                       std::span<std::array<K,6>> faces, std::span<uint8_t> deg,
                       batch::BatchStateView state, std::span<T> out_surface_area,
                       Workspace ws = {});
template <typename T, typename K>
SyclEvent volume(SyclQueue& Q, std::span<std::array<T,3>> xyz, int N, int capacity,
                 std::span<std::array<K,6>> faces, std::span<uint8_t> deg,
                 batch::BatchStateView state, std::span<T> out_volume,
                 Workspace ws = {});
