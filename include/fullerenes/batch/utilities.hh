#pragma once

// batch/utilities.hh — convenience helpers that live on top of the core
// Batch / BatchView / BatchQueue layer. These replace the old
// FullereneBatch / FullereneQueue QoL utilities (per-row assignment,
// whole-batch copy, Graph/Polyhedron conversions) without reintroducing
// those legacy types.
//
// All helpers are free functions in namespace `batch`.

#include "fullerenes/batch/batch.hh"
#include "fullerenes/batch/batch_queue.hh"
#include "fullerenes/dense_graph.hh"
#include "fullerenes/graph.hh"

#include <array>
#include <cassert>
#include <cstdint>
#include <span>
#include <type_traits>
#include <vector>

namespace batch {

// ---------------------------------------------------------------------------
// Per-slot assignment: copy a host graph into a single `RSRAdjacencyView<Kd>`.
// Source can be any RSRAdjacencyView<Ks> (including Graph, which inherits
// from RSRAdjacencyView<node_t>). Index types may differ (narrowing cast).
// Empty slots in the destination row are filled with sentinel (-1) so that
// `deg[u] < dmax` rows have a well-defined tail.
// ---------------------------------------------------------------------------
template <class Kd, class Ks>
inline void assign(Spanify::RSRAdjacencyView<Kd> dst,
                   const Spanify::RSRAdjacencyView<Ks>& src)
{
    assert(int(dst.N) == int(src.N));
    const int N    = int(src.N);
    const int sdm  = src.dmax;
    const int ddm  = dst.dmax;
    for (int u = 0; u < N; ++u) {
        const int du = int(src.deg[u]);
        assert(du <= ddm);
        dst.deg[u]   = uint8_t(du);
        for (int k = 0; k < ddm; ++k) {
            if (k < du && k < sdm)
                dst.neighbours[u * ddm + k] = Kd(src.neighbours[u * sdm + k]);
            else
                dst.neighbours[u * ddm + k] = Kd(-1);
        }
    }
}

// Polyhedron variant: also copies per-vertex xyz from an external span.
template <class T, class Kd, class Ks>
inline void assign(Spanify::RSRPolyhedronView<T, Kd> dst,
                   const Spanify::RSRAdjacencyView<Ks>& src_graph,
                   std::span<const std::array<T, 3>> src_xyz)
{
    assign(static_cast<Spanify::RSRAdjacencyView<Kd>>(dst), src_graph);
    assert(src_xyz.size() >= std::size_t(src_graph.N));
    for (int u = 0; u < int(src_graph.N); ++u) dst.points[u] = src_xyz[u];
}

// ---------------------------------------------------------------------------
// push_back overloads: append one host graph (+ optional xyz) to a Batch or
// BatchQueue, filling any metadata. Returns the slot index (for Batch:
// logical/physical, identical pre-pop; for BatchQueue: physical slot).
// ---------------------------------------------------------------------------
template <class Kd, class Ks>
inline int push_back(Batch<Spanify::RSRAdjacencyView<Kd>>& dst,
                     const Spanify::RSRAdjacencyView<Ks>& src)
{
    assert(dst.size() < dst.capacity());
    int slot = dst.size();
    dst.resize(slot + 1);
    assign(dst.view_capacity()[std::size_t(slot)], src);
    return slot;
}

template <class Kd, class Ks>
inline int push_back(BatchQueue<Spanify::RSRAdjacencyView<Kd>>& q,
                     const Spanify::RSRAdjacencyView<Ks>& src,
                     uint64_t id = 0, StatusFlag flag = StatusFlag{}, int32_t iter = 0)
{
    // Allocate a slot using push_back(default-constructed V) would fail
    // batchability contract; instead we reserve a slot via a temporary V
    // whose spans we populate from `src`. Simpler: materialise a one-slot
    // Batch<V>, assign, then push by value.
    using V = Spanify::RSRAdjacencyView<Kd>;
    Batch<V> tmp(int(src.N), 1, q.dmax());
    push_back(tmp, src);
    return q.push_back(tmp.view()[0], id, flag, iter);
}

template <class T, class Kd, class Ks>
inline int push_back(Batch<Spanify::RSRPolyhedronView<T, Kd>>& dst,
                     const Spanify::RSRAdjacencyView<Ks>& src,
                     std::span<const std::array<T, 3>> xyz)
{
    assert(dst.size() < dst.capacity());
    int slot = dst.size();
    dst.resize(slot + 1);
    assign(dst.view_capacity()[std::size_t(slot)], src, xyz);
    return slot;
}

template <class T, class Kd, class Ks>
inline int push_back(BatchQueue<Spanify::RSRPolyhedronView<T, Kd>>& q,
                     const Spanify::RSRAdjacencyView<Ks>& src,
                     std::span<const std::array<T, 3>> xyz,
                     uint64_t id = 0, StatusFlag flag = StatusFlag{}, int32_t iter = 0)
{
    using V = Spanify::RSRPolyhedronView<T, Kd>;
    Batch<V> tmp(int(src.N), 1, q.dmax());
    push_back(tmp, src, xyz);
    return q.push_back(tmp.view()[0], id, flag, iter);
}

// ---------------------------------------------------------------------------
// Whole-batch copy: element-wise copy of each batchable field for every
// slot in `src` into the corresponding slot in `dst`. Both views must have
// identical N, dmax and size.
// ---------------------------------------------------------------------------
template <class V>
    requires batchable_view<V>
inline void copy(BatchView<V> dst, BatchView<V> src)
{
    assert(dst.N()    == src.N());
    assert(dst.dmax() == src.dmax());
    assert(dst.size() == src.size());
    auto ds = dst.spans();
    auto ss = src.spans();
    detail::for_each_field<V>([&](auto Ic) {
        constexpr std::size_t k = Ic;
        auto& d = std::get<k>(ds);
        const auto& s = std::get<k>(ss);
        const std::size_t n = std::min<std::size_t>(d.size(), s.size());
        for (std::size_t i = 0; i < n; ++i) d[i] = s[i];
    });
}

// ---------------------------------------------------------------------------
// Broadcast a single-isomer view `src` into every slot of `dst`.
// `src`'s field spans must have size N * size_factor[k] (one isomer).
// ---------------------------------------------------------------------------
template <class V>
    requires batchable_view<V>
inline void fill(BatchView<V> dst, const V& src)
{
    assert(dst.N()    == int(src.N));
    assert(dst.dmax() == src.dmax);
    for (int i = 0; i < dst.size(); ++i) {
        V slot = dst[std::size_t(i)];
        auto dt = slot.to_tuple();
        auto st = src.to_tuple();
        detail::for_each_field<V>([&](auto Ic) {
            constexpr std::size_t k = Ic;
            auto& d = std::get<k>(dt);
            const auto& s = std::get<k>(st);
            const std::size_t n = std::min<std::size_t>(d.size(), s.size());
            for (std::size_t j = 0; j < n; ++j) d[j] = s[j];
        });
    }
}

// ---------------------------------------------------------------------------
// Conversions back to host-side types.
// ---------------------------------------------------------------------------

// Build a host Graph from a single-isomer RSRAdjacencyView (any K).
template <class K>
inline Graph to_graph(const Spanify::RSRAdjacencyView<K>& v)
{
    Graph g(std::size_t(v.N), int(v.dmax));
    for (int u = 0; u < int(v.N); ++u) {
        g.deg[u] = v.deg[u];
        for (int k = 0; k < v.dmax; ++k)
            g.neighbours[u * v.dmax + k] = node_t(v.neighbours[u * v.dmax + k]);
    }
    return g;
}

// Convert an xyz span of std::array<T,3> to a host vector<coord3d>.
template <class T>
inline std::vector<coord3d> to_points(std::span<const std::array<T, 3>> xyz)
{
    std::vector<coord3d> out(xyz.size());
    for (std::size_t u = 0; u < xyz.size(); ++u)
        out[u] = coord3d(double(xyz[u][0]), double(xyz[u][1]), double(xyz[u][2]));
    return out;
}

} // namespace batch
