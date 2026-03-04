#pragma once

#include <vector>
#include <span>
#include <cstdint>
#include <cassert>
#include <algorithm>

namespace Spanify {

// Storage aliases
template<typename T> using Owned   = std::vector<T>;
template<typename T> using View    = std::span<const T>;
template<typename T> using MutView = std::span<T>;

// ---------------------------------------------------------------------------
// PlanarCSR: Compressed Sparse Row representation of an oriented planar graph.
//
// Template parameters:
//   Array<T> — storage policy: Owned<T> = vector<T>, View<T> = span<const T>
//   K        — index type: int32_t (CPU), uint16_t (GPU)
// ---------------------------------------------------------------------------
template<template<typename> class Array, typename K = int32_t>
struct PlanarCSR {
    K N       = 0;    // Number of vertices
    K n_arcs  = 0;    // Number of directed arcs (= 2 * edges)
    Array<K> offsets;  // [N+1]    — row pointers (offsets[v]..offsets[v+1])
    Array<K> values;   // [n_arcs] — target vertex of each arc
    Array<K> twin;     // [n_arcs] — position of the reverse arc
};

// Convenience aliases
using CSRGraph     = PlanarCSR<Owned,  int32_t>;
using CSRGraphView = PlanarCSR<View,   int32_t>;

// ---------------------------------------------------------------------------
// Conversion: Owned → View
// ---------------------------------------------------------------------------
template<typename K>
PlanarCSR<View, K> as_view(const PlanarCSR<Owned, K>& g) {
    return {g.N, g.n_arcs, g.offsets, g.values, g.twin};
}

// ---------------------------------------------------------------------------
// Core navigation (O(1) per call)
// ---------------------------------------------------------------------------

template<template<typename> class A, typename K>
K degree(const PlanarCSR<A,K>& g, K v) {
    return g.offsets[v+1] - g.offsets[v];
}

template<template<typename> class A, typename K>
std::span<const K> neighbours(const PlanarCSR<A,K>& g, K v) {
    return {&g.values[g.offsets[v]], std::size_t(degree(g, v))};
}

// Next arc CCW around the same source vertex.
// Caller provides source vertex u (avoids reverse lookup).
template<template<typename> class A, typename K>
K next_arc(const PlanarCSR<A,K>& g, K p, K u) {
    K base = g.offsets[u];
    K deg  = g.offsets[u+1] - base;
    return base + (p - base + 1) % deg;
}

// Previous arc CCW around the same source vertex.
template<template<typename> class A, typename K>
K prev_arc(const PlanarCSR<A,K>& g, K p, K u) {
    K base = g.offsets[u];
    K deg  = g.offsets[u+1] - base;
    return base + (p - base + deg - 1) % deg;
}

// Face traversal: cross edge via twin, then step CW around target.
// Returns the next arc on the same face (face to the left of arc p).
// Equivalent to: next_on_face(u,v) = prev(v,u) in the vertex-based API.
template<template<typename> class A, typename K>
K next_on_face_arc(const PlanarCSR<A,K>& g, K p) {
    K twin_p = g.twin[p];
    K v = g.values[p];  // twin_p originates at v
    return prev_arc(g, twin_p, v);
}

// Previous arc on the same face.
template<template<typename> class A, typename K>
K prev_on_face_arc(const PlanarCSR<A,K>& g, K p) {
    K u = source(g, p);
    K q = next_arc(g, p, u);
    return g.twin[q];
}

// Source vertex of arc p (O(log N) binary search in offsets).
// For cubic graphs where offsets[v]=3v, caller should use p/3 directly.
template<template<typename> class A, typename K>
K source(const PlanarCSR<A,K>& g, K p) {
    // Find v such that offsets[v] <= p < offsets[v+1]
    // offsets is sorted, so use upper_bound on offsets[0..N]
    auto it = std::upper_bound(&g.offsets[0], &g.offsets[g.N+1], p);
    return K(std::distance(&g.offsets[0], it) - 1);
}

// ---------------------------------------------------------------------------
// Backward-compatible vertex-based API (O(degree) for arc lookup)
// ---------------------------------------------------------------------------

// Find position of arc u→v. Returns -1 if not found.
template<template<typename> class A, typename K>
K find_arc(const PlanarCSR<A,K>& g, K u, K v) {
    for (K p = g.offsets[u]; p < g.offsets[u+1]; ++p)
        if (g.values[p] == v) return p;
    return K(-1);
}

// Next CCW neighbour of u after v (vertex-based, O(degree)).
template<template<typename> class A, typename K>
K next(const PlanarCSR<A,K>& g, K u, K v) {
    K p = find_arc(g, u, v);
    assert(p != K(-1));
    return g.values[next_arc(g, p, u)];
}

// Previous CCW neighbour of u before v (vertex-based, O(degree)).
template<template<typename> class A, typename K>
K prev(const PlanarCSR<A,K>& g, K u, K v) {
    K p = find_arc(g, u, v);
    assert(p != K(-1));
    return g.values[prev_arc(g, p, u)];
}

// Next vertex on face containing arc u→v (vertex-based, O(degree)).
// This is the vertex after v when walking the face to the left of u→v.
// Equivalent to old Graph::next_on_face(u, v) = prev(v, u).
template<template<typename> class A, typename K>
K next_on_face(const PlanarCSR<A,K>& g, K u, K v) {
    return prev(g, v, u);
}

// ---------------------------------------------------------------------------
// FaceData: optional per-arc face assignment (computed lazily)
// ---------------------------------------------------------------------------
template<template<typename> class Array, typename K = int32_t>
struct FaceData {
    K n_faces = 0;
    Array<K> face_id;      // [n_arcs] — face index per arc
    Array<K> face_start;   // [n_faces] — one representative arc per face
};

// Compute face assignments from topology. O(n_arcs).
template<typename K>
FaceData<Owned, K> compute_faces(const PlanarCSR<View, K>& g) {
    FaceData<Owned, K> fd;
    fd.face_id.resize(g.n_arcs, K(-1));

    K f = 0;
    for (K p = 0; p < g.n_arcs; ++p) {
        if (fd.face_id[p] != K(-1)) continue;

        fd.face_start.push_back(p);
        K q = p;
        do {
            fd.face_id[q] = f;
            q = next_on_face_arc(g, q);
        } while (q != p);
        ++f;
    }
    fd.n_faces = f;
    return fd;
}

// Size of face f: count arcs in the face cycle.
template<template<typename> class AG, template<typename> class AF, typename K>
K face_size(const PlanarCSR<AG,K>& g, const FaceData<AF,K>& fd, K f) {
    K start = fd.face_start[f];
    K count = 0;
    K p = start;
    do {
        ++count;
        p = next_on_face_arc(g, p);
    } while (p != start);
    return count;
}

// Collect vertices of face f into a vector.
template<template<typename> class AG, template<typename> class AF, typename K>
std::vector<K> face_vertices(const PlanarCSR<AG,K>& g, const FaceData<AF,K>& fd, K f) {
    std::vector<K> verts;
    K start = fd.face_start[f];
    K p = start;
    do {
        verts.push_back(source(g, p));
        p = next_on_face_arc(g, p);
    } while (p != start);
    return verts;
}

// ---------------------------------------------------------------------------
// Construction: compute_twin
// ---------------------------------------------------------------------------

// Compute twin array from offsets + values. O(sum of degree^2).
// For cubic graphs this is O(9N) = O(N).
template<typename K>
std::vector<K> compute_twin(K N, const std::vector<K>& offsets,
                            const std::vector<K>& values)
{
    K n_arcs = K(values.size());
    std::vector<K> twin(n_arcs, K(-1));
    for (K u = 0; u < N; ++u) {
        for (K p = offsets[u]; p < offsets[u+1]; ++p) {
            K v = values[p];
            // Find arc v→u
            for (K q = offsets[v]; q < offsets[v+1]; ++q) {
                if (values[q] == u) {
                    twin[p] = q;
                    break;
                }
            }
        }
    }
    return twin;
}

// ---------------------------------------------------------------------------
// CSRBuilder: incremental graph construction
// ---------------------------------------------------------------------------
template<typename K = int32_t>
struct CSRBuilder {
    K N;
    std::vector<std::vector<K>> adj;

    explicit CSRBuilder(K N) : N(N), adj(N) {}

    void add_arc(K u, K v) { adj[u].push_back(v); }
    void add_edge(K u, K v) { add_arc(u, v); add_arc(v, u); }

    PlanarCSR<Owned, K> freeze() const {
        PlanarCSR<Owned, K> g;
        g.N = N;
        g.offsets.resize(N + 1);
        g.offsets[0] = 0;
        for (K v = 0; v < N; ++v)
            g.offsets[v+1] = g.offsets[v] + K(adj[v].size());
        g.n_arcs = g.offsets[N];

        g.values.reserve(g.n_arcs);
        for (K v = 0; v < N; ++v)
            g.values.insert(g.values.end(), adj[v].begin(), adj[v].end());

        g.twin = compute_twin(g.N, g.offsets, g.values);
        return g;
    }
};

// ---------------------------------------------------------------------------
// Conversion from old neighbours_t (vector<vector<int>>) — for testing
// ---------------------------------------------------------------------------
template<typename K = int32_t>
PlanarCSR<Owned, K> from_neighbours(int N_verts,
                                     const std::vector<std::vector<int>>& neighbours)
{
    PlanarCSR<Owned, K> g;
    g.N = K(N_verts);
    g.offsets.resize(N_verts + 1);
    g.offsets[0] = 0;
    for (int v = 0; v < N_verts; ++v)
        g.offsets[v+1] = g.offsets[v] + K(neighbours[v].size());
    g.n_arcs = g.offsets[N_verts];

    g.values.reserve(g.n_arcs);
    for (int v = 0; v < N_verts; ++v)
        for (int u : neighbours[v])
            g.values.push_back(K(u));

    g.twin = compute_twin(g.N, g.offsets, g.values);
    return g;
}

// Validate CSR invariants. Returns true if all checks pass.
template<template<typename> class A, typename K>
bool validate(const PlanarCSR<A,K>& g) {
    if (g.offsets[0] != 0) return false;
    if (g.offsets[g.N] != g.n_arcs) return false;

    // Check offsets are non-decreasing
    for (K v = 0; v < g.N; ++v)
        if (g.offsets[v+1] < g.offsets[v]) return false;

    // Check twin is an involution: twin[twin[p]] == p
    for (K p = 0; p < g.n_arcs; ++p) {
        if (g.twin[p] < 0 || g.twin[p] >= g.n_arcs) return false;
        if (g.twin[g.twin[p]] != p) return false;
    }

    // Check twin targets match: values[twin[p]] == source(p)
    for (K v = 0; v < g.N; ++v) {
        for (K p = g.offsets[v]; p < g.offsets[v+1]; ++p) {
            K tp = g.twin[p];
            if (g.values[tp] != v) return false;
        }
    }

    return true;
}

} // namespace Spanify
