#pragma once

#include <vector>
#include <span>
#include <cstdint>
#include <cassert>
#include <algorithm>
#include <iostream>
#include <initializer_list>

#include "planar_csr.hh"

namespace Spanify {

// ---------------------------------------------------------------------------
// DenseGraph: fixed-stride flat array representation for oriented planar graphs.
//
// Each vertex has Dmax slots; a degree counter tracks how many are filled.
// Supports O(Dmax) insert/remove at arbitrary cyclic positions.
//
// Also serves as a backward-compatible replacement for vector<vector<node_t>>
// via RowProxy: operator[] returns a proxy with push_back, size, begin/end,
// operator=, and implicit conversion to vector.
//
// Template parameters:
//   Dmax — maximum degree per vertex (3 for cubic, 6 for fullerene duals, etc.)
//   K    — index type: int32_t (CPU), uint16_t (GPU)
// ---------------------------------------------------------------------------
template<int Dmax, typename K = int32_t>
struct DenseGraph {
    K Nv = 0;
    std::vector<K> values;         // [Nv * Dmax] — neighbor IDs, row-major
    std::vector<uint8_t> deg;      // [Nv] — current degree of each vertex

    // --- Row proxy types for operator[] backward compat ---

    struct ConstRowProxy {
        const DenseGraph& g;
        K u;

        const K& operator[](int i) const { return g.values[u * Dmax + i]; }
        int size() const { return g.deg[u]; }
        bool empty() const { return g.deg[u] == 0; }
        const K* begin() const { return g.values.data() + u * Dmax; }
        const K* end() const { return begin() + g.deg[u]; }
        const K* data() const { return begin(); }

        operator std::vector<K>() const { return std::vector<K>(begin(), end()); }

        friend std::ostream& operator<<(std::ostream& s, const ConstRowProxy& r) {
            s << '[';
            for (int i = 0; i < r.size(); ++i) { if (i) s << ','; s << r[i]; }
            s << ']';
            return s;
        }
    };

    struct RowProxy {
        DenseGraph& g;
        K u;

        K& operator[](int i) { return g.values[u * Dmax + i]; }
        const K& operator[](int i) const { return g.values[u * Dmax + i]; }
        int size() const { return g.deg[u]; }
        bool empty() const { return g.deg[u] == 0; }
        K* begin() { return g.values.data() + u * Dmax; }
        K* end() { return begin() + g.deg[u]; }
        const K* begin() const { return g.values.data() + u * Dmax; }
        const K* end() const { return begin() + g.deg[u]; }
        K* data() { return begin(); }
        const K* data() const { return begin(); }

        void push_back(K v) {
            assert(g.deg[u] < Dmax);
            g.values[u * Dmax + g.deg[u]] = v;
            g.deg[u]++;
        }

        void clear() { g.deg[u] = 0; }

        void resize(int n) {
            assert(n <= Dmax);
            for (int i = g.deg[u]; i < n; ++i)
                g.values[u * Dmax + i] = K(-1);
            g.deg[u] = uint8_t(n);
        }

        RowProxy& operator=(const std::vector<K>& v) {
            assert(int(v.size()) <= Dmax);
            g.deg[u] = uint8_t(v.size());
            for (int i = 0; i < int(v.size()); ++i)
                g.values[u * Dmax + i] = v[i];
            return *this;
        }

        RowProxy& operator=(std::initializer_list<K> v) {
            assert(int(v.size()) <= Dmax);
            g.deg[u] = uint8_t(v.size());
            int i = 0;
            for (K x : v) g.values[u * Dmax + i++] = x;
            return *this;
        }

        operator std::vector<K>() const { return std::vector<K>(begin(), end()); }

        friend std::ostream& operator<<(std::ostream& s, const RowProxy& r) {
            s << '[';
            for (int i = 0; i < r.size(); ++i) { if (i) s << ','; s << r[i]; }
            s << ']';
            return s;
        }
    };

    // --- Constructors ---

    DenseGraph() = default;

    explicit DenseGraph(K Nv) : Nv(Nv), values(Nv * Dmax, K(-1)), deg(Nv, 0) {}

    // Fill constructor: create Nv vertices, each initialized with the given row
    DenseGraph(K Nv, const std::vector<K>& initial_row)
        : Nv(Nv), values(Nv * Dmax, K(-1)), deg(Nv, uint8_t(initial_row.size())) {
        assert(int(initial_row.size()) <= Dmax);
        for (K v = 0; v < Nv; ++v)
            for (int i = 0; i < int(initial_row.size()); ++i)
                values[v * Dmax + i] = initial_row[i];
    }

    // Brace-init constructor: neighbours_t{{1,2,3},{4,5,6},...}
    DenseGraph(std::initializer_list<std::initializer_list<K>> adj)
        : Nv(K(adj.size())), values(adj.size() * Dmax, K(-1)), deg(adj.size(), 0) {
        int v = 0;
        for (auto& row : adj) {
            assert(int(row.size()) <= Dmax);
            deg[v] = uint8_t(row.size());
            int i = 0;
            for (K x : row) values[v * Dmax + i++] = x;
            v++;
        }
    }

    // Converting constructor from vector<vector<K>> (backward compat with neighbours_t)
    DenseGraph(const std::vector<std::vector<K>>& adj)
        : Nv(K(adj.size())), values(adj.size() * Dmax, K(-1)), deg(adj.size(), 0) {
        for (int v = 0; v < int(adj.size()); ++v) {
            assert(int(adj[v].size()) <= Dmax);
            deg[v] = uint8_t(adj[v].size());
            for (int i = 0; i < int(adj[v].size()); ++i)
                values[v * Dmax + i] = adj[v][i];
        }
    }

    // --- operator[] returning row proxies ---

    RowProxy operator[](K u) { return {*this, u}; }
    ConstRowProxy operator[](K u) const { return {*this, u}; }

    // --- Span-based accessors ---

    std::span<const K> nbrs(K u) const {
        return {values.data() + u * Dmax, deg[u]};
    }

    std::span<K> nbrs_mut(K u) {
        return {values.data() + u * Dmax, deg[u]};
    }

    int degree(K u) const { return deg[u]; }

    // --- Mutation (vertex-targeted) ---

    // Append v to end of u's neighbour list. O(1).
    void push_back(K u, K v) {
        assert(deg[u] < Dmax);
        values[u * Dmax + deg[u]] = v;
        deg[u]++;
    }

    // Insert v at position pos in u's neighbour list (shifting later entries). O(Dmax).
    void insert_at(K u, K v, int pos) {
        assert(deg[u] < Dmax);
        K* row = values.data() + u * Dmax;
        int d = deg[u];
        for (int i = d; i > pos; --i)
            row[i] = row[i-1];
        row[pos] = v;
        deg[u]++;
    }

    // Remove entry at position pos in u's neighbour list (shifting later entries). O(Dmax).
    void erase_at(K u, int pos) {
        K* row = values.data() + u * Dmax;
        int d = deg[u];
        for (int i = pos; i < d - 1; ++i)
            row[i] = row[i+1];
        row[d-1] = K(-1);
        deg[u]--;
    }

    // Find position of v in u's neighbour list. Returns -1 if not found.
    int find(K u, K v) const {
        const K* row = values.data() + u * Dmax;
        for (int i = 0; i < deg[u]; ++i)
            if (row[i] == v) return i;
        return -1;
    }

    // --- Backward compat with vector<vector<>> interface ---

    // Number of vertices (like vector::size())
    int size() const { return Nv; }

    // Resize to new_N vertices (like vector::resize())
    void resize(int new_N) {
        values.resize(new_N * Dmax, K(-1));
        deg.resize(new_N, 0);
        Nv = K(new_N);
    }

    // Add a vertex with given neighbors (like vector::push_back(vector))
    void push_back(const std::vector<K>& row) {
        assert(int(row.size()) <= Dmax);
        values.resize((Nv + 1) * Dmax, K(-1));
        deg.push_back(uint8_t(row.size()));
        for (int i = 0; i < int(row.size()); ++i)
            values[Nv * Dmax + i] = row[i];
        Nv++;
    }

    // Remove last vertex
    void pop_back() {
        assert(Nv > 0);
        Nv--;
        values.resize(Nv * Dmax);
        deg.resize(Nv);
    }

    // Explicit conversion back to vector<vector<K>>
    std::vector<std::vector<K>> to_vectors() const {
        std::vector<std::vector<K>> adj(Nv);
        for (K v = 0; v < Nv; ++v) {
            adj[v].resize(deg[v]);
            for (int i = 0; i < deg[v]; ++i)
                adj[v][i] = values[v * Dmax + i];
        }
        return adj;
    }
};

// Convenience aliases
using DenseCubic         = DenseGraph<3>;
using DenseTriangulation = DenseGraph<6>;
using DenseFulleroid     = DenseGraph<10>;

// --- ostream support ---

template<int Dmax, typename K>
std::ostream& operator<<(std::ostream& s, const DenseGraph<Dmax,K>& g) {
    s << '[';
    for (K v = 0; v < g.Nv; ++v) {
        if (v) s << ',';
        s << '[';
        for (int i = 0; i < g.deg[v]; ++i) {
            if (i) s << ',';
            s << g.values[v * Dmax + i];
        }
        s << ']';
    }
    s << ']';
    return s;
}

// ---------------------------------------------------------------------------
// Conversion: Dense -> CSR (freeze). O(N).
// ---------------------------------------------------------------------------
template<int Dmax, typename K>
PlanarCSR<Owned, K> freeze(const DenseGraph<Dmax, K>& g) {
    PlanarCSR<Owned, K> csr;
    csr.N = g.Nv;
    csr.offsets.resize(g.Nv + 1);
    csr.offsets[0] = 0;
    for (K v = 0; v < g.Nv; ++v)
        csr.offsets[v+1] = csr.offsets[v] + K(g.deg[v]);
    csr.n_arcs = csr.offsets[g.Nv];

    csr.values.resize(csr.n_arcs);
    for (K v = 0; v < g.Nv; ++v) {
        K base = csr.offsets[v];
        for (int i = 0; i < g.deg[v]; ++i)
            csr.values[base + i] = g.values[v * Dmax + i];
    }

    csr.twin = compute_twin(csr.N, csr.offsets, csr.values);
    return csr;
}

// ---------------------------------------------------------------------------
// Conversion: CSR -> Dense (thaw). O(N).
// Caller must choose Dmax large enough for the graph's maximum degree.
// ---------------------------------------------------------------------------
template<int Dmax, typename K>
DenseGraph<Dmax, K> thaw(const PlanarCSR<Owned, K>& csr) {
    DenseGraph<Dmax, K> g(csr.N);
    for (K v = 0; v < csr.N; ++v) {
        K d = csr.offsets[v+1] - csr.offsets[v];
        assert(d <= Dmax);
        g.deg[v] = uint8_t(d);
        for (K i = 0; i < d; ++i)
            g.values[v * Dmax + i] = csr.values[csr.offsets[v] + i];
    }
    return g;
}

// ---------------------------------------------------------------------------
// Conversion: vector<vector<int>> -> Dense. O(N).
// ---------------------------------------------------------------------------
template<int Dmax, typename K = int32_t>
DenseGraph<Dmax, K> dense_from_neighbours(int N_verts,
                                           const std::vector<std::vector<int>>& neighbours) {
    DenseGraph<Dmax, K> g{K(N_verts)};
    for (int v = 0; v < N_verts; ++v) {
        assert(int(neighbours[v].size()) <= Dmax);
        g.deg[v] = uint8_t(neighbours[v].size());
        for (int i = 0; i < int(neighbours[v].size()); ++i)
            g.values[v * Dmax + i] = K(neighbours[v][i]);
    }
    return g;
}

// ---------------------------------------------------------------------------
// Conversion: Dense -> vector<vector<int>>. O(N).
// ---------------------------------------------------------------------------
template<int Dmax, typename K>
std::vector<std::vector<int>> to_neighbours(const DenseGraph<Dmax, K>& g) {
    std::vector<std::vector<int>> adj(g.Nv);
    for (K v = 0; v < g.Nv; ++v) {
        adj[v].resize(g.deg[v]);
        for (int i = 0; i < g.deg[v]; ++i)
            adj[v][i] = int(g.values[v * Dmax + i]);
    }
    return adj;
}

} // namespace Spanify
