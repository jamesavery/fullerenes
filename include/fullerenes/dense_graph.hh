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
// Each vertex has dmax slots; a degree counter tracks how many are filled.
// Supports O(dmax) insert/remove at arbitrary cyclic positions.
//
// Also serves as a backward-compatible replacement for vector<vector<node_t>>
// via RowProxy: operator[] returns a proxy with push_back, size, begin/end,
// operator=, and implicit conversion to vector.
//
// dmax is a runtime parameter set at construction:
//   - Auto-detected from data (brace-init, vector<vector>, fill constructors)
//   - Explicit for empty graphs: DenseGraph(Nv, dmax)
//   - CubicGraph uses dmax=3, Triangulation dmax=10, etc.
//
// Template parameter:
//   K — index type: int32_t (CPU), uint16_t (GPU)
// ---------------------------------------------------------------------------
template<typename K = int32_t>
struct DenseGraph {
    K Nv = 0;
    int dmax = 0;
    std::vector<K> values;         // [Nv * dmax] — neighbor IDs, row-major
    std::vector<uint8_t> deg;      // [Nv] — current degree of each vertex

    // --- Row proxy types for operator[] backward compat ---

    struct ConstRowProxy {
        const DenseGraph& g;
        K u;

        const K& operator[](int i) const { return g.values[u * g.dmax + i]; }
        int size() const { return g.deg[u]; }
        bool empty() const { return g.deg[u] == 0; }
        const K* begin() const { return g.values.data() + u * g.dmax; }
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

        K& operator[](int i) { return g.values[u * g.dmax + i]; }
        const K& operator[](int i) const { return g.values[u * g.dmax + i]; }
        int size() const { return g.deg[u]; }
        bool empty() const { return g.deg[u] == 0; }
        K* begin() { return g.values.data() + u * g.dmax; }
        K* end() { return begin() + g.deg[u]; }
        const K* begin() const { return g.values.data() + u * g.dmax; }
        const K* end() const { return begin() + g.deg[u]; }
        K* data() { return begin(); }
        const K* data() const { return begin(); }

        void push_back(K v) {
            assert(g.deg[u] < g.dmax);
            g.values[u * g.dmax + g.deg[u]] = v;
            g.deg[u]++;
        }

        void clear() { g.deg[u] = 0; }

        void resize(int n) {
            assert(n <= g.dmax);
            for (int i = g.deg[u]; i < n; ++i)
                g.values[u * g.dmax + i] = K(-1);
            g.deg[u] = uint8_t(n);
        }

        RowProxy& operator=(const std::vector<K>& v) {
            assert(int(v.size()) <= g.dmax);
            g.deg[u] = uint8_t(v.size());
            for (int i = 0; i < int(v.size()); ++i)
                g.values[u * g.dmax + i] = v[i];
            return *this;
        }

        RowProxy& operator=(std::initializer_list<K> v) {
            assert(int(v.size()) <= g.dmax);
            g.deg[u] = uint8_t(v.size());
            int i = 0;
            for (K x : v) g.values[u * g.dmax + i++] = x;
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

    // Explicit: Nv vertices, given stride, all empty.
    explicit DenseGraph(K Nv, int dmax) : Nv(Nv), dmax(dmax), values(Nv * dmax, K(-1)), deg(Nv, 0) {}

    // Fill constructor: Nv vertices, each initialized with the given row. dmax = row size.
    DenseGraph(K Nv, const std::vector<K>& initial_row)
        : Nv(Nv), dmax(int(initial_row.size())), values(Nv * dmax, K(-1)),
          deg(Nv, uint8_t(initial_row.size())) {
        for (K v = 0; v < Nv; ++v)
            for (int i = 0; i < dmax; ++i)
                values[v * dmax + i] = initial_row[i];
    }

    // Brace-init constructor: DenseGraph{{1,2,3},{4,5,6},...}. dmax = max row size.
    DenseGraph(std::initializer_list<std::initializer_list<K>> adj)
        : Nv(K(adj.size())), dmax(0) {
        for (auto& row : adj) dmax = std::max(dmax, int(row.size()));
        values.resize(Nv * dmax, K(-1));
        deg.resize(Nv, 0);
        int v = 0;
        for (auto& row : adj) {
            deg[v] = uint8_t(row.size());
            int i = 0;
            for (K x : row) values[v * dmax + i++] = x;
            v++;
        }
    }

    // Converting constructor from vector<vector<K>>. dmax = max row size.
    DenseGraph(const std::vector<std::vector<K>>& adj)
        : Nv(K(adj.size())), dmax(0) {
        for (auto& row : adj) dmax = std::max(dmax, int(row.size()));
        values.resize(Nv * dmax, K(-1));
        deg.resize(Nv, 0);
        for (int v = 0; v < int(adj.size()); ++v) {
            deg[v] = uint8_t(adj[v].size());
            for (int i = 0; i < int(adj[v].size()); ++i)
                values[v * dmax + i] = adj[v][i];
        }
    }

    // --- operator[] returning row proxies ---

    RowProxy operator[](K u) { return {*this, u}; }
    ConstRowProxy operator[](K u) const { return {*this, u}; }

    // --- Span-based accessors ---

    std::span<const K> nbrs(K u) const {
        return {values.data() + u * dmax, (size_t)deg[u]};
    }

    std::span<K> nbrs_mut(K u) {
        return {values.data() + u * dmax, (size_t)deg[u]};
    }

    int degree(K u) const { return deg[u]; }

    // --- Mutation (vertex-targeted) ---

    void push_back(K u, K v) {
        assert(deg[u] < dmax);
        values[u * dmax + deg[u]] = v;
        deg[u]++;
    }

    void insert_at(K u, K v, int pos) {
        assert(deg[u] < dmax);
        K* row = values.data() + u * dmax;
        int d = deg[u];
        for (int i = d; i > pos; --i)
            row[i] = row[i-1];
        row[pos] = v;
        deg[u]++;
    }

    void erase_at(K u, int pos) {
        K* row = values.data() + u * dmax;
        int d = deg[u];
        for (int i = pos; i < d - 1; ++i)
            row[i] = row[i+1];
        row[d-1] = K(-1);
        deg[u]--;
    }

    int find(K u, K v) const {
        const K* row = values.data() + u * dmax;
        for (int i = 0; i < deg[u]; ++i)
            if (row[i] == v) return i;
        return -1;
    }

    // --- Backward compat with vector<vector<>> interface ---

    int size() const { return Nv; }

    void resize(int new_N) {
        values.resize(new_N * dmax, K(-1));
        deg.resize(new_N, 0);
        Nv = K(new_N);
    }

    void push_back(const std::vector<K>& row) {
        assert(int(row.size()) <= dmax);
        values.resize((Nv + 1) * dmax, K(-1));
        deg.push_back(uint8_t(row.size()));
        for (int i = 0; i < int(row.size()); ++i)
            values[Nv * dmax + i] = row[i];
        Nv++;
    }

    void pop_back() {
        assert(Nv > 0);
        Nv--;
        values.resize(Nv * dmax);
        deg.resize(Nv);
    }

    std::vector<std::vector<K>> to_vectors() const {
        std::vector<std::vector<K>> adj(Nv);
        for (K v = 0; v < Nv; ++v) {
            adj[v].resize(deg[v]);
            for (int i = 0; i < deg[v]; ++i)
                adj[v][i] = values[v * dmax + i];
        }
        return adj;
    }

    // --- Stride change ---

    // Create a copy with a different dmax (narrowing or widening).
    DenseGraph restride(int new_dmax) const {
        DenseGraph g(Nv, new_dmax);
        for (K v = 0; v < Nv; ++v) {
            assert(deg[v] <= new_dmax);
            g.deg[v] = deg[v];
            for (int i = 0; i < deg[v]; ++i)
                g.values[v * new_dmax + i] = values[v * dmax + i];
        }
        return g;
    }
};

// --- ostream support ---

template<typename K>
std::ostream& operator<<(std::ostream& s, const DenseGraph<K>& g) {
    s << '[';
    for (K v = 0; v < g.Nv; ++v) {
        if (v) s << ',';
        s << '[';
        for (int i = 0; i < g.deg[v]; ++i) {
            if (i) s << ',';
            s << g.values[v * g.dmax + i];
        }
        s << ']';
    }
    s << ']';
    return s;
}

// ---------------------------------------------------------------------------
// Conversion: Dense -> CSR (freeze). O(N).
// ---------------------------------------------------------------------------
template<typename K>
PlanarCSR<Owned, K> freeze(const DenseGraph<K>& g) {
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
            csr.values[base + i] = g.values[v * g.dmax + i];
    }

    csr.twin = compute_twin(csr.N, csr.offsets, csr.values);
    return csr;
}

// ---------------------------------------------------------------------------
// Conversion: CSR -> Dense (thaw). O(N).
// Caller must choose dmax large enough for the graph's maximum degree.
// ---------------------------------------------------------------------------
template<typename K>
DenseGraph<K> thaw(const PlanarCSR<Owned, K>& csr, int dmax) {
    DenseGraph<K> g(csr.N, dmax);
    for (K v = 0; v < csr.N; ++v) {
        K d = csr.offsets[v+1] - csr.offsets[v];
        assert(d <= dmax);
        g.deg[v] = uint8_t(d);
        for (K i = 0; i < d; ++i)
            g.values[v * dmax + i] = csr.values[csr.offsets[v] + i];
    }
    return g;
}

// ---------------------------------------------------------------------------
// Conversion: vector<vector<int>> -> Dense. O(N).
// ---------------------------------------------------------------------------
template<typename K = int32_t>
DenseGraph<K> dense_from_neighbours(int N_verts,
                                     const std::vector<std::vector<int>>& neighbours,
                                     int dmax) {
    DenseGraph<K> g(K(N_verts), dmax);
    for (int v = 0; v < N_verts; ++v) {
        assert(int(neighbours[v].size()) <= dmax);
        g.deg[v] = uint8_t(neighbours[v].size());
        for (int i = 0; i < int(neighbours[v].size()); ++i)
            g.values[v * dmax + i] = K(neighbours[v][i]);
    }
    return g;
}

// ---------------------------------------------------------------------------
// Conversion: Dense -> vector<vector<int>>. O(N).
// ---------------------------------------------------------------------------
template<typename K>
std::vector<std::vector<int>> to_neighbours(const DenseGraph<K>& g) {
    std::vector<std::vector<int>> adj(g.Nv);
    for (K v = 0; v < g.Nv; ++v) {
        adj[v].resize(g.deg[v]);
        for (int i = 0; i < g.deg[v]; ++i)
            adj[v][i] = int(g.values[v * g.dmax + i]);
    }
    return adj;
}

} // namespace Spanify
