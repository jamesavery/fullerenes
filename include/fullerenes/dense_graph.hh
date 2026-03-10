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

template<typename K> struct OwnedDenseGraph;  // forward declaration

// ---------------------------------------------------------------------------
// DenseGraph: non-owning view over fixed-stride flat adjacency data.
//
// Each vertex has dmax slots; a degree counter tracks how many are filled.
// Supports O(dmax) insert/remove at arbitrary cyclic positions.
//
// operator[] returns a span over the row's active entries.
//
// Template parameter:
//   K -- index type: int32_t (CPU), uint16_t (GPU)
// ---------------------------------------------------------------------------
template<typename K = int32_t>
struct DenseGraph {
    K N = 0;
    int dmax = 0;
    std::span<K> values;         // [N * dmax] -- neighbor IDs, row-major
    std::span<uint8_t> deg;      // [N] -- current degree of each vertex

    // --- Constructors ---

    DenseGraph() = default;

    // View constructor: wrap existing spans.
    DenseGraph(K N, int dmax, std::span<K> values, std::span<uint8_t> deg)
        : N(N), dmax(dmax), values(values), deg(deg) {}

    // --- operator[] returns span over active entries ---

    std::span<K> operator[](K u) {
        return {values.data() + u * dmax, (size_t)deg[u]};
    }
    std::span<const K> operator[](K u) const {
        return {values.data() + u * dmax, (size_t)deg[u]};
    }

    // --- Span-based accessors ---

    std::span<const K> nbrs(K u) const {
        return {values.data() + u * dmax, (size_t)deg[u]};
    }

    std::span<K> nbrs_mut(K u) {
        return {values.data() + u * dmax, (size_t)deg[u]};
    }

    int degree(K u) const { return deg[u]; }

    // --- Per-row mutation (no reallocation) ---

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

    void clear_row(K u) { deg[u] = 0; }

    void resize_row(K u, int n) {
        assert(n <= dmax);
        for (int i = deg[u]; i < n; ++i)
            values[u * dmax + i] = K(-1);
        deg[u] = uint8_t(n);
    }

    void assign_row(K u, std::initializer_list<K> v) {
        assert(int(v.size()) <= dmax);
        deg[u] = uint8_t(v.size());
        int i = 0;
        for (K x : v) values[u * dmax + i++] = x;
    }

    void assign_row(K u, const std::vector<K>& v) {
        assert(int(v.size()) <= dmax);
        deg[u] = uint8_t(v.size());
        for (int i = 0; i < int(v.size()); ++i)
            values[u * dmax + i] = v[i];
    }

    void assign_row(K u, std::span<const K> v) {
        assert(int(v.size()) <= dmax);
        deg[u] = uint8_t(v.size());
        for (int i = 0; i < int(v.size()); ++i)
            values[u * dmax + i] = v[i];
    }

    // --- Backward compat with vector<vector<>> interface ---

    int size() const { return N; }

    std::vector<std::vector<K>> to_vectors() const {
        std::vector<std::vector<K>> adj(N);
        for (K v = 0; v < N; ++v) {
            adj[v].resize(deg[v]);
            for (int i = 0; i < deg[v]; ++i)
                adj[v][i] = values[v * dmax + i];
        }
        return adj;
    }
};

// ---------------------------------------------------------------------------
// OwnedDenseGraph: owning version with vector storage.
// Inherits all view/mutation methods from DenseGraph.
// ---------------------------------------------------------------------------
template<typename K = int32_t>
struct OwnedDenseGraph : DenseGraph<K> {
    std::vector<K> owned_values;
    std::vector<uint8_t> owned_deg;

    void repoint() {
        this->values = std::span<K>(owned_values.data(), owned_values.size());
        this->deg = std::span<uint8_t>(owned_deg.data(), owned_deg.size());
    }

    // Bring base push_back(K,K) into scope (otherwise hidden by push_back(vector))
    using DenseGraph<K>::push_back;

    // --- Constructors ---

    OwnedDenseGraph() = default;

    // Explicit: N vertices, given stride, all empty.
    explicit OwnedDenseGraph(K N, int dmax)
        : owned_values(N * dmax, K(-1)), owned_deg(N, 0) {
        this->N = N; this->dmax = dmax; repoint();
    }

    // Fill constructor: N vertices, each initialized with the given row.
    OwnedDenseGraph(K N, const std::vector<K>& initial_row)
        : owned_values(N * int(initial_row.size()), K(-1)),
          owned_deg(N, uint8_t(initial_row.size())) {
        this->N = N; this->dmax = int(initial_row.size()); repoint();
        for (K v = 0; v < N; ++v)
            for (int i = 0; i < this->dmax; ++i)
                this->values[v * this->dmax + i] = initial_row[i];
    }

    // Brace-init constructor: OwnedDenseGraph{{1,2,3},{4,5,6},...}.
    OwnedDenseGraph(std::initializer_list<std::initializer_list<K>> adj) {
        this->N = K(adj.size());
        this->dmax = 0;
        for (auto& row : adj) this->dmax = std::max(this->dmax, int(row.size()));
        owned_values.resize(this->N * this->dmax, K(-1));
        owned_deg.resize(this->N, 0);
        repoint();
        int v = 0;
        for (auto& row : adj) {
            this->deg[v] = uint8_t(row.size());
            int i = 0;
            for (K x : row) this->values[v * this->dmax + i++] = x;
            v++;
        }
    }

    // Converting constructor from vector<vector<K>>.
    OwnedDenseGraph(const std::vector<std::vector<K>>& adj) {
        this->N = K(adj.size());
        this->dmax = 0;
        for (auto& row : adj) this->dmax = std::max(this->dmax, int(row.size()));
        owned_values.resize(this->N * this->dmax, K(-1));
        owned_deg.resize(this->N, 0);
        repoint();
        for (int v = 0; v < int(adj.size()); ++v) {
            this->deg[v] = uint8_t(adj[v].size());
            for (int i = 0; i < int(adj[v].size()); ++i)
                this->values[v * this->dmax + i] = adj[v][i];
        }
    }

    // Converting constructor from DenseGraph view (copies data).
    OwnedDenseGraph(const DenseGraph<K>& v) {
        this->N = v.N; this->dmax = v.dmax;
        if (v.N > 0 && v.values.data()) {
            owned_values.assign(v.values.begin(), v.values.end());
            owned_deg.assign(v.deg.begin(), v.deg.end());
        }
        repoint();
    }

    // --- Rule of 5 ---

    OwnedDenseGraph(const OwnedDenseGraph& o)
        : owned_values(o.owned_values), owned_deg(o.owned_deg) {
        this->N = o.N; this->dmax = o.dmax; repoint();
    }

    OwnedDenseGraph(OwnedDenseGraph&& o) noexcept
        : owned_values(std::move(o.owned_values)), owned_deg(std::move(o.owned_deg)) {
        this->N = o.N; this->dmax = o.dmax; repoint();
        o.values = {}; o.deg = {}; o.N = 0;
    }

    OwnedDenseGraph& operator=(const OwnedDenseGraph& o) {
        if (this != &o) {
            this->N = o.N; this->dmax = o.dmax;
            owned_values = o.owned_values;
            owned_deg = o.owned_deg;
            repoint();
        }
        return *this;
    }

    OwnedDenseGraph& operator=(OwnedDenseGraph&& o) noexcept {
        this->N = o.N; this->dmax = o.dmax;
        owned_values = std::move(o.owned_values);
        owned_deg = std::move(o.owned_deg);
        repoint();
        o.values = {}; o.deg = {}; o.N = 0;
        return *this;
    }

    // --- Resize operations (require reallocation) ---

    void resize(int new_N) {
        owned_values.resize(new_N * this->dmax, K(-1));
        owned_deg.resize(new_N, 0);
        this->N = K(new_N);
        repoint();
    }

    void push_back(const std::vector<K>& row) {
        assert(int(row.size()) <= this->dmax);
        owned_values.resize((this->N + 1) * this->dmax, K(-1));
        owned_deg.push_back(uint8_t(row.size()));
        repoint();
        for (int i = 0; i < int(row.size()); ++i)
            this->values[this->N * this->dmax + i] = row[i];
        this->N++;
    }

    void pop_back() {
        assert(this->N > 0);
        this->N--;
        owned_values.resize(this->N * this->dmax);
        owned_deg.resize(this->N);
        repoint();
    }

    // --- Stride change ---

    OwnedDenseGraph restride(int new_dmax) const {
        OwnedDenseGraph g(this->N, new_dmax);
        for (K v = 0; v < this->N; ++v) {
            assert(this->deg[v] <= new_dmax);
            g.deg[v] = this->deg[v];
            for (int i = 0; i < this->deg[v]; ++i)
                g.values[v * new_dmax + i] = this->values[v * this->dmax + i];
        }
        return g;
    }
};

// --- ostream support ---

template<typename K>
std::ostream& operator<<(std::ostream& s, const DenseGraph<K>& g) {
    s << '[';
    for (K v = 0; v < g.N; ++v) {
        if (v) s << ',';
        s << '[';
        auto row = g[v];
        for (int i = 0; i < int(row.size()); ++i) {
            if (i) s << ',';
            s << row[i];
        }
        s << ']';
    }
    s << ']';
    return s;
}

// ---------------------------------------------------------------------------
// Free-function restride: works on any DenseGraph view. O(N*dmax).
// ---------------------------------------------------------------------------
template<typename K>
OwnedDenseGraph<K> restride(const DenseGraph<K>& src, int new_dmax) {
    OwnedDenseGraph<K> g(src.N, new_dmax);
    for (K v = 0; v < src.N; ++v) {
        assert(src.deg[v] <= new_dmax);
        g.deg[v] = src.deg[v];
        for (int i = 0; i < src.deg[v]; ++i)
            g.values[v * new_dmax + i] = src.values[v * src.dmax + i];
    }
    return g;
}

// ---------------------------------------------------------------------------
// Conversion: Dense -> CSR (freeze). O(N).
// ---------------------------------------------------------------------------
template<typename K>
PlanarCSR<Owned, K> freeze(const DenseGraph<K>& g) {
    PlanarCSR<Owned, K> csr;
    csr.N = g.N;
    csr.offsets.resize(g.N + 1);
    csr.offsets[0] = 0;
    for (K v = 0; v < g.N; ++v)
        csr.offsets[v+1] = csr.offsets[v] + K(g.deg[v]);
    csr.n_arcs = csr.offsets[g.N];

    csr.values.resize(csr.n_arcs);
    for (K v = 0; v < g.N; ++v) {
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
OwnedDenseGraph<K> thaw(const PlanarCSR<Owned, K>& csr, int dmax) {
    OwnedDenseGraph<K> g(csr.N, dmax);
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
OwnedDenseGraph<K> dense_from_neighbours(int N_verts,
                                     const std::vector<std::vector<int>>& neighbours,
                                     int dmax) {
    OwnedDenseGraph<K> g(K(N_verts), dmax);
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
    std::vector<std::vector<int>> adj(g.N);
    for (K v = 0; v < g.N; ++v) {
        adj[v].resize(g.deg[v]);
        for (int i = 0; i < g.deg[v]; ++i)
            adj[v][i] = int(g.values[v * g.dmax + i]);
    }
    return adj;
}

} // namespace Spanify
