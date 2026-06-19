#pragma once

#include <vector>
#include <span>
#include <array>
#include <tuple>
#include <cstddef>
#include <cstdint>
#include <cassert>
#include <algorithm>
#include <iostream>
#include <initializer_list>
#include <type_traits>
#include <utility>

#include "planar_csr.hh"
#include "fullerenes/batch/batch.hh"

namespace Spanify {

template<typename K> struct OwnedDenseGraph;  // forward declaration

// ---------------------------------------------------------------------------
// RSRAdjacencyView: non-owning view over fixed-stride flat adjacency data.
//
// Each vertex has dmax slots; a degree counter tracks how many are filled.
// Supports O(dmax) insert/remove at arbitrary cyclic positions.
// Optionally carries a twin array for O(1) reverse-arc lookup.
//
// operator[] returns a span over the row's active entries.
//
// Template parameter:
//   K -- index type: int32_t (CPU), uint16_t (GPU)
// ---------------------------------------------------------------------------
template<typename K = int32_t>
struct RSRAdjacencyView {
    using node_type = K;
    using arc_t   = std::pair<K, K>;       // {source, target}
    using arcix_t = std::pair<K, uint8_t>; // {source, position_in_row}

    static constexpr uint8_t default_dmax = 3;  // cubic

    K N = 0;
    int dmax = 0;
    std::span<K> neighbours;     // [N * dmax] -- neighbor IDs, row-major
    std::span<uint8_t> deg;      // [N] -- current degree of each vertex
    std::span<uint8_t> twin;     // [N * dmax] -- twin position per arc (optional)
    //   twin[u * dmax + i] = j  means: neighbours[v][j] == u,
    //   where v = neighbours[u][i]. Empty when not computed.

    // --- Constructors ---

    RSRAdjacencyView() = default;

    // View constructor: wrap existing spans.
    RSRAdjacencyView(K N, int dmax, std::span<K> neighbours, std::span<uint8_t> deg,
                     std::span<uint8_t> twin = {})
        : N(N), dmax(dmax), neighbours(neighbours), deg(deg), twin(twin) {}

    // -----------------------------------------------------------------------
    // Batchability contract (see include/fullerenes/batch/batchable.hh).
    //
    // Tuple of span-typed fields in canonical order: {neighbours, deg, twin}.
    // Per-field size factor: element count per vertex = factor, so the total
    // number of elements in field k for one batch entry with N vertices is
    // N * size_factor[k].
    //   neighbours : dmax per vertex  (N * dmax total)
    //   deg        : 1    per vertex  (N      total)
    //   twin       : dmax per vertex  (N * dmax total; empty until computed)
    //
    // Derived views that add span fields (e.g. PolyhedronView adds `points`)
    // override n_fields / to_tuple / get_size_factors to extend the tuple.
    // -----------------------------------------------------------------------
    static constexpr std::size_t n_fields = 3;

    auto to_tuple() {
        return std::forward_as_tuple(neighbours, deg, twin);
    }
    auto to_tuple() const {
        return std::forward_as_tuple(neighbours, deg, twin);
    }

    static constexpr std::array<std::size_t, n_fields>
    get_size_factors(int /*N*/, int dmax) {
        return { (std::size_t)dmax, (std::size_t)1, (std::size_t)dmax };
    }

    // --- operator[] returns span over active entries ---

    std::span<K> operator[](K u) {
        return {neighbours.data() + u * dmax, (size_t)deg[u]};
    }
    std::span<const K> operator[](K u) const {
        return {neighbours.data() + u * dmax, (size_t)deg[u]};
    }

    // --- Span-based accessors ---

    std::span<const K> nbrs(K u) const {
        return {neighbours.data() + u * dmax, (size_t)deg[u]};
    }

    std::span<K> nbrs_mut(K u) {
        return {neighbours.data() + u * dmax, (size_t)deg[u]};
    }

    int degree(K u) const { return deg[u]; }

    // --- Twin query ---
    bool has_twin() const { return !twin.empty(); }

    // --- Vertex-pair navigation (O(degree) per call) ---

    int find(K u, K v) const {
        const K* row = neighbours.data() + u * dmax;
        for (int i = 0; i < deg[u]; ++i)
            if (row[i] == v) return i;
        return -1;
    }

    // --- Arc-index navigation (O(1) per call, except find_arc which is O(degree)) ---

    arcix_t find_arc(K u, K v) const {
        return {u, uint8_t(find(u, v))};
    }

    // Flat arc id in [0, N*dmax): the dense canonical index of the directed arc at
    // position i in u's row (the arc u -> neighbours[u*dmax+i]). It indexes both
    // `neighbours` and any arc-attribute vector<T>(N*dmax) identically, so an arc-flag
    // array replaces a hash map keyed by arc_t -- denser, allocation-free, and parallel
    // arcs / self-loops stay distinct (an (u,v) key conflates them). The arc VALUE is
    // recoverable O(1): source = id/dmax (arc_of), target = neighbours[id]. size_t so
    // the index can never overflow before the (size_t-sized) backing buffer does.
    size_t  arcid(K u, int i) const { return size_t(u) * dmax + i; }
    size_t  arcid(arcix_t a)  const { return arcid(a.first, a.second); }
    arcix_t arc_of(size_t id) const { return {K(id / dmax), uint8_t(id % dmax)}; }

    K target(arcix_t a) const {
        return neighbours[a.first * dmax + a.second];
    }

    arcix_t next(arcix_t a) const {
        auto [u, i] = a;
        uint8_t d = deg[u];
        return {u, uint8_t((i + 1) % d)};
    }

    arcix_t prev(arcix_t a) const {
        auto [u, i] = a;
        uint8_t d = deg[u];
        return {u, uint8_t((i + d - 1) % d)};
    }

    arcix_t reverse_arc(arcix_t a) const {
        assert(has_twin());
        auto [u, i] = a;
        K v = neighbours[u * dmax + i];
        uint8_t j = twin[u * dmax + i];
        return {v, j};
    }

    arcix_t next_on_face(arcix_t a) const {
        return next(reverse_arc(a));
    }

    // --- Per-row mutation (no reallocation) ---

    void push_back(K u, K v) {
        if (deg[u] >= dmax) {
            std::cerr << "RSRAdjacencyView::push_back: degree overflow at vertex " << u
                      << " (deg=" << int(deg[u]) << ", dmax=" << dmax << ")\n";
            std::abort();
        }
        neighbours[u * dmax + deg[u]] = v;
        deg[u]++;
    }

    void insert_at(K u, K v, int pos) {
        if (deg[u] >= dmax) {
            std::cerr << "RSRAdjacencyView::insert_at: degree overflow at vertex " << u
                      << " (deg=" << int(deg[u]) << ", dmax=" << dmax << ")\n";
            std::abort();
        }
        K* row = neighbours.data() + u * dmax;
        int d = deg[u];
        for (int i = d; i > pos; --i)
            row[i] = row[i-1];
        row[pos] = v;
        deg[u]++;
    }

    void erase_at(K u, int pos) {
        K* row = neighbours.data() + u * dmax;
        int d = deg[u];
        for (int i = pos; i < d - 1; ++i)
            row[i] = row[i+1];
        row[d-1] = K(-1);
        deg[u]--;
    }

    void clear_row(K u) { deg[u] = 0; }

    void resize_row(K u, int n) {
        assert(n <= dmax);
        for (int i = deg[u]; i < n; ++i)
            neighbours[u * dmax + i] = K(-1);
        deg[u] = uint8_t(n);
    }

    void assign_row(K u, std::initializer_list<K> v) {
        assert(int(v.size()) <= dmax);
        deg[u] = uint8_t(v.size());
        int i = 0;
        for (K x : v) neighbours[u * dmax + i++] = x;
    }

    void assign_row(K u, const std::vector<K>& v) {
        assert(int(v.size()) <= dmax);
        deg[u] = uint8_t(v.size());
        for (int i = 0; i < int(v.size()); ++i)
            neighbours[u * dmax + i] = v[i];
    }

    void assign_row(K u, std::span<const K> v) {
        assert(int(v.size()) <= dmax);
        deg[u] = uint8_t(v.size());
        for (int i = 0; i < int(v.size()); ++i)
            neighbours[u * dmax + i] = v[i];
    }

    // --- Backward compat with vector<vector<>> interface ---

    int size() const { return N; }

    bool operator==(const RSRAdjacencyView& other) const {
        if (N != other.N || dmax != other.dmax) return false;
        for (K u = 0; u < N; ++u) {
            auto a = (*this)[u], b = other[u];
            if (a.size() != b.size()) return false;
            if (!std::equal(a.begin(), a.end(), b.begin())) return false;
        }
        return true;
    }

    std::vector<std::vector<K>> to_vectors() const {
        std::vector<std::vector<K>> adj(N);
        for (K v = 0; v < N; ++v) {
            adj[v].resize(deg[v]);
            for (int i = 0; i < deg[v]; ++i)
                adj[v][i] = neighbours[v * dmax + i];
        }
        return adj;
    }
};

// Backward-compat alias
template<typename K = int32_t>
using DenseGraph = RSRAdjacencyView<K>;

// TC verification: RSRAdjacencyView with default K is trivially copyable.
static_assert(std::is_trivially_copyable_v<RSRAdjacencyView<int32_t>>,
    "RSRAdjacencyView<int32_t> must be trivially copyable");
static_assert(std::is_trivially_copyable_v<RSRAdjacencyView<uint16_t>>,
    "RSRAdjacencyView<uint16_t> must be trivially copyable");

// ---------------------------------------------------------------------------
// RSRPolyhedronView<T,K>: RSR adjacency extended with per-vertex 3D coordinates.
//
// Carries the same three adjacency fields (neighbours, deg, twin) as
// RSRAdjacencyView plus a points span of std::array<T,3> (xyz).  This makes
// (graph + geometry) atomically transferable through Batch<V>/BatchQueue<V>.
//
// Canonical tuple order: {neighbours, deg, twin, points}.  The first three
// fields' size_factors match RSRAdjacencyView, so a BatchView<RSRPolyhedronView>
// can be sliced into a BatchView<RSRAdjacencyView> via as_adjacency_view().
// ---------------------------------------------------------------------------
template<typename T, typename K = int32_t>
struct RSRPolyhedronView : RSRAdjacencyView<K> {
    using typename RSRAdjacencyView<K>::node_type;

    static constexpr uint8_t default_dmax = 3;  // cubic

    std::span<std::array<T,3>> points;          // [N] -- per-vertex xyz

    RSRPolyhedronView() = default;

    RSRPolyhedronView(K N, int dmax,
                      std::span<K> neighbours, std::span<uint8_t> deg,
                      std::span<std::array<T,3>> pts,
                      std::span<uint8_t> twin = {})
        : RSRAdjacencyView<K>(N, dmax, neighbours, deg, twin), points(pts) {}

    // Construct from an existing adjacency view + points span.
    RSRPolyhedronView(const RSRAdjacencyView<K>& g,
                      std::span<std::array<T,3>> pts)
        : RSRAdjacencyView<K>(g), points(pts) {}

    // -- Batchability contract (extends adjacency with `points`) ------------
    static constexpr std::size_t n_fields = 4;

    auto to_tuple() {
        return std::forward_as_tuple(this->neighbours, this->deg, this->twin, points);
    }
    auto to_tuple() const {
        return std::forward_as_tuple(this->neighbours, this->deg, this->twin, points);
    }

    static constexpr std::array<std::size_t, n_fields>
    get_size_factors(int /*N*/, int dmax) {
        return { (std::size_t)dmax, (std::size_t)1, (std::size_t)dmax, (std::size_t)1 };
    }
};

static_assert(std::is_trivially_copyable_v<RSRPolyhedronView<float,  uint16_t>>,
    "RSRPolyhedronView<float,uint16_t> must be trivially copyable");
static_assert(std::is_trivially_copyable_v<RSRPolyhedronView<double, uint16_t>>,
    "RSRPolyhedronView<double,uint16_t> must be trivially copyable");

// ---------------------------------------------------------------------------
// Graph-prefix slicing: view a BatchView<RSRPolyhedronView<T,K>> as a
// BatchView<RSRAdjacencyView<K>> over the shared adjacency fields.
// The points span is available via points_span().
// ---------------------------------------------------------------------------
template<typename T, typename K>
batch::BatchView<RSRAdjacencyView<K>>
as_adjacency_view(batch::BatchView<RSRPolyhedronView<T,K>> pv) {
    const auto& s = pv.spans();
    typename batch::detail::span_tuple_t<RSRAdjacencyView<K>> adj_spans{
        std::get<0>(s), std::get<1>(s), std::get<2>(s)
    };
    return batch::BatchView<RSRAdjacencyView<K>>(pv.N(), pv.dmax(), pv.size(), adj_spans);
}

template<typename T, typename K>
std::span<std::array<T,3>>
points_span(batch::BatchView<RSRPolyhedronView<T,K>> pv) {
    return std::get<3>(pv.spans());
}

// ---------------------------------------------------------------------------
// OwnedDenseGraph: owning version with vector storage.
// Inherits all view/mutation methods from RSRAdjacencyView.
// ---------------------------------------------------------------------------
template<typename K = int32_t>
struct OwnedDenseGraph : RSRAdjacencyView<K> {
    std::vector<K> owned_neighbours;
    std::vector<uint8_t> owned_deg;

    void repoint() {
        this->neighbours = std::span<K>(owned_neighbours.data(), owned_neighbours.size());
        this->deg = std::span<uint8_t>(owned_deg.data(), owned_deg.size());
    }

    // Bring base push_back(K,K) into scope (otherwise hidden by push_back(vector))
    using RSRAdjacencyView<K>::push_back;

    // --- Constructors ---

    OwnedDenseGraph() = default;

    // Explicit: N vertices, given stride, all empty.
    explicit OwnedDenseGraph(K N, int dmax)
        : owned_neighbours(N * dmax, K(-1)), owned_deg(N, 0) {
        this->N = N; this->dmax = dmax; repoint();
    }

    // Fill constructor: N vertices, each initialized with the given row.
    OwnedDenseGraph(K N, const std::vector<K>& initial_row)
        : owned_neighbours(N * int(initial_row.size()), K(-1)),
          owned_deg(N, uint8_t(initial_row.size())) {
        this->N = N; this->dmax = int(initial_row.size()); repoint();
        for (K v = 0; v < N; ++v)
            for (int i = 0; i < this->dmax; ++i)
                this->neighbours[v * this->dmax + i] = initial_row[i];
    }

    // Brace-init constructor: OwnedDenseGraph{{1,2,3},{4,5,6},...}.
    OwnedDenseGraph(std::initializer_list<std::initializer_list<K>> adj) {
        this->N = K(adj.size());
        this->dmax = 0;
        for (auto& row : adj) this->dmax = std::max(this->dmax, int(row.size()));
        owned_neighbours.resize(this->N * this->dmax, K(-1));
        owned_deg.resize(this->N, 0);
        repoint();
        int v = 0;
        for (auto& row : adj) {
            this->deg[v] = uint8_t(row.size());
            int i = 0;
            for (K x : row) this->neighbours[v * this->dmax + i++] = x;
            v++;
        }
    }

    // Converting constructor from vector<vector<K>>.
    OwnedDenseGraph(const std::vector<std::vector<K>>& adj) {
        this->N = K(adj.size());
        this->dmax = 0;
        for (auto& row : adj) this->dmax = std::max(this->dmax, int(row.size()));
        owned_neighbours.resize(this->N * this->dmax, K(-1));
        owned_deg.resize(this->N, 0);
        repoint();
        for (int v = 0; v < int(adj.size()); ++v) {
            this->deg[v] = uint8_t(adj[v].size());
            for (int i = 0; i < int(adj[v].size()); ++i)
                this->neighbours[v * this->dmax + i] = adj[v][i];
        }
    }

    // Converting constructor from RSRAdjacencyView (copies data).
    OwnedDenseGraph(const RSRAdjacencyView<K>& v) {
        this->N = v.N; this->dmax = v.dmax;
        if (v.N > 0 && v.neighbours.data()) {
            owned_neighbours.assign(v.neighbours.begin(), v.neighbours.end());
            owned_deg.assign(v.deg.begin(), v.deg.end());
        }
        repoint();
    }

    // --- Rule of 5 ---

    OwnedDenseGraph(const OwnedDenseGraph& o)
        : owned_neighbours(o.owned_neighbours), owned_deg(o.owned_deg) {
        this->N = o.N; this->dmax = o.dmax; repoint();
    }

    OwnedDenseGraph(OwnedDenseGraph&& o) noexcept
        : owned_neighbours(std::move(o.owned_neighbours)), owned_deg(std::move(o.owned_deg)) {
        this->N = o.N; this->dmax = o.dmax; repoint();
        o.neighbours = {}; o.deg = {}; o.N = 0;
    }

    OwnedDenseGraph& operator=(const OwnedDenseGraph& o) {
        if (this != &o) {
            this->N = o.N; this->dmax = o.dmax;
            owned_neighbours = o.owned_neighbours;
            owned_deg = o.owned_deg;
            repoint();
        }
        return *this;
    }

    OwnedDenseGraph& operator=(OwnedDenseGraph&& o) noexcept {
        this->N = o.N; this->dmax = o.dmax;
        owned_neighbours = std::move(o.owned_neighbours);
        owned_deg = std::move(o.owned_deg);
        repoint();
        o.neighbours = {}; o.deg = {}; o.N = 0;
        return *this;
    }

    // --- Resize operations (require reallocation) ---

    void resize(int new_N) {
        owned_neighbours.resize(new_N * this->dmax, K(-1));
        owned_deg.resize(new_N, 0);
        this->N = K(new_N);
        repoint();
    }

    void push_back(const std::vector<K>& row) {
        assert(int(row.size()) <= this->dmax);
        owned_neighbours.resize((this->N + 1) * this->dmax, K(-1));
        owned_deg.push_back(uint8_t(row.size()));
        repoint();
        for (int i = 0; i < int(row.size()); ++i)
            this->neighbours[this->N * this->dmax + i] = row[i];
        this->N++;
    }

    void pop_back() {
        assert(this->N > 0);
        this->N--;
        owned_neighbours.resize(this->N * this->dmax);
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
                g.neighbours[v * new_dmax + i] = this->neighbours[v * this->dmax + i];
        }
        return g;
    }
};

// --- ostream support ---

template<typename K>
std::ostream& operator<<(std::ostream& s, const RSRAdjacencyView<K>& g) {
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
// Free-function restride: works on any RSRAdjacencyView. O(N*dmax).
// ---------------------------------------------------------------------------
template<typename K>
OwnedDenseGraph<K> restride(const RSRAdjacencyView<K>& src, int new_dmax) {
    OwnedDenseGraph<K> g(src.N, new_dmax);
    for (K v = 0; v < src.N; ++v) {
        assert(src.deg[v] <= new_dmax);
        g.deg[v] = src.deg[v];
        for (int i = 0; i < src.deg[v]; ++i)
            g.neighbours[v * new_dmax + i] = src.neighbours[v * src.dmax + i];
    }
    return g;
}

// ---------------------------------------------------------------------------
// Conversion: RSR -> CSR (freeze). O(N).
// ---------------------------------------------------------------------------
template<typename K>
PlanarCSR<Owned, K> freeze(const RSRAdjacencyView<K>& g) {
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
            csr.values[base + i] = g.neighbours[v * g.dmax + i];
    }

    csr.twin = compute_twin(csr.N, csr.offsets, csr.values);
    return csr;
}

// ---------------------------------------------------------------------------
// Conversion: CSR -> RSR (thaw). O(N).
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
            g.neighbours[v * dmax + i] = csr.values[csr.offsets[v] + i];
    }
    return g;
}

// ---------------------------------------------------------------------------
// Conversion: vector<vector<int>> -> RSR. O(N).
// ---------------------------------------------------------------------------
template<typename K = int32_t>
OwnedDenseGraph<K> dense_from_neighbours(int N_verts,
                                     const std::vector<std::vector<int>>& nbrs,
                                     int dmax) {
    OwnedDenseGraph<K> g(K(N_verts), dmax);
    for (int v = 0; v < N_verts; ++v) {
        assert(int(nbrs[v].size()) <= dmax);
        g.deg[v] = uint8_t(nbrs[v].size());
        for (int i = 0; i < int(nbrs[v].size()); ++i)
            g.neighbours[v * dmax + i] = K(nbrs[v][i]);
    }
    return g;
}

// ---------------------------------------------------------------------------
// Conversion: RSR -> vector<vector<int>>. O(N).
// ---------------------------------------------------------------------------
template<typename K>
std::vector<std::vector<int>> to_neighbours(const RSRAdjacencyView<K>& g) {
    std::vector<std::vector<int>> adj(g.N);
    for (K v = 0; v < g.N; ++v) {
        adj[v].resize(g.deg[v]);
        for (int i = 0; i < g.deg[v]; ++i)
            adj[v][i] = int(g.neighbours[v * g.dmax + i]);
    }
    return adj;
}

} // namespace Spanify
