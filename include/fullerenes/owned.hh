#pragma once

// Owned<View>: universal ownership template for the view hierarchy.
//
// Adds vector storage for any view type (GraphView, PolyhedronView, etc.).
// Points the base view's spans into its own vectors.
//
// Detects at compile time whether the view carries geometry (a `points`
// span) and conditionally owns coordinate storage. One template, one
// Rule-of-5, handles all owned types.

#include "fullerenes/graphview.hh"
#include <variant>
#include <type_traits>

template<typename View>
struct Owned : View {
    using node = typename View::node_type;

    // Compile-time detection: does View have a `points` member?
    static constexpr bool has_geometry = requires(View v) { v.points; };

    // --- Storage ---
    std::vector<node>    owned_neighbours;
    std::vector<uint8_t> owned_deg;
    std::vector<uint8_t> owned_twin;

    // For geometry views: vector<coord3d>.  For non-geometry: monostate (zero overhead).
    std::conditional_t<has_geometry, std::vector<coord3d>,
                       std::monostate> owned_points{};

    // --- Repoint all spans to owned storage ---
    void repoint() {
        this->neighbours = std::span<node>(owned_neighbours);
        this->deg = std::span<uint8_t>(owned_deg);
        this->twin = owned_twin.empty()
            ? std::span<uint8_t>()
            : std::span<uint8_t>(owned_twin);
        if constexpr (has_geometry)
            this->points = std::span<coord3d>(owned_points);
    }

    bool owns_memory() const { return !owned_neighbours.empty(); }

    // --- Allocating constructor ---
    // default_dmax comes from the View type:
    //   Owned<GraphView>(60)          -> dmax=10
    //   Owned<CubicGraphView>(60)     -> dmax=3
    //   Owned<TriangulationView>(32)  -> dmax=6
    //   Owned<PolyhedronView>(60)     -> dmax=10, points[60]
    explicit Owned(int N, uint8_t dmax = View::default_dmax)
        : owned_neighbours(N * dmax, node(-1)), owned_deg(N, 0) {
        this->N = node(N);
        this->dmax = dmax;
        if constexpr (has_geometry) owned_points.resize(N);
        repoint();
    }

    // --- Deep copy from any compatible view (adjacency only) ---
    // For geometry types, allocates default-initialized coordinates.
    explicit Owned(const GraphView& src)
        : owned_neighbours(src.neighbours.begin(), src.neighbours.end()),
          owned_deg(src.deg.begin(), src.deg.end()),
          owned_twin(src.twin.begin(), src.twin.end()) {
        this->N = src.N;
        this->dmax = src.dmax;
        if constexpr (has_geometry) owned_points.resize(src.N);
        repoint();
    }

    // --- Deep copy adjacency + move coordinates (geometry types only) ---
    Owned(const GraphView& src, std::vector<coord3d> pts) requires has_geometry
        : owned_neighbours(src.neighbours.begin(), src.neighbours.end()),
          owned_deg(src.deg.begin(), src.deg.end()),
          owned_twin(src.twin.begin(), src.twin.end()),
          owned_points(std::move(pts)) {
        this->N = src.N;
        this->dmax = src.dmax;
        repoint();
    }

    // --- Deep copy from geometry view (adjacency + coordinates) ---
    explicit Owned(const View& src) requires has_geometry
        : owned_neighbours(src.neighbours.begin(), src.neighbours.end()),
          owned_deg(src.deg.begin(), src.deg.end()),
          owned_twin(src.twin.begin(), src.twin.end()),
          owned_points(src.points.begin(), src.points.end()) {
        this->N = src.N;
        this->dmax = src.dmax;
        repoint();
    }

    // --- Default + Rule-of-5 ---
    Owned() = default;

    Owned(const Owned& o)
        : View(o),
          owned_neighbours(o.owned_neighbours),
          owned_deg(o.owned_deg),
          owned_twin(o.owned_twin),
          owned_points(o.owned_points) {
        if (!owned_neighbours.empty()) repoint();
    }

    Owned(Owned&& o) noexcept
        : View(o),
          owned_neighbours(std::move(o.owned_neighbours)),
          owned_deg(std::move(o.owned_deg)),
          owned_twin(std::move(o.owned_twin)),
          owned_points(std::move(o.owned_points)) {
        if (!owned_neighbours.empty()) repoint();
        o.neighbours = {}; o.deg = {}; o.N = 0;
        if constexpr (has_geometry) o.points = {};
    }

    Owned& operator=(const Owned& o) {
        if (this != &o) {
            View::operator=(o);
            owned_neighbours = o.owned_neighbours;
            owned_deg = o.owned_deg;
            owned_twin = o.owned_twin;
            owned_points = o.owned_points;
            if (!owned_neighbours.empty()) repoint();
            else {
                this->neighbours = o.neighbours;
                this->deg = o.deg;
                this->twin = o.twin;
                if constexpr (has_geometry) this->points = o.points;
            }
        }
        return *this;
    }

    Owned& operator=(Owned&& o) noexcept {
        if (this != &o) {
            View::operator=(o);
            owned_neighbours = std::move(o.owned_neighbours);
            owned_deg = std::move(o.owned_deg);
            owned_twin = std::move(o.owned_twin);
            owned_points = std::move(o.owned_points);
            if (!owned_neighbours.empty()) repoint();
            else {
                this->neighbours = o.neighbours;
                this->deg = o.deg;
                this->twin = o.twin;
                if constexpr (has_geometry) this->points = o.points;
            }
            o.neighbours = {}; o.deg = {}; o.N = 0;
            if constexpr (has_geometry) o.points = {};
        }
        return *this;
    }

    // --- Assignment from geometry view (deep copy adjacency + coordinates) ---
    Owned& operator=(const View& src) requires has_geometry {
        if (src.neighbours.data() == owned_neighbours.data()) return *this;
        this->N = src.N;
        this->dmax = src.dmax;
        owned_neighbours.assign(src.neighbours.begin(), src.neighbours.end());
        owned_deg.assign(src.deg.begin(), src.deg.end());
        owned_twin.assign(src.twin.begin(), src.twin.end());
        owned_points.assign(src.points.begin(), src.points.end());
        repoint();
        return *this;
    }

    // --- Assignment from any view (deep copy adjacency only) ---
    Owned& operator=(const GraphView& src) {
        if (src.neighbours.data() == owned_neighbours.data()) return *this;
        this->N = src.N;
        this->dmax = src.dmax;
        owned_neighbours.assign(src.neighbours.begin(), src.neighbours.end());
        owned_deg.assign(src.deg.begin(), src.deg.end());
        owned_twin.assign(src.twin.begin(), src.twin.end());
        if constexpr (has_geometry) owned_points.resize(src.N);
        repoint();
        return *this;
    }

    // --- Twin computation (allocates owned_twin) ---
    void compute_twin() {
        owned_twin.resize(this->N * this->dmax, 0);
        this->twin = std::span<uint8_t>(owned_twin);
        for (node u = 0; u < this->N; ++u) {
            for (int i = 0; i < this->deg[u]; ++i) {
                node v = this->neighbours[u * this->dmax + i];
                int j = this->find(v, u);
                assert(j >= 0);
                owned_twin[u * this->dmax + i] = uint8_t(j);
            }
        }
    }

    // --- Restride (owned only -- reallocates) ---
    void restride_inplace(uint8_t new_dmax) {
        std::vector<node> new_neighbours(this->N * new_dmax, node(-1));
        std::vector<uint8_t> new_deg(this->N);
        for (int u = 0; u < this->N; ++u) {
            int d = std::min((int)this->degree(u), (int)new_dmax);
            new_deg[u] = d;
            for (int j = 0; j < d; ++j)
                new_neighbours[u * new_dmax + j] = (*this)[u][j];
        }
        owned_neighbours = std::move(new_neighbours);
        owned_deg = std::move(new_deg);
        owned_twin.clear();
        this->dmax = new_dmax;
        repoint();
    }

    // --- Resize (owned only -- reallocates) ---
    void resize(int new_N) {
        owned_neighbours.resize(new_N * this->dmax, node(-1));
        owned_deg.resize(new_N, 0);
        this->N = node(new_N);
        if constexpr (has_geometry) owned_points.resize(new_N);
        repoint();
    }

    // --- Push/pop vertex rows (owned only) ---
    void push_back(const std::vector<node>& row) {
        assert(int(row.size()) <= this->dmax);
        owned_neighbours.resize((this->N + 1) * this->dmax, node(-1));
        owned_deg.push_back(uint8_t(row.size()));
        repoint();
        for (int i = 0; i < int(row.size()); ++i)
            this->neighbours[this->N * this->dmax + i] = row[i];
        this->N++;
        if constexpr (has_geometry) {
            owned_points.push_back(coord3d());
            this->points = std::span<coord3d>(owned_points);
        }
    }

    void pop_back() {
        assert(this->N > 0);
        this->N--;
        owned_neighbours.resize(this->N * this->dmax);
        owned_deg.resize(this->N);
        if constexpr (has_geometry) {
            owned_points.resize(this->N);
        }
        repoint();
    }

    // Bring base push_back(node, node) into scope
    using View::push_back;

    // --- Remove isolated vertices (compacts graph) ---
    void remove_isolated_vertices() {
        std::vector<int> new_id(this->N);
        int u_new = 0;
        for (int u = 0; u < this->N; u++)
            if (!(*this)[u].empty())
                new_id[u] = u_new++;

        Owned g(u_new, this->dmax);
        for (int u = 0; u < this->N; u++)
            for (node v : (*this)[u])
                g.push_back(node(new_id[u]), node(new_id[v]));

        *this = std::move(g);
    }

    // --- Remove specific vertices and compact ---
    void remove_vertices(std::set<int>& sv) {
        const int N_naught = this->N;
        for (int u : sv) {
            while (!(*this)[u].empty()) {
                node v = (*this)[u][0];
                this->remove_edge({node(u), v});
            }
        }
        remove_isolated_vertices();
        if (N_naught != int(sv.size()) + this->N)
            std::cerr << "removed more vertices than intended" << std::endl;
        assert(this->is_connected());
    }
};
