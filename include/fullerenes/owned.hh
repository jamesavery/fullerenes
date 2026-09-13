#pragma once

// Owned<View>: the storage behind a view, and nothing else.
//
// A view (graphview.hh) is a bundle of spans over fields it does not own.
// This owner allocates ONE std::vector per field of the view's batchability
// contract -- to_tuple / n_fields / get_element_counts, the contract
// Batch<V> slices (batch/batchable.hh) -- sized for a vertex CAPACITY that
// may exceed the live vertex count N, keeps the view's spans pointing at
// the WHOLE storage (capacity vertices' worth of each buffer), frees it
// when it dies, and hands out its view.  A reader of a view is bounded by
// N, so the rows [N, capacity()) are storage and not vertices of the
// graph; the view's capacity() is the rows its spans hold (dense_graph.hh)
// and here that is the storage's size -- one number, derived, not a field
// (Batch<V>'s size and capacity are the same pair one dimension over).  So:
//
//   - every field of every view rides one law: adjacency, degrees, twin,
//     coordinates, the 12-entry pentagon list -- no per-field members and no
//     "does this view carry geometry" branches; a view that adds a field to
//     its tuple is owned without a line here changing;
//   - a graph edited in place for a whole run (an enumerator's working
//     graph) is a VIEW copied from an owner sized for the run's bound: it
//     moves its own N within capacity(), grows nothing and re-forms no
//     span, and the owner's storage outlives it; within capacity the
//     owner's own resize changes its count and pads the rows it exposes,
//     and allocates nothing;
//   - the owner IS the view (it derives from it), so a reader takes
//     `const View&` and an algorithm is written once against the view.
//
// What lives here: allocation and capacity, repointing, the deep copies
// (batch::copy_entry into freshly sized storage) and the row-count changes
// those imply, and the owner's half of the twin table (its storage; the
// derivation is the view's compute_twin).  The whole-graph relayouts --
// restride, relabel, compaction -- are the view's words (dense_graph.hh:
// pad_rows, restride_into, relabel_into) composed here with one allocation
// each.
//
// THE TWIN IS OPTIONAL and its buffer is not: the table is allocated with
// the other fields but the view's span is switched on only by
// compute_twin(), so has_twin() means "computed" exactly as it does for a
// view over caller arrays, and a graph that never asks for a table never
// pays its maintenance.
//
// STORAGE BACKEND: std::vector on host.  (Batch<View> uses BatchAlloc<T>,
// which aliases to SyclVector<T> in SYCL builds; the owner stays host-vector
// so the SYCL translation units that include it compile unchanged.)

#include "fullerenes/graphview.hh"
#include "fullerenes/batch/batchable.hh"

#include <algorithm>
#include <set>
#include <string>
#include <tuple>
#include <type_traits>
#include <utility>
#include <vector>

// An owning graph as a filler sees it (BuckyGen::dual_slot): storage its
// view spans point into, reshaped on request to an empty N-by-dmax graph
// with any twin table dropped, exposing the view's fields the fill writes.
// Owned<View> and the legacy Graph both satisfy it, so a filler is written
// once against any owner rather than once per concrete graph class.
template<typename G>
concept owning_graph = requires(G& g, typename G::node_type n, int d) {
    g.reshape(n, d);
    g.repoint();
    g.N;
    g.dmax;
    g.neighbours;
    g.deg;
};

namespace owned_detail {
    // tuple<vector<T0>, vector<T1>, ...>: one host buffer per field of V.
    template<class V, std::size_t... Is>
    auto make_buffer_tuple(std::index_sequence<Is...>)
        -> std::tuple<std::vector<std::remove_const_t<batch::field_element_t<V, Is>>>...>;
    template<class V>
    using buffer_tuple_t =
        decltype(make_buffer_tuple<V>(std::make_index_sequence<V::n_fields>{}));
}

template<typename View>
struct Owned : View {
    using node      = typename View::node_type;
    using buffers_t = owned_detail::buffer_tuple_t<View>;

    // The graph triple {neighbours, deg, twin} is the base adjacency
    // contract's field set, in its canonical order (batchable.hh); the twin
    // is its last field, and the one whose span is switched on by
    // computation rather than by allocation.  Fields past the triple are the
    // view's own (coordinates, the pentagon list).
    static constexpr std::size_t n_graph_fields = Spanify::RSRAdjacencyView<node>::n_fields;
    static constexpr std::size_t twin_field     = n_graph_fields - 1;
    static_assert(std::is_same_v<batch::field_element_t<View, twin_field>, uint8_t>,
                  "the last field of the graph triple must be the twin table");
    // The degree field holds one element per vertex, so its buffer's length
    // is the vertex capacity the buffers are sized for: the one number,
    // read off the storage before the spans are formed (storage_capacity)
    // and off the view's capacity() -- deg.size() -- after.
    static constexpr std::size_t deg_field = 1;
    static_assert(std::is_same_v<batch::field_element_t<View, deg_field>, uint8_t>,
                  "the middle field of the graph triple must be the degree table");

    buffers_t buffers;            // one vector per field, holding capacity() vertices' worth
    bool      twin_computed = false;

    // The k-th field's buffer -- what a deep copy or a test reads; the view's
    // span covers all of it, and the first N vertices' worth is the graph.
    template<std::size_t k> auto&       buffer()       { return std::get<k>(buffers); }
    template<std::size_t k> const auto& buffer() const { return std::get<k>(buffers); }

    bool owns_memory() const { return storage_capacity() > 0; }

    // --- Repoint: every span over the whole storage -- capacity vertices'
    //     worth of its buffer -- the twin span at nothing until computed.
    //     Never past a buffer's end: that clamp is what makes an owner
    //     without storage total (a moved-from one, or a default-constructed
    //     one before reserve) -- its spans are empty. ---
    void repoint() {
        const auto counts = View::get_element_counts(storage_capacity(), this->dmax);
        auto spans = this->to_tuple();
        batch::for_each_field<View>([&](auto Ic) {
            constexpr std::size_t k = Ic;
            using span_t = batch::field_span_t<View, k>;
            auto& buf = std::get<k>(buffers);
            const bool on = (k != twin_field) || twin_computed;
            std::get<k>(spans) = on ? span_t(buf.data(), std::min(counts[k], buf.size())) : span_t{};
        });
    }

    // --- Capacity: size every buffer for at least `cap` vertices at the
    //     current stride, keeping contents, and repoint.  Buffers only ever
    //     grow. ---
    // @anchor owned-reserve
    // @post capacity() >= cap
    // @post sized: every buffer holds at least get_element_counts(capacity(), dmax) elements
    // @post kept:  the first N vertices' worth of every field is unchanged
    void reserve(int cap) {
        grow_buffers(cap);
        repoint();
    }

    // --- Live vertex count.  Within capacity nothing is allocated and no
    //     buffer moves.  Rows exposed by growth are padded (the view's
    //     pad_rows), so they are empty and a computed twin stays valid (an
    //     empty row has no arcs); the view's own fields are value-initialised
    //     over the exposed slice.  Shrinking abandons rows; they are padded
    //     again if ever re-exposed by resize. ---
    // @anchor owned-resize
    // @post N == new_N && capacity() >= new_N
    // @post empty:   all_of(indices(old_N, new_N), [&](node u){ return degree(u) == 0; })
    // @post padded:  all_of(indices(old_N, new_N), [&](node u){ return row_is_padded(u); })
    // @post twin:    implies(has_twin(), twin_is_valid() == twin_is_valid_before)
    // @post storage: implies(new_N <= capacity_before, every buffer's data() is unchanged)
    void resize(int new_N) {
        const int old_N = this->N;
        grow_buffers(new_N);
        this->N = node(new_N);
        repoint();
        if (new_N > old_N) {
            const auto was = View::get_element_counts(old_N, this->dmax);
            const auto now = View::get_element_counts(new_N, this->dmax);
            batch::for_each_field<View>([&](auto Ic) {
                constexpr std::size_t k = Ic;
                if constexpr (k >= n_graph_fields) {
                    auto& buf = std::get<k>(buffers);
                    using elem_t = typename std::remove_reference_t<decltype(buf)>::value_type;
                    std::fill(buf.begin() + was[k], buf.begin() + now[k], elem_t{});
                }
            });
            this->pad_rows(node(old_N), node(new_N));
        }
    }

    // --- An EMPTY N-by-dmax graph in this storage: the shape a filler wants
    //     (BuckyGen::dual_slot), reusing the buffers when they are large
    //     enough and dropping any twin table (it described the previous
    //     graph).  Every row is padded on every call -- a filler that skips
    //     a row can no longer inherit the previous graph's. ---
    // @anchor owned-reshape
    // @post N == N_ && dmax == dmax_ && !has_twin()
    // @post empty:   all_of(indices(N), [&](node u){ return degree(u) == 0; })
    // @post storage: implies(N_ <= capacity_before && dmax_ <= dmax_before,
    //                        every buffer's data() is unchanged)
    void reshape(node N_, int dmax_) {
        this->N = 0;
        this->dmax = dmax_;
        twin_computed = false;
        resize(int(N_));
    }

    // --- Construction ---

    // No storage and no vertices -- except that the constant-size fields (a
    // dual's pentagon list) get theirs at once, so a view invariant such as
    // "the pentagon span has 12 slots" holds from the first moment.
    Owned() { reserve(0); }

    // An empty graph: N vertices with dmax slots each, sized for
    // max(N, cap) vertices (an enumerator's storage is sized for its run's
    // bound here).  Every row empty, every slot padding.
    explicit Owned(int N, uint8_t dmax = View::default_dmax, int cap = 0) {
        this->dmax = dmax;
        grow_buffers(std::max(N, cap));
        resize(N);
    }

    // Deep copy from any view of the hierarchy (batch::copy_entry): the
    // fields the two contracts share are copied (the source's twin comes
    // along iff it has one), the rest value-initialised -- a FullereneDual
    // derives its pentagon list at its boundary, a Polyhedron sets its
    // points.  The buffers are sized for src.N; a larger existing capacity
    // is kept (copy-assignment from another owner, below, takes that
    // owner's capacity instead).
    template<class Src> requires std::is_base_of_v<GraphView, Src>
    explicit Owned(const Src& src) { assign(src); }

    template<class Src> requires std::is_base_of_v<GraphView, Src>
    Owned& operator=(const Src& src) { assign(src); return *this; }

    // --- Rule of 5: the buffers move or copy, the spans follow. ---

    Owned(const Owned& o)
        : View(o), buffers(o.buffers), twin_computed(o.twin_computed) {
        repoint();
    }
    Owned(Owned&& o) noexcept
        : View(o), buffers(std::move(o.buffers)), twin_computed(o.twin_computed) {
        repoint();
        o.release();
    }
    Owned& operator=(const Owned& o) {
        if (this != &o) {
            View::operator=(o);
            buffers = o.buffers;
            twin_computed = o.twin_computed;
            repoint();
        }
        return *this;
    }
    Owned& operator=(Owned&& o) noexcept {
        if (this != &o) {
            View::operator=(o);
            buffers = std::move(o.buffers);
            twin_computed = o.twin_computed;
            repoint();
            o.release();
        }
        return *this;
    }

    // --- Twin: the owner's half is the storage.  Switch the span on and
    //     derive through the view's compute_twin (the locator, one body for
    //     owner and view alike); an arc the locator could not reverse
    //     (first_twinless_arc) means an asymmetric graph, refused by name
    //     with the table switched off again. ---
    // @anchor owned-compute-twin
    // @pre  symmetric: adjacency_is_symmetric() -- violation throws
    //       graph_surgery_error{AsymmetricAdjacency} naming the arc's
    //       endpoints (this is the from-scratch oracle the surgery tests
    //       compare against, so the guard holds in every build configuration)
    // @post twin:   has_twin() && twin_is_valid()
    // @post atomic: a throwing call leaves !has_twin()
    void compute_twin() {
        twin_computed = true;
        repoint();
        this->View::compute_twin();
        const auto a = this->first_twinless_arc();
        if (View::source(a) != node(-1)) {
            const node u = View::source(a);
            const node v = this->neighbours[this->arcid(u, View::slot(a))];
            twin_computed = false;
            repoint();
            this->surgery_fail("compute_twin", View::Code::AsymmetricAdjacency, u, v);
        }
    }

    // --- Row-count changes: resize plus a row write. ---

    // @anchor owned-push-back
    // @pre  fits: row.size() <= size_t(dmax) -- violation throws
    //       graph_surgery_error{RowFull}, naming the new row's index
    // @post N == N_before + 1 && equal(row, nbrs(N - 1)); a computed twin is
    //       stale on the new row until recomputed (its arcs have no reverses)
    void push_back(const std::vector<node>& row) {
        if (row.size() > size_t(this->dmax))
            this->surgery_fail("push_back", View::Code::RowFull, this->N, node(-1),
                               " (" + std::to_string(row.size()) + " entries for a stride of "
                               + std::to_string(this->dmax) + ")");
        resize(this->N + 1);
        const node u = this->N - 1;
        this->deg[u] = uint8_t(row.size());
        std::copy(row.begin(), row.end(), this->neighbours.begin() + u * this->dmax);
    }

    // @anchor owned-pop-back
    // @pre  nonempty: N > 0 -- violation throws graph_surgery_error{VertexOutOfRange}
    // @post N == N_before - 1
    void pop_back() {
        if (this->N == 0)
            this->surgery_fail("pop_back", View::Code::VertexOutOfRange, node(-1), node(-1));
        resize(this->N - 1);
    }

    using View::push_back;

    // --- Whole-graph relayouts: a view word into a fresh owner of the same
    //     capacity, whose storage this one then takes. ---

    // The same rotation system at another stride (restride_into); the
    // fields past the graph triple ride along unchanged; the twin is
    // dropped (stale at the new stride).
    // @anchor owned-restride
    // @pre  fits: all_of(indices(N), [&](node u){ return degree(u) <= new_dmax; })
    //       -- violation throws graph_surgery_error{RowFull} and leaves this owner unchanged
    // @post dmax == new_dmax && !has_twin() && capacity() == capacity_before
    // @post rows: all_of(indices(N), [&](node u){ return equal(nbrs(u), nbrs_before(u)); })
    void restride_inplace(uint8_t new_dmax) {
        Owned g(this->N, new_dmax, this->capacity());
        this->restride_into(g);
        g.take_tail_fields(*this, {});
        *this = std::move(g);
    }

    // Relabel by pi (pi[u_old] = u_new): row u_old lands in row pi[u_old]
    // with every target relabelled (relabel_into); a per-vertex field
    // follows its vertex; a computed twin travels with its rows and stays
    // valid.  A constant-size field is copied as it stands: a pentagon
    // LIST names vertices, so it is STALE afterwards and its producer
    // re-derives it at the boundary (graphview.hh's pentagon contract).
    // @anchor owned-apply-permutation
    // @pre  permutation: pi.size() == size_t(N) && is_permutation(pi, identity(N))
    //       -- violation throws graph_surgery_error{BadPermutation}, this owner unchanged
    // @post rows: relabel_into's @post rows, over all of [0, N)
    // @post twin: implies(has_twin_before, has_twin() && twin_is_valid() == twin_is_valid_before)
    // @post per_vertex_fields: field[pi[u]] == field_before[u]
    // @post constant_fields: unchanged (stale)
    // @post capacity() == capacity_before
    void apply_permutation(const Permutation& pi) {
        require_permutation("apply_permutation", pi);
        Owned g(this->N, uint8_t(this->dmax), this->capacity());
        g.twin_computed = twin_computed;
        g.repoint();
        this->relabel_into(pi, g);
        g.take_tail_fields(*this, pi);
        *this = std::move(g);
    }

    // Drop every vertex of degree 0, relabelling the survivors in order --
    // relabel_into's partial relabelling, with per-vertex fields and a
    // computed twin following their vertices.
    // @anchor owned-remove-isolated
    // @pre  symmetric: adjacency_is_symmetric() (so no kept row points at a dropped vertex)
    // @post N == count_if(indices(N_before), [&](node u){ return degree_before(u) > 0; })
    // @post order: survivors keep their relative order
    // @post capacity() == capacity_before
    void remove_isolated_vertices() {
        std::vector<int> pi(this->N, -1);
        int kept = 0;
        for (node u = 0; u < this->N; ++u)
            if (this->deg[u] > 0) pi[u] = kept++;
        Owned g(kept, uint8_t(this->dmax), this->capacity());
        g.twin_computed = twin_computed;
        g.repoint();
        this->relabel_into(pi, g);
        g.take_tail_fields(*this, pi);
        *this = std::move(g);
    }

    // Remove the named vertices: every edge at each of them (the view's
    // surgery), then the compaction above -- which also drops any vertex
    // that was ALREADY isolated.  Connectivity is not preserved in general.
    // @anchor owned-remove-vertices
    // @pre  vertices: all_of(sv, [&](int u){ return size_t(u) < size_t(N); })
    //       -- violation throws graph_surgery_error{VertexOutOfRange} before any write
    // @post absent: no vertex of sv survives; every other vertex of positive
    //       degree survives, relabelled in order
    void remove_vertices(const std::set<int>& sv) {
        for (int u : sv) this->require_vertices("remove_vertices", node(u), node(u));
        for (int u : sv)
            while (!(*this)[u].empty())
                this->remove_edge({node(u), (*this)[u][0]});
        remove_isolated_vertices();
    }

  private:
    // The vertex capacity the buffers are sized for, read off the storage
    // itself (the degree buffer's length; see deg_field) -- what repoint and
    // grow_buffers consult before the spans say it.
    int storage_capacity() const {
        return static_cast<int>(std::get<deg_field>(buffers).size());
    }

    // The buffers alone, sized for at least `cap` vertices at the current
    // stride, keeping contents; the spans are NOT repointed.  Its callers
    // here set N and repoint once afterwards (resize, the sized constructor,
    // the deep copy), so a growth is followed by exactly one repoint;
    // reserve, the public form, is this plus that repoint.
    void grow_buffers(int cap) {
        const auto counts = View::get_element_counts(std::max(storage_capacity(), cap), this->dmax);
        batch::for_each_field<View>([&](auto Ic) {
            constexpr std::size_t k = Ic;
            auto& buf = std::get<k>(buffers);
            if (buf.size() < counts[k]) buf.resize(counts[k]);
        });
    }

    // A moved-from owner: no storage, no live vertices, every span empty.
    void release() {
        buffers = buffers_t{};
        twin_computed = false;
        this->N = 0;
        repoint();
    }

    void require_permutation(const char* op, std::span<const int> pi) const {
        std::vector<bool> hit(this->N, false);
        bool ok = pi.size() == size_t(this->N);
        for (std::size_t u = 0; ok && u < pi.size(); ++u) {
            ok = size_t(pi[u]) < size_t(this->N) && !hit[pi[u]];
            if (ok) hit[pi[u]] = true;
        }
        if (!ok)
            this->surgery_fail(op, View::Code::BadPermutation, node(-1), node(-1),
                               " (" + std::to_string(pi.size()) + " entries for "
                               + std::to_string(this->N) + " vertices)");
    }

    // The deep copy behind the converting constructor and assignment: size
    // the buffers for the source, then the library's one field-wise copy of
    // an entry (batch::copy_entry) into this owner's view.  Assigning a view
    // of this owner's own storage is the identity (the aliasing rule the
    // legacy Graph::operator= keeps as well).
    template<class Src>
    void assign(const Src& src) {
        if (owns_memory() && static_cast<const void*>(src.neighbours.data())
                             == static_cast<const void*>(std::get<0>(buffers).data()))
            return;
        this->N = 0;
        this->dmax = src.dmax;
        twin_computed = src.has_twin();
        grow_buffers(src.N);
        repoint();
        batch::copy_entry(*this, src);
    }

    // The fields past the graph triple, taken from `o` into this (freshly
    // sized) owner: a PER-VERTEX field -- one element per vertex, read off
    // the contract (batch::elements_per_vertex) -- follows its vertex through
    // pi, a dropped vertex (pi[u] < 0) taking its element with it; with pi
    // empty, and for a constant-size field, the contents are copied as they
    // stand.
    void take_tail_fields(const Owned& o, std::span<const int> pi) {
        const auto ours   = View::get_element_counts(this->N, this->dmax);
        const auto theirs = View::get_element_counts(o.N, o.dmax);
        const auto per_vertex = batch::elements_per_vertex<View>(o.dmax);
        auto dst = this->to_tuple();
        const auto src = o.to_tuple();
        batch::for_each_field<View>([&](auto Ic) {
            constexpr std::size_t k = Ic;
            if constexpr (k >= n_graph_fields) {
                auto& d = std::get<k>(dst);
                const auto& s = std::get<k>(src);
                if (per_vertex[k] == 1 && !pi.empty()) {
                    for (node u = 0; u < o.N; ++u)
                        if (pi[u] >= 0) d[pi[u]] = s[u];
                } else {
                    std::copy_n(s.data(), std::min(ours[k], theirs[k]), d.data());
                }
            }
        });
    }
};
