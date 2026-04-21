#pragma once

// Batch<V> / BatchView<V>: generic owning / non-owning container for any
// type satisfying batch::batchable_view<V>.
//
// Memory layout
// -------------
//   For each of V's `n_fields` span-typed fields, Batch<V> owns one
//   contiguous BatchAlloc<T> buffer of length
//
//       capacity * N * get_size_factors(N, dmax)[k].
//
//   Per-entry views are stamped by slicing each field buffer with
//   `(i * N * factor[k], N * factor[k])` and repointing a value-type V
//   through its to_tuple() reference pack.
//
// Scope (Phase 4)
// ---------------
//   - construction/destruction, resize, push_back (host copy),
//     operator[], begin/end iteration, slice(offset, count).
//   - No BatchState, no BatchQueue (separate phases).
//   - Layout-compatible assignment between Batch<A> and Batch<B> (e.g.
//     Batch<PolyhedronView<double>> -> BatchView<PlanarGraphView> over
//     the shared graph-field prefix) is expressed via explicit
//     rebind<To>() helpers when the graph-field prefix matches.

#include "fullerenes/batch/batchable.hh"
#include "fullerenes/batch/storage-policy.hh"

#include <algorithm>
#include <cassert>
#include <cstddef>
#include <span>
#include <tuple>
#include <type_traits>
#include <utility>

namespace batch {

namespace detail {

// tuple_t = std::tuple<std::span<T0>&, std::span<T1>&, ...> -- what V::to_tuple()
// returns.  Strip cvref to get the underlying tuple of lvalue-reference spans.
template<class V>
using field_tuple_t = std::remove_cvref_t<
    decltype(std::declval<V&>().to_tuple())>;

// span<Ti> (value, no reference) for field I.
template<class V, std::size_t I>
using field_span_t = std::remove_reference_t<
    std::tuple_element_t<I, field_tuple_t<V>>>;

// element_type of that span (e.g. int32_t, uint8_t, coord3<double>).
template<class V, std::size_t I>
using field_element_t = typename field_span_t<V, I>::element_type;

// tuple<span<T0>, span<T1>, ...> -- storage-view tuple used by BatchView<V>.
template<class V, std::size_t... Is>
auto make_span_tuple_t(std::index_sequence<Is...>)
    -> std::tuple<field_span_t<V, Is>...>;
template<class V>
using span_tuple_t = decltype(make_span_tuple_t<V>(
    std::make_index_sequence<V::n_fields>{}));

// tuple<BatchAlloc<T0>, BatchAlloc<T1>, ...> -- owning tuple used by Batch<V>.
template<class V, std::size_t... Is>
auto make_buffer_tuple_t(std::index_sequence<Is...>)
    -> std::tuple<BatchAlloc<std::remove_const_t<field_element_t<V, Is>>>...>;
template<class V>
using buffer_tuple_t = decltype(make_buffer_tuple_t<V>(
    std::make_index_sequence<V::n_fields>{}));

// Apply a callable at each field index 0..V::n_fields-1.
template<class V, class F>
constexpr void for_each_field(F&& f) {
    [&]<std::size_t... Is>(std::index_sequence<Is...>) {
        (f(std::integral_constant<std::size_t, Is>{}), ...);
    }(std::make_index_sequence<V::n_fields>{});
}

} // namespace detail

// ---------------------------------------------------------------------------
// BatchView<V>: non-owning batch slice.
// ---------------------------------------------------------------------------
template<class V>
    requires batchable_view<V>
class BatchView {
public:
    using view_type  = V;
    using spans_type = detail::span_tuple_t<V>;

    BatchView() = default;

    BatchView(int N, int dmax, int size, spans_type spans)
        : N_(N), dmax_(dmax), size_(size), spans_(spans) {}

    // Stamp a per-entry view by slicing each field span.
    V operator[](std::size_t i) const {
        assert(int(i) < size_);
        V v;
        v.N    = typename V::node_type(N_);
        v.dmax = dmax_;
        const auto factors = V::get_size_factors(N_, dmax_);
        auto vtup = v.to_tuple();
        detail::for_each_field<V>([&](auto Ic) {
            constexpr std::size_t k = Ic;
            const std::size_t per = std::size_t(N_) * factors[k];
            std::get<k>(vtup) = std::get<k>(spans_).subspan(i * per, per);
        });
        return v;
    }

    // Sub-range with shared storage.
    BatchView slice(std::size_t offset, int count) const {
        assert(int(offset) + count <= size_);
        const auto factors = V::get_size_factors(N_, dmax_);
        spans_type sub{};
        detail::for_each_field<V>([&](auto Ic) {
            constexpr std::size_t k = Ic;
            const std::size_t per = std::size_t(N_) * factors[k];
            std::get<k>(sub) = std::get<k>(spans_).subspan(offset * per, std::size_t(count) * per);
        });
        return BatchView(N_, dmax_, count, sub);
    }

    // Lightweight random-access iterator yielding V by value.
    struct iterator {
        using difference_type   = std::ptrdiff_t;
        using value_type        = V;
        using reference         = V;
        using pointer           = void;
        using iterator_category = std::random_access_iterator_tag;

        const BatchView* bv_ = nullptr;
        std::ptrdiff_t   i_  = 0;

        V operator*()  const { return (*bv_)[std::size_t(i_)]; }
        V operator[](difference_type n) const { return (*bv_)[std::size_t(i_ + n)]; }

        iterator& operator++()    { ++i_; return *this; }
        iterator  operator++(int) { auto t = *this; ++i_; return t; }
        iterator& operator--()    { --i_; return *this; }
        iterator  operator--(int) { auto t = *this; --i_; return t; }

        iterator& operator+=(difference_type n) { i_ += n; return *this; }
        iterator& operator-=(difference_type n) { i_ -= n; return *this; }
        iterator  operator+ (difference_type n) const { return {bv_, i_ + n}; }
        iterator  operator- (difference_type n) const { return {bv_, i_ - n}; }
        difference_type operator-(const iterator& o) const { return i_ - o.i_; }

        friend bool operator==(const iterator& a, const iterator& b) { return a.i_ == b.i_ && a.bv_ == b.bv_; }
        friend bool operator!=(const iterator& a, const iterator& b) { return !(a == b); }
        friend bool operator< (const iterator& a, const iterator& b) { return a.i_ <  b.i_; }
        friend bool operator<=(const iterator& a, const iterator& b) { return a.i_ <= b.i_; }
        friend bool operator> (const iterator& a, const iterator& b) { return a.i_ >  b.i_; }
        friend bool operator>=(const iterator& a, const iterator& b) { return a.i_ >= b.i_; }
    };

    iterator begin() const { return {this, 0}; }
    iterator end()   const { return {this, size_}; }

    // --- Accessors ---
    int  size()     const { return size_; }
    int  capacity() const { return size_; }   // views expose their extent only
    int  N()        const { return N_; }
    int  dmax()     const { return dmax_; }
    bool empty()    const { return size_ == 0; }

    const spans_type& spans() const { return spans_; }

private:
    int        N_    = 0;
    int        dmax_ = 0;
    int        size_ = 0;
    spans_type spans_{};
};

// ---------------------------------------------------------------------------
// Batch<V>: owning batch.
// ---------------------------------------------------------------------------
template<class V>
    requires batchable_view<V>
class Batch {
public:
    using view_type     = V;
    using buffers_type  = detail::buffer_tuple_t<V>;
    using spans_type    = detail::span_tuple_t<V>;

    Batch() = default;

    // Allocate storage for up to `capacity` entries, each with `N` vertices
    // and row stride `dmax`.  Size starts at 0.
    Batch(int N, int capacity, int dmax = V::default_dmax)
        : N_(N), dmax_(dmax), size_(0), capacity_(capacity) {
        allocate_();
    }

    // --- Capacity management ---
    void reserve(int new_capacity) {
        if (new_capacity <= capacity_) return;
        capacity_ = new_capacity;
        allocate_();
    }
    void resize(int new_size) {
        assert(new_size <= capacity_);
        size_ = new_size;
    }
    void clear() { size_ = 0; }

    // --- View accessors ---
    BatchView<V> view() const {
        return BatchView<V>(N_, dmax_, size_, build_spans_(capacity_));
    }
    BatchView<V> slice(std::size_t offset, int count) const {
        return view().slice(offset, count);
    }
    V operator[](std::size_t i) const { return view()[i]; }

    // --- Mutation ---
    // Append one entry by copying each field from `v`.  The source spans
    // in `v` must have size N * get_size_factors(N, dmax)[k].
    void push_back(const V& v) {
        assert(size_ < capacity_);
        assert(v.N == typename V::node_type(N_));
        assert(v.dmax == dmax_);
        const auto factors = V::get_size_factors(N_, dmax_);
        auto src = v.to_tuple();
        detail::for_each_field<V>([&](auto Ic) {
            constexpr std::size_t k = Ic;
            const std::size_t per = std::size_t(N_) * factors[k];
            auto* dst = std::get<k>(buffers_).data() + std::size_t(size_) * per;
            const auto& s = std::get<k>(src);
            if (s.size() >= per)
                std::copy_n(s.data(), per, dst);
            else
                std::copy_n(s.data(), s.size(), dst); // empty optional fields
        });
        ++size_;
    }

    // --- Accessors ---
    int  size()     const { return size_; }
    int  capacity() const { return capacity_; }
    int  N()        const { return N_; }
    int  dmax()     const { return dmax_; }
    bool empty()    const { return size_ == 0; }

private:
    void allocate_() {
        const auto factors = V::get_size_factors(N_, dmax_);
        detail::for_each_field<V>([&](auto Ic) {
            constexpr std::size_t k = Ic;
            const std::size_t total = std::size_t(capacity_) * std::size_t(N_) * factors[k];
            std::get<k>(buffers_).resize(total);
        });
    }

    // Build span_tuple covering the first `entries` entries of the buffers.
    spans_type build_spans_(int entries) const {
        const auto factors = V::get_size_factors(N_, dmax_);
        spans_type spans{};
        detail::for_each_field<V>([&](auto Ic) {
            constexpr std::size_t k = Ic;
            const std::size_t total = std::size_t(entries) * std::size_t(N_) * factors[k];
            using elem_t = detail::field_element_t<V, k>;
            // BatchAlloc<T> exposes T* data().  Construct the span explicitly.
            auto* p = const_cast<elem_t*>(std::get<k>(buffers_).data());
            std::get<k>(spans) = std::span<elem_t>(p, total);
        });
        return spans;
    }

    int           N_        = 0;
    int           dmax_     = 0;
    int           size_     = 0;
    int           capacity_ = 0;
    buffers_type  buffers_{};
};

} // namespace batch
