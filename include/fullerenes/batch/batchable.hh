#pragma once

// batchable_view concept.
//
// A view V is batchable if it provides an intrinsic layout description
// that Batch<V> and BatchView<V> can use to slice a single contiguous
// allocation per field into per-entry views:
//
//   static constexpr std::size_t V::n_fields
//     Tuple arity -- number of std::span fields in the view.
//
//   auto V::to_tuple()       -> std::tuple of references to the span fields
//   auto V::to_tuple() const -> same, but over const-view fields
//     Returned in canonical order:
//       graph-like:  { neighbours, deg, twin }
//       geometry:    { neighbours, deg, twin, points }
//     The tuple must hold std::span<T>& references so that Batch<V> can
//     repoint them when stamping per-entry views.
//
//   static std::array<std::size_t, V::n_fields>
//     V::get_element_counts(int N, int dmax)
//       Element count of each field for ONE batch entry -- absolute counts,
//       not per-vertex factors, so N-proportional fields (neighbours: N*dmax)
//       and constant-size fields (FullereneDualView's 12 pentagon ids) ride
//       the same law.
//
// Views that add fields (e.g. PolyhedronView<T> adds `points`) override
// n_fields / to_tuple / get_element_counts to extend the base graph tuple.
//
// No external trait table is needed: batchability is expressed entirely
// by the view type itself.

#include <algorithm>
#include <array>
#include <concepts>
#include <cstddef>
#include <stdexcept>
#include <string>
#include <tuple>
#include <type_traits>
#include <utility>

namespace batch {

// -- Concept ---------------------------------------------------------------

template<class V>
concept batchable_view =
    std::is_trivially_copyable_v<V> &&
    requires(V v, const V cv) {
        { V::n_fields } -> std::convertible_to<std::size_t>;
        { v.to_tuple() };
        { cv.to_tuple() };
        { V::get_element_counts(0, 0) }
            -> std::same_as<std::array<std::size_t, V::n_fields>>;
    };

// -- Derived types of the contract -----------------------------------------
//
// The field tuple as V::to_tuple() returns it, each field's span type and
// element type, and the compile-time loop over the field indices.  Whatever
// STORES or SLICES a batchable view's fields -- Batch<V> and BatchView<V>
// (batch.hh), Owned<V> (owned.hh) -- derives its per-field storage from
// these, so the fields are named in one place: the view's to_tuple().

// tuple<span<T0>&, span<T1>&, ...> with cvref stripped.
template<class V>
using field_tuple_t = std::remove_cvref_t<decltype(std::declval<V&>().to_tuple())>;

// span<Ti> (a value, not a reference) for field I.
template<class V, std::size_t I>
using field_span_t = std::remove_reference_t<std::tuple_element_t<I, field_tuple_t<V>>>;

// The element type of that span (int32_t, uint8_t, coord3<double>, ...).
template<class V, std::size_t I>
using field_element_t = typename field_span_t<V, I>::element_type;

// Apply f(std::integral_constant<std::size_t, k>) for k = 0 .. V::n_fields-1.
template<class V, class F>
constexpr void for_each_field(F&& f) {
    [&]<std::size_t... Is>(std::index_sequence<Is...>) {
        (f(std::integral_constant<std::size_t, Is>{}), ...);
    }(std::make_index_sequence<V::n_fields>{});
}

// -- Layout compatibility --------------------------------------------------

// Two batchable views share a layout iff their element counts agree
// field-wise for the given (N, dmax).  This is the prerequisite for
// batch-of-A to be reinterpretable as batch-of-B (e.g. slicing a
// PolyhedronView batch into its underlying graph layout).
template<class A, class B>
constexpr bool layout_compatible(int N, int dmax) {
    constexpr std::size_t K =
        A::n_fields < B::n_fields ? A::n_fields : B::n_fields;
    auto a = A::get_element_counts(N, dmax);
    auto b = B::get_element_counts(N, dmax);
    for (std::size_t k = 0; k < K; ++k)
        if (a[k] != b[k]) return false;
    return true;
}

// -- Field shape -----------------------------------------------------------

// How each field scales with the vertex count, read off the contract
// itself: its element count at N = 2 minus at N = 1 -- 0 for a
// constant-size field (a dual's pentagon list), 1 for one element per
// vertex (degrees, coordinates), dmax for one per arc (adjacency, twin).
// What a relabelling needs to move a field with its vertex (Owned<V>), and
// what a batch will need the first time it permutes or compacts.
template<class V>
constexpr std::array<std::size_t, V::n_fields> elements_per_vertex(int dmax) {
    const auto one = V::get_element_counts(1, dmax);
    const auto two = V::get_element_counts(2, dmax);
    std::array<std::size_t, V::n_fields> per{};
    for (std::size_t k = 0; k < V::n_fields; ++k) per[k] = two[k] - one[k];
    return per;
}

// -- The one field-wise copy of an entry ------------------------------------
//
// dst takes src's graph.  For every field the two contracts share, the
// source entry's elements -- get_element_counts(src.N, src.dmax) of them,
// what src holds for its N vertices -- are copied to the front of dst's
// span; a field dst has and src has not is value-initialised over the same
// extent; then dst.N = src.N.  A field whose SOURCE span is empty is absent
// (an uncomputed twin table) and dst's is left alone.  dst's spans may
// cover more rows than src.N -- a view over storage sized for growth: an
// owner's, a slot of a batch, an enumerator's working graph -- and the
// rows past src.N are not touched.  The two must agree on the stride
// (std::invalid_argument otherwise) and every destination span must hold
// the source's count (std::length_error otherwise); this header is generic
// and has no graph error type.  Owned<V>::assign is this after sizing its
// buffers; Batch<V>::push_back is a candidate.
//
// @pre  dst.dmax == src.dmax
// @pre  every shared field with a present source span: dst span holds
//       Src::get_element_counts(src.N, src.dmax)[k] elements
// @post dst.N == src.N and the first N vertices' worth of every shared
//       field equals src's
template<class Dst, class Src>
void copy_entry(Dst& dst, const Src& src) {
    if (dst.dmax != src.dmax)
        throw std::invalid_argument("batch::copy_entry: the strides differ ("
                                    + std::to_string(int(dst.dmax)) + " vs "
                                    + std::to_string(int(src.dmax)) + ")");
    const auto theirs = Src::get_element_counts(src.N, src.dmax);
    const auto ours   = Dst::get_element_counts(src.N, src.dmax);
    auto d = dst.to_tuple();
    const auto s = src.to_tuple();
    const auto require = [](std::size_t held, std::size_t needed, std::size_t k) {
        if (held < needed)
            throw std::length_error("batch::copy_entry: field " + std::to_string(k)
                                    + " of the destination holds " + std::to_string(held)
                                    + " elements, the source entry needs " + std::to_string(needed));
    };
    for_each_field<Dst>([&](auto Ic) {
        constexpr std::size_t k = Ic;
        auto& dk = std::get<k>(d);
        using elem_t = typename std::remove_reference_t<decltype(dk)>::element_type;
        if constexpr (k < Src::n_fields) {
            static_assert(std::is_same_v<field_element_t<Src, k>, elem_t>,
                          "a shared field must have one element type in both contracts");
            const auto& sk = std::get<k>(s);
            if (sk.empty()) return;
            require(dk.size(), theirs[k], k);
            std::copy_n(sk.data(), theirs[k], dk.data());
        } else {
            require(dk.size(), ours[k], k);
            std::fill_n(dk.data(), ours[k], elem_t{});
        }
    });
    dst.N = src.N;
}

} // namespace batch
