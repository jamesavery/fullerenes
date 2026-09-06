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

#include <array>
#include <concepts>
#include <cstddef>
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

} // namespace batch
