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

} // namespace batch
