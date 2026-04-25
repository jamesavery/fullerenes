#pragma once
// Phase 9: Span<T> is a thin wrapper around std::span<T> that keeps the
// extra `as_span<U>()` reinterpret helper used pervasively by the SYCL
// kernels. All other behaviour is inherited from std::span.

#include <span>
#include <ostream>
#include <algorithm>
#include <limits>
#include <cmath>
#include <type_traits>
#include <cstddef>

template <typename T>
struct Span : public std::span<T> {
    using base = std::span<T>;
    using base::base;                          // inherit std::span constructors
    using typename base::element_type;
    using typename base::value_type;
    using typename base::size_type;

    constexpr Span() noexcept = default;
    constexpr Span(const Span&) noexcept = default;
    constexpr Span& operator=(const Span&) noexcept = default;
    constexpr Span(base b) noexcept : base(b) {}

    // Convert any pointer-like object (e.g. sycl::multi_ptr) to T* and wrap
    // it as a span of `sz` elements. Needed because std::span's (It, size)
    // constructor requires `It` to satisfy std::contiguous_iterator, which
    // SYCL's multi_ptr does not. Excluded for raw pointers to avoid
    // ambiguity with the inherited std::span(T*, size_type) ctor.
    template <typename Ptr,
              typename = std::enable_if_t<!std::is_pointer_v<std::decay_t<Ptr>>>,
              typename = decltype(static_cast<T*>(std::declval<Ptr&>()))>
    constexpr Span(Ptr&& p, std::size_t sz)
        : base(static_cast<T*>(p), sz) {}

    // Reinterpret-cast view of the underlying storage as a span of U.
    template <typename U>
    constexpr Span<U> as_span() const {
        return Span<U>(reinterpret_cast<U*>(this->data()),
                       (this->size() * sizeof(T)) / sizeof(U));
    }

    // Override subspan so it returns a Span<T>, preserving as_span<U>()
    // chaining in existing call sites.
    constexpr Span<T> subspan(size_type offset) const {
        return Span<T>(this->data() + offset, this->size() - offset);
    }
    constexpr Span<T> subspan(size_type offset, size_type count) const {
        return Span<T>(this->data() + offset, count);
    }
};

// Deduction guides – cover raw pointers, iterator pairs, and iterator+size.
template <typename T> Span(T*, std::size_t) -> Span<T>;
template <typename It> Span(It, It) -> Span<typename std::iterator_traits<It>::value_type>;

template <typename T>
inline std::ostream& operator<<(std::ostream& os, std::span<T> v) {
    os << "[";
    for (std::size_t i = 0; i < v.size(); ++i) {
        os << v[i];
        if (i + 1 < v.size()) os << ", ";
    }
    os << "]";
    return os;
}

template <typename T>
inline bool span_fuzzy_equal(std::span<T> a, std::span<T> b) {
    if (a.size() != b.size()) return false;
    if (a.data() == b.data()) return true;
    if constexpr (std::is_floating_point_v<std::decay_t<T>>) {
        return std::equal(a.begin(), a.end(), b.begin(), [](const auto& x, const auto& y) {
            T eps = std::numeric_limits<T>::epsilon() * 20;
            T max_v = std::max(std::abs(x), std::abs(y));
            return std::abs(x - y) / (max_v > eps ? max_v : 1) < eps;
        });
    } else {
        return std::equal(a.begin(), a.end(), b.begin());
    }
}

// std::span does not define operator==; provide a fuzzy-element comparison
// so legacy code that compared SyclVector/FullereneData members through
// spans keeps working.
template <typename T>
inline bool operator==(std::span<T> a, std::span<T> b) {
    return span_fuzzy_equal(a, b);
}
