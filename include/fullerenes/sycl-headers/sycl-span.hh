#pragma once
#include <span>
#include <ostream>
#include <array>
#include <algorithm>
#include <limits>
#include <cmath>
#include <type_traits>
#include <cstddef>

// Reinterpret-cast a std::span<T> as a std::span<U>.
template <typename U, typename T>
constexpr std::span<U> as_span(std::span<T> s) {
    return std::span<U>(reinterpret_cast<U*>(s.data()),
                        (s.size() * sizeof(T)) / sizeof(U));
}

// Forward declaration for the std::array stream operator (defined generically over
// element type and size in sycl-util-impl.cc, with more-specialized <T,2>/<T,3>
// overloads in coord3d.cc). The generic span operator<< below streams each element
// (os << v[i]); when the element is a std::array, that resolves to this overload.
// Under the strict two-phase lookup of a conforming compiler, os<<element is bound at
// THIS template's definition point, and std::array's only associated namespace is std,
// so this global-namespace overload must be visible here -- otherwise the call is found
// by neither template-definition lookup nor ADL. (gcc accepted it; clang requires the
// prior declaration.) One generic declaration covers every array size used.
template <typename U, std::size_t N>
constexpr std::ostream& operator<<(std::ostream& os, const std::array<U,N>& a);

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
