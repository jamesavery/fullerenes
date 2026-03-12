#pragma once

#include <vector>
#include <span>
#include <cstddef>
#include <utility>

namespace Spanify {

// ---------------------------------------------------------------------------
// SpanVector<T>: a span<T> view with optional owned vector<T> storage.
//
// Unifies the view/owned pattern used across the codebase:
//   - Viewing external memory: SpanVector(span<T>)
//   - Owning storage:          SpanVector(vector<T>) or sv = vector<T>{...}
//
// The span view is always the active reference for element access.
// When owning, repoint() keeps the span in sync with the vector.
//
// Rule of 5 is correct: copies/moves propagate ownership correctly
// and repoint the span.
// ---------------------------------------------------------------------------
template<typename T>
struct SpanVector {
    std::span<T> view;              // always the active reference
    std::vector<T> owned;           // storage (empty when viewing external data)

    // Repoint view to owned storage.  Call after mutating owned.
    void repoint() { view = std::span<T>(owned); }

    // --- Constructors ---

    SpanVector() = default;

    // Owning: take a vector.
    SpanVector(std::vector<T> v) : owned(std::move(v)) { repoint(); }

    // Viewing: wrap external span (caller manages lifetime).
    SpanVector(std::span<T> s) : view(s) {}

    // --- Rule of 5 ---

    SpanVector(const SpanVector& o) : owned(o.owned) {
        if (!owned.empty()) repoint();
        else view = o.view;
    }

    SpanVector(SpanVector&& o) noexcept
        : owned(std::move(o.owned)) {
        if (!owned.empty()) repoint();
        else view = o.view;
        o.view = {};
    }

    SpanVector& operator=(const SpanVector& o) {
        if (this != &o) {
            owned = o.owned;
            if (!owned.empty()) repoint();
            else view = o.view;
        }
        return *this;
    }

    SpanVector& operator=(SpanVector&& o) noexcept {
        owned = std::move(o.owned);
        if (!owned.empty()) repoint();
        else view = o.view;
        o.view = {};
        return *this;
    }

    // Assign from vector (take ownership).
    // Replaces set_points / set_values / etc.
    SpanVector& operator=(std::vector<T> v) {
        owned = std::move(v);
        repoint();
        return *this;
    }

    // --- Element access ---

    T& operator[](size_t i) { return view[i]; }
    const T& operator[](size_t i) const { return view[i]; }

    size_t size() const { return view.size(); }
    bool empty() const { return view.empty(); }
    bool owns_memory() const { return !owned.empty(); }
    T* data() { return view.data(); }
    const T* data() const { return view.data(); }

    auto begin() { return view.begin(); }
    auto end() { return view.end(); }
    auto begin() const { return view.begin(); }
    auto end() const { return view.end(); }

    // --- Implicit conversion to span ---
    // Three overloads cover all calling conventions:
    //   non-const SpanVector → span<T>       (mutable view)
    //   non-const SpanVector → span<const T> (read-only view from mutable obj)
    //   const SpanVector     → span<const T> (read-only view from const obj)

    operator std::span<T>() { return view; }
    operator std::span<const T>() { return view; }
    operator std::span<const T>() const { return view; }
};

// --- ostream support (delegates to span<T> operator<< from auxiliary.hh) ---
template<typename T>
std::ostream& operator<<(std::ostream& s, const SpanVector<T>& sv) {
    return s << sv.view;
}

} // namespace Spanify
