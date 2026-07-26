#pragma once

// Span-backed containers: SpanVector<T> (a span view with optional owned
// vector storage), and two bounded, trivially-copyable workspace containers
// over caller-owned memory -- SpanStack<T> (vector-LIFO semantics, loud
// overflow) and BitSpan (vector<bool> flag-array semantics).

#include <cassert>
#include <cstdint>
#include <iosfwd>
#include <vector>
#include <span>
#include <cstddef>
#include <type_traits>
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

// ---------------------------------------------------------------------------
// SpanStack<T>: bounded LIFO stack over caller-owned storage.
//
// The no-allocation counterpart of a std::vector<T> used purely LIFO
// (back()/pop_back()/push_back()) or a std::stack<T>: a fresh stack starts
// EMPTY, matching a freshly constructed vector, and any op sequence that
// stays within capacity is indistinguishable from the vector's.  On
// exhaustion push_back writes nothing, returns false, and latches a sticky
// overflowed() flag -- overflow is loud, never silent growth or truncation.
//
// @inv  0 <= size() <= capacity()  and  peak() >= size()
// @pre  back() / pop_back(): !empty().  The empty end is the CALLER's
//       obligation (guard with empty(), the vector idiom); only overflow
//       is reported.
//
// peak() is the lifetime high-water mark -- the calibration input for
// sizing backing storage.  clear() re-establishes the fresh-empty state
// and resets the overflow latch; peak() survives it.  Neither pop_back()
// nor clear() touches the storage contents.  live() is the live prefix,
// bottom to top -- the sequence byte-equality gates compare.
//
// A default-constructed stack has capacity 0: every push is rejected and
// latches, so an unbound workspace member fails loudly, never by UB.
//
// Trivially copyable, so it can live inside device-capturable view
// structs.  A copy SHARES the storage but FORKS size/peak/latch: at most
// one live copy may mutate, or the copies corrupt each other's live
// prefix.
// ---------------------------------------------------------------------------
template<typename T>
class SpanStack {
    std::span<T> buf_;
    long n_ = 0, peak_ = 0;
    bool overflowed_ = false;

public:
    SpanStack() = default;
    explicit SpanStack(std::span<T> b) : buf_(b) {}

    // Re-bind the storage, keeping size/peak/latch: the owner repoint idiom
    // (an owner whose backing vector reallocated re-binds the stack to the
    // new storage; the live prefix was copied by the reallocation).
    // @pre b.size() >= size() -- a shrink would orphan the live prefix.
    void rebind(std::span<T> b) {
      assert((long)b.size() >= n_);
      buf_ = b;
    }

    bool empty()      const { return n_ == 0; }
    long size()       const { return n_; }
    long capacity()   const { return (long)buf_.size(); }   // element slots
    bool overflowed() const { return overflowed_; }
    long peak()       const { return peak_; }

    std::span<const T> live() const { return buf_.first((std::size_t)n_); }

    T    back() const { return buf_[n_ - 1]; }
    void pop_back()   { n_--; }

    bool push_back(T x) {
        // Unsigned compare: also contains an out-of-contract negative n_
        // (an unbalanced pop) -- the push then rejects and latches instead
        // of writing before the buffer.  Identical to == for in-contract n_.
        if ((std::size_t)n_ >= buf_.size()) { overflowed_ = true; return false; }
        buf_[n_++] = x;
        if (n_ > peak_) peak_ = n_;
        return true;
    }

    void clear() { n_ = 0; overflowed_ = false; }
};

static_assert(std::is_trivially_copyable_v<SpanStack<int>>,
    "SpanStack must be trivially copyable");

// ---------------------------------------------------------------------------
// BitSpan: bit array over caller-owned words, 64 bits per word.
//
// The no-allocation counterpart of a std::vector<bool> flag array:
// test/set/clear individual bits, test_and_set as the "did I claim this?"
// dedup primitive (returns the PREVIOUS value, then sets -- one
// read-modify-write, so a frontier mark and its test cannot interleave),
// clear_all to establish the all-false fresh state.
//
// Size the backing span as words_for(nbits) words.  capacity() is the
// ADDRESSABLE bit count 64 * words: it rounds the caller's logical bit
// count UP to a word boundary.  Indices in [0, capacity()) are legal;
// bits above the logical count are storage, not data -- scan the logical
// range, never capacity(), when looking for set bits.
//
// @pre  test / set / clear / test_and_set:  0 <= i < capacity()
//
// Construction does NOT zero the storage -- call clear_all() (or own
// zero-initialised storage) before first use.  A default-constructed
// BitSpan has capacity 0, making every index out of contract: bind the
// span before use.
//
// Trivially copyable; copies share the words (flags are shared state --
// there are no per-copy counters to fork).
// ---------------------------------------------------------------------------
class BitSpan {
    std::span<std::uint64_t> bits_;

public:
    // Words needed for nbits logical flags -- THE sizing formula.
    static constexpr long words_for(long nbits) { return (nbits + 63) / 64; }

    BitSpan() = default;
    explicit BitSpan(std::span<std::uint64_t> bits) : bits_(bits) {}

    bool test(long i) const { return (bits_[i >> 6] >> (i & 63)) & 1ull; }
    void set(long i)        { bits_[i >> 6] |=  (1ull << (i & 63)); }
    void clear(long i)      { bits_[i >> 6] &= ~(1ull << (i & 63)); }
    void clear_all()        { for (auto& w : bits_) w = 0; }

    bool test_and_set(long i) {
        const std::uint64_t m = 1ull << (i & 63);
        std::uint64_t& w = bits_[i >> 6];
        const bool prev = w & m;
        w |= m;
        return prev;
    }

    long capacity() const { return (long)bits_.size() * 64; }   // addressable bits
};

static_assert(std::is_trivially_copyable_v<BitSpan>,
    "BitSpan must be trivially copyable");

} // namespace Spanify
