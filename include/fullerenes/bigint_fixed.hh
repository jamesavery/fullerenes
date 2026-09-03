#pragma once
// ============================================================================
// bigint_fixed.hh -- FixedBigInt<Limbs>: a fixed-capacity signed big integer.
//
// Arithmetic substrate for exact algebraic-number work (consumer: the
// cyclotomic sign oracles' dyadic fallbacks, cyclotomic.hh and
// cyclotomic_ambient.hh).  Bounded and allocation-free by construction:
// capacity is the template parameter, chosen by the CONSUMER from its own
// precision policy.
//
// Overflow discipline -- the poison monoid: every operation that could
// exceed the capacity sets the sticky poison instead of wrapping, and the
// flag survives every operator (including the zero-operand early returns:
// a zero with a history is still poisoned).  The flag sits behind the
// overflowed() accessor, the house idiom for capacity poison (the
// SpanStack precedent, span_vector.hh); callers must consult it before
// trusting a sign or magnitude.  cmp_mag compares magnitudes only and
// does not consult it.
// ============================================================================

#include <array>
#include <cstdint>

template <int Limbs>
struct FixedBigInt {
  static_assert(Limbs >= 2, "from_i128 writes two limbs unconditionally");

  int sgn = 0;                        // -1, 0, +1
  int n = 0;                          // limbs in use
  std::array<uint64_t, Limbs> l{};    // little-endian magnitude

  bool overflowed() const { return overflowed_; }

  static FixedBigInt from_i128(__int128 v) {
    FixedBigInt b;
    // Negate in the unsigned domain: defined for INT128_MIN too.
    const unsigned __int128 m =
        v < 0 ? -(unsigned __int128)v : (unsigned __int128)v;
    b.sgn = v == 0 ? 0 : (v < 0 ? -1 : 1);
    b.l[0] = (uint64_t)m;
    b.l[1] = (uint64_t)(m >> 64);
    b.n = b.l[1] ? 2 : (b.l[0] ? 1 : 0);
    return b;
  }

  // The poison monoid, stated once.
  static bool poisoned(const FixedBigInt& x, const FixedBigInt& y) {
    return x.overflowed_ || y.overflowed_;
  }

  static int cmp_mag(const FixedBigInt& x, const FixedBigInt& y) {
    if (x.n != y.n) return x.n < y.n ? -1 : 1;
    for (int i = x.n - 1; i >= 0; i--)
      if (x.l[i] != y.l[i]) return x.l[i] < y.l[i] ? -1 : 1;
    return 0;
  }

  friend FixedBigInt operator+(const FixedBigInt& x, const FixedBigInt& y) {
    if (!x.sgn) { FixedBigInt r = y; r.overflowed_ |= x.overflowed_; return r; }
    if (!y.sgn) { FixedBigInt r = x; r.overflowed_ |= y.overflowed_; return r; }
    FixedBigInt r;
    if (x.sgn == y.sgn) {
      r = add_mag(x, y);
      r.sgn = r.n ? x.sgn : 0;
    } else {
      const int c = cmp_mag(x, y);
      if (c == 0) { r.overflowed_ = poisoned(x, y); return r; }
      r = c > 0 ? sub_mag(x, y) : sub_mag(y, x);
      r.sgn = r.n ? (c > 0 ? x.sgn : y.sgn) : 0;
    }
    r.overflowed_ = r.overflowed_ || poisoned(x, y);
    return r;
  }

  friend FixedBigInt operator*(const FixedBigInt& x, const FixedBigInt& y) {
    FixedBigInt r;
    r.overflowed_ = poisoned(x, y);
    if (!x.sgn || !y.sgn) return r;
    if (x.n + y.n > Limbs) { r.overflowed_ = true; return r; }
    for (int i = 0; i < x.n; i++) {
      unsigned __int128 carry = 0;
      for (int j = 0; j < y.n; j++) {
        const unsigned __int128 cur =
            (unsigned __int128)x.l[i] * y.l[j] + r.l[i + j] + carry;
        r.l[i + j] = (uint64_t)cur;
        carry = cur >> 64;
      }
      int k = i + y.n;
      while (carry) {
        if (k >= Limbs) { r.overflowed_ = true; return r; }
        const unsigned __int128 cur = (unsigned __int128)r.l[k] + carry;
        r.l[k] = (uint64_t)cur;
        carry = cur >> 64;
        k++;
      }
    }
    r.n = x.n + y.n;
    while (r.n > 0 && !r.l[r.n - 1]) r.n--;
    r.sgn = r.n ? x.sgn * y.sgn : 0;
    return r;
  }

  friend FixedBigInt operator<<(const FixedBigInt& x, int bits) {
    FixedBigInt r = x;
    if (bits < 0) { r.overflowed_ = true; return r; }   // out of contract
    if (!x.sgn || !bits) return r;
    const int words = bits / 64, rem = bits % 64;
    if (x.n + words + 1 > Limbs) { r.overflowed_ = true; return r; }
    r = FixedBigInt{};
    r.sgn = x.sgn;
    r.overflowed_ = x.overflowed_;
    for (int i = x.n - 1; i >= 0; i--) {
      const uint64_t v = x.l[i];
      if (rem) {
        r.l[i + words + 1] |= (uint64_t)(v >> (64 - rem));
        r.l[i + words] |= v << rem;
      } else {
        r.l[i + words] = v;
      }
    }
    r.n = x.n + words + 1;
    while (r.n > 0 && !r.l[r.n - 1]) r.n--;
    return r;
  }

 private:
  bool overflowed_ = false;           // sticky capacity poison

  static FixedBigInt add_mag(const FixedBigInt& x, const FixedBigInt& y) {
    FixedBigInt r;
    unsigned __int128 carry = 0;
    const int m = x.n > y.n ? x.n : y.n;
    for (int i = 0; i < m; i++) {
      const unsigned __int128 s =
          carry + (i < x.n ? x.l[i] : 0) + (i < y.n ? y.l[i] : 0);
      r.l[i] = (uint64_t)s;
      carry = s >> 64;
    }
    r.n = m;
    if (carry) {
      if (m >= Limbs) { r.overflowed_ = true; return r; }
      r.l[m] = (uint64_t)carry;
      r.n = m + 1;
    }
    return r;
  }
  // @pre |x| >= |y|
  static FixedBigInt sub_mag(const FixedBigInt& x, const FixedBigInt& y) {
    FixedBigInt r;
    unsigned __int128 borrow = 0;
    for (int i = 0; i < x.n; i++) {
      const unsigned __int128 d =
          (unsigned __int128)x.l[i] - (i < y.n ? y.l[i] : 0) - borrow;
      r.l[i] = (uint64_t)d;
      borrow = (d >> 64) ? 1 : 0;
    }
    r.n = x.n;
    while (r.n > 0 && !r.l[r.n - 1]) r.n--;
    return r;
  }
};
