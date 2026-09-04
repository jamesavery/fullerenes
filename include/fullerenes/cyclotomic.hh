#pragma once
// ============================================================================
// cyclotomic.hh -- exact arithmetic in the RANK-4 ring Z[gamma],
// gamma = 2 cos(pi/15): the working algebra of the FLATTENED kis surface
// (regular unit pentagons and hexagons; curvature k*pi/15 at the sixty
// pentagon corners).  Notation, here and below: gamma_k := two_cos(k)
// = 2 cos(k pi/15), so gamma = gamma_1 and golden = gamma_3.
//
// Mathematics: claude-projects/delaunay/cyclotomic-idt.tex (the whitepaper;
// its sigma-descent, @ref sec:rank4 there, supersedes the conductor-60
// account that kis-metrics.tex sec. 4 carried before its 2026-08-25
// revision).  Acceptance frame: CANONICAL-TESSELATION.md (fullerenes repo
// root); this module's role in it is the 2026-08-24 certificate-key ruling
// (claude-projects/parallel-primitives/benchmarks/README.md).
// Implementation debt entry: claude-projects/delaunay/refactor-debt.md
// (2026-08-24-cyclotomic-algebra-and-idt-extension), of which this module
// is layer 1; layer 2 (the CyclotomicMetric iDT policy) builds on the
// Diamond classifier at the bottom of this file.
//
// The structure implemented (whitepaper @ref sec:rank4, sec:ring): every
// vertex of the kis complex lies in (1/5) Z[zeta_30] (the sigma-descent --
// the pentagon spoke is sigma-even, witnessed by the walk identity
// @ref eq:walk), so every squared length AND every face wedge w -- the
// delta-normalized cross product, with 16 Area^2 = (2 - gamma_2) w^2
// (@ref eq:heron) and 2 - gamma_2 a UNIT (@ref eq:unit) -- lies in
// (1/25) Z[gamma].  One global scaling by kLsqScale = kPointScale^2 makes
// both integral, and the geometric input of the cubic chain is the four
// constants of Real30's table (plus the Heron factor and its inverse, a
// ring identity): no development coordinate is ever computed.  Every
// predicate is then the exact sign of ONE ring element (delaunay:
// s_u w_l + s_l w_u; convexity: Q w_u + P w_l), and flips TRANSPORT the
// carry by ring arithmetic plus two exact divisions by 2e
// (Diamond::flipped, @ref eq:flip).  The ambient rings Z[zeta_30] /
// Z[zeta_60] are proof scaffolding and live in cyclotomic_ambient.hh as
// construction and cross-check machinery -- nothing on this header's run
// path touches them.
//
// What is deliberately NOT here (the Eisenstein gifts this surface lacks):
// no discrete lattice, hence no llround-style float->exact snap -- and no
// integer square root, hence the wedge is CARRIED, never re-extracted from
// the Heron form.  Exactness enters at CONSTRUCTION (the constants table),
// never by rounding a double.
//
// Representation: Real30 = int64 coordinates over the power basis
// 1, gamma, gamma^2, gamma^3, reduced by
// psi(y) = y^4 + y^3 - 4 y^2 - 4 y + 1 (the minimal polynomial; conjugates
// 2 cos(k pi/15), k in {1,7,11,13}).
//
// OVERFLOW DISCIPLINE -- checked-and-poisoned, never silent: EVERY ring
// operation checks its own coefficient arithmetic (int64 overflow on
// add/sub/scale, the int64 store-back and the unsigned-magnitude operand
// guard on the __int128 product accumulation) and POISONS the value
// (ok = false) instead of wrapping; poison is sticky, and every consumer
// refuses a poisoned value with a NAMED reason (Refusal, carried by
// SignTrace / DivTrace / the Diamond factory) -- a wrong sign from
// overflow is therefore structurally impossible; large inputs refuse
// loudly.  The guaranteed-unpoisoned input envelopes are the NAMED
// constants kPredicateCoeffMax / kFlipCoeffMax below, each pinned by a
// static_assert spelling its binding term.  Today's kis constants have
// |c| <= 100.
//
// EXACT DIVISION (the flip divides by 2e, twice): the quotient is
// computed in F_p[y]/(psi) for a fixed prime p (extended Euclid),
// symmetrically lifted, and then VERIFIED by the ring's own checked
// product -- exact_div returns q only if divisor * q == numerator holds
// exactly in Z[gamma].  The mod-p step is a fast guess; the verification
// is the correctness argument, so a wrong lift, an unlucky prime, or
// genuine non-divisibility all degrade to a named refusal, never to a
// wrong quotient.  (The primes' primality is not load-bearing for
// correctness -- a composite would only lower the success rate; the
// verification gate is total.)  psi SPLITS COMPLETELY mod the first prime
// 2^61 - 1, so first-prime-singular divisors exist from coefficient
// height ~p^(1/4) ~ 3.3e4 (a pinned witness lives in the test suite's [X]
// group); psi is IRREDUCIBLE mod the second prime 2^63 - 25, so the
// fallback succeeds for every divisor nonzero mod it.  For real kis
// inputs (|c| <= 100) |N(2e)| < 2^61 - 1, and the first prime always
// succeeds.  A failed verification retries the other prime before
// refusing.  Divisibility itself is guaranteed by the module geometry
// (whitepaper @ref sec:flip); a refusal on module inputs falsifies the
// caller's premise and must trip loudly upstream.
//
// TIER NOTE: add/sub/mul and rungs (0)-(1) of sign_real are allocation-free
// fixed-array code with NO transcendental call anywhere -- gamma enters
// rung (1) as a correctly rounded stored constant (kGammaDouble), so the
// double rung is bit-identical across hosts and devices up to FMA
// contraction, which the error bound covers.  The exact dyadic fallback
// (rung 2) is a SEPARATE, host-tier function with a fixed multi-KB frame
// (FixedBigInt arrays); the fast path never pays it.  Nothing here is
// device-TESTED yet -- the device story is layer-2 work.
// ============================================================================

#include "bigint_fixed.hh"
#include "diamond_forms.hh"
#include "metric_forms.hh"

#include <array>
#include <cmath>
#include <cstdint>
#include <limits>
#include <optional>

namespace cyclotomic {

// The module scale: kis points live in (1/kPointScale) Z[zeta_30], so the
// bilinear quantities -- squared lengths and wedges -- carry denominator
// kPointScale^2, cleared once by kLsqScale.
inline constexpr long long kPointScale = 5;
inline constexpr long long kLsqScale = kPointScale * kPointScale;

// ---------------------------------------------------------------------------
// The verdict and refusal vocabulary.  Sign is the sign group alone; a
// refusal travels as the optional's nullopt (the house convention, as in
// DiamondSq's optional forms) PLUS a named Refusal through the observability
// channels (SignTrace, DivTrace, the Diamond factory's out-parameter) --
// so a caller can always distinguish WHY it was refused.
// ---------------------------------------------------------------------------

enum class Sign : int { Negative = -1, Zero = 0, Positive = 1 };
inline Sign operator-(Sign s) { return (Sign)(-(int)s); }
using SignOr = std::optional<Sign>;   // nullopt = refused (see Refusal)
inline SignOr operator-(SignOr s) { return s ? SignOr(-*s) : s; }

enum class Refusal : int {
  None = 0,
  Poisoned,           // a coefficient overflowed somewhere in the history
  NotReal,            // ambient tier: sign of a non-real element requested
  DyadicCap,          // rung 2 hit kMaxDyadicBits (fail-loud backstop)
  BigWidth,           // rung 2 overflowed a FixedBigInt capacity
  DivisorZero,        // exact_div by zero
  NotDivisible,       // exact_div: no quotient verified (either prime)
  NonPositiveWedge,   // Diamond: a face wedge not strictly positive
  InconsistentCarry,  // Diamond: Heron identity fails against the carry
  NotConvex,          // Diamond::flip on a non-convex diamond
};

// Optional observability of the sign ladder (the counting-metric pattern):
// which rung decided, how many bisections rung 2 spent, and -- on a
// refusal -- its name.
struct SignTrace {
  int rung = -1;
  int bisections = 0;
  Refusal refusal = Refusal::None;
};

namespace detail {

inline Sign sign_from_int(int s) {
  return s > 0 ? Sign::Positive : (s < 0 ? Sign::Negative : Sign::Zero);
}

// psi, the minimal polynomial of gamma (monic, degree 4), index = power --
// THE authority; the reduction row below is derived from it.
inline constexpr std::array<long long, 5> kPsi30 = {1, -4, -4, 1, 1};

// gamma^4 = -(psi's tail): the fold row every reduction uses.
inline constexpr auto make_reduce_row() {
  std::array<long long, 4> r{};
  for (int i = 0; i < 4; i++) r[i] = -kPsi30[i];
  return r;
}
inline constexpr auto kGamReduceRow = make_reduce_row();

// kGamPow[k] = the reduced coordinates of gamma^k, k = 0..6 (a product of
// two cubics has degree <= 6; indices 4..6 are the fold rows).
inline constexpr auto make_gam_pow() {
  std::array<std::array<long long, 4>, 7> P{};
  P[0][0] = 1;
  for (int k = 1; k < 7; k++) {
    const long long top = P[k - 1][3];
    for (int i = 3; i > 0; i--) P[k][i] = P[k - 1][i - 1];
    P[k][0] = 0;
    for (int i = 0; i < 4; i++) P[k][i] += top * kGamReduceRow[i];
  }
  return P;
}
inline constexpr auto kGamPow = make_gam_pow();
static_assert(kGamPow[5] == std::array<long long, 4>{1, -5, 0, 5});
static_assert(kGamPow[6] == std::array<long long, 4>{-5, 21, 15, -5});

// kTwoCos15[k] = 2 cos(k pi/15) in the power basis, k = 0..15, by the
// dilated-Chebyshev recurrence p_k = gamma p_{k-1} - p_{k-2} (the shift
// spill folded through kGamReduceRow, the same idiom as make_gam_pow).
inline constexpr auto make_two_cos15() {
  std::array<std::array<long long, 4>, 16> P{};
  P[0][0] = 2;
  P[1][1] = 1;
  for (int k = 2; k < 16; k++) {
    long long t[5] = {};
    for (int i = 0; i < 4; i++) t[i + 1] = P[k - 1][i];
    const long long top = t[4];
    for (int i = 0; i < 4; i++)
      P[k][i] = t[i] + top * kGamReduceRow[i] - P[k - 2][i];
  }
  return P;
}
inline constexpr auto kTwoCos15 = make_two_cos15();
static_assert(kTwoCos15[2] == std::array<long long, 4>{-2, 0, 1, 0});
static_assert(kTwoCos15[3] == std::array<long long, 4>{0, -3, 0, 1});
static_assert(kTwoCos15[4] == std::array<long long, 4>{1, 4, 0, -1});

// Operand bound under which a product's __int128 accumulation cannot
// overflow: the convolution slots fold into per-coordinate amplification
// <= 45 * B^2 for the rank-4 ring (exactly (11,45,30,22) per coordinate)
// and <= 58 * B^2 for the ambient rank-16 ring, and 58 * (2^60)^2 < 2^126
// < 2^127.  Compared UNSIGNED against the unsigned magnitude (an INT64_MIN
// coordinate is 2^63 and must trip; the 2026-08-25 review exhibited the
// signed comparison letting it through into undefined behaviour).
inline constexpr unsigned long long kMulOperandMax = 1ULL << 60;

// gamma as a correctly rounded double (verified against an exact rational
// Newton iteration on psi and by the test suite's long-double cross-check).
// 1.9562952014676114 round-trips this exact bit pattern.
inline constexpr double kGammaDouble = 0x1.f4cfc327a0080p+0;

}  // namespace detail

// ---------------------------------------------------------------------------
// CheckedCoeffRing: the checked-and-poisoned coefficient ring over the
// power basis 1, x, .., x^(Rank-1) -- the ONE spelling of the scaffold the
// working ring (Real30) and the ambient ring (Zeta60) share.  `ok == false`
// marks a value whose coefficients overflowed somewhere in its history;
// every consumer refuses it by name.  The derived ring supplies its
// reduced-power table via pow_row(k) (indices Rank .. 2 Rank - 2 fold a
// product; the ambient conjugation reads the full table).
// ---------------------------------------------------------------------------
template <class D, int Rank>
struct CheckedCoeffRing {
  static constexpr int kRank = Rank;

  std::array<long long, Rank> a{};
  bool ok = true;

  static D from_coords(const std::array<long long, Rank>& c) {
    D z;
    z.a = c;
    return z;
  }
  static D integer(long long n) {
    D z;
    z.a[0] = n;
    return z;
  }

  // Coordinate-zero AND unpoisoned: a poisoned value never claims zero.
  bool is_zero() const {
    if (!ok) return false;
    for (long long c : a)
      if (c) return false;
    return true;
  }
  // Coordinate equality (poison compared too: a poisoned value equals
  // nothing, itself included, on the safe side).
  friend bool operator==(const D& x, const D& y) {
    return x.ok && y.ok && x.a == y.a;
  }

  friend D operator+(const D& x, const D& y) {
    D r;
    r.ok = x.ok && y.ok;
    for (int i = 0; i < Rank; i++)
      r.ok &= !__builtin_add_overflow(x.a[i], y.a[i], &r.a[i]);
    return r;
  }
  friend D operator-(const D& x, const D& y) {
    D r;
    r.ok = x.ok && y.ok;
    for (int i = 0; i < Rank; i++)
      r.ok &= !__builtin_sub_overflow(x.a[i], y.a[i], &r.a[i]);
    return r;
  }
  D operator-() const {
    D r;
    r.ok = ok;
    for (int i = 0; i < Rank; i++)
      r.ok &= !__builtin_sub_overflow(0LL, a[i], &r.a[i]);
    return r;
  }
  friend D operator*(long long n, const D& x) {
    D r;
    r.ok = x.ok;
    for (int i = 0; i < Rank; i++)
      r.ok &= !__builtin_mul_overflow(n, x.a[i], &r.a[i]);
    return r;
  }

  // The coefficient magnitude, UNSIGNED: 2^63 (an INT64_MIN coordinate)
  // must compare above kMulOperandMax, so no signed cast may intervene.
  unsigned long long max_abs_coord() const {
    unsigned long long m = 0;
    for (long long c : a) {
      const unsigned long long u =
          c < 0 ? -(unsigned long long)c : (unsigned long long)c;
      if (u > m) m = u;
    }
    return m;
  }

  // Product with the ring's polynomial reduction.  Accumulation in
  // __int128 behind the kMulOperandMax guard; the int64 store-back is
  // checked -- either way out of range poisons, never wraps.
  friend D operator*(const D& x, const D& y) {
    D r;
    r.ok = x.ok && y.ok;
    if (!r.ok) return r;
    if (x.max_abs_coord() > detail::kMulOperandMax ||
        y.max_abs_coord() > detail::kMulOperandMax) {
      r.ok = false;
      return r;
    }
    __int128 conv[2 * Rank - 1] = {};
    for (int i = 0; i < Rank; i++) {
      if (!x.a[i]) continue;
      for (int j = 0; j < Rank; j++)
        conv[i + j] += (__int128)x.a[i] * y.a[j];
    }
    __int128 red[Rank] = {};
    for (int i = 0; i < Rank; i++) red[i] = conv[i];
    for (int k = Rank; k < 2 * Rank - 1; k++) {
      if (!conv[k]) continue;
      const auto& row = D::pow_row(k);
      for (int i = 0; i < Rank; i++) red[i] += conv[k] * row[i];
    }
    for (int i = 0; i < Rank; i++) {
      r.a[i] = (long long)red[i];
      r.ok &= ((__int128)r.a[i] == red[i]);
    }
    return r;
  }
};

// ---------------------------------------------------------------------------
// The working ring Z[gamma].  Every element is real by construction --
// there is no non-real refusal here.
// ---------------------------------------------------------------------------
struct Real30 : CheckedCoeffRing<Real30, 4> {
  static const std::array<long long, 4>& pow_row(int k) {
    return detail::kGamPow[k];
  }

  static Real30 gamma() { return two_cos(1); }
  // 2 cos(k pi/15), any k (folded by the symmetries of the cosine; the
  // exponent convention is zeta_30^k + zeta_30^-k -- each ring's two_cos
  // is in its OWN zeta, so Real30::two_cos(k) == Zeta60::two_cos(2k)).
  static Real30 two_cos(int k) {
    k = ((k % 30) + 30) % 30;
    if (k > 15) k = 30 - k;
    return from_coords(detail::kTwoCos15[k]);
  }
  // The golden ratio: 2 cos(pi/5) = gamma_3; golden^2 == golden + 1.
  static Real30 golden() { return two_cos(3); }

  // ---- the constants table (whitepaper @ref tab:constants): the four
  // ---- geometric inputs of the cubic chain (at the x kLsqScale
  // ---- convention), plus the Heron factor and its inverse certificate.
  // Cubic edges -- and equally hexagon spokes -- squared.
  static Real30 lsq_cubic_edge() { return integer(kLsqScale); }
  // Pentagon spokes, squared: 25 R5^2 = 10 + 5 gamma_3.
  static Real30 lsq_pentagon_spoke() {
    return from_coords({10, -15, 0, 5});
  }
  // Hexagon kis-triangle wedge: 25 (1 + gamma_2 + gamma_4)
  // [= 25 sin(60deg)/sin(12deg); the Heron identity ties it to H = 3*25^2].
  static Real30 wedge_hexagon_kis() {
    return from_coords({0, 100, 25, -25});
  }
  // Pentagon kis-triangle wedge: 5 (2 + gamma_3)(1 + gamma + gamma_3)
  // [= 25 R5^2 sin(72deg)/sin(12deg)].
  static Real30 wedge_pentagon_kis() {
    return from_coords({10, -30, 5, 15});
  }
  // The Heron factor 2 - gamma_2 = 4 - gamma^2: 16 Area^2 =
  // (2 - gamma_2) w^2 for every module triangle (@ref eq:heron).  A ring
  // UNIT -- the certificate (4 - gamma^2) gamma (1 + gamma) = 1
  // (@ref eq:unit) is heron_unit() * heron_unit_inv() == 1, stated exactly
  // by the test suite's [I] group (heron_unit_inv exists to state it).
  static Real30 heron_unit() { return from_coords({4, 0, -1, 0}); }
  static Real30 heron_unit_inv() { return from_coords({0, 1, 1, 0}); }

  // Numeric value (diagnostics and the sign oracle's double rung): Horner
  // at the correctly rounded stored constant -- no cos call.
  double value() const {
    double v, abs_sum;
    value_with_abs(v, abs_sum);
    return v;
  }
  void value_with_abs(double& value, double& abs_sum) const {
    const double g = detail::kGammaDouble;
    value = (((double)a[3] * g + (double)a[2]) * g + (double)a[1]) * g +
            (double)a[0];
    abs_sum = 0;
    for (long long c : a) abs_sum += std::fabs((double)c);
  }
};

// ---------------------------------------------------------------------------
// The guaranteed-unpoisoned input envelopes, derived by exact worst-case
// per-coordinate bound composition through the reduction rows (the
// amplification vector of one product is (11,45,30,22); the 2026-08-25
// review re-derived both bounds to the digit and located the binding
// terms).  Inputs within the envelope NEVER poison; larger inputs refuse
// by name, never lie.
//   - kPredicateCoeffMax binds one whole predicate chain (the validity
//     gate + delaunay/convexity); the binding term is the store-back of
//     heron_unit * w^2, per-coordinate <= 410 * C^2.
//   - kFlipCoeffMax binds the flip numerators (the division itself is
//     mod-p, growth-free; its verification product runs under the operand
//     guard); the binding term is P*Q - heron_unit * wu * wl,
//     per-coordinate <= 815 * C^2.
// ---------------------------------------------------------------------------
inline constexpr long long kPredicateCoeffMax = 149'986'763;
inline constexpr long long kFlipCoeffMax = 106'381'487;
static_assert((__int128)410 * kPredicateCoeffMax * kPredicateCoeffMax <=
                  std::numeric_limits<long long>::max(),
              "predicate binding term (heron_unit * w^2 store-back)");
static_assert((__int128)815 * kFlipCoeffMax * kFlipCoeffMax <=
                  std::numeric_limits<long long>::max(),
              "flip binding term (PQ - unit*wu*wl store-back)");
static_assert((__int128)45 * kPredicateCoeffMax * kPredicateCoeffMax <=
                  (__int128)detail::kMulOperandMax,
              "second-level operands (w*w) must clear the product guard");

// The double rung's relative error bound, derived: 4 int->double
// conversions and the 3+3 Horner operations contribute <= 7u per term path
// (u = 2^-53), against term magnitudes |a_k| gamma^k <= 8 |a_k|, so
// <= 56u * sum|a_k|; the argument error |kGammaDouble - gamma| <= u
// (correctly rounded, gamma in [1,2)) enters through |B'| <= 12 sum|a_k|
// as <= 12u * sum|a_k|.  Total <= 68u ~ 7.6e-15 * sum|a_k| under
// round-to-nearest; FMA contraction only removes rounding steps, and even
// a 4-ulp-perturbed constant stays under ~1.7e-14 (the review measured a
// 4.6M-vector worst of 3.95e-15, incl. LLL near-cancelling families).
// kRung1RelErr = 1e-13 keeps >= 5x margin over the worst reading;
// -ffast-math reassociation would void the analysis (do not build this
// file with it).
inline constexpr double kRung1RelErr = 1e-13;
static_assert(kRung1RelErr >= 5.1e-14,
              "must keep >= 3x margin over the 4-ulp worst reading 1.7e-14");

// ---------------------------------------------------------------------------
// The exact sign oracle.
//
// Ladder: (0) exact zero by coordinates (the power basis is a Z-basis, so
// zero has one representation); (1) double Horner against kRung1RelErr;
// (2) sign_at_isolated_root -- bisect the isolating interval [15/8, 2] of
// gamma (exactly: 2^12 psi(15/8) = -6599 < 0 < psi(2) = 1, and the
// nearest conjugate is 2 cos(7 pi/15) = 0.20905..) with EXACT integer
// arithmetic until |B(mid)| exceeds the derivative bound times the
// interval width.  No transcendental constant enters rung (2).
//
// Termination of rung (2), honestly bounded for FULL int64 coordinates
// (the poison discipline guarantees stored coordinates, nothing smaller):
// B is a nonzero integer polynomial in gamma of degree <= 3 and height
// H <= 2^63; its integer norm is >= 1 and the three conjugate values are
// bounded by (1+2+4+8) H = 15 H, so |B(gamma)| >= (15 H)^-3 > 2^-201.
// The derivative bound is M <= (1+4+12) H < 2^68, so the bisection
// concludes by s ~ 201 + 68 + 1 = 270 bits -- inside kMaxDyadicBits = 384.
// Refused-at-cap is a fail-loud backstop, never a tolerance.
// ---------------------------------------------------------------------------

namespace detail {

// Rung-2 capacities.  The binding intermediate is lhs = Bmid << smid at
// the deepest iteration, ~ (deg+1)*smid + 63 bits: provably 27 of the 28
// limbs here (the ambient tier runs 99 of 100) -- one limb of margin,
// checked at every operation, refusing (BigWidth) rather than wrapping.
// Review-measured extremes (2026-08-25 record, LLL-adversarial inputs at
// |B(gamma)| = 2^-190.3): 252 bisections, 16-limb high water.
inline constexpr int kMaxDyadicBits = 384;
inline constexpr int kLimbs = (4 * kMaxDyadicBits) / 64 + 4;
using Big = FixedBigInt<kLimbs>;

// Exact integer 2^(s*deg) * poly(m/2^s), coefficients B, by Horner.
template <class B>
inline B eval_scaled(const B* coeff, int deg, const B& m, int s) {
  B r = coeff[deg];
  for (int k = deg - 1; k >= 0; k--)
    r = r * m + (coeff[k] << (s * (deg - k)));
  return r;
}

// The shared rung-2 driver: the sign of B(alpha) for alpha the root of
// the monic Psi isolated in [m0, m0+1]/2^s0 -- ONE body for the rank-4
// and the ambient conductor-60 oracles (each supplies coefficients, its
// isolating seed, and its capacities).  In its own frame (the FixedBigInt
// arrays are multi-KB; the fast path must not pay them).  Host tier.
// @pre B not identically zero; Psi(lo) and Psi(hi) of opposite signs with
// exactly one root of Psi inside.  @variant max_bits - s.
template <class B>
[[gnu::noinline]] inline SignOr sign_at_isolated_root(
    const B* Bc, int deg, const B* Psi, int psi_deg, long long m0, int s0,
    int max_bits, SignTrace* tr) {
  const int hdeg = deg > 0 ? deg : 1;   // constant B still Horners at deg 1

  // |B'(x)| on [0,2] is bounded by M = sum k |b_k| 2^{k-1}: magnitudes
  // only, so each term enters with sgn forced positive (no cancellation).
  B M = B{};
  for (int k = 1; k <= deg; k++) {
    B t = Bc[k];
    t.sgn = t.n ? 1 : 0;
    for (int rep = 0; rep < k; rep++) M = M + (t << (k - 1));
  }

  B m = B::from_i128(m0);
  int s = s0;
  auto psi_sign_at = [&](const B& mm, int ss) -> std::optional<int> {
    const B val = eval_scaled(Psi, psi_deg, mm, ss);
    if (val.overflowed()) return std::nullopt;
    return val.sgn;
  };
  const auto lo0 = psi_sign_at(m, s);
  if (!lo0) {
    if (tr) tr->refusal = Refusal::BigWidth;
    return std::nullopt;
  }
  const int lo_sign = *lo0;

  while (s < max_bits) {
    if (tr) tr->bisections++;
    // Conclusive when |B(mid)| > M * width (mid is the parent interval's
    // midpoint, so |alpha - mid| <= 2^-smid): |Bmid|*2^{smid} vs
    // M*2^{smid*hdeg}, both exactly scaled.
    const B mid = (m << 1) + B::from_i128(1);
    const int smid = s + 1;
    const B Bmid = eval_scaled(Bc, hdeg, mid, smid);
    const B lhs = Bmid << smid;
    const B rhs = M << (smid * hdeg);
    if (Bmid.overflowed() || lhs.overflowed() || rhs.overflowed()) {
      if (tr) tr->refusal = Refusal::BigWidth;
      return std::nullopt;
    }
    if (B::cmp_mag(lhs, rhs) > 0) return sign_from_int(Bmid.sgn);

    const auto ms = psi_sign_at(mid, smid);
    if (!ms || *ms == 0) {   // *ms == 0: alpha rational -- corrupt input
      if (tr) tr->refusal = Refusal::BigWidth;
      return std::nullopt;
    }
    m = (*ms == lo_sign) ? mid : (m << 1);
    s = smid;
  }
  if (tr) tr->refusal = Refusal::DyadicCap;
  return std::nullopt;   // fail-loud backstop (see the oracle banner)
}

// Rung 1, the shared sentence: conclusive iff |val| clears the bound.
inline std::optional<Sign> rung1(double val, double abs_sum, double relerr,
                                 SignTrace* tr) {
  const double err = abs_sum * relerr;
  if (val > err) {
    if (tr) tr->rung = 1;
    return Sign::Positive;
  }
  if (val < -err) {
    if (tr) tr->rung = 1;
    return Sign::Negative;
  }
  return std::nullopt;
}

// Rung 2 for the working ring: the element IS an integer polynomial in
// gamma -- no rewrite step.
[[gnu::noinline]] inline SignOr sign_real_exact(const Real30& v,
                                                SignTrace* tr) {
  Big B[4], Psi[5];
  int deg = 0;
  for (int i = 0; i < 4; i++) {
    B[i] = Big::from_i128(v.a[i]);
    if (B[i].sgn) deg = i;
  }
  for (int i = 0; i <= 4; i++) Psi[i] = Big::from_i128(kPsi30[i]);
  return sign_at_isolated_root(B, deg, Psi, 4, /*m0=*/15, /*s0=*/3,
                               kMaxDyadicBits, tr);
}

}  // namespace detail

// Exact sign of a ring element; nullopt = refused by name (see Refusal;
// the trace carries which).  @anchor cyclotomic-sign-real
inline SignOr sign_real(const Real30& v, SignTrace* tr = nullptr) {
  if (!v.ok) {
    if (tr) tr->refusal = Refusal::Poisoned;
    return std::nullopt;
  }
  if (v.is_zero()) {
    if (tr) tr->rung = 0;
    return Sign::Zero;
  }
  double val, abs_sum;
  v.value_with_abs(val, abs_sum);
  if (const auto s = detail::rung1(val, abs_sum, kRung1RelErr, tr)) return *s;
  if (tr) tr->rung = 2;
  return detail::sign_real_exact(v, tr);
}

// ---------------------------------------------------------------------------
// Exact division (the flip divides by 2e).  See the file banner's EXACT
// DIVISION paragraph for the scheme and its guarantees; the trace reports
// the named refusal and, on success, which prime carried the quotient.
// ---------------------------------------------------------------------------

struct DivTrace {
  Refusal refusal = Refusal::None;
  int prime_index = -1;   // 0 or 1 on success
};

namespace detail {

inline constexpr unsigned long long kDivPrimes[2] = {
    (1ULL << 61) - 1,          // 2^61 - 1 (Mersenne; psi splits completely)
    (1ULL << 63) - 25,         // 2^63 - 25 (psi irreducible: never fails
                               // for a divisor nonzero mod it)
};

inline unsigned long long mulmod(unsigned long long x, unsigned long long y,
                                 unsigned long long p) {
  return (unsigned long long)((unsigned __int128)x * y % p);
}
inline unsigned long long powmod(unsigned long long x, unsigned long long e,
                                 unsigned long long p) {
  unsigned long long r = 1;
  while (e) {
    if (e & 1) r = mulmod(r, x, p);
    x = mulmod(x, x, p);
    e >>= 1;
  }
  return r;
}
inline unsigned long long tomod(long long v, unsigned long long p) {
  const long long m = v % (long long)p;
  return (unsigned long long)(m < 0 ? m + (long long)p : m);
}

// Inverse of d in F_p[y]/(psi) by extended Euclid on (psi, d); nullopt if
// gcd(d, psi) != 1 mod p (then p | N(d): the caller tries the other
// prime).  The degree invariant deg(t1) = 4 - deg(r0) makes the t-update
// bound exact: deg(t1) + shift = 4 - d1 <= 3, so no term is ever dropped
// (proven in the 2026-08-25 review).  @variant d1 (strictly decreases per
// swap).
inline std::optional<std::array<unsigned long long, 4>> ring_inverse_mod(
    const std::array<long long, 4>& d, unsigned long long p) {
  // r0 = psi (degree 4), r1 = d mod p; t0 = 0, t1 = 1; invariant
  // t_i * d == r_i  (mod psi, p).
  unsigned long long r0[5], r1[5] = {}, t0[5] = {}, t1[5] = {};
  for (int i = 0; i <= 4; i++) r0[i] = tomod(kPsi30[i], p);
  for (int i = 0; i < 4; i++) r1[i] = tomod(d[i], p);
  t1[0] = 1;
  auto degree = [](const unsigned long long* r) {
    for (int k = 4; k >= 0; k--)
      if (r[k]) return k;
    return -1;
  };
  int d0 = 4, d1 = degree(r1);
  if (d1 < 0) return std::nullopt;   // d == 0 mod p
  while (d1 > 0) {
    const unsigned long long inv_lead = powmod(r1[d1], p - 2, p);
    while (d0 >= d1) {
      const unsigned long long f = mulmod(r0[d0], inv_lead, p);
      const int shift = d0 - d1;
      for (int i = 0; i <= d1; i++)
        r0[i + shift] = (r0[i + shift] + p - mulmod(f, r1[i], p)) % p;
      for (int i = 0; i <= 4 - shift; i++)
        t0[i + shift] = (t0[i + shift] + p - mulmod(f, t1[i], p)) % p;
      const int nd0 = degree(r0);
      if (nd0 == d0) return std::nullopt;   // dead guard: cancellation exact
      d0 = nd0;
      if (d0 < 0) break;
    }
    for (int i = 0; i <= 4; i++) {
      const unsigned long long tr = r0[i], tt = t0[i];
      r0[i] = r1[i]; t0[i] = t1[i];
      r1[i] = tr;    t1[i] = tt;
    }
    const int td = d0;
    d0 = d1;
    d1 = td;
    if (d1 < 0) return std::nullopt;   // gcd has positive degree
  }
  // r1 is a nonzero constant: inverse = t1 / r1[0].
  const unsigned long long ic = powmod(r1[0], p - 2, p);
  std::array<unsigned long long, 4> inv{};
  for (int i = 0; i < 4; i++) inv[i] = mulmod(t1[i], ic, p);
  return inv;
}

// (x * y) mod (psi, p): the SAME conv + kGamPow fold as the checked
// product, over F_p scalars (a guess path -- its correctness is not
// load-bearing, the verification product is).
inline std::array<unsigned long long, 4> mulmod_ring(
    const std::array<unsigned long long, 4>& x,
    const std::array<unsigned long long, 4>& y, unsigned long long p) {
  unsigned long long conv[7] = {};
  for (int i = 0; i < 4; i++)
    for (int j = 0; j < 4; j++)
      conv[i + j] = (conv[i + j] + mulmod(x[i], y[j], p)) % p;
  std::array<unsigned long long, 4> r{};
  for (int i = 0; i < 4; i++) r[i] = conv[i];
  for (int k = 4; k < 7; k++) {
    if (!conv[k]) continue;
    const auto& row = kGamPow[k];
    for (int i = 0; i < 4; i++)
      r[i] = (r[i] + mulmod(conv[k], tomod(row[i], p), p)) % p;
  }
  return r;
}

}  // namespace detail

// numerator / divisor in Z[gamma], exact or refused by name.
// @anchor cyclotomic-exact-div
inline std::optional<Real30> exact_div(const Real30& numerator,
                                       const Real30& divisor,
                                       DivTrace* tr = nullptr) {
  if (!numerator.ok || !divisor.ok) {
    if (tr) tr->refusal = Refusal::Poisoned;
    return std::nullopt;
  }
  if (divisor.is_zero()) {
    if (tr) tr->refusal = Refusal::DivisorZero;
    return std::nullopt;
  }
  for (int pi = 0; pi < 2; pi++) {
    const unsigned long long p = detail::kDivPrimes[pi];
    const auto inv = detail::ring_inverse_mod(divisor.a, p);
    if (!inv) continue;
    std::array<unsigned long long, 4> nm{};
    for (int i = 0; i < 4; i++) nm[i] = detail::tomod(numerator.a[i], p);
    const auto qm = detail::mulmod_ring(nm, *inv, p);
    Real30 q;
    for (int i = 0; i < 4; i++)
      q.a[i] = qm[i] > p / 2 ? (long long)(qm[i] - p) : (long long)qm[i];
    if (divisor * q == numerator) {   // the exact verification
      if (tr) tr->prime_index = pi;
      return q;
    }
  }
  if (tr) tr->refusal = Refusal::NotDivisible;
  return std::nullopt;
}

// ---------------------------------------------------------------------------
// Diamond: the exact wedge-carrying diamond classifier -- the
// CyclotomicMetric's predicate core (the DiamondForms skeleton is shared
// with DiamondSq and Diamond60, diamond_forms.hh; delaunay_geometry.hh
// carries the picture).  Five squared lengths PLUS the two face wedges,
// all in the module's x kLsqScale convention; sqrt(H) = sqrt(2 - gamma_2) w
// makes every verdict the sign of ONE ring element, and the flip
// transports the carry (whitepaper @ref sec:pred, sec:flip).
//
// THE INVARIANT IS THE CONSTRUCTOR'S: Diamond::make refuses (by name)
// unless both faces are positively oriented (wu, wl strictly positive)
// and the carry is CONSISTENT -- H_upper == (2 - gamma_2) wu^2 and
// H_lower == (2 - gamma_2) wl^2, exactly.  The check is TIGHT: the factor
// is a unit in an integral domain, so the two identities pin each wedge
// up to sign and the positivity pins the sign -- a wrong carry cannot
// exist in a constructed Diamond, and the predicates never re-check.
// @inv H_upper() == heron_unit * wu^2  &&  H_lower() == heron_unit * wl^2
//      && wu, wl > 0   (for the life of the object; fields immutable)
//
// The bool convenience forms fold a refusal CONSERVATIVELY -- is_convex /
// is_cocircular to false ("do not act"), and is_delaunay ALSO to false,
// which for a bool-only caller reads as "must flip": exactly DiamondSq's
// convention, kept for consistency, but a mutating caller must consume
// delaunay_form_sign() and trip loudly on nullopt, never the bool.
// (Refusals cannot arise from a constructed Diamond except by later
// poison, which the signs surface.)  @anchor cyclotomic-diamond
// ---------------------------------------------------------------------------
struct Diamond {
 public:
  // The validated constructor: nullopt (with the named reason) on a
  // poisoned input, a non-positive wedge, or an inconsistent carry.
  static std::optional<Diamond> make(const Real30& e, const Real30& a,
                                     const Real30& b, const Real30& c,
                                     const Real30& d, const Real30& wu,
                                     const Real30& wl,
                                     Refusal* why = nullptr) {
    Diamond D{{e, a, b, c, d}, wu, wl};
    if (!(e.ok && a.ok && b.ok && c.ok && d.ok && wu.ok && wl.ok)) {
      if (why) *why = Refusal::Poisoned;
      return std::nullopt;
    }
    const SignOr su = sign_real(wu), sl = sign_real(wl);
    if (!su || !sl || *su != Sign::Positive || *sl != Sign::Positive) {
      if (why) *why = Refusal::NonPositiveWedge;
      return std::nullopt;
    }
    const Real30 u = Real30::heron_unit();
    if (!(D.H_upper() == u * (wu * wu) && D.H_lower() == u * (wl * wl))) {
      if (why) *why = Refusal::InconsistentCarry;
      return std::nullopt;
    }
    return D;
  }

  const Real30& e() const { return f_.e; }
  const Real30& a() const { return f_.a; }
  const Real30& b() const { return f_.b; }
  const Real30& c() const { return f_.c; }
  const Real30& d() const { return f_.d; }
  const Real30& wu() const { return wu_; }
  const Real30& wl() const { return wl_; }

  Real30 s_upper() const { return f_.s_upper(); }
  Real30 s_lower() const { return f_.s_lower(); }
  Real30 P() const { return f_.P(); }
  Real30 Q() const { return f_.Q(); }
  Real30 H_upper() const {
    return metric_forms::heron_product_sq(f_.e, f_.a, f_.b);
  }
  Real30 H_lower() const {
    return metric_forms::heron_product_sq(f_.e, f_.c, f_.d);
  }

  // THE ring elements whose signs are the verdicts (one element each --
  // the module's thesis; certificate-able as such).
  Real30 delaunay_form() const { return s_upper() * wl_ + s_lower() * wu_; }
  Real30 convexity_form_at_origin() const { return Q() * wu_ + P() * wl_; }

  // sign(F), F = s_upper sqrt(H_lower) + s_lower sqrt(H_upper)
  //           = sqrt(2 - gamma_2) * delaunay_form():
  // Positive = strictly locally Delaunay, Zero = cocircular, Negative =
  // must flip; nullopt = refused (poison surfaced after construction).
  SignOr delaunay_form_sign(SignTrace* tr = nullptr) const {
    return sign_real(delaunay_form(), tr);
  }
  bool is_delaunay() const {
    const SignOr s = delaunay_form_sign();
    return s && *s != Sign::Negative;
  }
  bool is_cocircular() const {
    const SignOr s = delaunay_form_sign();
    return s && *s == Sign::Zero;
  }

  // Strict convexity at the diagonal's origin endpoint: sign(Q wu + P wl)
  // (Q pairs with the upper wedge, exactly as DiamondSq pairs Q with
  // tau_upper); the other endpoint is the reversed diamond's test.
  // 2e times the flip's new origin-side wedge IS this form -- the
  // legality test and the transport share one expression.
  SignOr convex_at_origin_sign(SignTrace* tr = nullptr) const {
    return sign_real(convexity_form_at_origin(), tr);
  }
  // The reversal involution (diamond_forms.hh): the same two faces read
  // from the other endpoint, so the carry rides along unchanged and the
  // invariant is preserved verbatim (heron_product_sq is symmetric in the
  // swapped pair) -- a trusted construction.
  Diamond reversed() const { return Diamond{f_.reversed(), wu_, wl_}; }
  bool is_convex() const {
    const SignOr u = convex_at_origin_sign();
    const SignOr v = reversed().convex_at_origin_sign();
    return u && v && *u == Sign::Positive && *v == Sign::Positive;
  }

  // The raw flip carry: the squared diagonal and the two transported
  // wedges,
  //   f^2       = a + c - (P Q - (2 - gamma_2) wu wl) / (2e),
  //   w_origin  = (Q wu + P wl) / (2e),        [the origin-side new face]
  //   w_far     = wu + wl - w_origin,          [the far-side new face]
  // all exact ring arithmetic (whitepaper @ref eq:flip); the divisions
  // are exact by the module geometry, and a divisibility failure refuses
  // by name -- on module inputs that falsifies the caller's premise (a
  // bug, not a border case).  The formulas hold UNCONDITIONALLY (also on
  // non-convex diamonds, where a transported wedge is <= 0); the wedge
  // signs are the two convexity certificates.
  struct Flipped {
    Real30 f2, w_origin, w_far;
  };
  std::optional<Flipped> flipped(DivTrace* tr = nullptr) const {
    const Real30 two_e = 2 * f_.e;
    const auto q1 = exact_div(
        P() * Q() - Real30::heron_unit() * (wu_ * wl_), two_e, tr);
    if (!q1) return std::nullopt;
    const auto wo = exact_div(Q() * wu_ + P() * wl_, two_e, tr);
    if (!wo) return std::nullopt;
    const Real30 f2 = f_.a + f_.c - *q1;
    const Real30 wf = (wu_ + wl_) - *wo;
    if (!f2.ok || !wf.ok) {
      if (tr) tr->refusal = Refusal::Poisoned;
      return std::nullopt;
    }
    return Flipped{f2, *wo, wf};
  }

  // The flip as a WORD: the diamond of the new diagonal A-C, apexes the
  // old endpoints (v upper, u lower), fields by the relabelling theorem
  //   (e,a,b,c,d,wu,wl) -> (f2, b, d, a, c, w_far, w_origin),
  // valid only when this diamond is strictly convex (both transported
  // wedges positive) -- else the "flipped faces" are not faces and the
  // word refuses (NotConvex).  Routed through make(): the transported
  // carry is consistent by theorem, so a make() refusal here is a loud
  // implementation-bug trap, not a border case.
  std::optional<Diamond> flip(Refusal* why = nullptr) const {
    DivTrace dt;
    const auto g = flipped(&dt);
    if (!g) {
      if (why) *why = dt.refusal;
      return std::nullopt;
    }
    const SignOr so = sign_real(g->w_origin), sf = sign_real(g->w_far);
    if (!so || !sf || *so != Sign::Positive || *sf != Sign::Positive) {
      if (why) *why = Refusal::NotConvex;
      return std::nullopt;
    }
    return make(g->f2, f_.b, f_.d, f_.a, f_.c, g->w_far, g->w_origin, why);
  }

 private:
  Diamond(const DiamondForms<Real30>& f, const Real30& wu, const Real30& wl)
      : f_(f), wu_(wu), wl_(wl) {}

  DiamondForms<Real30> f_;
  Real30 wu_, wl_;   // the carried face wedges (upper, lower)
};

}  // namespace cyclotomic
