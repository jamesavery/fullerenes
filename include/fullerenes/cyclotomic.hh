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
// WIDTHS.  The ring is a template on its coefficient type: Real30T<Coef>
// holds coordinates over the power basis 1, gamma, gamma^2, gamma^3
// (reduced by psi(y) = y^4 + y^3 - 4 y^2 - 4 y + 1, the minimal polynomial;
// conjugates 2 cos(k pi/15), k in {1,7,11,13}), and does its arithmetic in
// that width.  Two tiers are used:
//   32/64  -- Real30 = Real30T<long long> arithmetic over a 32-bit STORAGE
//             of the carry (Stored30<int32_t>): the batch tier, every
//             fullerene through ~C15,000 (the envelope below), the same
//             code on the host and in a GPU work-item.  This header spells
//             no integer wider than 64 bits anywhere: the sign oracle
//             works in 32-bit limbs and the modular division in 64-bit
//             folds, so a device compile of this tier sees no 128-bit
//             type at all.
//   64/128 -- Real30Wide, Real30T at 128-bit coefficients, over a 64-bit
//             storage (cyclotomic_wide.hh, host only): single huge isomers
//             beyond the 64-bit envelope, through N of order 10^9.
// widen / narrow convert between a storage and a ring; a value the
// storage width cannot hold refuses (nullopt), and the policy trips by
// name.
//
// OVERFLOW DISCIPLINE -- checked-and-poisoned, never silent: EVERY ring
// operation checks its own coefficient arithmetic and POISONS the value
// (ok = false) instead of wrapping; poison is sticky, and every consumer
// refuses a poisoned value with a NAMED reason (Refusal, carried by
// SignTrace / DivTrace / the Diamond factory) -- a wrong sign from
// overflow is therefore structurally impossible; large inputs refuse
// loudly.  Add, subtract and scale are builtin-checked.  A PRODUCT is
// guarded on its operands' magnitudes: a coordinate of x*y is at most
// kAmplification * max|x| * max|y| (the fold rows' worst case, derived at
// compile time), so the product is accumulated in the coefficient width
// exactly whenever that bound fits.  The guaranteed-unpoisoned input
// envelope per width is carry_coeff_max<Coef>() (kCarryCoeffMax for the
// 64-bit ring), pinned by static_asserts spelling the binding terms.
// Today's kis constants have |c| <= 100.
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
// ONE TIER OF CODE: everything on the run path -- add/sub/mul, the sign
// oracle (sign_real: exact zero by coordinates, else ONE fixed-point
// evaluation at a precision derived from the Liouville bound; its banner
// below), and the exact division (modular products reduced by shifts and
// adds) -- is integer arithmetic over fixed-size arrays with no loop of
// data-dependent length, no allocation, no exception and no floating
// point (the double value() is a SHADOW for the mesh's edge lengths,
// never a decision).  The same code runs in a GPU work-item and on the
// host; it is the cubic chain's device reduction
// (claude-projects/parallel-primitives).
// ============================================================================

#include "diamond_forms.hh"
#include "metric_forms.hh"

#include <array>
#include <cmath>
#include <cstdint>
#include <limits>
#include <optional>
#include <type_traits>

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
  DyadicCap,          // ambient tier: the bisection hit its dyadic cap
  BigWidth,           // ambient tier: the bisection overflowed a limb capacity
  Undecided,          // working ring: the fixed-point evaluation of a NONZERO
                      // element returned zero -- impossible by the Liouville
                      // bound the precision is derived from; the constants
                      // are corrupt (a fail-loud backstop, never a tolerance)
  DivisorZero,        // exact_div by zero
  NotDivisible,       // exact_div: no quotient verified (either prime)
  NonPositiveWedge,   // Diamond: a face wedge not strictly positive
  InconsistentCarry,  // Diamond: Heron identity fails against the carry
  NotConvex,          // Diamond::flip on a non-convex diamond
};

// Optional observability of the sign oracle (the counting-metric pattern):
// which step decided (0 = zero by coordinates, 1 = the fixed-point
// evaluation; 2 = the ambient ring's bisection), how many bisections the
// ambient ring spent, and -- on a refusal -- its name.
struct SignTrace {
  int rung = -1;
  int bisections = 0;
  Refusal refusal = Refusal::None;
#ifdef FULLERENES_SIGN_FILTER_CHECK
  // The differential's per-isomer counters, when one is being run (see
  // SignFilterCounters below).  A raw pointer into the caller's storage:
  // device-legal, nothing owned, null when no check is running.
  struct SignFilterCounters* ctr = nullptr;
#endif
};

#ifdef FULLERENES_SIGN_FILTER_CHECK
// ── THE FLOATING-POINT FILTER DIFFERENTIAL (instrument, compiled in only
//    under this macro; production and device builds carry none of it).
//
//    A ladder in front of the exact oracle would evaluate B(gamma) in
//    single precision and return the sign whenever |value| exceeds a
//    CERTIFIED bound on that evaluation's error, falling through to double
//    and then to the oracle otherwise -- never a tolerance, always a
//    proof.  This counts what each rung would claim AND checks every claim
//    against the sign the oracle actually returned: a rung may decline
//    (the fall-through the ladder is built around), but a rung that
//    answers must answer correctly.  disagree32 / disagree64 must stay 0;
//    one of them non-zero falsifies the bound.
//
//    THE BOUND.  With H = max|b_k| (so H >= 1 for B != 0), gamma < 2 and
//    the power basis of degree 3:
//      - argument:    |B'| <= 17 H on [0,2] and |gamma~ - gamma| <= 2.01 u
//                     (the stored constant is rounded to double and then to
//                     the working precision -- a double rounding), so
//                     <= 34.2 H u;
//      - coefficients: fl(b_k) = b_k (1 + delta), so <= 15 H u;
//      - Horner:      Higham Thm 5.1, gamma_6 * sum |b_k| gamma^k
//                     <= 6.01 u * 15 H <= 90.2 H u.
//    Total <= 139.4 H u = 69.7 H eps; the constant below is 256, a factor
//    of ~3.7 of slack.  No cancellation is assumed anywhere -- every term
//    is coefficient-wise -- which is what makes the bound sound however
//    badly B(gamma) cancels.
//
//    SIDE CONDITION.  The derivation needs 15 H representable in the
//    working precision.  At the batch tier (Stored30<int32_t> over Real30,
//    H <= 2^63) 15 H ~ 1.4e20 against float's 3.4e38 -- safe by eighteen
//    orders.  At the 128-bit tier H can reach 2^127 and 15 H OVERFLOWS, so
//    a non-finite evaluation is never trusted: it falls through instead.
//
//    INCOMPLETE BY CONSTRUCTION.  A nonzero B has |B(gamma)| >= (15H)^-3
//    (Liouville), which is far below either float threshold, so ring
//    elements provably exist that no float rung can decide.  The oracle is
//    not vestigial and the fall-through is not padding. ──
struct SignFilterCounters {
  long long zero = 0;          // decided by coordinates, before any rung
  long long claim32 = 0, agree32 = 0, disagree32 = 0;
  long long claim64 = 0, agree64 = 0, disagree64 = 0;
  long long fell_through = 0;  // reached the oracle with no rung deciding
  long long nonfinite = 0;     // the side condition refused a rung
  // The first disagreement, kept whole so it can be reproduced.
  int  witness_rung = 0, witness_filter = 0, witness_oracle = 0;
  double witness_value = 0, witness_bound = 0, witness_coef[4] = {0, 0, 0, 0};
};
#define FULLERENES_BIND_SIGN_CTR(trace, counters) ((trace).ctr = (counters))
#define FULLERENES_BIND_SIGN_CTR_STORE(policy, counters) ((policy).sign_ctr = (counters))
#else
#define FULLERENES_BIND_SIGN_CTR(trace, counters) ((void)0)
#define FULLERENES_BIND_SIGN_CTR_STORE(policy, counters) ((void)0)
#endif

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

// The per-coordinate AMPLIFICATION of a product, from the fold rows:
// |(x*y)_i| <= A_i max|x| max|y| with A_i the sum, over the convolution
// slots s = j + k, of (the number of pairs (j,k) landing in s) times (the
// magnitude with which slot s folds into coordinate i: 1 for s = i, the
// row kGamPow[s][i] for s >= Rank).  For a FIXED left factor c the pair
// count becomes sum_{j+k=s} |c_j|.  Both are compile-time integers.
constexpr int amplification_of_product() {
  int worst = 0;
  for (int i = 0; i < 4; i++) {
    long long tot = 0;
    for (int s = 0; s < 7; s++) {
      const int pairs = (s < 4 ? s : 6 - s) + 1;
      const long long fold = s < 4 ? (s == i ? 1 : 0) : (kGamPow[s][i] < 0 ? -kGamPow[s][i] : kGamPow[s][i]);
      tot += pairs * fold;
    }
    if (tot > worst) worst = (int)tot;
  }
  return worst;
}
constexpr int amplification_of_factor(const std::array<long long, 4>& c) {
  int worst = 0;
  for (int i = 0; i < 4; i++) {
    long long tot = 0;
    for (int s = 0; s < 7; s++) {
      long long pairs = 0;
      for (int j = 0; j < 4; j++)
        if (s - j >= 0 && s - j < 4) pairs += c[j] < 0 ? -c[j] : c[j];
      const long long fold = s < 4 ? (s == i ? 1 : 0) : (kGamPow[s][i] < 0 ? -kGamPow[s][i] : kGamPow[s][i]);
      tot += pairs * fold;
    }
    if (tot > worst) worst = (int)tot;
  }
  return worst;
}
inline constexpr int kProductAmplification = amplification_of_product();
inline constexpr int kHeronUnitAmplification = amplification_of_factor({4, 0, -1, 0});
inline constexpr int kGammaAmplification = amplification_of_factor({0, 1, 0, 0});
static_assert(kProductAmplification == 45, "per-coordinate amplification (11,45,30,22)");
static_assert(kHeronUnitAmplification == 13 && kGammaAmplification == 5);

// gamma as a correctly rounded double (verified against an exact rational
// Newton iteration on psi and by the test suite's long-double cross-check).
// 1.9562952014676114 round-trips this exact bit pattern.
inline constexpr double kGammaDouble = 0x1.f4cfc327a0080p+0;

// Integer square root (floor) of an unsigned value of any width, by bits:
// the envelope constants are floor(sqrt(MAX / amplification)).
template <class U>
constexpr U isqrt(U n) {
  U r = 0;
  for (int bit = (int)(4 * sizeof(U)) - 1; bit >= 0; bit--) {
    const U t = r | ((U)1 << bit);
    if (t * t <= n) r = t;
  }
  return r;
}

}  // namespace detail

// ---------------------------------------------------------------------------
// CheckedCoeffRing: the checked-and-poisoned coefficient ring over the
// power basis 1, x, .., x^(Rank-1) at coefficient type Coef -- the ONE
// spelling of the scaffold the working rings (Real30T) and the ambient
// ring (Zeta60) share.  `ok == false` marks a value whose coefficients
// overflowed somewhere in its history; every consumer refuses it by name.
// The derived ring supplies its reduced-power table via pow_row(k)
// (indices Rank .. 2 Rank - 2 fold a product; the ambient conjugation
// reads the full table) and its per-coordinate product amplification
// kAmplification (the fold rows' worst case).
// ---------------------------------------------------------------------------
template <class D, int Rank, class Coef>
struct CheckedCoeffRing {
  static constexpr int kRank = Rank;
  using coef_type = Coef;
  using UCoef = std::make_unsigned_t<Coef>;
  static constexpr int kCoefBits = 8 * (int)sizeof(Coef);

  std::array<Coef, Rank> a{};
  bool ok = true;

  static D from_coords(const std::array<Coef, Rank>& c) {
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
    for (const Coef& c : a)
      if (c) return false;
    return true;
  }
  // Coordinate equality (poison compared too: a poisoned value equals
  // nothing, itself included, on the safe side).  Spelled as a fold without
  // an early exit: an element-wise compare loop with early exit is a memory
  // comparison idiom the optimizer replaces by a library call (bcmp), which
  // a GPU kernel has no library to resolve.
  friend bool operator==(const D& x, const D& y) {
    bool eq = x.ok && y.ok;
    for (int i = 0; i < Rank; i++) eq &= (x.a[i] == y.a[i]);
    return eq;
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
      r.ok &= !__builtin_sub_overflow((Coef)0, a[i], &r.a[i]);
    return r;
  }
  friend D operator*(long long n, const D& x) {
    D r;
    r.ok = x.ok;
    for (int i = 0; i < Rank; i++)
      r.ok &= !__builtin_mul_overflow((Coef)n, x.a[i], &r.a[i]);
    return r;
  }

  // The coefficient magnitude, UNSIGNED: the most negative coordinate's
  // magnitude is 2^(bits-1) and must enter the product guard as such, so
  // no signed cast may intervene (the 2026-08-25 review exhibited a
  // signed comparison letting it through into undefined behaviour).
  UCoef max_abs_coord() const {
    UCoef m = 0;
    for (const Coef& c : a) {
      const UCoef u = c < 0 ? (UCoef)0 - (UCoef)c : (UCoef)c;
      if (u > m) m = u;
    }
    return m;
  }

  // The product guard: every coordinate of x*y is bounded by
  // kAmplification * max|x| * max|y|, so the accumulation below is exact
  // in Coef iff that bound fits it -- decided on the magnitudes BEFORE any
  // arithmetic, the magnitude product itself overflow-checked.
  static bool product_fits(UCoef mx, UCoef my) {
    constexpr UCoef kCap = ((((UCoef)1) << (kCoefBits - 1)) - 1) / (UCoef)D::kAmplification;
    UCoef p;
    return !__builtin_mul_overflow(mx, my, &p) && p <= kCap;
  }

  // Product with the ring's polynomial reduction, accumulated in Coef
  // behind the guard (every partial sum is bounded by the sum of the
  // magnitudes it folds, which the guard bounds) -- out of range poisons,
  // never wraps.
  friend D operator*(const D& x, const D& y) {
    D r;
    r.ok = x.ok && y.ok && product_fits(x.max_abs_coord(), y.max_abs_coord());
    if (!r.ok) return r;
    Coef conv[2 * Rank - 1] = {};
    for (int i = 0; i < Rank; i++) {
      if (!x.a[i]) continue;
      for (int j = 0; j < Rank; j++)
        conv[i + j] += x.a[i] * y.a[j];
    }
    for (int i = 0; i < Rank; i++) r.a[i] = conv[i];
    for (int k = Rank; k < 2 * Rank - 1; k++) {
      if (!conv[k]) continue;
      const auto& row = D::pow_row(k);
      for (int i = 0; i < Rank; i++) r.a[i] += conv[k] * (Coef)row[i];
    }
    return r;
  }
};

// ---------------------------------------------------------------------------
// The working ring Z[gamma] at coefficient type Coef.  Every element is
// real by construction -- there is no non-real refusal here.
// ---------------------------------------------------------------------------
template <class Coef>
struct Real30T : CheckedCoeffRing<Real30T<Coef>, 4, Coef> {
  using Base = CheckedCoeffRing<Real30T<Coef>, 4, Coef>;
  using Base::from_coords;
  using Base::integer;
  static constexpr int kAmplification = detail::kProductAmplification;
  static const std::array<long long, 4>& pow_row(int k) {
    return detail::kGamPow[k];
  }
  static Real30T from_rows(const std::array<long long, 4>& c) {
    Real30T z;
    for (int i = 0; i < 4; i++) z.a[i] = (Coef)c[i];
    return z;
  }

  static Real30T gamma() { return two_cos(1); }
  // 2 cos(k pi/15), any k (folded by the symmetries of the cosine; the
  // exponent convention is zeta_30^k + zeta_30^-k -- each ring's two_cos
  // is in its OWN zeta, so Real30::two_cos(k) == Zeta60::two_cos(2k)).
  static Real30T two_cos(int k) {
    k = ((k % 30) + 30) % 30;
    if (k > 15) k = 30 - k;
    return from_rows(detail::kTwoCos15[k]);
  }
  // The golden ratio: 2 cos(pi/5) = gamma_3; golden^2 == golden + 1.
  static Real30T golden() { return two_cos(3); }

  // ---- the constants table (whitepaper @ref tab:constants): the four
  // ---- geometric inputs of the cubic chain (at the x kLsqScale
  // ---- convention), plus the Heron factor and its inverse certificate.
  // Cubic edges -- and equally hexagon spokes -- squared.
  static Real30T lsq_cubic_edge() { return integer(kLsqScale); }
  // Pentagon spokes, squared: 25 R5^2 = 10 + 5 gamma_3.
  static Real30T lsq_pentagon_spoke() { return from_rows({10, -15, 0, 5}); }
  // Hexagon kis-triangle wedge: 25 (1 + gamma_2 + gamma_4)
  // [= 25 sin(60deg)/sin(12deg); the Heron identity ties it to H = 3*25^2].
  static Real30T wedge_hexagon_kis() { return from_rows({0, 100, 25, -25}); }
  // Pentagon kis-triangle wedge: 5 (2 + gamma_3)(1 + gamma + gamma_3)
  // [= 25 R5^2 sin(72deg)/sin(12deg)].
  static Real30T wedge_pentagon_kis() { return from_rows({10, -30, 5, 15}); }
  // The Heron factor 2 - gamma_2 = 4 - gamma^2: 16 Area^2 =
  // (2 - gamma_2) w^2 for every module triangle (@ref eq:heron).  A ring
  // UNIT -- the certificate (4 - gamma^2) gamma (1 + gamma) = 1
  // (@ref eq:unit) is heron_unit() * heron_unit_inv() == 1, stated exactly
  // by the test suite's [I] group (heron_unit_inv exists to state it).
  static Real30T heron_unit() { return from_rows({4, 0, -1, 0}); }
  static Real30T heron_unit_inv() { return from_rows({0, 1, 1, 0}); }

  // Numeric value (the float SHADOW of an exact quantity -- the mesh's
  // double edge lengths -- and diagnostics): Horner at the correctly
  // rounded stored constant, no cos call.  Never a decision: every sign is
  // the exact oracle's (sign_real).
  double value() const {
    const double g = detail::kGammaDouble;
    return (((double)this->a[3] * g + (double)this->a[2]) * g + (double)this->a[1]) * g +
           (double)this->a[0];
  }
};

using Real30 = Real30T<long long>;

// ---------------------------------------------------------------------------
// STORAGE.  A carry may be stored at a width other than its ring's:
// Stored30<Coef> holds the power-basis coordinates alone (a stored value
// is never poisoned: narrow refuses poison), and a ring value can itself
// serve as the storage of a wider ring (the 64-bit Real30 under the
// 128-bit arithmetic of cyclotomic_wide.hh).  widen<R>() is exact;
// narrow<S>() refuses (nullopt) a coordinate the width cannot hold and a
// poisoned value, and the policy turns that refusal into a named trip: an
// overflow of the storage width refuses, never wraps.  Every consumer is
// written once over the storage type S and the ring R.
// ---------------------------------------------------------------------------
template <class Coef>
struct Stored30 {
  using coef_type = Coef;
  std::array<Coef, 4> a{};
  friend bool operator==(const Stored30& x, const Stored30& y) {
    bool eq = true;
    for (int i = 0; i < 4; i++) eq &= (x.a[i] == y.a[i]);
    return eq;
  }
};

template <class R, class Coef>
inline R widen(const Stored30<Coef>& s) {
  R r;
  for (int i = 0; i < 4; i++) r.a[i] = (typename R::coef_type)s.a[i];
  return r;
}
template <class R, class Coef>
inline R widen(const Real30T<Coef>& s) {
  R r;
  r.ok = s.ok;
  for (int i = 0; i < 4; i++) r.a[i] = (typename R::coef_type)s.a[i];
  return r;
}

template <class S, class Coef>
inline std::optional<S> narrow(const Real30T<Coef>& v) {
  if (!v.ok) return std::nullopt;
  using SC = typename S::coef_type;
  S s;
  for (int i = 0; i < 4; i++) {
    if constexpr (sizeof(SC) < sizeof(Coef)) {
      if (v.a[i] < (Coef)std::numeric_limits<SC>::min() ||
          v.a[i] > (Coef)std::numeric_limits<SC>::max())
        return std::nullopt;
    }
    s.a[i] = (SC)v.a[i];
  }
  return s;
}

// ---------------------------------------------------------------------------
// The guaranteed-unpoisoned input envelope of a width: inputs -- squared
// lengths and wedges -- of coordinate magnitude <= carry_coeff_max<Coef>()
// pass every product guard and every store of the predicate chain (the
// validity gate, delaunay / convexity) and of the flip.  Along the chain
// the operands are C-sized (sides, wedges), 3C-sized (the law-of-cosines
// forms) or kHeronUnitAmplification*C-sized (the Heron unit times ONE
// wedge: the consistency check and the flip's cross term are associated
// as (unit * w) * w on purpose, so that no operand grows quadratically),
// and the largest stored sum is the Heron form, nine products.  Binding:
//   the guard of (heron_unit * w) * w,   45 * 13 * C^2 <= MAX,
// i.e. C <= 125 564 516 (2^26.9) at 64 bits; the Heron store and the
// guard of P * Q both allow 405 C^2.  Inputs within the envelope NEVER
// poison; larger inputs refuse by name, never lie.  (Exact division is
// mod-p, growth-free; its verification product is guarded like any
// product, so a wrong lift refuses instead of verifying.)
// ---------------------------------------------------------------------------
template <class Coef>
constexpr Coef carry_coeff_max() {
  using U = std::make_unsigned_t<Coef>;
  constexpr U kMax = (((U)1) << (8 * sizeof(Coef) - 1)) - 1;
  return (Coef)detail::isqrt<U>(kMax / (U)(detail::kProductAmplification *
                                           detail::kHeronUnitAmplification));
}
template <class Coef>
constexpr bool envelope_holds() {
  using U = std::make_unsigned_t<Coef>;
  constexpr U kMax = (((U)1) << (8 * sizeof(Coef) - 1)) - 1;
  constexpr U C = (U)carry_coeff_max<Coef>();
  constexpr U bind = (U)(detail::kProductAmplification * detail::kHeronUnitAmplification);
  constexpr U sums = (U)(9 * detail::kProductAmplification);
  return C * C <= kMax / bind &&              // the binding guard
         (C + 1) * (C + 1) > kMax / bind &&   // ...and C is the largest it admits
         C * C <= kMax / sums;                // the Heron store, and P * Q at 3C
}
inline constexpr long long kCarryCoeffMax = carry_coeff_max<long long>();
static_assert(kCarryCoeffMax == 125'564'516, "the 64-bit envelope");
static_assert(envelope_holds<long long>(), "the 64-bit envelope's binding terms");

// ---------------------------------------------------------------------------
// The exact sign oracle of the working ring: ONE fixed-point evaluation.
//
// Every verdict is the sign of B(gamma) for an integer polynomial B of
// degree <= 3 in gamma = 2 cos(pi/15) with coefficients b_k of height
// H = max |b_k| <= 2^(bits-1), bits the coefficient width (the poison
// discipline guarantees nothing smaller).  Step (0): exact zero by
// coordinates (the power basis is a Z-basis, so zero has one
// representation).  Step (1), for B != 0: with m the integer
// floor(gamma 2^p) the sign of the EXACT integer
//     S = 2^{3p} B(m / 2^p) = sum_k b_k * m^k * 2^{p (3 - k)}
// is the sign of B(gamma), because the two are closer than any nonzero
// value can be to zero:
//   * |B(gamma) - B(m/2^p)| <= max |B'| on [0, 2] * |gamma - m/2^p|
//                            <= (1 + 4 + 12) H * 2^-p = 17 H 2^-p,
//     H = max |b_k| (a coefficient-wise bound, no cancellation assumed);
//   * a nonzero B has integer norm >= 1 and three conjugate values each
//     bounded by (1 + 2 + 4 + 8) H = 15 H, so |B(gamma)| >= (15 H)^-3
//     (Liouville);
//   * hence sign S = sign B(gamma) whenever 2^p > 2 * 17 * 15^3 * H^4
//     = 114 750 H^4 < 2^16.81 H^4, i.e. p >= 4 (bits - 1) + 17: p >= 269
//     at 64 bits (the module takes 320), p >= 525 at 128 bits (576).
// The constants m^k 2^{p(3-k)} are compile-time integers; m itself is
// derived by exact bisection of psi offline
// (claude-projects/delaunay/tools/derive_gamma_fixed.py) and VERIFIED
// here at compile time by psi(m/2^p) < 0 < psi((m+1)/2^p), evaluated
// exactly.  The evaluation is a fixed number of multiply-accumulates of
// 32-bit coefficient chunks into a fixed-width integer of 32-bit limbs
// -- no loop with a data-dependent trip count, no branch on the data
// beyond the final compare, no floating point, no allocation, no
// integer wider than 64 bits: the same code on a GPU work-item and on
// the host.  Each width supplies its parameters through SignParams<Coef>
// (the 128-bit instantiation lives in cyclotomic_wide.hh).
// ---------------------------------------------------------------------------

namespace detail {

// Fixed-width unsigned integers, little-endian 32-bit limbs with 64-bit
// carries, for the evaluation and its compile-time verification.  Sizes
// are chosen so that no operation below can overflow (each is stated at
// its use).
template <int N>
struct FixedU {
  uint32_t l[N] = {};

  constexpr bool is_zero() const {
    bool z = true;
    for (int i = 0; i < N; i++) z &= (l[i] == 0);
    return z;
  }
  static constexpr FixedU one() { FixedU r; r.l[0] = 1; return r; }
  // Three-way order of magnitudes: -1, 0, +1.
  static constexpr int cmp(const FixedU& x, const FixedU& y) {
    for (int i = N - 1; i >= 0; i--)
      if (x.l[i] != y.l[i]) return x.l[i] < y.l[i] ? -1 : 1;
    return 0;
  }
  friend constexpr FixedU operator+(const FixedU& x, const FixedU& y) {
    FixedU r;
    uint64_t c = 0;
    for (int i = 0; i < N; i++) {
      c += (uint64_t)x.l[i] + y.l[i];
      r.l[i] = (uint32_t)c;
      c >>= 32;
    }
    return r;
  }
  // x << bits, bits a multiple of 32 or not; limbs shifted out are lost
  // (callers size N so none are).
  constexpr FixedU shl(int bits) const {
    FixedU r;
    const int w = bits / 32, s = bits % 32;
    for (int i = N - 1; i >= 0; i--) {
      const int j = i - w;
      if (j < 0) break;
      uint32_t v = l[j] << s;
      if (s && j > 0) v |= l[j - 1] >> (32 - s);
      r.l[i] = v;
    }
    return r;
  }
  // Schoolbook product truncated at N limbs (callers size N so the true
  // product fits).  The carry chain never exceeds 64 bits: a 32x32
  // product plus two 32-bit addends is at most 2^64 - 1.
  friend constexpr FixedU operator*(const FixedU& x, const FixedU& y) {
    FixedU r;
    for (int i = 0; i < N; i++) {
      if (!x.l[i]) continue;
      uint64_t c = 0;
      for (int j = 0; i + j < N; j++) {
        c += (uint64_t)x.l[i] * y.l[j] + r.l[i + j];
        r.l[i + j] = (uint32_t)c;
        c >>= 32;
      }
    }
    return r;
  }
  constexpr FixedU times_small(uint32_t k) const {
    FixedU r;
    uint64_t c = 0;
    for (int i = 0; i < N; i++) {
      c += (uint64_t)l[i] * k;
      r.l[i] = (uint32_t)c;
      c >>= 32;
    }
    return r;
  }
  // this += (v * k) << (32 offset), for a constant v of M limbs and a
  // 32-bit k: the RUN-PATH operation of the oracle (one fixed pass per
  // coefficient chunk).  Callers size N >= M + offset + 1.
  template <int M>
  constexpr void mul_add(const FixedU<M>& v, uint32_t k, int offset) {
    uint64_t c = 0;
    for (int i = 0; i < M; i++) {
      c += (uint64_t)v.l[i] * k + l[i + offset];
      l[i + offset] = (uint32_t)c;
      c >>= 32;
    }
    for (int i = M + offset; i < N && c; i++) {
      c += l[i];
      l[i] = (uint32_t)c;
      c >>= 32;
    }
  }
};

// The oracle's parameters per coefficient width: the precision p
// (kBits), the limb counts it implies, and the fixed-point gamma as
// 64-bit words (little-endian; the derivation script prints them).
template <class Coef>
struct SignParams;

template <>
struct SignParams<long long> {
  static constexpr int kBits = 320;                  // p
  static constexpr int kConstLimbs = 31;             // m^k 2^{p(3-k)} < 2^{3p+3} = 2^963
  static constexpr int kAccLimbs = 33;               // + 64-bit coefficients, four terms: < 2^1029
  static constexpr int kVerifyLimbs = 41;            // psi at scale 2^{4p}: < 2^1285
  static constexpr std::array<uint64_t, 6> kGammaWords = {
      0x45effbeef3b55e15ULL, 0x10c1517b2ab04fe0ULL, 0x5af6b5d8ea29e11eULL,
      0xf5be43c6e4340270ULL, 0xf4cfc327a007f8a9ULL, 0x0000000000000001ULL};
};

// The precision suffices for every representable coefficient, and the
// limb counts hold every intermediate (with room for the chunked
// coefficient in the accumulator).
template <class Coef>
constexpr bool sign_params_sound() {
  using P = SignParams<Coef>;
  constexpr int bits = 8 * (int)sizeof(Coef);
  return P::kBits >= 4 * (bits - 1) + 17 &&
         P::kConstLimbs * 32 >= 3 * P::kBits + 3 &&
         P::kAccLimbs * 32 >= 3 * P::kBits + 3 + bits + 2 &&
         P::kAccLimbs >= P::kConstLimbs + bits / 32 &&   // mul_add at offsets 0 .. bits/32 - 1
         P::kVerifyLimbs * 32 >= 4 * P::kBits + 5 &&
         P::kConstLimbs >= 2 * (int)P::kGammaWords.size();
}

template <class Coef, int N>
constexpr FixedU<N> gamma_fixed() {
  FixedU<N> m;
  const auto& w = SignParams<Coef>::kGammaWords;
  for (std::size_t i = 0; i < w.size(); i++) {
    m.l[2 * i] = (uint32_t)w[i];
    m.l[2 * i + 1] = (uint32_t)(w[i] >> 32);
  }
  return m;
}

// sign of psi(M / 2^p) at scale 2^{4p}, exactly: psi(y) = y^4 + y^3 - 4y^2
// - 4y + 1 (kPsi30), so 2^{4p} psi(M/2^p) = M^4 + M^3 2^p + 2^{4p}
// - 4 M^2 2^{2p} - 4 M 2^{3p}; compared as positive part vs negative part.
template <class Coef>
constexpr int psi_sign_at(const FixedU<SignParams<Coef>::kVerifyLimbs>& M) {
  using P = SignParams<Coef>;
  using F = FixedU<P::kVerifyLimbs>;
  const F M2 = M * M, M3 = M2 * M, M4 = M3 * M;
  const F pos = M4 + M3.shl(P::kBits) + F::one().shl(4 * P::kBits);
  const F neg = M2.shl(2 * P::kBits).times_small(4) + M.shl(3 * P::kBits).times_small(4);
  return F::cmp(pos, neg);
}
// The isolating property of m, verified at compile time (exact).
template <class Coef>
constexpr bool gamma_fixed_isolates() {
  using P = SignParams<Coef>;
  using F = FixedU<P::kVerifyLimbs>;
  const F m = gamma_fixed<Coef, P::kVerifyLimbs>();
  return psi_sign_at<Coef>(m) < 0 && psi_sign_at<Coef>(m + F::one()) > 0;
}

// The evaluation constants m^k 2^{p(3-k)}, k = 0..3.
template <class Coef>
constexpr std::array<FixedU<SignParams<Coef>::kConstLimbs>, 4> make_sign_powers() {
  using P = SignParams<Coef>;
  using F = FixedU<P::kConstLimbs>;
  const F m = gamma_fixed<Coef, P::kConstLimbs>();
  std::array<F, 4> Pw{};
  Pw[0] = F::one().shl(3 * P::kBits);
  Pw[1] = m.shl(2 * P::kBits);
  Pw[2] = (m * m).shl(P::kBits);
  Pw[3] = m * m * m;
  return Pw;
}
template <class Coef>
inline constexpr std::array<FixedU<SignParams<Coef>::kConstLimbs>, 4> kSignPowers =
    make_sign_powers<Coef>();

static_assert(sign_params_sound<long long>(), "sign oracle, 64-bit: parameters");
static_assert(gamma_fixed_isolates<long long>(),
              "sign oracle, 64-bit: psi(m / 2^p) < 0 < psi((m + 1) / 2^p) must hold");

// sign of S = sum_k b_k m^k 2^{p(3-k)}: the positive and the negative terms
// accumulated apart (magnitudes masked by the coefficient's sign, so both
// passes run on every lane, one per 32-bit chunk of the coefficient),
// then compared.  ONE out-of-line body per width: every predicate calls
// it, and inlining its multiply-accumulates at each of the ~30 call sites
// of a reduction kernel multiplied the device compile time, not the run
// time.
template <class Coef>
[[gnu::noinline]] inline Sign sign_of_fixed_eval(const std::array<Coef, 4>& b) {
  using P = SignParams<Coef>;
  using U = std::make_unsigned_t<Coef>;
  constexpr int kChunks = (int)sizeof(Coef) / 4;
  FixedU<P::kAccLimbs> pos, neg;
  for (int k = 0; k < 4; k++) {
    const uint32_t neg_mask = (uint32_t)0 - (uint32_t)(b[k] < 0);
    const U mag = b[k] < 0 ? (U)0 - (U)b[k] : (U)b[k];
    for (int j = 0; j < kChunks; j++) {
      const uint32_t chunk = (uint32_t)(mag >> (32 * j));
      pos.mul_add(kSignPowers<Coef>[k], chunk & ~neg_mask, j);
      neg.mul_add(kSignPowers<Coef>[k], chunk & neg_mask, j);
    }
  }
  return sign_from_int(FixedU<P::kAccLimbs>::cmp(pos, neg));
}

#ifdef FULLERENES_SIGN_FILTER_CHECK
// One rung: evaluate B at the working precision by Horner and claim the
// sign iff |value| clears the certified bound (SignFilterCounters' banner
// carries the derivation).  A non-finite value means the side condition
// (15 H representable) failed -- claim nothing.
template <class F, class Coef>
inline bool filter_rung(const Real30T<Coef>& v, double H, double& value, double& bound) {
  const F g = (F)kGammaDouble;
  const F x = (((F)v.a[3] * g + (F)v.a[2]) * g + (F)v.a[1]) * g + (F)v.a[0];
  value = (double)x;
  bound = 256.0 * H * (double)std::numeric_limits<F>::epsilon();
  return std::isfinite(value) && std::isfinite(bound) && std::fabs(value) > bound;
}

inline void keep_filter_witness(SignFilterCounters& c, int rung, int filter, int oracle,
                                double value, double bound, const double coef[4]) {
  if (c.witness_rung) return;
  c.witness_rung = rung; c.witness_filter = filter; c.witness_oracle = oracle;
  c.witness_value = value; c.witness_bound = bound;
  for (int k = 0; k < 4; k++) c.witness_coef[k] = coef[k];
}

// The ladder as a ladder: single first, double only on what single
// declined, and every claim checked against the oracle's sign.
template <class Coef>
inline void check_sign_filter(const Real30T<Coef>& v, Sign oracle, SignFilterCounters& c) {
  double H = 0.0, coef[4];
  for (int k = 0; k < 4; k++) {
    coef[k] = (double)v.a[k];
    const double m = std::fabs(coef[k]);
    if (m > H) H = m;
  }
  const int o = (int)oracle;
  double value = 0, bound = 0;
  if (filter_rung<float>(v, H, value, bound)) {
    c.claim32++;
    const int f = value > 0 ? 1 : -1;
    if (f == o) c.agree32++;
    else { c.disagree32++; keep_filter_witness(c, 1, f, o, value, bound, coef); }
    return;
  }
  if (!std::isfinite(value)) c.nonfinite++;
  if (filter_rung<double>(v, H, value, bound)) {
    c.claim64++;
    const int f = value > 0 ? 1 : -1;
    if (f == o) c.agree64++;
    else { c.disagree64++; keep_filter_witness(c, 2, f, o, value, bound, coef); }
    return;
  }
  c.fell_through++;
}
#endif

}  // namespace detail

// Exact sign of a ring element; nullopt = refused by name (see Refusal;
// the trace carries which, and the step that decided: 0 = zero by
// coordinates, 1 = the fixed-point evaluation).  @anchor cyclotomic-sign-real
template <class Coef>
inline SignOr sign_real(const Real30T<Coef>& v, SignTrace* tr = nullptr) {
  if (!v.ok) {
    if (tr) tr->refusal = Refusal::Poisoned;
    return std::nullopt;
  }
  if (v.is_zero()) {
    if (tr) tr->rung = 0;
#ifdef FULLERENES_SIGN_FILTER_CHECK
    if (tr && tr->ctr) tr->ctr->zero++;
#endif
    return Sign::Zero;
  }
  const Sign s = detail::sign_of_fixed_eval<Coef>(v.a);
  if (tr) tr->rung = 1;
#ifdef FULLERENES_SIGN_FILTER_CHECK
  // Checked against the sign this call RETURNS, so it runs after the
  // oracle.  The oracle still decides: this only observes.
  if (tr && tr->ctr && s != Sign::Zero) detail::check_sign_filter(v, s, *tr->ctr);
#endif
  if (s == Sign::Zero) {   // impossible for a nonzero element (oracle banner)
    if (tr) tr->refusal = Refusal::Undecided;
    return std::nullopt;
  }
  return s;
}

// The ring's order, three-way: sign(a - b).  Z[gamma] is a subring of R,
// so this IS the real order; nullopt = refused as sign_real refuses.
template <class Coef>
inline SignOr compare(const Real30T<Coef>& a, const Real30T<Coef>& b, SignTrace* tr = nullptr) {
  return sign_real(a - b, tr);
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

// The 128-bit product of two 64-bit words as (hi, lo), by 32-bit halves:
// the only wide multiply the division needs, and no wide type.
inline void mul64(uint64_t x, uint64_t y, uint64_t& hi, uint64_t& lo) {
  const uint64_t x0 = (uint32_t)x, x1 = x >> 32, y0 = (uint32_t)y, y1 = y >> 32;
  const uint64_t p00 = x0 * y0, p01 = x0 * y1, p10 = x1 * y0, p11 = x1 * y1;
  const uint64_t mid = (p00 >> 32) + (uint32_t)p01 + (uint32_t)p10;   // < 3 * 2^32
  lo = (mid << 32) | (uint32_t)p00;
  hi = p11 + (p01 >> 32) + (p10 >> 32) + (mid >> 32);
}

// Reduction modulo the two primes by shifts and adds, of a value given as
// hi 2^64 + lo -- no wide division, no wide type (a GPU emulates a
// 128-bit division in hundreds of instructions).
//   2^61 - 1:  a 64-bit x == (x & M) + (x >> 61) folded once, and
//              2^64 == 8, so hi 2^64 + lo == 8 hi + lo;
//   2^63 - 25: a 64-bit x == (x & L) + 25 (x >> 63) folded once, and
//              2^64 == 50 = 32 + 16 + 2, each 2^k hi folded through
//              2^63 == 25.
inline uint64_t fold61(uint64_t x) {
  const uint64_t M = (1ULL << 61) - 1;
  const uint64_t r = (x & M) + (x >> 61);   // < 2^61 + 8
  return r >= M ? r - M : r;
}
inline uint64_t fold63(uint64_t x) {
  const uint64_t P = (1ULL << 63) - 25, L = (1ULL << 63) - 1;
  const uint64_t r = (x & L) + (x >> 63) * 25;   // < 2^63 + 25
  return r >= P ? r - P : r;
}
// (a + b) mod p for residues a, b < p < 2^63.
inline uint64_t addmod(uint64_t a, uint64_t b, uint64_t p) {
  const uint64_t s = a + b;
  return s >= p ? s - p : s;
}
// h 2^k mod (2^63 - 25) for a residue h and k <= 6: split at 2^63.
inline uint64_t pow2mod63(uint64_t h, int k) {
  const uint64_t L = (1ULL << 63) - 1;
  const uint64_t high = h >> (63 - k), low = (h << k) & L;
  return fold63(low + 25 * high);   // low < 2^63, 25 high < 2^11
}
inline uint64_t reduce_p61(uint64_t hi, uint64_t lo) {
  const uint64_t M = (1ULL << 61) - 1;
  const uint64_t a = fold61(hi);            // hi == a, a < M
  const uint64_t b = fold61(a << 3);        // 8 a < 2^64
  return addmod(b, fold61(lo), M);
}
inline uint64_t reduce_p63(uint64_t hi, uint64_t lo) {
  const uint64_t P = (1ULL << 63) - 25;
  const uint64_t h = fold63(hi);
  const uint64_t fifty_h = addmod(addmod(pow2mod63(h, 5), pow2mod63(h, 4), P),
                                  pow2mod63(h, 1), P);
  return addmod(fifty_h, fold63(lo), P);
}
inline uint64_t reduce_mod(uint64_t hi, uint64_t lo, uint64_t p) {
  return p == kDivPrimes[0] ? reduce_p61(hi, lo) : reduce_p63(hi, lo);
}
inline uint64_t mulmod(uint64_t x, uint64_t y, uint64_t p) {
  uint64_t hi, lo;
  mul64(x, y, hi, lo);
  return reduce_mod(hi, lo, p);
}
inline uint64_t powmod(uint64_t x, uint64_t e, uint64_t p) {
  uint64_t r = 1;
  while (e) {
    if (e & 1) r = mulmod(r, x, p);
    x = mulmod(x, x, p);
    e >>= 1;
  }
  return r;
}
// The residue of a signed coefficient of any width: its magnitude as
// (hi, lo) words, reduced, then negated in the field.
template <class Coef>
inline uint64_t tomod(Coef v, uint64_t p) {
  using U = std::make_unsigned_t<Coef>;
  const U mag = v < 0 ? (U)0 - (U)v : (U)v;
  uint64_t hi = 0, lo = (uint64_t)mag;
  if constexpr (sizeof(U) > 8) hi = (uint64_t)(mag >> 64);
  const uint64_t r = reduce_mod(hi, lo, p);
  return v < 0 && r ? p - r : r;
}

// Inverse of d in F_p[y]/(psi) by extended Euclid on (psi, d); nullopt if
// gcd(d, psi) != 1 mod p (then p | N(d): the caller tries the other
// prime).  The degree invariant deg(t1) = 4 - deg(r0) makes the t-update
// bound exact: deg(t1) + shift = 4 - d1 <= 3, so no term is ever dropped
// (proven in the 2026-08-25 review).  @variant d1 (strictly decreases per
// swap).
template <class Coef>
inline std::optional<std::array<uint64_t, 4>> ring_inverse_mod(
    const std::array<Coef, 4>& d, uint64_t p) {
  // r0 = psi (degree 4), r1 = d mod p; t0 = 0, t1 = 1; invariant
  // t_i * d == r_i  (mod psi, p).
  uint64_t r0[5], r1[5] = {}, t0[5] = {}, t1[5] = {};
  for (int i = 0; i <= 4; i++) r0[i] = tomod<long long>(kPsi30[i], p);
  for (int i = 0; i < 4; i++) r1[i] = tomod<Coef>(d[i], p);
  t1[0] = 1;
  auto degree = [](const uint64_t* r) {
    for (int k = 4; k >= 0; k--)
      if (r[k]) return k;
    return -1;
  };
  int d0 = 4, d1 = degree(r1);
  if (d1 < 0) return std::nullopt;   // d == 0 mod p
  while (d1 > 0) {
    const uint64_t inv_lead = powmod(r1[d1], p - 2, p);
    while (d0 >= d1) {
      const uint64_t f = mulmod(r0[d0], inv_lead, p);
      const int shift = d0 - d1;
      for (int i = 0; i <= d1; i++)
        r0[i + shift] = addmod(r0[i + shift], p - mulmod(f, r1[i], p), p);
      for (int i = 0; i <= 4 - shift; i++)
        t0[i + shift] = addmod(t0[i + shift], p - mulmod(f, t1[i], p), p);
      const int nd0 = degree(r0);
      if (nd0 == d0) return std::nullopt;   // dead guard: cancellation exact
      d0 = nd0;
      if (d0 < 0) break;
    }
    for (int i = 0; i <= 4; i++) {
      const uint64_t tr = r0[i], tt = t0[i];
      r0[i] = r1[i]; t0[i] = t1[i];
      r1[i] = tr;    t1[i] = tt;
    }
    const int td = d0;
    d0 = d1;
    d1 = td;
    if (d1 < 0) return std::nullopt;   // gcd has positive degree
  }
  // r1 is a nonzero constant: inverse = t1 / r1[0].
  const uint64_t ic = powmod(r1[0], p - 2, p);
  std::array<uint64_t, 4> inv{};
  for (int i = 0; i < 4; i++) inv[i] = mulmod(t1[i], ic, p);
  return inv;
}

// (x * y) mod (psi, p): the SAME conv + kGamPow fold as the checked
// product, over F_p scalars (a guess path -- its correctness is not
// load-bearing, the verification product is).
inline std::array<uint64_t, 4> mulmod_ring(const std::array<uint64_t, 4>& x,
                                           const std::array<uint64_t, 4>& y,
                                           uint64_t p) {
  uint64_t conv[7] = {};
  for (int i = 0; i < 4; i++)
    for (int j = 0; j < 4; j++)
      conv[i + j] = addmod(conv[i + j], mulmod(x[i], y[j], p), p);
  std::array<uint64_t, 4> r{};
  for (int i = 0; i < 4; i++) r[i] = conv[i];
  for (int k = 4; k < 7; k++) {
    if (!conv[k]) continue;
    const auto& row = kGamPow[k];
    for (int i = 0; i < 4; i++)
      r[i] = addmod(r[i], mulmod(conv[k], tomod<long long>(row[i], p), p), p);
  }
  return r;
}

}  // namespace detail

// numerator / divisor in Z[gamma], exact or refused by name.  The quotient
// is lifted symmetrically from its residue, so a correct quotient must
// have coordinates below p/2 (2^60 for the first prime, 2^62 for the
// second) -- true of every module quantity inside the envelope; a larger
// true quotient is refused (NotDivisible), never mis-lifted.
// @anchor cyclotomic-exact-div
template <class Coef>
inline std::optional<Real30T<Coef>> exact_div(const Real30T<Coef>& numerator,
                                              const Real30T<Coef>& divisor,
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
    const uint64_t p = detail::kDivPrimes[pi];
    const auto inv = detail::ring_inverse_mod<Coef>(divisor.a, p);
    if (!inv) continue;
    std::array<uint64_t, 4> nm{};
    for (int i = 0; i < 4; i++) nm[i] = detail::tomod<Coef>(numerator.a[i], p);
    const auto qm = detail::mulmod_ring(nm, *inv, p);
    Real30T<Coef> q;
    for (int i = 0; i < 4; i++)
      q.a[i] = qm[i] > p / 2 ? -(Coef)(p - qm[i]) : (Coef)qm[i];
    if (divisor * q == numerator) {   // the exact verification
      if (tr) tr->prime_index = pi;
      return q;
    }
  }
  if (tr) tr->refusal = Refusal::NotDivisible;
  return std::nullopt;
}

// ---------------------------------------------------------------------------
// Zeta30T: Z[zeta_30] as the rank-2 module over Z[gamma] with basis
// {1, zeta}, zeta = zeta_30 = e^{i pi/15}, zeta^2 = gamma zeta - 1.  The
// POINT ring: kis module points (x kPointScale) are Zeta30 values, and
// the carried quantities fall out exactly --
//   lsq(u)      = |u|^2       = u.x^2 + gamma u.x u.y + u.y^2   (Real30)
//   wedge(u, v) = the delta-normalized cross product = x y' - y x'
// (with 16 Area^2 = (2 - gamma_2) wedge^2, whitepaper @ref eq:heron).
// RUN-PATH since layer 2: the CyclotomicMetric's star development
// (delaunay_cyclotomic.hh) develops flat stars in this module -- the
// tenth policy word the flat-vertex removal needs; construction and
// cross-check uses stay too (the walk identity as five_pentagon_spoke).
// KINSHIP: this is the same quadratic-tower body as Eisenstein
// (eisenstein.hh) at the trace gamma = 1 = 2 cos(pi/3) -- product,
// conjugation (zeta -> gamma - zeta), norm x^2 + gamma x y + y^2, and the
// wedge x y' - y x' all specialize term-for-term; the two stay separate
// because Eisenstein is the unchecked, device-legal, int-based hot path
// and this is the checked tier.
// ---------------------------------------------------------------------------
template <class R>
struct Zeta30T {
  using ring = R;
  R x, y;   // x + y * zeta

  bool ok() const { return x.ok && y.ok; }

  static Zeta30T integer(long long n) { return {R::integer(n), {}}; }
  static Zeta30T zeta_pow(int k) {
    k = ((k % 30) + 30) % 30;
    Zeta30T r = integer(1);
    const Zeta30T z{{}, R::integer(1)};
    for (int i = 0; i < k; i++) r = r * z;
    return r;
  }
  // The walk identity (whitepaper @ref eq:walk): the kPointScale-scaled
  // pentagon spoke 5 R5 zeta_60^9 = 2 + z^-3 + z^3 + z^6 + 2 z^9 + z^12
  // (z = zeta_30) -- sigma-even, hence HERE, one conductor down from its
  // two sigma-odd factors.
  static Zeta30T five_pentagon_spoke() {
    return integer(2) + zeta_pow(-3) + zeta_pow(3) + zeta_pow(6) +
           2 * zeta_pow(9) + zeta_pow(12);
  }

  bool is_zero() const { return x.is_zero() && y.is_zero(); }
  friend bool operator==(const Zeta30T& u, const Zeta30T& v) {
    return (u.x == v.x) & (u.y == v.y);   // no early exit (see the ring's operator==)
  }
  friend Zeta30T operator+(const Zeta30T& u, const Zeta30T& v) {
    return {u.x + v.x, u.y + v.y};
  }
  friend Zeta30T operator-(const Zeta30T& u, const Zeta30T& v) {
    return {u.x - v.x, u.y - v.y};
  }
  Zeta30T operator-() const { return {-x, -y}; }
  friend Zeta30T operator*(long long n, const Zeta30T& u) {
    return {n * u.x, n * u.y};
  }
  // (x + y zeta)(x' + y' zeta) with zeta^2 = gamma zeta - 1.  The gamma
  // term is associated (gamma * y) * y' so that both operands of every
  // product stay input-sized (the envelope's discipline).
  friend Zeta30T operator*(const Zeta30T& u, const Zeta30T& v) {
    return {u.x * v.x - u.y * v.y,
            u.x * v.y + u.y * v.x + (R::gamma() * u.y) * v.y};
  }
  // Complex conjugation: zeta -> gamma - zeta.
  Zeta30T conj() const { return {x + R::gamma() * y, -y}; }

  // |u|^2, a ring element (the identity (conj(u) u).y == 0 holds by algebra).
  R lsq() const { return x * x + (R::gamma() * x) * y + y * y; }
};
using Zeta30 = Zeta30T<Real30>;

// The delta-normalized wedge of two vectors: Im(conj(u) v)/sin(pi/15).
template <class R>
inline R wedge(const Zeta30T<R>& u, const Zeta30T<R>& v) {
  return u.x * v.y - u.y * v.x;
}

// ---------------------------------------------------------------------------
// DiamondT: the exact wedge-carrying diamond classifier -- the
// CyclotomicMetric's predicate core (the DiamondForms skeleton is shared
// with DiamondSq and Diamond60, diamond_forms.hh; delaunay_geometry.hh
// carries the picture).  Five squared lengths PLUS the two face wedges,
// all in the module's x kLsqScale convention; sqrt(H) = sqrt(2 - gamma_2) w
// makes every verdict the sign of ONE ring element, and the flip
// transports the carry (whitepaper @ref sec:pred, sec:flip).
//
// THE INVARIANT IS THE CONSTRUCTOR'S: DiamondT::make refuses (by name)
// unless both faces are positively oriented (wu, wl strictly positive)
// and the carry is CONSISTENT -- H_upper == (2 - gamma_2) wu^2 and
// H_lower == (2 - gamma_2) wl^2, exactly.  The check is TIGHT: the factor
// is a unit in an integral domain, so the two identities pin each wedge
// up to sign and the positivity pins the sign -- a wrong carry cannot
// exist in a constructed diamond, and the predicates never re-check.
// @inv H_upper() == heron_unit * wu^2  &&  H_lower() == heron_unit * wl^2
//      && wu, wl > 0   (for the life of the object; fields immutable)
//
// The bool convenience forms fold a refusal CONSERVATIVELY -- is_convex /
// is_cocircular to false ("do not act"), and is_delaunay ALSO to false,
// which for a bool-only caller reads as "must flip": exactly DiamondSq's
// convention, kept for consistency, but a mutating caller must consume
// delaunay_form_sign() and trip loudly on nullopt, never the bool.
// (Refusals cannot arise from a constructed diamond except by later
// poison, which the signs surface.)  @anchor cyclotomic-diamond
// ---------------------------------------------------------------------------
template <class R>
struct DiamondT {
 public:
  using ring = R;

  // The validated constructor: nullopt (with the named reason) on a
  // poisoned input, a non-positive wedge, or an inconsistent carry.
  static std::optional<DiamondT> make(const R& e, const R& a, const R& b,
                                      const R& c, const R& d, const R& wu,
                                      const R& wl, Refusal* why = nullptr) {
    DiamondT D{{e, a, b, c, d}, wu, wl};
    if (!(e.ok && a.ok && b.ok && c.ok && d.ok && wu.ok && wl.ok)) {
      if (why) *why = Refusal::Poisoned;
      return std::nullopt;
    }
    const SignOr su = sign_real(wu), sl = sign_real(wl);
    if (!su || !sl || *su != Sign::Positive || *sl != Sign::Positive) {
      if (why) *why = Refusal::NonPositiveWedge;
      return std::nullopt;
    }
    // (unit * w) * w, not unit * (w * w): both operands of each product
    // stay input-sized (the envelope's binding term, see carry_coeff_max).
    // A product that poisons here is the input leaving the envelope --
    // reported as such, not as an inconsistency.
    const R u = R::heron_unit();
    const R Hu = D.H_upper(), Hl = D.H_lower();
    const R Uu = (u * wu) * wu, Ul = (u * wl) * wl;
    if (!(Hu.ok && Hl.ok && Uu.ok && Ul.ok)) {
      if (why) *why = Refusal::Poisoned;
      return std::nullopt;
    }
    if (!(Hu == Uu && Hl == Ul)) {
      if (why) *why = Refusal::InconsistentCarry;
      return std::nullopt;
    }
    return D;
  }

  const R& e() const { return f_.e; }
  const R& a() const { return f_.a; }
  const R& b() const { return f_.b; }
  const R& c() const { return f_.c; }
  const R& d() const { return f_.d; }
  const R& wu() const { return wu_; }
  const R& wl() const { return wl_; }

  R s_upper() const { return f_.s_upper(); }
  R s_lower() const { return f_.s_lower(); }
  R P() const { return f_.P(); }
  R Q() const { return f_.Q(); }
  R H_upper() const { return metric_forms::heron_product_sq(f_.e, f_.a, f_.b); }
  R H_lower() const { return metric_forms::heron_product_sq(f_.e, f_.c, f_.d); }

  // THE ring elements whose signs are the verdicts (one element each --
  // the module's thesis; certificate-able as such).
  R delaunay_form() const { return s_upper() * wl_ + s_lower() * wu_; }
  R convexity_form_at_origin() const { return Q() * wu_ + P() * wl_; }

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
  DiamondT reversed() const { return DiamondT{f_.reversed(), wu_, wl_}; }
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
    R f2, w_origin, w_far;
  };
  std::optional<Flipped> flipped(DivTrace* tr = nullptr) const {
    const R two_e = 2 * f_.e;
    const auto q1 = exact_div(P() * Q() - (R::heron_unit() * wu_) * wl_, two_e, tr);
    if (!q1) return std::nullopt;
    const auto wo = exact_div(Q() * wu_ + P() * wl_, two_e, tr);
    if (!wo) return std::nullopt;
    const R f2 = f_.a + f_.c - *q1;
    const R wf = (wu_ + wl_) - *wo;
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
  std::optional<DiamondT> flip(Refusal* why = nullptr) const {
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
  DiamondT(const DiamondForms<R>& f, const R& wu, const R& wl)
      : f_(f), wu_(wu), wl_(wl) {}

  DiamondForms<R> f_;
  R wu_, wl_;   // the carried face wedges (upper, lower)
};
using Diamond = DiamondT<Real30>;

}  // namespace cyclotomic
