#pragma once
// ============================================================================
// cyclotomic_ambient.hh -- the AMBIENT rings above Z[gamma]: construction
// and cross-check machinery for cyclotomic.hh, NOT run-path code.
//
//   sigma7   The generator of Gal(Q(gamma)/Q) (cyclic of order 4,
//            gamma -> 2 cos(7 pi/15)), as an integer matrix on the power
//            basis: lets tests state Galois facts (norm rationality, the
//            Heron unit's norm 1, the divisibility premise behind
//            exact_div) exactly.
//
//   Zeta60   Z[zeta_60] in the rank-16 power basis, with its own sign
//            oracle and the two-square-root scheme (Diamond60) -- the
//            module's first landing, kept in full as the independent
//            cross-check: the sigma-descent says the rank-4 verdicts must
//            agree with it wherever both apply, and the embeddings below
//            make that a testable ring statement (the suite compares both
//            pinned and randomized diamonds).  Superseded on the run path
//            by cyclotomic.hh's wedge-carrying Diamond (whitepaper
//            @ref sec:pred, sec:flip); do not consume it from new
//            algorithm code.
//
// Everything here inherits the checked-and-poisoned overflow discipline
// through CheckedCoeffRing (cyclotomic.hh states it once); host tier
// throughout (real_part_with_abs calls std::cos).
// ============================================================================

#include "bigint_fixed.hh"
#include "cyclotomic.hh"

namespace cyclotomic {

inline constexpr int kOrder60 = 60;

namespace detail {

// Phi_60(x) = x^16 + x^14 - x^10 - x^8 - x^6 + x^2 + 1 -- THE authority;
// the reduction row is its negated tail (the same derivation idiom as
// kGamReduceRow from kPsi30).
inline constexpr std::array<long long, 17> kPhi60 = {
    1, 0, 1, 0, 0, 0, -1, 0, -1, 0, -1, 0, 0, 0, 1, 0, 1};

inline constexpr auto make_reduce_row60() {
  std::array<long long, 16> r{};
  for (int i = 0; i < 16; i++) r[i] = -kPhi60[i];
  return r;
}
inline constexpr auto kReduceRow60 = make_reduce_row60();

// kPow60[k] = the reduced coordinates of zeta^k, k = 0..59.
inline constexpr auto make_pow_table60() {
  std::array<std::array<long long, 16>, kOrder60> P{};
  P[0][0] = 1;
  for (int k = 1; k < kOrder60; k++) {
    const long long top = P[k - 1][15];
    for (int i = 15; i > 0; i--) P[k][i] = P[k - 1][i - 1];
    P[k][0] = 0;
    for (int i = 0; i < 16; i++) P[k][i] += top * kReduceRow60[i];
  }
  return P;
}
inline constexpr auto kPow60 = make_pow_table60();

}  // namespace detail

// ---------------------------------------------------------------------------
// The ring Z[zeta_60] (its rank-16 product folds through kPow60 with
// per-coordinate amplification <= 58 B^2, the guard's constant; the
// accumulation is 128-bit -- this is host-tier scaffolding).
// ---------------------------------------------------------------------------
struct Zeta60 : CheckedCoeffRing<Zeta60, 16, __int128> {
  static constexpr int kAmplification = 58;
  static const std::array<long long, 16>& pow_row(int k) {
    return detail::kPow60[k];
  }

  // zeta^k, any k (mod 60).
  static Zeta60 zeta(int k) {
    return from_coords(detail::kPow60[((k % kOrder60) + kOrder60) % kOrder60]);
  }
  // 2 cos(k*pi/30) = zeta^k + zeta^-k -- in THIS ring's own zeta, so
  // Zeta60::two_cos(2k) == the embedded Real30::two_cos(k).
  static Zeta60 two_cos(int k) { return zeta(k) + zeta(-k); }

  // Complex conjugate: zeta^k -> zeta^{-k}.  Checked like the product
  // (provably safe even unchecked: |kPow60 entries| <= 1, so the
  // accumulation is <= 16 * 2^63 = 2^67 << 2^127).
  Zeta60 conj() const {
    Zeta60 r;
    r.ok = ok;
    __int128 acc[16] = {};
    for (int k = 0; k < 16; k++) {
      if (!a[k]) continue;
      const auto& row = detail::kPow60[(kOrder60 - k) % kOrder60];
      for (int i = 0; i < 16; i++) acc[i] += (__int128)a[k] * row[i];
    }
    for (int i = 0; i < 16; i++) {
      r.a[i] = (long long)acc[i];
      r.ok &= ((__int128)r.a[i] == acc[i]);
    }
    return r;
  }
  bool is_real() const {
    const Zeta60 c = conj();
    return ok && c.ok && a == c.a;
  }

  // |z|^2 = z * conj(z)  (a real element).
  Zeta60 norm() const { return *this * conj(); }

  // Re(z) as a double, with the abs-sum the rung-1 bound scales by.  (The
  // name says what it computes for ANY element; sign_real guards realness
  // separately.)
  double real_part() const {
    double v, abs_sum;
    real_part_with_abs(v, abs_sum);
    return v;
  }
  void real_part_with_abs(double& value, double& abs_sum) const {
    value = 0;
    abs_sum = 0;
    for (int k = 0; k < 16; k++) {
      if (!a[k]) continue;
      value += (double)a[k] * std::cos(k * 3.14159265358979323846 / 30.0);
      abs_sum += std::fabs((double)a[k]);
    }
  }
};

// ---------------------------------------------------------------------------
// Embedding Real30 -> Z[zeta_60]: gamma = zeta_60^2 + zeta_60^-2.  (The
// Zeta30 embedding follows its definition below.)  These make every
// rank-4 / rank-8 statement a conductor-60 ring statement -- the
// cross-check bridge.
// ---------------------------------------------------------------------------
inline Zeta60 to_zeta60(const Real30& v) {
  const Zeta60 g = Zeta60::two_cos(2);
  Zeta60 r = Zeta60::integer(v.a[Real30::kRank - 1]);
  for (int k = Real30::kRank - 2; k >= 0; k--)
    r = r * g + Zeta60::integer(v.a[k]);
  r.ok &= v.ok;
  return r;
}

// The named constants of the flattened kis surface at conductor 60.  The
// first three are sigma-ODD -- genuinely outside the real subring, which
// is the descent's whole content -- and are spelled natively; everything
// sigma-even is DERIVED from the rank-4 authority through the embedding
// (one spelling, no test-enforced twins).
// s = 2 sin(pi/5) and t = 2 sin(2*pi/5); (s t)^2 == 5 exactly.
inline Zeta60 two_sin_pi5() { return Zeta60::two_cos(9); }
inline Zeta60 two_sin_2pi5() { return Zeta60::two_cos(3); }
// 5 * R5 (an algebraic integer): R5 = 1/(2 sin(pi/5)) = s t^2 / 5.
inline Zeta60 five_R5() {
  const Zeta60 t = two_sin_2pi5();
  return two_sin_pi5() * t * t;
}
// The golden ratio; golden^2 == golden + 1 exactly.
inline Zeta60 golden() { return to_zeta60(Real30::golden()); }
// The kLsqScale-scaled squared kis lengths, embedded from the authority.
inline Zeta60 lsq_unit() { return to_zeta60(Real30::lsq_cubic_edge()); }
inline Zeta60 lsq_pentagon_spoke() {
  return to_zeta60(Real30::lsq_pentagon_spoke());
}

// ---------------------------------------------------------------------------
// The conductor-60 sign oracle: the same three-rung ladder as Real30's,
// over the degree-8 real subfield Q(2 cos(pi/30)).  Rung (1) evaluates
// per-term cosines: no Horner amplification (every term is bounded by
// |a_k| outright), but a 16-term accumulation and a libm cos -- a
// DIFFERENT error structure from Real30's, with a SMALLER derived bound
// and a correspondingly tighter constant (see kRung1RelErr60).  Rung (2)
// rewrites 2v as an integer polynomial in c = 2 cos(pi/30) via the
// Chebyshev table and hands it to the shared bisection driver.
// Termination at full int64 heights: |B(c)| >= 2^-469 (degree <= 8;
// conjugate bound 16 * 2^63 = 2^67 across the 8 real embeddings, and
// 7 * 67 = 469), derivative bound M <= 769 max|b_k| <= 2^80, so
// s ~ 469 + 81 + 1 = 551 bits inside kMaxDyadicBits60 = 768 (the first
// landing's review measured an LLL-adversarial extreme of s = 509).
// ---------------------------------------------------------------------------

// Rung-1 bound for the conductor-60 evaluation (per-term std::cos):
// int->double conversion <= u|a_k|, cos argument+libm error <= 6.4e-16
// absolute (|pi_double - pi| plus <= 1 ulp cos), product rounding
// <= u|a_k|, accumulation over <= 16 terms gives 15u * abs_sum; total
// <= 2.6e-15 * abs_sum under round-to-nearest and a <= 1-ulp cos, ~3.2e-15
// at SYCL's 4-ulp cos allowance.  kRung1RelErr60 = 1e-14 keeps >= 3x
// margin over either.
inline constexpr double kRung1RelErr60 = 1e-14;
static_assert(kRung1RelErr60 >= 9.6e-15,
              "must keep >= 3x margin over the 4-ulp worst reading 3.2e-15");

namespace detail {

// psi_60, the minimal polynomial of c = 2 cos(pi/30) (monic, degree 8):
// c^8 - 7 c^6 + 14 c^4 - 8 c^2 + 1  (index = power of c), and its
// negated-tail reduction row (the shared derivation idiom).
inline constexpr std::array<long long, 9> kPsi60 = {1, 0, -8, 0,  14,
                                                    0, -7, 0, 1};
inline constexpr auto make_psi60_reduce_row() {
  std::array<long long, 8> r{};
  for (int i = 0; i < 8; i++) r[i] = -kPsi60[i];
  return r;
}
inline constexpr auto kPsi60ReduceRow = make_psi60_reduce_row();

// kTwoCosInC[k] = the coefficients of two_cos(k) = zeta^k + zeta^-k as an
// integer polynomial in c, reduced mod psi_60 (degree <= 7): p_0 = 2,
// p_1 = c, p_k = c p_{k-1} - p_{k-2} (the shift spill folded through
// kPsi60ReduceRow -- the same idiom as make_two_cos15).
inline constexpr auto make_two_cos_table60() {
  std::array<std::array<long long, 8>, 16> P{};
  P[0][0] = 2;
  P[1][1] = 1;
  for (int k = 2; k < 16; k++) {
    long long t[9] = {};
    for (int i = 0; i < 8; i++) t[i + 1] = P[k - 1][i];
    const long long top = t[8];
    for (int i = 0; i < 8; i++)
      P[k][i] = t[i] + top * kPsi60ReduceRow[i] - P[k - 2][i];
  }
  return P;
}
inline constexpr auto kTwoCosInC = make_two_cos_table60();

// ---- The bisection driver (host tier): the sign of B(alpha) for alpha
// the root of a monic Psi isolated in a dyadic interval, by exact
// fixed-capacity integer arithmetic.  Once the working ring's rung 2; the
// working ring now decides every sign by one fixed-point evaluation
// (cyclotomic.hh), and this driver serves the conductor-60 oracle below
// and its cross-checks only.  Capacities: the binding intermediate is
// lhs = Bmid << smid at the deepest iteration, ~ (deg+1)*smid + 63 bits;
// the rank-4 sizing (kMaxDyadicBits, Big) stays as the driver's smallest
// instantiation, exercised by the property suite.
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

// @pre B not identically zero; Psi(lo) and Psi(hi) of opposite signs with
// exactly one root of Psi inside.  @variant max_bits - s.  In its own
// frame (the FixedBigInt arrays are multi-KB).
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

// Rung 1, the double evaluation against its derived bound: conclusive iff
// |val| clears it.
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

// Rung-2 capacities (termination bound in the oracle banner; the binding
// intermediate is again lhs = Bmid << smid, provably 99 of the 100 limbs).
inline constexpr int kMaxDyadicBits60 = 768;
inline constexpr int kLimbs60 = (8 * kMaxDyadicBits60) / 64 + 4;
using Big60 = FixedBigInt<kLimbs60>;

// Rung (2) for the conductor-60 ring: rewrite 2v as an integer polynomial
// B in c through the Chebyshev table, then the shared driver on the
// isolating interval [63/32, 2] (exactly: 2^40 psi_60(63/32) =
// -1752146664959 < 0 < psi_60(2) = 1, and the next root below
// c = 1.98904.. is 2 cos(7 pi/30) = 1.48629..).
// @pre v real, unpoisoned, nonzero.
[[gnu::noinline]] inline SignOr sign_real_exact60(const Zeta60& v,
                                                  SignTrace* tr) {
  __int128 b[8] = {};
  b[0] = 2 * (__int128)v.a[0];
  for (int k = 1; k < 16; k++) {
    if (!v.a[k]) continue;
    for (int i = 0; i < 8; i++)
      b[i] += (__int128)v.a[k] * kTwoCosInC[k][i];
  }
  Big60 B[8], Psi[9];
  int deg = 0;
  for (int i = 0; i < 8; i++) {
    B[i] = Big60::from_i128(b[i]);
    if (B[i].sgn) deg = i;
  }
  for (int i = 0; i <= 8; i++) Psi[i] = Big60::from_i128(kPsi60[i]);
  return sign_at_isolated_root(B, deg, Psi, 8, /*m0=*/63, /*s0=*/5,
                               kMaxDyadicBits60, tr);
}

}  // namespace detail

// Exact sign of a real conductor-60 element; nullopt = refused by name
// (Poisoned / NotReal / the rung-2 caps, via the trace).
// @pre v.is_real() -- checked, a violation refuses rather than mis-signs.
inline SignOr sign_real(const Zeta60& v, SignTrace* tr = nullptr) {
  if (!v.ok) {
    if (tr) tr->refusal = Refusal::Poisoned;
    return std::nullopt;
  }
  if (v.is_zero()) {
    if (tr) tr->rung = 0;
    return Sign::Zero;
  }
  if (!v.is_real()) {
    if (tr) tr->refusal = Refusal::NotReal;
    return std::nullopt;
  }
  double val, abs_sum;
  v.real_part_with_abs(val, abs_sum);
  if (const auto s = detail::rung1(val, abs_sum, kRung1RelErr60, tr))
    return *s;
  if (tr) tr->rung = 2;
  return detail::sign_real_exact60(v, tr);
}

// ---------------------------------------------------------------------------
// sign(x * sqrt(A) + y * sqrt(B)) for real ring elements with A, B >= 0:
// the two-square-root scheme the first landing reduced every Delaunay-form
// case to (superseded by the wedge carry, kept as the cross-check).  Total
// on its stated domain: a zero radicand drops its term, mixed live terms
// compare their squares (G = x^2 A - y^2 B).
// ---------------------------------------------------------------------------
inline SignOr sign_x_sqrtA_plus_y_sqrtB(const Zeta60& x, const Zeta60& A,
                                        const Zeta60& y, const Zeta60& B) {
  const SignOr sA = sign_real(A), sB = sign_real(B);
  if (!sA || !sB || *sA == Sign::Negative || *sB == Sign::Negative)
    return std::nullopt;   // outside the domain
  SignOr sx = sign_real(x), sy = sign_real(y);
  if (!sx || !sy) return std::nullopt;
  const Sign ex = (*sA == Sign::Zero) ? Sign::Zero : *sx;   // effective terms
  const Sign ey = (*sB == Sign::Zero) ? Sign::Zero : *sy;
  const int ix = (int)ex, iy = (int)ey;
  if (ix >= 0 && iy >= 0)
    return (ix == 0 && iy == 0) ? Sign::Zero : Sign::Positive;
  if (ix <= 0 && iy <= 0)
    return (ix == 0 && iy == 0) ? Sign::Zero : Sign::Negative;
  const SignOr g = sign_real(x * x * A - y * y * B);
  if (!g) return std::nullopt;
  return ix > 0 ? *g : -*g;
}

// ---------------------------------------------------------------------------
// Diamond60: the first landing's diamond classifier -- five squared
// lengths over Z[zeta_60] on the DiamondForms skeleton, verdicts assembled
// by the two-square-root scheme.  Superseded on the run path by
// cyclotomic::Diamond (the wedge carry); kept as the independent
// cross-check the tests compare against.  Domain conventions exactly as
// cyclotomic::Diamond / DiamondSq (refusal folds conservatively in the
// bool forms).
// ---------------------------------------------------------------------------
struct Diamond60 : DiamondForms<Zeta60> {
  Diamond60() = default;
  Diamond60(const Zeta60& e_, const Zeta60& a_, const Zeta60& b_,
            const Zeta60& c_, const Zeta60& d_)
      : DiamondForms<Zeta60>{e_, a_, b_, c_, d_} {}
  Diamond60(const DiamondForms<Zeta60>& f) : DiamondForms<Zeta60>(f) {}

  Zeta60 H_upper() const { return metric_forms::heron_product_sq(e, a, b); }
  Zeta60 H_lower() const { return metric_forms::heron_product_sq(e, c, d); }

  SignOr delaunay_form_sign() const {
    const Zeta60 Hu = H_upper(), Hl = H_lower();
    if (!strictly_positive(Hu) || !strictly_positive(Hl)) return std::nullopt;
    return sign_x_sqrtA_plus_y_sqrtB(s_upper(), Hl, s_lower(), Hu);
  }
  bool is_delaunay() const {
    const SignOr s = delaunay_form_sign();
    return s && *s != Sign::Negative;
  }
  bool is_cocircular() const {
    const SignOr s = delaunay_form_sign();
    return s && *s == Sign::Zero;
  }

  SignOr convex_at_origin_sign() const {
    const Zeta60 Hu = H_upper(), Hl = H_lower();
    if (!strictly_positive(Hu) || !strictly_positive(Hl)) return std::nullopt;
    return sign_x_sqrtA_plus_y_sqrtB(Q(), Hu, P(), Hl);
  }
  Diamond60 reversed() const {
    return Diamond60{DiamondForms<Zeta60>::reversed()};
  }
  bool is_convex() const {
    const SignOr u = convex_at_origin_sign();
    const SignOr v = reversed().convex_at_origin_sign();
    return u && v && *u == Sign::Positive && *v == Sign::Positive;
  }

 private:
  static bool strictly_positive(const Zeta60& H) {
    const SignOr s = sign_real(H);
    return s && *s == Sign::Positive;
  }
};

// ---------------------------------------------------------------------------
// sigma7: the generator of Gal(Q(gamma)/Q) = <sigma>, sigma(gamma) =
// 2 cos(7 pi/15), cyclic of order 4 (7^2 = 49 == +-11, 7^3 == +-13,
// 7^4 == +-1 mod 30).  Integer matrix on the power basis; checked like
// every ring operation.  Verification vocabulary: norms x sigma(x)
// sigma^2(x) sigma^3(x) are rational, N(2 - gamma_2) = 1, and the
// norm identity behind exact division's divisibility premise is its
// Galois statement.
// ---------------------------------------------------------------------------
namespace detail {
// Row i = the coordinates of sigma(gamma^i)  (= (2 cos(7 pi/15))^i mod psi).
inline constexpr long long kSigma7[4][4] = {
    {1, 0, 0, 0}, {-2, 3, 1, -1}, {2, -1, 0, 0}, {-5, 12, 3, -4}};
}  // namespace detail

inline Real30 sigma7(const Real30& x) {
  Real30 r;
  r.ok = x.ok;
  __int128 acc[Real30::kRank] = {};
  for (int i = 0; i < Real30::kRank; i++)
    for (int j = 0; j < Real30::kRank; j++)
      acc[j] += (__int128)x.a[i] * detail::kSigma7[i][j];
  for (int i = 0; i < Real30::kRank; i++) {
    r.a[i] = (long long)acc[i];
    r.ok &= ((__int128)r.a[i] == acc[i]);
  }
  return r;
}

// (Zeta30, the Z[zeta_30] point ring, was PROMOTED to cyclotomic.hh when
// the CyclotomicMetric's star development made it run-path -- layer 2;
// its embedding below stays here with the conductor-60 machinery.)

// Embedding Z[zeta_30] -> Z[zeta_60]: zeta_30 = zeta_60^2.
inline Zeta60 to_zeta60(const Zeta30& u) {
  return to_zeta60(u.x) + to_zeta60(u.y) * Zeta60::zeta(2);
}

}  // namespace cyclotomic
