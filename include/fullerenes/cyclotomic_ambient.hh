#pragma once
// ============================================================================
// cyclotomic_ambient.hh -- the AMBIENT rings above Z[gamma]: construction
// and cross-check machinery for cyclotomic.hh, NOT run-path code.
//
//   Zeta30   Z[zeta_30] as the rank-2 module Z[gamma] + Z[gamma] zeta
//            (zeta^2 = gamma zeta - 1): the ring the kis module's POINTS
//            live in (x kPointScale).  Squared lengths and wedges of
//            coordinate constructions fall out as Real30 elements (lsq,
//            wedge), so every carried quantity of the working surface can
//            be re-derived from explicit points and compared exactly --
//            the verification role.  The walk identity (whitepaper
//            @ref eq:walk) is five_pentagon_spoke().  KINSHIP: this is
//            the same quadratic-tower body as Eisenstein (eisenstein.hh)
//            at the trace gamma = 1 = 2 cos(pi/3) -- product, conjugation
//            (zeta -> gamma - zeta), norm x^2 + gamma x y + y^2, and the
//            wedge x y' - y x' all specialize term-for-term; the two stay
//            separate because Eisenstein is the unchecked, device-legal,
//            int-based hot path and this is checked host-tier scaffolding.
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
// per-coordinate amplification <= 58 B^2, inside the shared operand
// guard).
// ---------------------------------------------------------------------------
struct Zeta60 : CheckedCoeffRing<Zeta60, 16> {
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

// ---------------------------------------------------------------------------
// Zeta30: Z[zeta_30] as the rank-2 module over Z[gamma] with basis
// {1, zeta}, zeta = zeta_30 = e^{i pi/15}, zeta^2 = gamma zeta - 1.  The
// construction ring: kis module points (x kPointScale) are Zeta30 values,
// and the working surface's carried quantities fall out exactly:
//   lsq(u)      = |u|^2       = u.x^2 + gamma u.x u.y + u.y^2   (Real30)
//   wedge(u, v) = the delta-normalized cross product = x y' - y x'
// (with 16 Area^2 = (2 - gamma_2) wedge^2, whitepaper @ref eq:heron).
// See the file banner for the Eisenstein kinship (gamma = 1 instance).
// ---------------------------------------------------------------------------
struct Zeta30 {
  Real30 x, y;   // x + y * zeta

  bool ok() const { return x.ok && y.ok; }

  static Zeta30 integer(long long n) { return {Real30::integer(n), {}}; }
  static Zeta30 zeta_pow(int k) {
    k = ((k % 30) + 30) % 30;
    Zeta30 r = integer(1);
    const Zeta30 z{{}, Real30::integer(1)};
    for (int i = 0; i < k; i++) r = r * z;
    return r;
  }
  // The walk identity (whitepaper @ref eq:walk): the kPointScale-scaled
  // pentagon spoke 5 R5 zeta_60^9 = 2 + z^-3 + z^3 + z^6 + 2 z^9 + z^12
  // (z = zeta_30) -- sigma-even, hence HERE, one conductor down from its
  // two sigma-odd factors.
  static Zeta30 five_pentagon_spoke() {
    return integer(2) + zeta_pow(-3) + zeta_pow(3) + zeta_pow(6) +
           2 * zeta_pow(9) + zeta_pow(12);
  }

  bool is_zero() const { return x.is_zero() && y.is_zero(); }
  friend bool operator==(const Zeta30& u, const Zeta30& v) {
    return u.x == v.x && u.y == v.y;
  }
  friend Zeta30 operator+(const Zeta30& u, const Zeta30& v) {
    return {u.x + v.x, u.y + v.y};
  }
  friend Zeta30 operator-(const Zeta30& u, const Zeta30& v) {
    return {u.x - v.x, u.y - v.y};
  }
  Zeta30 operator-() const { return {-x, -y}; }
  friend Zeta30 operator*(long long n, const Zeta30& u) {
    return {n * u.x, n * u.y};
  }
  // (x + y zeta)(x' + y' zeta) with zeta^2 = gamma zeta - 1.
  friend Zeta30 operator*(const Zeta30& u, const Zeta30& v) {
    return {u.x * v.x - u.y * v.y,
            u.x * v.y + u.y * v.x + Real30::gamma() * (u.y * v.y)};
  }
  // Complex conjugation: zeta -> gamma - zeta.
  Zeta30 conj() const { return {x + Real30::gamma() * y, -y}; }

  // |u|^2, a Real30 (the identity (conj(u) u).y == 0 holds by algebra).
  Real30 lsq() const { return x * x + Real30::gamma() * (x * y) + y * y; }
};

// The delta-normalized wedge of two vectors: Im(conj(u) v)/sin(pi/15).
inline Real30 wedge(const Zeta30& u, const Zeta30& v) {
  return u.x * v.y - u.y * v.x;
}

// Embedding Z[zeta_30] -> Z[zeta_60]: zeta_30 = zeta_60^2.
inline Zeta60 to_zeta60(const Zeta30& u) {
  return to_zeta60(u.x) + to_zeta60(u.y) * Zeta60::zeta(2);
}

}  // namespace cyclotomic
