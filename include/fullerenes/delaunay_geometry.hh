#pragma once

// The intrinsic-geometry vocabulary of the metric DCEL: the numeric
// primitives and tolerance bands (delaunay_detail), the Diamond -- the
// five-length local geometry around an edge -- in both coordinate systems
// (Diamond over real lengths, DiamondSq over integer squared lengths), and
// the fan development of a flat vertex's star (FanPolygon / FanTriangulation
// + the ear acceptance tests of both regimes).  Pure mathematics over
// lengths: no DCEL, no workspace, no I/O -- provable in isolation and
// reusable by any metric-surface code.
//
// The exact regime IS Eisenstein-lattice geometry (every face a lattice
// triangle, the Heron-Eisenstein identity H = 3*tau^2), so this header
// includes eisenstein.hh for the lattice primitives (lattice_tau,
// apex_factor, wedge, norm2_ll).

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <numbers>
#include <optional>
#include <span>
#include <type_traits>

#include "fullerenes/eisenstein.hh"

// ============================================================================
// Intrinsic geometry primitives.
// ============================================================================
namespace delaunay_detail {

inline constexpr double two_pi = 2 * std::numbers::pi_v<double>;

// Tolerance bands of the floating-point predicates (the general-metric
// regime).  These absorb FP noise in the cotangent/Heron/development
// evaluations, nothing else: the mathematics is is_delaunay <=> cot-sum
// >= 0, strict convexity <=> both endpoint forms > 0, and the ear
// acceptance criteria of paper sec:ear-clip.  The exact integer regime
// (DiamondSq) has no bands at all.  Never widen a band to make a case pass.
inline constexpr double delaunay_band   = -1e-10;  // accept near-tight as tight
inline constexpr double convexity_band  =  1e-12;  // strictness margin
inline constexpr double cot_degenerate  =  1e15;   // sentinel for a degenerate triangle
inline constexpr double ear_area_band   =  1e-10;  // CCW ear area must exceed
inline constexpr double ear_angle_band  =  1e-10;  // reflex slack past pi
inline constexpr double ear_length_floor = 1e-15;  // degenerate-diagonal floor
inline constexpr double tie_cocircular_tol = 1e-12; // tie-break spoke test
inline constexpr double default_flat_tol   = 1e-6;  // equilateral-metric flatness band
inline constexpr double lsq_integrality_band = 1e-9; // |len^2 - round| relative band
                                                     // for entering the exact regime

// Heron product: H(a,b,c) = (a+b+c)(-a+b+c)(a-b+c)(a+b-c) = 16*Area^2.
// Returns 0 if the triangle inequality is violated.
// (The squared-length integer form heron_product_sq lives in eisenstein.hh,
// beside the lattice primitives that consume it.)
inline double heron_product(double a, double b, double c) {
  double s1 = -a + b + c;
  double s2 =  a - b + c;
  double s3 =  a + b - c;
  if (s1 < 0 || s2 < 0 || s3 < 0) return 0;
  return (a + b + c) * s1 * s2 * s3;
}

// Cotangent of the angle opposite side `opp` in a triangle with sides
// (opp, b, c): cot(alpha) = (b^2 + c^2 - opp^2) / sqrt(H).
// Returns +/-cot_degenerate on a degenerate triangle (H <= 0) -- a value no
// sane tolerance ever brackets, so degenerate diamonds classify as
// non-tight rather than tripping the tight test.
inline double cot_opposite(double opp, double b, double c) {
  double H = heron_product(opp, b, c);
  double num = b*b + c*c - opp*opp;
  if (H <= 0) return (num >= 0) ? cot_degenerate : -cot_degenerate;
  return num / std::sqrt(H);
}

// Angle of the corner adjacent to sides `adj1`, `adj2` in a triangle whose
// opposite side is `opp`.  Law of cosines, clamped for floating-point safety
// at triangle-inequality boundaries.
inline double triangle_angle(double adj1, double adj2, double opp) {
  double c = (adj1*adj1 + adj2*adj2 - opp*opp) / (2 * adj1 * adj2);
  return std::acos(std::clamp(c, -1.0, 1.0));
}

}  // namespace delaunay_detail

// ============================================================================
// The exact regime's width envelope.
//
// Squared lengths up to exact_lsq_max keep every integer intermediate in
// range.  An INDUCTIVE invariant: the initial carry is verified at entry
// (remove_flat_vertices_exact), and the two sites that write new squared
// lengths (the flip diagonal, the ear diagonal) trip loudly past the
// bound, so every predicate and development afterwards runs on in-envelope
// inputs.  The THREE binding constraints (M = exact_lsq_max, from
// tau <= M, |K| <= 2M, factor components |m| <= 3M/2, and lattice
// positions of length <= (2/sqrt 3)*sqrt(M)):
//   - apex_factor's factor components (<= 3M/2) must fit Eisenstein's int32;
//   - the apex placement's factor*base product components (<= 4*M^1.5,
//     Eisenstein::operator* is int32) -- THE binding constraint, measured
//     to first diverge near Lsq ~ 1.7e6 when violated;
//   - the flipped diagonal's factor-difference norm (<= 19*M^2), int64.
// The corpus-measured maxima are 27 (final iDT lengths) and 64
// (intermediates, instrumented probe) -- four orders below this envelope;
// a trip is loud, never silent.
// ============================================================================
inline constexpr long long exact_lsq_max = 500'000;
static_assert(3 * exact_lsq_max / 2 <= INT32_MAX,
              "apex_factor's m + n*w components must fit Eisenstein's int32");
static_assert(16 * exact_lsq_max * exact_lsq_max * exact_lsq_max <=
              (long long)INT32_MAX * INT32_MAX,
              "the apex placement's int32 product (4*M^1.5) must fit: 16*M^3 <= INT32_MAX^2");
static_assert(exact_lsq_max <=
              std::numeric_limits<long long>::max() / 19 / exact_lsq_max,
              "the flipped diagonal's factor-difference norm (19*M^2) must fit int64");

// An edge length in the two coordinates the pipeline carries: the real
// length and -- in the exact regime -- its integer square (a lattice norm).
// The banded regime carries len alone (sq stays 0 and is never read).
struct Length {
  double    len = 0;
  long long sq  = 0;
};

// ============================================================================
// DiamondSq: the diamond in SQUARED-length coordinates -- the exact regime.
//
// Same picture as Diamond below, but each field is the SQUARE of the edge
// length, an integer whenever the metric is equilateral-and-flips: every
// face is then an Eisenstein-lattice triangle, so by the Heron-Eisenstein
// identity (paper prop:heron-eisenstein) H = 16*Area^2 = 3*tau^2 with
// tau the integer lattice area number (lattice_tau).  The Delaunay form
//
//     F = s_upper*sqrt(H_lower) + s_lower*sqrt(H_upper)
//       = sqrt(3) * (s_upper*tau_lower + s_lower*tau_upper)
//
// is therefore sqrt(3) times an INTEGER, and the cotangent sum satisfies
// cot(angle_B) + cot(angle_D) = F / sqrt(H_upper * H_lower): sign(F)
// classifies the edge -- > 0 Delaunay, == 0 cocircular (tight), < 0 must
// flip -- as the sign of one integer linear form, tolerance-free.
// Convexity at an endpoint is the same shape (Q*tau_upper + P*tau_lower),
// with the second endpoint the SAME predicate on the reversed diamond (the
// reversal involution sigma(e,a,b,c,d) = (e,b,a,d,c); whole-diamond
// quantities are sigma-invariant, per-endpoint ones transform by sigma).
//
// Degeneracy: a diamond with a degenerate or non-lattice side (tau <= 0)
// is OUTSIDE the predicates' domain and classifies as none of
// Delaunay / cocircular / convex -- the conservative refusal (the float
// path's cot_degenerate sentinel agrees on the negative-numerator side and
// may still accept the positive-numerator side; the refusal funnels a
// corrupt Lsq carry into the sweep's loud trip instead).  A valid exact
// complex has strictly positive lattice face areas, so live diamonds are
// never refused (0 refusals measured over 3.4M live evaluations).
// @pre fields within [0, exact_lsq_max] for totality of the integer
// arithmetic (the inductive envelope; hostile larger values are out of
// contract -- heron_product_sq overflows int64 near 1.2e9).
// ============================================================================
struct DiamondSq {
  long long e, a, b, c, d;   // squared lengths

  // sigma, the reversal involution: the diamond read from the diagonal's
  // other endpoint.  sigma(sigma(D)) = D.
  DiamondSq reversed() const { return {e, b, a, d, c}; }

  // The law-of-cosines numerators (squared coordinates): s_* are the
  // cotangent numerators at the two apexes, P/Q the endpoint-convexity
  // numerators at the origin endpoint u.
  long long s_upper() const { return a + b - e; }
  long long s_lower() const { return c + d - e; }
  long long P() const { return e + a - b; }
  long long Q() const { return e + c - d; }

  // The lattice area numbers of the two triangles (H = 3*tau^2); -1 when a
  // side bounds no lattice triangle, 0 when degenerate.
  long long tau_upper() const { return lattice_tau(e, a, b); }
  long long tau_lower() const { return lattice_tau(e, c, d); }

  // F / sqrt(3): THE integer Delaunay form.  nullopt when either triangle
  // is degenerate or non-lattice (tau <= 0) -- outside the domain.
  std::optional<long long> delaunay_form() const {
    const long long tu = tau_upper(), tl = tau_lower();
    if (tu <= 0 || tl <= 0) return std::nullopt;
    return s_upper() * tl + s_lower() * tu;
  }
  bool is_delaunay()   const { auto F = delaunay_form(); return F && *F >= 0; }
  bool is_cocircular() const { auto F = delaunay_form(); return F && *F == 0; }

  // Strict convexity of the diamond at the diagonal's ORIGIN endpoint u:
  // angle(u) < pi <=> Q*sqrt(H_upper) + P*sqrt(H_lower) > 0
  //               <=> Q*tau_upper + P*tau_lower > 0.
  // The v-endpoint is convex_at_origin() of the reversed diamond.
  bool convex_at_origin() const {
    const long long tu = tau_upper(), tl = tau_lower();
    return tu > 0 && tl > 0 && Q() * tu + P() * tl > 0;
  }
  bool is_convex() const {
    return convex_at_origin() && reversed().convex_at_origin();
  }

  // Squared length of the flipped diagonal BD, exact: develop the diamond
  // on the lattice over its shared base u->v.  The two apex displacement
  // factors (apex_factor; B at chirality +1 = left of u->v, D at -1) differ
  // by a lattice vector with norm2(B - D) = norm2(fB - fD) / e, and the
  // division is exact whenever the five squared lengths bound a lattice
  // diamond -- nullopt otherwise (degenerate base and out-of-envelope
  // inputs included; never on a valid exact complex, so the caller trips).
  std::optional<long long> flipped_length_sq() const {
    if (e <= 0) return std::nullopt;                 // degenerate base: no diamond
    if (e > exact_lsq_max || a > exact_lsq_max || b > exact_lsq_max ||
        c > exact_lsq_max || d > exact_lsq_max) return std::nullopt;
    auto fB = apex_factor((int)e, (int)a, (int)b, +1);
    auto fD = apex_factor((int)e, (int)c, (int)d, -1);
    if (!fB || !fD) return std::nullopt;
    const long long n2 = norm2_ll(*fB - *fD);
    if (n2 % e) return std::nullopt;
    return n2 / e;
  }
};

// ============================================================================
// Diamond: the local geometry around an edge in a metrized triangulation.
//
//      B
//     / \          upper triangle: sides (e, a, b)
//    a   b         a = side adjacent to u, b = side adjacent to v
//   /     \  .
//  u---e---v       e = diagonal being tested/flipped
//   \     /
//    c   d         lower triangle: sides (e, c, d)
//     \ /          c = side adjacent to u, d = side adjacent to v
//      D
//
// All geometric predicates (Delaunay, convexity, flipped length) depend only
// on these five edge lengths.  No vertex IDs or topology needed.
// ============================================================================
struct Diamond {
  double e, a, b, c, d;

  // cot(angle_B) + cot(angle_D) >= 0, within delaunay_band.
  bool is_delaunay() const {
    using namespace delaunay_detail;
    return cot_opposite(e, a, b) + cot_opposite(e, c, d) >= delaunay_band;
  }

  // Angle sum < pi at both u and v, strictly (convexity_band margin).
  // The two endpoints are the SAME half-predicate applied to the two ends of
  // the diagonal (the diamond read from v swaps a<->b, c<->d): sin(angle) is
  // proportional to sqrt(Ha)*Q + P*sqrt(Hd), with (P, Q) the endpoint's two
  // law-of-cosines numerators.  Both tests share one evaluation of the two
  // Heron roots, so the arithmetic per endpoint is identical.  (The exact
  // sibling spells this as the sigma-involution pair -- DiamondSq -- where
  // recomputation is free of FP order effects; here the shared-root spelling
  // is deliberate, per the R-B decision.)
  bool is_convex() const {
    using namespace delaunay_detail;
    double e2 = e*e;
    double Ha = heron_product(e, a, b), Hd = heron_product(e, c, d);
    double sHa = (Ha > 0) ? std::sqrt(Ha) : 0;
    double sHd = (Hd > 0) ? std::sqrt(Hd) : 0;
    auto convex_at = [&](double P, double Q) {
      return sHa * Q + P * sHd > convexity_band;
    };
    return convex_at(e2 + a*a - b*b, e2 + c*c - d*d)      // at u
        && convex_at(e2 + b*b - a*a, e2 + d*d - c*c);     // at v
  }

  // Length of BD, the other diagonal.
  double flipped_length() const {
    using delaunay_detail::heron_product;
    // f^2 = a^2 + c^2 - (PQ - sqrt(Ha*Hd)) / (2e^2)
    double e2 = e*e, a2 = a*a, b2 = b*b, c2 = c*c, d2 = d*d;
    double P = e2 + a2 - b2;
    double Q = e2 + c2 - d2;
    double Ha = heron_product(e, a, b), Hd = heron_product(e, c, d);
    double sqrtHH = (Ha > 0 && Hd > 0) ? std::sqrt(Ha * Hd) : 0;
    double f2 = a2 + c2 - (P * Q - sqrtHH) / (2 * e2);
    return (f2 > 0) ? std::sqrt(f2) : 0;
  }

  // The squared-length diamond, entered through the integrality trust
  // boundary: nullopt unless every length squares to an integer within
  // lsq_integrality_band (relative).  Post-reduction utilities that only
  // hold doubles enter the exact regime here; the reduction itself carries
  // integer squared lengths and never rounds.
  std::optional<DiamondSq> squared() const {
    using delaunay_detail::lsq_integrality_band;
    long long L[5]; const double v[5] = {e, a, b, c, d};
    for (int i = 0; i < 5; i++) {
      const double sq = v[i] * v[i];
      L[i] = (long long)std::llround(sq);
      if (std::abs(sq - (double)L[i]) > lsq_integrality_band * std::max(1.0, sq))
        return std::nullopt;
    }
    return DiamondSq{L[0], L[1], L[2], L[3], L[4]};
  }

  // Cocircular ("tight") test: cot(angle_B) + cot(angle_D) == 0 exactly,
  // i.e. the four points u, v, B, D are concyclic on the surface.  In this
  // case both triangulations of the diamond are equally-valid Delaunay
  // refinements of the same cell.  F == 0 in exact integer arithmetic
  // (DiamondSq; see CANONICAL-TESSELATION.md for the derivation).  False
  // when the lengths do not square to integers (not the exact regime) or
  // the diamond is degenerate.
  bool is_cocircular_exact() const {
    auto ds = squared();
    return ds && ds->is_cocircular();
  }

  // Floating-point cocircular test for general (non-equilateral) metrics,
  // where length-squared is not integer so the exact predicate above does not
  // apply: tight iff |cot(angle_B) + cot(angle_D)| < tol.  Scale-invariant
  // (cotangents are dimensionless), so tol is a pure angle threshold.
  // cot_opposite returns +/-1e15 on a degenerate triangle, which never lands
  // within a sane tol, so degenerate diamonds are correctly reported non-tight.
  bool is_cocircular(double tol) const {
    using delaunay_detail::cot_opposite;
    double cotB = cot_opposite(e, a, b);
    double cotD = cot_opposite(e, c, d);
    return std::abs(cotB + cotD) < tol;
  }
};

// ============================================================================
// Fan polygon: isometric 2D development of a flat vertex's star.
//
// A flat vertex (cone angle = 2pi) has a star that unfolds without overlap
// into a planar disk.  The cumulative angle parameterization gives polar
// coordinates (spokes[i], cum[i]) for each boundary vertex; the exact
// regime additionally develops the boundary onto Z[w] (P[i], apex at the
// origin -- develop_fan_lattice).
//
//         nb[1]
//        / | \
//       /  |  \     spokes[i] = |v - nb[i]|
//      /   |   \    rims[i]   = |nb[i] - nb[(i+1)%k]|
//     v----+----    cum[i]    = sum of fan angles 0..i-1
//      \   |   /
//       \  |  /
//        \ | /
//         nb[0]
//
// Span-backed over caller workspace, capacity k_max (the declared max star
// degree); k <= k_max is the active size for this call.  Pure aggregate.
// ============================================================================
struct FanPolygon {
  int k = 0;                  // number of fan vertices (= star degree)
  std::span<int>    nb;       // [k_max] neighbor vertex IDs, CCW order
  std::span<int>    spoke_he; // [k_max] spoke half-edges: v -> nb[i]
  std::span<int>    inner_rim;// [k_max] inner rim half-edges: nb[i] -> nb[i+1]
  std::span<double> spokes;   // [k_max] spoke lengths
  std::span<double> rims;     // [k_max] rim edge lengths
  std::span<double> cum;      // [k_max+1] cumulative fan angles [0 .. 2pi]
  std::span<Eisenstein> P;    // [k_max] exact lattice development (exact
                              //         regime only, else untouched)

  // 2D fan coordinates of boundary vertex i.
  double x(int i) const { return spokes[i] * std::cos(cum[i]); }
  double y(int i) const { return spokes[i] * std::sin(cum[i]); }

  // Diagonal length between fan boundary vertices, as Euclidean distance in
  // the isometric development.
  double diag_length(int from, int to) const {
    double angle = (to > from) ? cum[to] - cum[from]
                               : (cum[k] - cum[from]) + cum[to];
    double sf = spokes[from], st = spokes[to];
    double len2 = sf*sf + st*st - 2*sf*st*std::cos(angle);
    return (len2 > 0) ? std::sqrt(len2) : 0;
  }

  // Signed area of triangle (pp, pi, pn) in fan coordinates.  Positive means
  // CCW orientation (valid ear).
  double ear_area(int pp, int pi, int pn) const {
    double rp = spokes[pp], ri = spokes[pi], rn = spokes[pn];
    double tp = cum[pp], ti = cum[pi], tn = cum[pn];
    return rp*ri*std::sin(ti - tp) + ri*rn*std::sin(tn - ti) + rn*rp*std::sin(tp - tn);
  }
};

// Fan triangulation: the result of ear-clipping a fan polygon.  Capacity
// k_max spans with explicit counts (n_diagonals <= k-3, n_triangles <= k-2,
// both < k <= k_max, so once extract_fan validates k these cannot overflow).
// len carries both length coordinates; len.sq is a lattice norm in the
// exact regime, 0 (never read) in the banded one.
struct FanTriangulation {
  struct Diagonal { int from, ear, to; Length len; };
  struct Triangle { int v0, v1, v2; };

  std::span<Diagonal> diagonals; int n_diagonals = 0;   // k-3 ear diagonals
  std::span<Triangle> triangles; int n_triangles = 0;   // k-2 ear triangles
};

// The outcome of one lattice development attempt (paper sec:exactness).
enum class Development { Closed, OffLattice, Unclosed };

namespace delaunay_detail {

// Ear acceptance test for candidate ear (pp, pi, pn): Meisters-style
// area/convexity/length.  Self-loop and multi-edge diagonals are legal
// delta-complex edges; splice_fan wires them by polygon POSITION, so
// repeated ids are unambiguous.  Returns the diagonal length if acceptable,
// else 0.
inline double ear_length_if_acceptable(const FanPolygon& fan, int pp, int pi, int pn) {
  if (fan.ear_area(pp, pi, pn) <= ear_area_band) return 0;
  double sub = (pn > pp) ? fan.cum[pn] - fan.cum[pp]
                         : (fan.cum[fan.k] - fan.cum[pp]) + fan.cum[pn];
  if (sub > std::numbers::pi_v<double> + ear_angle_band) return 0;
  double len = fan.diag_length(pp, pn);
  return (len > ear_length_floor) ? len : 0;
}

// Does the CCW sector from direction `from` to direction `to` around the
// origin subtend at most pi?  Exact, valid whenever the true sector lies in
// (0, 2pi) -- which holds for any sub-arc of a developed fan (every fan
// angle is strictly positive and the total is exactly 2pi): the plane wedge
// then decides the fan sub-arc.  wedge > 0 is (0, pi); wedge == 0 with
// negative dot is exactly pi (the accepted tie); wedge == 0 with positive
// dot would be a 0 / 2pi sector, excluded by the premise (defensive
// reject).  Shared by the exact ear test and the tie-break side order.
inline bool lattice_sector_at_most_pi(Eisenstein from, Eisenstein to) {
  const long w = wedge(from, to);
  if (w > 0) return true;
  return w == 0 && dot2(from, to) < 0;
}

// Exact ear acceptance over the lattice development P of a flat star
// (develop_fan_lattice): candidate ear (pp, pi, pn) is acceptable iff the
// ear triangle is strictly CCW (wedge sign) and the fan sector from pp to
// pn subtends at most pi at the apex.  Returns the squared diagonal length
// (a lattice norm, > 0 whenever accepted -- a strict CCW ear has distinct
// endpoints), 0 = reject.
inline long long ear_diag_sq_if_acceptable(std::span<const Eisenstein> P,
                                           int pp, int pi, int pn) {
  if (wedge(P[pi] - P[pp], P[pn] - P[pp]) <= 0) return 0;  // not a strict CCW ear
  if (!lattice_sector_at_most_pi(P[pp], P[pn])) return 0;  // sector > pi
  return norm2_ll(P[pn] - P[pp]);
}

}  // namespace delaunay_detail
