#pragma once

// The intrinsic-geometry vocabulary of the metric DCEL: the numeric
// primitives and tolerance bands (delaunay_detail), the Diamond -- the
// five-length local geometry around an edge -- with its predicates, and the
// fan development of a flat vertex's star (FanPolygon / FanTriangulation +
// the ear acceptance test).  Pure mathematics over lengths: no DCEL, no
// workspace, no I/O -- provable in isolation and reusable by any
// metric-surface code.  Device-includable.

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <numbers>
#include <span>
#include <type_traits>

// ============================================================================
// Intrinsic geometry primitives.
// ============================================================================
namespace delaunay_detail {

inline constexpr double two_pi = 2 * std::numbers::pi_v<double>;

// Tolerance bands of the floating-point predicates.  These absorb FP noise
// in the cotangent/Heron/development evaluations, nothing else: the
// mathematics is is_delaunay <=> cot-sum >= 0, strict convexity <=> both
// endpoint forms > 0, and the ear acceptance criteria of paper sec:ear-clip.
// Never widen a band to make a case pass.
inline constexpr double delaunay_band   = -1e-10;  // accept near-tight as tight
inline constexpr double convexity_band  =  1e-12;  // strictness margin
inline constexpr double cot_degenerate  =  1e15;   // sentinel for a degenerate triangle
inline constexpr double ear_area_band   =  1e-10;  // CCW ear area must exceed
inline constexpr double ear_angle_band  =  1e-10;  // reflex slack past pi
inline constexpr double ear_length_floor = 1e-15;  // degenerate-diagonal floor
inline constexpr double tie_cocircular_tol = 1e-12; // tie-break spoke test
inline constexpr double default_flat_tol   = 1e-6;  // exact-metric flatness

// Heron product: H(a,b,c) = (a+b+c)(-a+b+c)(a-b+c)(a+b-c) = 16*Area^2.
// Returns 0 if the triangle inequality is violated.
inline double heron_product(double a, double b, double c) {
  double s1 = -a + b + c;
  double s2 =  a - b + c;
  double s3 =  a + b - c;
  if (s1 < 0 || s2 < 0 || s3 < 0) return 0;
  return (a + b + c) * s1 * s2 * s3;
}

// The Heron product in SQUARED-length coordinates, exact over the integers:
// for squared sides x = a^2, y = b^2, z = c^2,
//   H(x,y,z) = 2(xy + yz + zx) - (x^2 + y^2 + z^2) = 16*Area^2.
// The integer form the exact cocircularity predicate rationalizes over.
inline long long heron_product_sq(long long x, long long y, long long z) {
  return 2*(x*y + y*z + x*z) - (x*x + y*y + z*z);
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
  // Heron roots, so the arithmetic per endpoint is identical.
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

  // Cocircular ("tight") test: cot(angle_B) + cot(angle_D) == 0 exactly,
  // i.e. the four points u, v, B, D are concyclic on the surface.  In this
  // case both triangulations of the diamond are equally-valid Delaunay
  // refinements of the same cell.  Uses exact integer arithmetic on
  // length-squared (valid when all five lengths square to non-negative
  // integers, e.g. equilateral triangulations and their flips).
  // See CANONICAL-TESSELATION.md for the derivation.
  bool is_cocircular() const {
    // Tight Delaunay: cot(angle_B) + cot(angle_D) == 0 exactly.  Equivalent
    // to s1 * area_2 + s2 * area_1 = 0 where s1 = a^2+b^2-e^2, s2 = c^2+d^2-e^2.
    // Squaring after sign-check: tight iff sign(s1) != sign(s2) and
    // s1^2 * H2 == s2^2 * H1 (with H = heron_product_sq = 16*area^2).  Or:
    // both s1, s2 == 0.  Done in integer length-squared arithmetic so the
    // predicate is exact for equilateral triangulations and any sequence of
    // flips.
    using delaunay_detail::heron_product_sq;
    long long Le = (long long)std::llround(e * e);
    long long La = (long long)std::llround(a * a);
    long long Lb = (long long)std::llround(b * b);
    long long Lc = (long long)std::llround(c * c);
    long long Ld = (long long)std::llround(d * d);
    long long s1 = La + Lb - Le;
    long long s2 = Lc + Ld - Le;
    if (s1 == 0 && s2 == 0) return true;
    if (s1 == 0 || s2 == 0) return false;
    if ((s1 > 0) == (s2 > 0)) return false;          // same sign: not tight
    return s1 * s1 * heron_product_sq(Le, Lc, Ld)
        == s2 * s2 * heron_product_sq(Le, La, Lb);
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
// coordinates (spokes[i], cum[i]) for each boundary vertex.
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
struct FanTriangulation {
  struct Diagonal { int from, ear, to; double length; };
  struct Triangle { int v0, v1, v2; };

  std::span<Diagonal> diagonals; int n_diagonals = 0;   // k-3 ear diagonals
  std::span<Triangle> triangles; int n_triangles = 0;   // k-2 ear triangles
};

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

}  // namespace delaunay_detail

