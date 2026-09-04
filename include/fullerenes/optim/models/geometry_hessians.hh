#pragma once

// =====================================================================
// Second-derivative geometry for the ExtWu force field: everything the
// matrix-free HVP  (grad^2 E) . v  needs beyond the library's first-
// derivative primitives.
//
// Strategy: FORWARD-MODE differentiation along the direction v.  Each
// term's Hessian-vector contribution is
//
//     (grad^2 E_t) . v = k/2 [ h''(d) (grad q . v) grad q
//                              + h'(d) (grad^2 q . v) ]
//
// and  grad^2 q . v  is the directional derivative of grad q -- so a
// dual-number pass through the GRADIENT formula yields it exactly (to
// roundoff), with no hand-derived second-derivative expressions to get
// wrong.  The gradient formulas here are:
//
//   * angle_grad<S>  -- a scalar-templated transcription of the library
//     coord3<double>::dangle (same expressions, same theta in {0, pi}
//     subgradient-0 guard);
//   * dihedral_grad<S> -- the compact Blondel-Karplus form of the
//     dihedral gradient, SIGN-MATCHED to the library convention
//     (coord3d::dihedral's atan2((F^ x bc^).G^, F^.G^), which is the
//     negative of the standard sign);
//
// both verified against coord3d::dangle / coord3d::ddihedral to
// roundoff in test_extwu_hvp (the pin that keeps these transcriptions
// tied to the library's central functions).  On graduation these
// templated forms can replace the double-only lib specializations.
//
// The bond term's second derivative is the textbook projector
// (I - u u^T)/r and is written directly.
//
// Deviation profile: h(d) = d^2 (hard) or d^2/sqrt(4+d^2) (soft), as
// wu_forcefield.cc's profile(); here with the second derivative
//   hard: h'' = 2,     soft: h'' = (32 - 4 d^2) (4 + d^2)^{-5/2}.
// =====================================================================

#include <cmath>

namespace optim {
namespace geom2 {

// --- Forward-mode dual scalar: value + directional derivative --------
struct Dual {
  double v = 0, d = 0;

  Dual() = default;
  Dual(double value) : v(value), d(0) {}
  Dual(double value, double deriv) : v(value), d(deriv) {}
};

inline Dual operator+(Dual a, Dual b) { return {a.v + b.v, a.d + b.d}; }
inline Dual operator-(Dual a, Dual b) { return {a.v - b.v, a.d - b.d}; }
inline Dual operator-(Dual a)         { return {-a.v, -a.d}; }
inline Dual operator*(Dual a, Dual b) { return {a.v * b.v, a.d * b.v + a.v * b.d}; }
inline Dual operator/(Dual a, Dual b) {
  const double inv = 1.0 / b.v;
  return {a.v * inv, (a.d - a.v * b.d * inv) * inv};
}
inline Dual sqrt(Dual a) {
  const double s = std::sqrt(a.v);
  return {s, a.d / (2 * s)};
}
inline bool operator>=(Dual a, double b) { return a.v >= b; }
inline bool operator<=(Dual a, double b) { return a.v <= b; }
inline bool operator==(Dual a, double b) { return a.v == b; }

// double passes through the same code paths.
inline double sqrt(double a) { return std::sqrt(a); }

// --- Scalar-templated 3-vector ---------------------------------------
template <class S>
struct Vec3 {
  S x{}, y{}, z{};

  Vec3 operator+(const Vec3& b) const { return {x + b.x, y + b.y, z + b.z}; }
  Vec3 operator-(const Vec3& b) const { return {x - b.x, y - b.y, z - b.z}; }
  Vec3 operator-() const { return {-x, -y, -z}; }
  Vec3 operator*(S s) const { return {x * s, y * s, z * s}; }
  S dot(const Vec3& b) const { return x * b.x + y * b.y + z * b.z; }
  Vec3 cross(const Vec3& b) const {
    return {y * b.z - z * b.y, z * b.x - x * b.z, x * b.y - y * b.x};
  }
  S norm() const { return geom2::sqrt(dot(*this)); }
};

using DVec = Vec3<Dual>;

// --- Angle gradient: transcription of coord3<double>::dangle ---------
// theta at b(0,0,0) between a and c; da, dc are grad_a theta, grad_c
// theta (grad_b = -(da + dc)).  The theta in {0, pi} guard takes the
// subgradient 0, exactly as the library.
template <class S>
inline void angle_grad(const Vec3<S>& a, const Vec3<S>& c,
                       Vec3<S>& da, Vec3<S>& dc) {
  const S L2 = a.dot(a);
  const S R2 = c.dot(c);
  const S M2 = (c - a).dot(c - a);
  const S den = S(2.0) * geom2::sqrt(L2 * R2);
  const S arg = (L2 + R2 - M2) / den;

  if (arg >= 1.0 || arg <= -1.0) { da = dc = Vec3<S>{}; return; }

  const Vec3<S> dM2__da = (a - c) * S(2.0);
  const Vec3<S> dL2__da = a * S(2.0);
  const Vec3<S> dden__da = dL2__da * (R2 / geom2::sqrt(L2 * R2));
  const Vec3<S> darg__da = dL2__da * (S(1.0) / den) - dM2__da * (S(1.0) / den)
                         - dden__da * ((L2 + R2 - M2) / (den * den));

  const Vec3<S> dM2__dc = (c - a) * S(2.0);
  const Vec3<S> dR2__dc = c * S(2.0);
  const Vec3<S> dden__dc = dR2__dc * (L2 / geom2::sqrt(L2 * R2));
  const Vec3<S> darg__dc = dR2__dc * (S(1.0) / den) - dM2__dc * (S(1.0) / den)
                         - dden__dc * ((L2 + R2 - M2) / (den * den));

  const S w = S(-1.0) / geom2::sqrt(S(1.0) - arg * arg);
  da = darg__da * w;
  dc = darg__dc * w;
}

// --- Dihedral gradient: compact form, library sign convention --------
// Inputs are the library's call shape (dihedral at a(0,0,0)):
//   B = b - a, C = c - a, D = d - a.
// Outputs db, dc, dd = grad_b phi, grad_c phi, grad_d phi
// (grad_a = -(db + dc + dd)), matching coord3d::ddihedral to roundoff.
//
// Derivation: with bond vectors b1 = B, b2 = C - B, b3 = D - C and
// normals F = b1 x b2, G = b2 x b3, the compact (Blondel-Karplus form)
// gradient of phi_std with p = (b1.b2)/|b2|^2, q = (b3.b2)/|b2|^2 is
//   grad_a = -( |b2| / |F|^2 ) F
//   grad_d = +( |b2| / |G|^2 ) G
//   grad_b = -grad_a - p grad_a + q grad_d
//   grad_c = -grad_d + p grad_a - q grad_d          (sum == 0)
// and the library's y = (F^ x b2^).G^ = -(F x G).b2^ makes
// phi_lib = -phi_std, so every gradient is negated.  The endpoint
// gradients and both cross-term signs are verified against
// coord3d::ddihedral at machine precision (DdihedralParity).  The
// degenerate-quadruple guard (|F| = 0 or |G| = 0 or |b2| = 0) takes the
// subgradient 0, exactly as the library.
template <class S>
inline void dihedral_grad(const Vec3<S>& B, const Vec3<S>& C, const Vec3<S>& D,
                          Vec3<S>& db, Vec3<S>& dc, Vec3<S>& dd) {
  const Vec3<S> b1 = B;
  const Vec3<S> b2 = C - B;
  const Vec3<S> b3 = D - C;
  const Vec3<S> F = b1.cross(b2);
  const Vec3<S> G = b2.cross(b3);

  const S F2 = F.dot(F), G2 = G.dot(G), n2sq = b2.dot(b2);
  if (F2 == 0.0 || G2 == 0.0 || n2sq == 0.0) {
    db = dc = dd = Vec3<S>{};
    return;
  }
  const S n2 = geom2::sqrt(n2sq);

  const Vec3<S> ga_std = F * (-(n2 / F2));         // grad_a phi_std
  const Vec3<S> gd_std = G * (n2 / G2);            // grad_d phi_std
  const S p = b1.dot(b2) / n2sq;
  const S q = b3.dot(b2) / n2sq;
  const Vec3<S> gb_std = -ga_std - ga_std * p + gd_std * q;
  const Vec3<S> gc_std = -gd_std + ga_std * p - gd_std * q;

  db = -gb_std;                                    // library sign
  dc = -gc_std;
  dd = -gd_std;
}

// --- Deviation profile with second derivative ------------------------
struct Profile2 {
  double h, dh, d2h;
};

inline Profile2 profile2(double d, bool soft) {
  if (!soft) return {d * d, 2 * d, 2.0};
  const double s = 1.0 / std::sqrt(4 + d * d);
  const double s3 = s * s * s;
  return {d * d * s, d * (8 + d * d) * s3, (32 - 4 * d * d) * s3 * s * s};
}

}  // namespace geom2
}  // namespace optim
