#pragma once
#include <assert.h>
#include <stdlib.h>
#include <inttypes.h>
#include <cmath>
#include <cstdlib>
#include <algorithm>
#include <limits>
#include <vector>
#include <optional>
#include <stdexcept>
#include <string>
#include <numeric>

using namespace std;

/*
  Multiplication:
  (a,b) = a*1 + b*w -> (a,b)*(c,d) = (a*1+b*w)(c*1+d*w) = (ac*1 + (bc+ad)*w + b*d*(w-1) = (ac-bd)*1 + (bc+ad+bd)*w
  
  Inverse: (a+b,-b)/(a^2+ab+b^2) = complex_conj(a,b)/norm2(a,b)
  
*/
class Eisenstein: public pair<int,int> {
public:
  Eisenstein(int a=0, int b=0) : pair<int,int>(a,b) {}
  Eisenstein(const pair<int,int>& x) : pair<int,int>(x.first,x.second) {}
  Eisenstein(const pair<double,double>& x) : pair<int,int>(round(x.first-x.second/sqrt(3)), round(2*x.second/sqrt(3)))
  { }
  Eisenstein operator*(const int y) const { return Eisenstein(first*y,second*y); }
  Eisenstein operator/(const int y) const { return Eisenstein(first/y,second/y); }
  Eisenstein operator*(const Eisenstein& y) const { 
    int a = first, b = second, c = y.first, d = y.second;
    return Eisenstein(a*c-b*d, b*c+(a+b)*d); 
  }
  vector<Eisenstein> operator*(const vector<Eisenstein>& v) const {
    vector<Eisenstein> result(v.size());
    for(size_t i=0;i<v.size();i++) result[i] = (*this)*v[i];
    return result;
  }
  Eisenstein operator+(const Eisenstein& y) const { return Eisenstein(first+y.first,second+y.second); }
  Eisenstein operator-(const Eisenstein& y) const { return Eisenstein(first-y.first,second-y.second); }
  Eisenstein operator-()                    const { return Eisenstein(-first,-second); }   
  pair<double,double>    operator-(const pair<double,double>& y)    const { return pair<double,double>(first-y.first,second-y.second); } 
  pair<double,double>    operator+(const pair<double,double>& y)    const { return pair<double,double>(first+y.first,second+y.second); } 
  Eisenstein& operator*=(const Eisenstein& y) { return (*this = *this * y); }
  Eisenstein& operator+=(const Eisenstein& y) { first += y.first; second += y.second; return *this; }
  Eisenstein& operator-=(const Eisenstein& y) { first -= y.first; second -= y.second; return *this; }

  static Eisenstein unit[7];
  
  bool isUnit() const { return norm2() == 1;  }

  // w-reflection: conjugation in Z[w] in its native (1, w) basis,
  Eisenstein eis_conj()  const { return Eisenstein(first,-second);  }
  // i-reflection: complex conjugation in C, restricted to Z[w].
  Eisenstein complex_conj() const { return Eisenstein((first+second), -second); }
  Eisenstein GCtransform(int k, int l) const {  return Eisenstein(k,l) * (*this);  }
  Eisenstein affine(const Eisenstein& x0, const Eisenstein w) const {
    Eisenstein x(*this);
    //    cout << x0 << " + " << w << " * " << *this << " = " << (x0+(w*(*this))) <<endl;
    return x0+w*x;
  }
  //  
  // (-1,1)   \ /  (0,1)
  // (-1,0)  --x-- (1,0)
  // ( 0,-1)  / \  (1,-1)
  // 
  Eisenstein nextCW()    const { return (*this) * Eisenstein(1,-1); }
  Eisenstein nextCCW()   const { return (*this) * Eisenstein(0,1);  }
  Eisenstein transpose() const { return Eisenstein(second,first);   }

  int unit_angle() const {
    assert(norm2() == 1);
    switch(first*10+second){
    case  1*10 + 0: return 0;	
    case  0*10 + 1: return 1;	
    case -1*10 + 1: return 2;	
    case -1*10 + 0: return 3;  	
    case  0*10 - 1: return 4;
    case  1*10 - 1: return 5;	
    default:
      abort();
    }
  }

  int nearest_unit_angle() const {
    int l = first, k = second;
    if(l>=0 && k<=0 && -k < l) return 0;
    if(l>=0 && k>0)            return 1;
    if(l<0  && k>0 && -l <= k) return 2;
    if(l<0  && k>=0 && -l > k) return 3;
    if(l<=0 && k<0)            return 4;
    if(l>0 && k<=0 && -k >= l) return 5;
    return -1;
  }


  pair<double,double> coord() const { 
    return make_pair(first + second/2., sqrt(3/4.)*second);
  }
  int norm2() const { return first*first + first*second + second*second; }
  double norm() const { return sqrt(norm2()); }

  Eisenstein abs() const { return Eisenstein(::abs(first),::abs(second)); }

  Eisenstein div(const Eisenstein& y) const {
    // Naive, possibly non-robust algorithm
    Eisenstein z(*this * y.complex_conj());

    double re = z.first, im = z.second;
    re /= y.norm2(); im /= y.norm2();

    return Eisenstein(round(re),round(im));
  }

  Eisenstein mod(const Eisenstein& y) const {
    Eisenstein q(this->div(y));
    return (*this) - q*y;
  }

  static pair<double,double> average(const vector<Eisenstein>& xs)
  {
    double re = 0, im = 0;
    for(size_t i=0;i<xs.size();i++){ re += xs[i].first; im += xs[i].second; }
    re /= xs.size(); im /= xs.size();
    
    return {re,im};
  }

  static Eisenstein gcd(Eisenstein a, Eisenstein b)  {
    Eisenstein c;
    while ( a != Eisenstein(0,0) ) {
      c = a; a = b.mod(a);  b = c;
    }
    return b;
  }

  static Eisenstein gcd(const vector<Eisenstein>& xs) {
    if(xs.empty()) return Eisenstein(1,0);

    Eisenstein d(xs[0]);
    for(size_t i=1;i<xs.size();i++) d = gcd(xs[i],d);
    return d;
  }


  static int turn(const Eisenstein& a, const Eisenstein& b, const Eisenstein& c) {
    int x = (b.first-a.first)*(c.second-a.second) - (b.second-a.second)*(c.first-a.first);
    return x < 0 ? -1 : (x == 0 ? 0 : 1);
  }

};

namespace std {
  template<> struct hash<Eisenstein> { // Vectors of integers smaller than 32 bit
    size_t operator()(Eisenstein const &f) const {
      uint64_t combined_int = (uint64_t(f.first)<<32) + f.second;
      return std::hash<uint64_t>()(combined_int);
    }
  };
}


// =====================================================================
// Free-function helpers on Z[w].
//
// Things that are about Z[w] but read more naturally as binary or
// query operations than as methods on a single Eisenstein.
// =====================================================================

// d-th sixth-root of unity, with mod-6 wrap on d (handles d < 0 too).
// Eisenstein::unit[i] is the same value when 0 <= i <= 6 (where i=6 is
// a duplicate of i=0); use unit_direction when d may be out of range.
inline Eisenstein unit_direction(int d) {
  return Eisenstein::unit[((d % 6) + 6) % 6];
}

// Integer signed area in shear coords:  u.first*v.second - u.second*v.first.
// Sign matches the Cartesian wedge (Cartesian-wedge = (sqrt(3)/2) * this),
// so trigonometry-free orientation tests are exact.  Companion to
// Eisenstein::turn(a, b, c), which returns the SIGN of wedge(b-a, c-a).
inline long wedge(Eisenstein u, Eisenstein v) {
  return (long)u.first * v.second - (long)u.second * v.first;
}

// Sign-exact 2*Re(a * conj(b)): positive iff the angle between a and b is
// acute.  Companion to wedge (the sine sign); together they totally order
// directions around the circle.
inline long dot2(Eisenstein a, Eisenstein b) {
  const Eisenstein z = a * b.complex_conj();
  return 2L * z.first + z.second;
}

// a and b point the same way (parallel, same sense).
inline bool same_dir(Eisenstein a, Eisenstein b) { return wedge(a, b) == 0 && dot2(a, b) > 0; }

// v reduced to its primitive lattice direction (divide out the gcd).
inline Eisenstein primitive_dir(Eisenstein v) {
  const int g = std::gcd(std::abs(v.first), std::abs(v.second));
  return g ? Eisenstein(v.first / g, v.second / g) : v;
}

// First unit strictly CCW (resp. CW) of a non-zero direction v.  Always
// exists: the gap between consecutive units is 60 degrees, so exactly one
// unit is the nearest strictly on that side.
inline Eisenstein first_ccw(Eisenstein v) {
  Eisenstein best; bool have = false;
  for (int i = 0; i < 6; i++) {
    const Eisenstein u = Eisenstein::unit[i];
    if (wedge(v, u) <= 0) continue;                     // not strictly CCW of v
    if (!have || wedge(u, best) > 0) { best = u; have = true; }
  }
  if (!have) throw std::logic_error("first_ccw: null direction");
  return best;
}
inline Eisenstein first_cw(Eisenstein v) {
  Eisenstein best; bool have = false;
  for (int i = 0; i < 6; i++) {
    const Eisenstein u = Eisenstein::unit[i];
    if (wedge(v, u) >= 0) continue;                     // not strictly CW of v
    if (!have || wedge(u, best) < 0) { best = u; have = true; }
  }
  if (!have) throw std::logic_error("first_cw: null direction");
  return best;
}

inline long long isqrt_exact(long long n);   // defined below (numeric section)

// The norm-form scan's per-b leg: the a >= 0 with a^2 + a*b + b^2 == N, if
// any.  THE solver every norm-representative entry point runs.  Header-
// inline and allocation-free (run-path callers include the exact-metric fan
// development).
inline std::optional<int> a_of_norm_leg(int N, int b) {
  const long long s = isqrt_exact(4LL * N - 3LL * b * b);
  if (s < 0) return std::nullopt;
  const long long two_a = s - b;
  if (two_a < 0 || two_a % 2) return std::nullopt;
  return (int)(two_a / 2);
}

// Allocation-free range over the sector-0 representatives (a >= 0, b >= 0)
// of norm N, in b-order.  THE one spelling of the norm-form scan; empty for
// N < 0 (not a norm).  Distinct unit-orbit representatives of the same norm
// are NOT lattice-equivalent (split primes), so orbit-sensitive callers --
// e.g. the exact fan development's anchor search -- iterate all of them.
class Sector0Reps {
  int N;
public:
  explicit Sector0Reps(int N_) : N(N_) {}
  struct sentinel {};
  struct iterator {
    int N, b;
    Eisenstein rep;
    bool live;
    explicit iterator(int N_) : N(N_), b(0), live(false) { seek(); }
    void seek() {
      for (live = false; 3LL * b * b <= 4LL * N; b++)
        if (auto a = a_of_norm_leg(N, b)) { rep = {*a, b}; live = true; return; }
    }
    Eisenstein operator*() const { return rep; }
    iterator& operator++() { b++; seek(); return *this; }
    bool operator!=(sentinel) const { return live; }
  };
  iterator begin() const { return iterator(N); }
  sentinel  end()   const { return {}; }
};

// First sector-0 representative of norm N, or nullopt when N is not an
// Eisenstein norm.
inline std::optional<Eisenstein> first_rep_of_norm(int N) {
  for (Eisenstein z : Sector0Reps(N)) return z;
  return std::nullopt;
}

// Some Eisenstein integer (a, b) with a >= 0, b >= 0 and
// a^2 + a*b + b^2 == N.  Returns the first sector-0 representative
// found by scanning b = 0, 1, ... (= first_rep_of_norm); throws
// std::logic_error if no solution exists.
Eisenstein eisenstein_of_norm(int N);

// Enumerate ALL sector-0 Eisenstein reps (a >= 0, b >= 0) of norm N.
// Generic norms return 1 entry; split-prime norms return 2 entries in
// distinct rotation orbits.
std::vector<Eisenstein> sector0_reps_of_norm(int N);


// =====================================================================
// Rational points of the Eisenstein plane.
// =====================================================================

// A point num/den (den > 0) of Q(w), the field of fractions of Z[w]:
// the rational closure that EXACT lattice geometry lives in.  Midpoints
// of lattice segments, barycentric samples with integer weights, and
// crossing points of lattice segments all lie in Q(w), never further.
// Construction canonicalizes -- sign moved into num, the gcd of
// {|num.first|, |num.second|, den} divided out -- so component-wise
// equality IS field equality.  Trivially copyable.
struct EisensteinRational {
  Eisenstein num{0, 0};
  long       den = 1;

  EisensteinRational(Eisenstein n = {}, long d = 1) : num(n), den(d) {
    if (den == 0) throw std::logic_error("EisensteinRational: zero denominator");
    if (den < 0) { den = -den; num = -num; }
    if (den == 1) return;                       // already canonical
    const long g = std::gcd(std::gcd(std::labs((long)num.first),
                                     std::labs((long)num.second)), den);
    if (g > 1) { num = Eisenstein(num.first / (int)g, num.second / (int)g); den /= g; }
  }

  bool is_integral() const { return den == 1; }

  bool operator==(const EisensteinRational& o) const { return num == o.num && den == o.den; }
  bool operator!=(const EisensteinRational& o) const { return !(*this == o); }

  // Cartesian coordinates, for display/diagnostics only (exact code never
  // leaves Q(w)).
  pair<double, double> coord() const {
    const auto c = num.coord();
    return { c.first / den, c.second / den };
  }
};


// =====================================================================
// Lattice self-isometries of Z[w].
// =====================================================================

// An isometry of the plane mapping the Eisenstein lattice onto itself:
//
//   T(z) = t + u * z                    (reflect == false)
//        = t + u * z.complex_conj()     (reflect == true)
//
// with u one of the 6 Eisenstein units and t in Z[w].  These are ALL of
// them (the p6m crystallographic restriction): every plane isometry is
// affine; one carrying the lattice onto itself sends 0 to a lattice
// point, so t is in Z[w], and its linear part preserves the lattice and
// lengths, hence lies in the point group D6 = {6 units} x {conj}.
struct LatticeIsometry {
  Eisenstein u{1, 0}, t{0, 0};
  bool       reflect = false;

  Eisenstein apply(Eisenstein z) const {
    return t + u * (reflect ? z.complex_conj() : z);
  }

  // Exact on rational points, and PRESERVES CANONICAL FORM: unit
  // multiplication and conjugation act on the (a, b) components as
  // unimodular integer matrices, so they preserve gcd(|a|, |b|, den);
  // adding t*den changes num only by multiples of den.  Hence any common
  // divisor of the image's components and den already divided the
  // preimage's -- 1 on canonical input -- and the constructor's gcd pass
  // is a no-op.
  // Width contract: Eisenstein components are int (see the sector-geometry
  // banner); t*den must fit, so den is checked against INT_MAX -- a
  // violated width bound throws instead of silently truncating.
  EisensteinRational apply(const EisensteinRational& p) const {
    if (p.den > std::numeric_limits<int>::max())
      throw std::logic_error("LatticeIsometry::apply: denominator exceeds int width");
    const Eisenstein n = reflect ? p.num.complex_conj() : p.num;
    return EisensteinRational(t * (int)p.den + u * n, p.den);
  }

  // Composition in FUNCTION order: (A*B)(z) == A(B(z)).  Closed form,
  // with s1 = conj^{r1} the first factor's reflection:
  //   r = r1 ^ r2,  u = u1 * s1(u2),  t = t1 + u1 * s1(t2).
  LatticeIsometry operator*(const LatticeIsometry& o) const {
    const Eisenstein su = reflect ? o.u.complex_conj() : o.u;
    const Eisenstein st = reflect ? o.t.complex_conj() : o.t;
    return { u * su, t + u * st, reflect != o.reflect };
  }

  // Group inverse.  Rotation case: T^-1(z) = u^-1 z - u^-1 t with
  // u^-1 = complex_conj(u) (a unit's inverse is its conjugate).
  // Reflection case (T is an involution twist): w = conj(u^-1) conj(z)
  // - conj(u^-1 t), and conj(u^-1) = u, so T^-1 = { u, -(u * conj(t)),
  // reflect } -- reflect is preserved.
  LatticeIsometry inverse() const {
    if (!reflect) {
      const Eisenstein ui = u.complex_conj();
      return { ui, -(ui * t), false };
    }
    return { u, -(u * t.complex_conj()), true };
  }

  bool operator==(const LatticeIsometry& o) const {
    return u == o.u && t == o.t && reflect == o.reflect;
  }
  bool operator!=(const LatticeIsometry& o) const { return !(*this == o); }
};

// The unique orientation-preserving lattice isometry carrying the
// directed segment a_from -> b_from onto a_to -> b_to.  The segments
// must have equal nonzero norm2, and the aligning rotation must be a
// unit (both hold whenever the segments are honest copies of one
// lattice segment); violations are deep invariants -> std::logic_error.
LatticeIsometry isometry_from_segments(Eisenstein a_from, Eisenstein b_from,
                                       Eisenstein a_to,   Eisenstein b_to);

// Find the LINEAR (t = 0) lattice isometry T with T(z_from) == z_to.
// Both inputs must have the same norm.  Exactly one branch (rotation,
// reflection) gives a unit by the D6 symmetry of Z[w]; align returns
// that branch (rotation preferred when both work, e.g. on the mirror
// axes).  Throws std::logic_error on norm mismatch or non-divisibility
// (deep invariants; same channel as isometry_from_segments).
LatticeIsometry align(Eisenstein z_from, Eisenstein z_to);


// =====================================================================
// Triangle / sector geometry on Z[w].
//
// Exact-integer primitives for unfolding a metric delta-complex into the
// Eisenstein plane: place a triangle's apex from its three squared edge
// lengths, and the closed angular sectors that test whether a developed
// straight line stays inside a face strip.  Promoted from the delta-complex
// surface-metric project (validated on every fullerene dual C20-C160,
// 211,203,353 isomers, 0 failures); see the project's
// DELTA-COMPLEX-SURFACE-METRIC.md for the algorithm and proofs.
//
// int components suffice while unfolded-plane positions stay small (fullerene
// iDTs); widen Eisenstein to 64-bit if a component can exceed ~26000.
// =====================================================================

// Exact integer square root: the s >= 0 with s*s == n, or -1 if n is not a
// perfect square (n < 0 yields -1).
inline long long isqrt_exact(long long n) {
  if (n < 0) return -1;
  // round(sqrt(n)) lands within 1 of the true integer root; probe +-2 as a
  // safe margin against double rounding, and confirm exactness with s*s == n.
  long long t = (long long)std::round(std::sqrt((double)n));
  for (long long d = -2; d <= 2; d++) {
    long long s = t + d;
    if (s >= 0 && s*s == n) return s;
  }
  return -1;
}

// |z|^2 = a^2 + a*b + b^2 in int64: the width contract for components that
// exceed the int-product envelope of Eisenstein::norm2().
inline long long norm2_ll(Eisenstein z) {
  const long long a = z.first, b = z.second;
  return a*a + a*b + b*b;
}

// The Heron product in SQUARED-length coordinates, exact over the integers:
// for squared sides x = a^2, y = b^2, z = c^2,
//   H(x,y,z) = 2(xy + yz + zx) - (x^2 + y^2 + z^2) = 16*Area^2,
// equivalently 4xy - (x + y - z)^2.  Negative iff the triangle inequality
// fails.
inline long long heron_product_sq(long long x, long long y, long long z) {
  return 2*(x*y + y*z + x*z) - (x*x + y*y + z*z);
}

// The lattice area number tau of a triangle with integer squared sides: on
// Z[w] every triangle satisfies H = 16*Area^2 = 3*tau^2 with tau = |wedge|
// of its edge vectors, an integer (the Heron-Eisenstein identity, paper
// prop:heron-eisenstein).  Returns tau >= 0 when H = 3*(perfect square),
// else -1.  NOTE this is the NECESSARY condition for a lattice triangle,
// not sufficient (lattice_tau(5,5,5) = 5 although 5 is not a Loeschian
// norm); sufficiency is decided downstream by place_third_eis_total's
// divisibility test / an empty Sector0Reps.  tau == 0 iff degenerate.
inline long long lattice_tau(long long x, long long y, long long z) {
  const long long H = heron_product_sq(x, y, z);
  if (H < 0 || H % 3) return -1;
  return isqrt_exact(H / 3);
}

// True iff P lies strictly between the origin and C (collinear, 0 < s < 1 in
// P = s*C).  Integer-exact.
inline bool is_on_open_segment(Eisenstein P, Eisenstein C) {
  if (wedge(C, P) != 0) return false;                // collinear with 0 -> C
  if (P.first == 0 && P.second == 0) return false;
  if (P == C)                        return false;
  if (C.first == 0 && C.second == 0) return false;
  // Collinear, so P = s*C; read s off C's dominant ("major") component to
  // avoid division: same sign as C and strictly smaller magnitude there.
  bool major_a = (std::llabs(C.first) >= std::llabs(C.second));
  long long pm = major_a ? P.first : P.second;
  long long cm = major_a ? C.first : C.second;
  if (cm == 0) return false;
  if ((pm > 0) != (cm > 0)) return false;            // 0 < s
  return std::llabs(pm) < std::llabs(cm);            // s < 1
}

// A Sector (R, L) is the closed CCW angular sector from R to L at the origin.
//   - is_empty()       iff wedge(R, L) <= 0
//   - contains(d)      iff wedge(R, d) >= 0 && wedge(d, L) >= 0
//   - narrow_with(r,l) tightens each side independently against an inner
//                      candidate.
struct Sector { Eisenstein R, L;
  bool is_empty()      const { return wedge(R, L) <= 0; }
  bool contains(Eisenstein d) const {
    return wedge(R, d) >= 0 && wedge(d, L) >= 0;
  }
  // d strictly inside the OPEN sector (R, L) -- excludes both boundary rays.
  // Reflex-safe (sectors up to 300 degrees); throws on a degenerate sector
  // where R, L are parallel and same-sense (there is no interior).
  bool strictly_inside(Eisenstein d) const {
    const long rl = wedge(R, L);
    if (rl > 0) return wedge(R, d) > 0 && wedge(d, L) > 0;      // convex (< 180)
    if (rl < 0) return !(wedge(L, d) > 0 && wedge(d, R) > 0);   // reflex: complement of (L, R)
    if (dot2(R, L) < 0) return wedge(R, d) > 0;                 // exactly 180 (R, L antiparallel)
    throw std::logic_error("Sector::strictly_inside: degenerate sector (R, L parallel, same sense)");
  }
  // Reflex-safe HALF-OPEN membership (unlike contains(), which assumes a convex
  // <= 180-degree sector): the shared boundary ray belongs to exactly one side.
  bool contains_oc(Eisenstein d) const {   // (R, L] -- open at R, closed at L
    if (same_dir(d, L)) return true;
    if (same_dir(d, R)) return false;
    return strictly_inside(d);
  }
  bool contains_co(Eisenstein d) const {   // [R, L) -- closed at R, open at L
    if (same_dir(d, R)) return true;
    if (same_dir(d, L)) return false;
    return strictly_inside(d);
  }
  Sector narrow_with(Eisenstein r_new, Eisenstein l_new) const {
    return { wedge(R, r_new) > 0 ? r_new : R,
             wedge(l_new, L) > 0 ? l_new : L };
  }
};

// Order two non-zero directions into a CCW pair (R, L).  Colinear inputs
// (wedge == 0) return (B, A): the sector is then empty, but R, L keep their
// directions so contains() still answers the colinear half-plane correctly.
inline Sector ccw_order(Eisenstein A, Eisenstein B) {
  return wedge(A, B) > 0 ? Sector{A, B} : Sector{B, A};
}

// Entry sector for a child face's two adjacent corner positions.  nullopt iff
// there is no entry direction (both corners at origin, or two non-origin
// corners colinear).  One corner at origin collapses to {X, X}.
inline std::optional<Sector> entry_sector(Eisenstein A, Eisenstein B) {
  bool A_o = (A.first == 0 && A.second == 0);
  bool B_o = (B.first == 0 && B.second == 0);
  if (A_o && B_o)       return std::nullopt;
  if (A_o)              return Sector{B, B};
  if (B_o)              return Sector{A, A};
  if (wedge(A, B) == 0) return std::nullopt;
  return ccw_order(A, B);
}

// Apex of the CCW unit triangle on base a -> b: the +60-degree rotation of
// the base direction around its origin.  This is place_third_eis(a, b,
// 1, 1, 1, +1) specialized to the unit triangle, spelled directly because
// the unit apex always lies on the lattice (no optional, no throw).
inline Eisenstein unit_apex(Eisenstein a, Eisenstein b) { return a + (b - a).nextCCW(); }

// Thrown by place_third_eis for a triangle with no Z[w] apex -- a degenerate /
// non-Loeschian face, which cannot occur on a valid delta-complex.  Carries
// its inputs so the failing placement can be replayed from the caught
// exception.
struct DegenerateTriangle : std::logic_error {
  Eisenstein A, B;
  int abs2, acs2, bcs2, chirality;
  DegenerateTriangle(Eisenstein A_, Eisenstein B_,
                     int abs2_, int acs2_, int bcs2_, int chirality_)
    : std::logic_error(
        "place_third_eis: degenerate triangle  A=(" + std::to_string(A_.first)
        + "," + std::to_string(A_.second) + ") B=(" + std::to_string(B_.first)
        + "," + std::to_string(B_.second) + ")  |AB|2=" + std::to_string(abs2_)
        + " |AC|2=" + std::to_string(acs2_) + " |BC|2=" + std::to_string(bcs2_)
        + " chirality=" + std::to_string(chirality_)),
      A(A_), B(B_), abs2(abs2_), acs2(acs2_), bcs2(bcs2_), chirality(chirality_) {}
};

// The Z[w] displacement factor m+n*w of the apex over base A->B:
//   K = |AB|2+|AC|2-|BC|2  (law of cosines),  tau = lattice_tau of the
//   sides (16*Area^2 = 3*tau^2),  m = (K - chi*tau)/2,  n = chi*tau.
// nullopt iff these squared sides admit no Loeschian triangle.
//   @pre  chirality == +-1, abs2 > 0, acs2 >= 0, bcs2 >= 0.
inline std::optional<Eisenstein>
apex_factor(int abs2, int acs2, int bcs2, int chirality) {
  const long long K   = (long long)acs2 + abs2 - bcs2;   // law of cosines
  const long long tau = lattice_tau(abs2, acs2, bcs2);   // 16*Area^2 = 3*tau^2
  const long long n   = chirality * tau;
  if (tau < 0 || (K - n) % 2) return std::nullopt;
  return Eisenstein((int)((K - n) / 2), (int)n);          // m + n*w
}

// Total (non-throwing) apex placement: as place_third_eis below, with ALL
// failure modes -- degenerate base (abs2 <= 0), no Z[w] triangle for these
// sides, or apex off-lattice for this chirality -- collapsed into nullopt.
// The run-path form for no-throw callers (the DCEL exact-metric machinery
// converts nullopt to a Status trip); place_third_eis keeps the loud
// degenerate / soft off-lattice distinction for callers that prune on
// chirality.  TOTAL means total: a zero base must yield nullopt here, not
// reach the divisions below.
inline std::optional<Eisenstein>
place_third_eis_total(Eisenstein A, Eisenstein B, int abs2, int acs2, int bcs2, int chirality) {
  if (abs2 <= 0) return std::nullopt;                             // degenerate base
  const auto m_nw = apex_factor(abs2, acs2, bcs2, chirality);     // m + n*w, or empty
  if (!m_nw) return std::nullopt;                                 // no Z[w] triangle
  const Eisenstein num = (*m_nw) * (B - A);                       // (m+n*w)*(B-A)
  if (num.first % abs2 || num.second % abs2) return std::nullopt; // off-lattice for this chi
  return A + num / abs2;                                          // C = A + .../|AB|2
}

// Apex C of the triangle on base A->B, exact in Z[w]:
//   C = A + (m + n*w)*(B - A) / |AB|2.
// nullopt = the apex misses the lattice for this chirality (a soft prune);
// throws DegenerateTriangle when the sides admit no Z[w] triangle.
//   @pre  chirality == +-1, abs2 > 0, acs2 >= 0, bcs2 >= 0,
//         abs2 == (B - A).norm2().
inline std::optional<Eisenstein>
place_third_eis(Eisenstein A, Eisenstein B, int abs2, int acs2, int bcs2, int chirality) {
  if (!apex_factor(abs2, acs2, bcs2, chirality))
    throw DegenerateTriangle(A, B, abs2, acs2, bcs2, chirality);
  return place_third_eis_total(A, B, abs2, acs2, bcs2, chirality);
}




