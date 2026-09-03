#pragma once

// View-level dense linear algebra: the allocation-free bodies of the LinAlg
// solver family, written ONCE over caller storage, shared by
//   - the owner API in dense_linalg.cc (thin allocating wrappers),
//   - matrix<double>::operator* (matrix.hh delegates its product here), and
//   - device-legal ports (parallel-primitives' Alexandrov solver), whose
//     per-isomer workspaces bind these same bodies over batch arenas.
// BLAS/LAPACK-free by design -- dense_linalg.hh's banner records why.
//
// Every body is the historical dense_linalg.cc / matrix.hh loop, moved --
// with two stated exceptions: negate() is NEW at this level (the
// device-legal span form of auxiliary.hh's vector unary minus), and
// matvec() generalized its loop bounds from the square A.m to the honest
// m x n (identical on every square input).  The independent pins are
// dense-linalg-test's frozen historical oracle (the pre-promotion bodies,
// byte-compared), its strided-vs-packed leg, and -- for the batch port's
// BINDING layer only -- the Alexandrov pilot's 5,771 + 100 isomer corpus.
//
// Matrix scratch parameters are raw spans holding a PACKED n x n block
// (stride n); each body views them internally, so callers stay span-shaped
// (the batch arenas' native currency) while the arithmetic reads as A(i,j).
// Nothing here allocates, throws, or does I/O.  @pre for every function:
// no argument aliases another (the copy orders assume disjoint storage).
//
// Relation to the other matrix vocabularies: sycl-headers' MDSpan is the
// SYCL kernel container family; BatchLAS's DenseMatView (DEVICE-READINESS
// 2.2) is the future batched-engine surface.  This view is the lib-internal
// solver vocabulary; unification is the 2.2 engine's decision to make.
//
// (The symmetric-eigen family -- jacobi_eig, SymEigen -- stays owner-level
// in dense_linalg.cc: all its consumers are host-only.)

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <span>
#include <type_traits>
#include <utility>

namespace LinAlg {

// Row-major matrix view over caller storage: the live m x n block of a
// buffer with row stride lda.  Trivially copyable; capturable by value.
// @inv lda >= n >= 0, m >= 0, data.size() >= (m ? (m-1)*lda + n : 0)
template <class T>
struct MatView {
  std::span<T> data;
  int m = 0, n = 0, lda = 0;

  constexpr MatView() = default;
  constexpr MatView(std::span<T> d, int m_, int n_, int lda_)
      : data(d), m(m_), n(n_), lda(lda_) {}
  // const-adding conversion (MatView<double> -> MatView<const double>).
  template <class U>
    requires std::is_convertible_v<U (*)[], T (*)[]>
  constexpr MatView(const MatView<U>& o)
      : data(o.data), m(o.m), n(o.n), lda(o.lda) {}

  T& operator()(int i, int j) const {
    return data[(std::size_t)i * lda + j];
  }
};
using MatConstView = MatView<const double>;

// --- Vector reductions (dense_linalg.cc's loops; spans are exact-sized:
//     the span IS the vector, no separate count) ---

// @pre a.size() == b.size()
inline double dot(std::span<const double> a, std::span<const double> b) {
  double s = 0;
  for (std::size_t i = 0; i < a.size(); i++) s += a[i] * b[i];
  return s;
}
inline double sum_sq(std::span<const double> v) { return dot(v, v); }
inline double norm(std::span<const double> v) { return std::sqrt(sum_sq(v)); }
// max |v_i|; a NaN entry poisons the result to +inf, so a NaN residual can
// never pass a "< tol" convergence test.
inline double max_abs(std::span<const double> v) {
  double m = 0;
  for (double x : v) {
    if (std::isnan(x)) return HUGE_VAL;
    m = std::max(m, std::fabs(x));
  }
  return m;
}
// v := -v elementwise (exact sign flip; the span form of the vector unary
// minus the owner expressions use).  New at this level -- see the banner.
inline void negate(std::span<double> v) {
  for (double& x : v) x = -x;
}

// --- Vector assignments (the span forms of the owner V expressions the
//     solver family writes: r = r - dr, r_trial = r + delta,
//     F = kappa - t1*kappa1, result = result + r_j*basis).  Elementwise
//     with one rounding per element, so a call site composed of these is
//     bit-identical to the owner expression it restates.  New at this
//     level, like negate() -- see the banner.
//     @pre for each: src.size() == dst.size() (and no aliasing, per the
//     file banner). ---

// dst := src
inline void copy_into(std::span<double> dst, std::span<const double> src) {
  for (std::size_t i = 0; i < dst.size(); i++) dst[i] = src[i];
}
// dst := -src (exact sign flip, like negate)
inline void neg_into(std::span<double> dst, std::span<const double> src) {
  for (std::size_t i = 0; i < dst.size(); i++) dst[i] = -src[i];
}
// dst += src
inline void add_into(std::span<double> dst, std::span<const double> src) {
  for (std::size_t i = 0; i < dst.size(); i++) dst[i] += src[i];
}
// dst -= src
inline void sub_into(std::span<double> dst, std::span<const double> src) {
  for (std::size_t i = 0; i < dst.size(); i++) dst[i] -= src[i];
}
// dst += s * src
inline void add_scaled(std::span<double> dst, double s, std::span<const double> src) {
  for (std::size_t i = 0; i < dst.size(); i++) dst[i] += s * src[i];
}
// dst -= s * src
inline void sub_scaled(std::span<double> dst, double s, std::span<const double> src) {
  for (std::size_t i = 0; i < dst.size(); i++) dst[i] -= s * src[i];
}
// Sum of entries -- term-for-term identical to the owner idiom
// dot(v, ones) (v[i] * 1.0 == v[i] exactly, same accumulation order).
inline double sum(std::span<const double> v) {
  double s = 0;
  for (double x : v) s += x;
  return s;
}

// --- Solver policy (calibrated constants; NOT neutral linear algebra) ---

// Residual energy E = 1/2 ||v||^2 (the Gauss-Newton objective's 1/2).
inline double energy(std::span<const double> v) { return 0.5 * sum_sq(v); }

// Step-acceptance floor: ||v||^2 > 1e-30, i.e. ||v|| > 1e-15.  On the
// SQUARED norm deliberately -- a clip-capped trust-radius step with
// 0 < ||v|| <= 1e-15 must still count as a step (a max|v| > 0 test takes
// the opposite branch there, and the Gauss-Newton bisection routes ~31
// solves per fallback through this predicate).
inline constexpr double STEP_SQ_FLOOR = 1e-30;
inline bool is_usable_step(std::span<const double> v) {
  const double s = sum_sq(v);
  return std::isfinite(s) && s > STEP_SQ_FLOOR;
}

// --- Matrix products ---

// out := A * v.
// @pre v.size() >= A.n, out.size() >= A.m
inline void matvec(MatConstView A, std::span<const double> v,
                   std::span<double> out) {
  for (int i = 0; i < A.m; i++) {
    double s = 0;
    for (int j = 0; j < A.n; j++) s += A(i, j) * v[j];
    out[i] = s;
  }
}

// @anchor matmul-ijk-order
// out := A * B, packed row-major A.m x B.n.  The i-j-k loop with a
// k-ascending scalar accumulator.  LOAD-BEARING: the Alexandrov solver's
// byte gates are validated against this association, and
// matrix<double>::operator* delegates here.  Do not reassociate, block,
// or vectorize this loop.
// @pre A.n == B.m, out.size() >= A.m * B.n
inline void matmul(MatConstView A, MatConstView B, std::span<double> out) {
  const MatView<double> C{out, A.m, B.n, B.n};
  for (int i = 0; i < A.m; i++)
    for (int j = 0; j < B.n; j++) {
      double x = 0;
      for (int k = 0; k < A.n; k++) x += A(i, k) * B(k, j);
      C(i, j) = x;
    }
}

// --- LU with partial pivoting: THE shared core ---
//
// Singular (an exact zero pivot) is ONE event; each entry point presents
// it as its caller needs:
//   row_reduce        -> LuStatus::Singular; M/sign/b left partially reduced
//   solve (view)      -> LuStatus::Singular; x := 0
//   solve (owner)     -> x := 0  (status dropped; callers guard with
//                        is_usable_step)
//   solve_with_sign   -> unexpected(LuFail::Singular)
//   det               -> 0.0  (the TRUE value of det, not a failure)
// The enum is the two-outcome encoding, per jacobi_eig's bool
// (dense_linalg.hh).
enum class LuStatus : int { Ok, Singular };
struct LuReduction {
  LuStatus status = LuStatus::Ok;
  int      sign   = 1;   // (-1)^{#row swaps} x prod sign(diag U)
};

// Partial-pivot row reduction of a COPY of A staged in the packed scratch
// M (capacity n*n): M holds U on and above the diagonal; the returned sign
// is sign(det A).  When b is non-empty the same swaps and forward
// elimination are applied to it, leaving the triangular system U x = b
// ready for back_substitute.  L is deliberately NOT formed and the
// permutation is not recorded -- only U, sign det, and the eliminated RHS
// are outputs (hence not "lu_decompose": a second RHS cannot be solved
// from the result).  On Singular, M, sign, and b are only partially
// reduced.
// @pre A.m == A.n (square); M.size() >= n*n; b empty or b.size() >= n
inline LuReduction row_reduce(MatConstView A, std::span<double> Mbuf,
                              std::span<double> b) {
  const int n = A.n;
  const MatView<double> M{Mbuf, n, n, n};
  for (int i = 0; i < n; i++)
    for (int j = 0; j < n; j++) M(i, j) = A(i, j);
  LuReduction lu;
  for (int c = 0; c < n; c++) {
    int p = c;
    for (int q = c + 1; q < n; q++)
      if (std::fabs(M(q, c)) > std::fabs(M(p, c))) p = q;
    if (M(p, c) == 0) { lu.status = LuStatus::Singular; return lu; }
    if (p != c) {
      for (int j = c; j < n; j++) std::swap(M(c, j), M(p, j));
      if (!b.empty()) std::swap(b[c], b[p]);
      lu.sign = -lu.sign;                        // row swap parity
    }
    if (M(c, c) < 0) lu.sign = -lu.sign;         // sign of diag(U)
    for (int q = c + 1; q < n; q++) {
      const double mult = M(q, c) / M(c, c);
      if (mult == 0) continue;
      for (int j = c + 1; j < n; j++) M(q, j) -= mult * M(c, j);
      if (!b.empty()) b[q] -= mult * b[c];
    }
  }
  return lu;
}

// Back-substitution on the reduced system (Mbuf = packed U from
// row_reduce, x the forward-eliminated RHS).
inline void back_substitute(std::span<const double> Mbuf, int n,
                            std::span<double> x) {
  const MatView<const double> M{Mbuf, n, n, n};
  for (int c = n - 1; c >= 0; c--) {
    double s = x[c];
    for (int j = c + 1; j < n; j++) s -= M(c, j) * x[j];
    x[c] = s / M(c, c);
  }
}

// Solve A x = b over caller scratch M.  On Singular, x is the ZERO vector
// -- the owner solve()'s failure value, so callers that ignore the status
// still compute with the reference result (a zero step: guaranteed
// reject, radius shrink).
// @pre A.m == A.n; M.size() >= n*n; b.size() >= n; x.size() >= n
inline LuStatus solve(MatConstView A, std::span<const double> b,
                      std::span<double> M, std::span<double> x) {
  const int n = A.n;
  for (int i = 0; i < n; i++) x[i] = b[i];
  const LuReduction lu = row_reduce(A, M, x.first(n));
  if (lu.status != LuStatus::Ok) {
    for (int i = 0; i < n; i++) x[i] = 0.0;
    return lu.status;
  }
  back_substitute(M, n, x);
  return LuStatus::Ok;
}

// Solve (A + lambda I) x = b: the shifted copy staged in M2, then solve()
// on it with LU scratch M.  Copy-then-add-diagonal, the owner's exact
// arithmetic: adding 0.0 to off-diagonals instead would flip the sign of
// negative zeros.
// @pre A.m == A.n; M2.size() >= n*n; M.size() >= n*n
inline LuStatus solve_shifted(MatConstView A, std::span<const double> b,
                              double lambda, std::span<double> M2,
                              std::span<double> M, std::span<double> x) {
  const int n = A.n;
  const MatView<double> S{M2, n, n, n};
  for (int i = 0; i < n; i++)
    for (int j = 0; j < n; j++) S(i, j) = A(i, j);
  for (int i = 0; i < n; i++) S(i, i) += lambda;
  return solve(MatConstView{S}, b, M, x);
}

}  // namespace LinAlg
