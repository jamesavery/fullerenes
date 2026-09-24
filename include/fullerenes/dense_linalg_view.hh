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
// with three stated exceptions: negate() is NEW at this level (the
// device-legal span form of auxiliary.hh's vector unary minus), matvec()
// generalized its loop bounds from the square A.m to the honest m x n
// (identical on every square input), and the transposed products
// matmul_nt / matmul_tn are new (pinned against the plain product on
// materialised transposes).  The independent pins are
// dense-linalg-test's frozen historical oracle (the pre-promotion bodies,
// byte-compared), its strided-vs-packed leg, and -- for the batch port's
// BINDING layer only -- the Alexandrov pilot's 5,771 + 100 isomer corpus.
//
// SCALAR PARAMETERIZATION (2026-08-20).  Every body is a template on its
// scalar T, so a solver can run in float as well as double (the port's
// fp32 tier).  The double instantiation is the historical body VERBATIM --
// same loop order, same accumulator type, same constants -- so nothing on
// the double path moved; the frozen-oracle byte gates in dense-linalg-test
// are what say so.  (One store was added since: row_reduce keeps each
// multiplier in the entry it eliminates, a cell no output reads.)
// Deduction convention, chosen so that EVERY existing double call site
// compiles unchanged: T is deduced where an argument carries it
// unambiguously (a MatView operand, or the scalar of a scaled update) and
// is otherwise a DEFAULTED parameter (T = double) with the
// span arguments in a non-deduced context -- a std::vector<double> then
// binds exactly as it always did, and a float caller names its scalar
// (max_abs<float>(...)).
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
// The cyclic Jacobi eigensolver's rotation family is at this level since
// 2026-09-24 (jacobi_negligible, jacobi_rotation, rotate_pair, rotate_columns,
// rotate_rows, off_diagonal_max, jacobi_sweep, jacobi_diagonalize): the owner
// jacobi_eig in dense_linalg.cc composes them in its historical order, and
// the double instantiation is that loop VERBATIM (dense-linalg-test's
// frozen jacobi_eig oracle, byte-compared).  SymEigen stays owner-level.

#include <algorithm>
#include <bit>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>
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

// A span argument whose scalar is fixed by ANOTHER argument (or by the
// default T = double), never deduced from itself -- the one spelling of
// the deduction convention stated in the banner.
template <class T> using in_  = std::type_identity_t<std::span<const T>>;
template <class T> using out_ = std::type_identity_t<std::span<T>>;

// --- Vector reductions (dense_linalg.cc's loops; spans are exact-sized:
//     the span IS the vector, no separate count) ---

// @pre a.size() == b.size()
template <class T = double>
inline T dot(in_<T> a, in_<T> b) {
  T s = 0;
  for (std::size_t i = 0; i < a.size(); i++) s += a[i] * b[i];
  return s;
}
template <class T = double>
inline T sum_sq(in_<T> v) { return dot<T>(v, v); }
template <class T = double>
inline T norm(in_<T> v) { return std::sqrt(sum_sq<T>(v)); }
// max |v_i|; a NaN entry poisons the result to +inf, so a NaN residual can
// never pass a "< tol" convergence test.
template <class T = double>
inline T max_abs(in_<T> v) {
  T m = 0;
  for (T x : v) {
    if (std::isnan(x)) return std::numeric_limits<T>::infinity();
    m = std::max(m, std::fabs(x));
  }
  return m;
}
// v := -v elementwise (exact sign flip; the span form of the vector unary
// minus the owner expressions use).  New at this level -- see the banner.
template <class T>
inline void negate(std::span<T> v) {
  for (T& x : v) x = -x;
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
template <class T = double>
inline void copy_into(out_<T> dst, in_<T> src) {
  for (std::size_t i = 0; i < dst.size(); i++) dst[i] = src[i];
}
// dst := -src (exact sign flip, like negate)
template <class T = double>
inline void neg_into(out_<T> dst, in_<T> src) {
  for (std::size_t i = 0; i < dst.size(); i++) dst[i] = -src[i];
}
// dst += src
template <class T = double>
inline void add_into(out_<T> dst, in_<T> src) {
  for (std::size_t i = 0; i < dst.size(); i++) dst[i] += src[i];
}
// dst -= src
template <class T = double>
inline void sub_into(out_<T> dst, in_<T> src) {
  for (std::size_t i = 0; i < dst.size(); i++) dst[i] -= src[i];
}
// dst += s * src   (T deduced from the scale s)
template <class T>
inline void add_scaled(out_<T> dst, T s, in_<T> src) {
  for (std::size_t i = 0; i < dst.size(); i++) dst[i] += s * src[i];
}
// dst -= s * src   (T deduced from the scale s)
template <class T>
inline void sub_scaled(out_<T> dst, T s, in_<T> src) {
  for (std::size_t i = 0; i < dst.size(); i++) dst[i] -= s * src[i];
}
// v *= s   (T deduced from the scale s; the span form of the owner's
// feasibility clip v = s * v)
template <class T>
inline void scale(std::span<T> v, T s) {
  for (T& x : v) x = s * x;
}
// Sum of entries -- term-for-term identical to the owner idiom
// dot(v, ones) (v[i] * 1.0 == v[i] exactly, same accumulation order).
template <class T = double>
inline T sum(in_<T> v) {
  T s = 0;
  for (T x : v) s += x;
  return s;
}

// --- Solver policy (calibrated constants; NOT neutral linear algebra) ---

// Residual energy E = 1/2 ||v||^2 (the Gauss-Newton objective's 1/2).
template <class T = double>
inline T energy(in_<T> v) { return T(0.5) * sum_sq<T>(v); }

// Step-acceptance floor: ||v||^2 > 1e-30, i.e. ||v|| > 1e-15.  On the
// SQUARED norm deliberately -- a clip-capped trust-radius step with
// 0 < ||v|| <= 1e-15 must still count as a step (a max|v| > 0 test takes
// the opposite branch there, and the Gauss-Newton bisection routes ~31
// solves per fallback through this predicate).
// PRECISION NOTE: the value is calibrated on doubles and is used unscaled
// in every T (1e-30 is a normal float, so the float comparison is
// well-defined but far below float's own ~1e-7 relative resolution -- a
// float step of norm 1e-15 is numerical noise that this floor still
// accepts).  Whether the fp32 tier wants its own floor is an open pin.
inline constexpr double STEP_SQ_FLOOR = 1e-30;
template <class T = double>
inline bool is_usable_step(in_<T> v) {
  const T s = sum_sq<T>(v);
  return std::isfinite(s) && s > T(STEP_SQ_FLOOR);
}

// --- Matrix products ---

// out := A * v.
// @pre v.size() >= A.n, out.size() >= A.m
template <class T = double>
inline void matvec(MatView<const T> A, in_<T> v, out_<T> out) {
  for (int i = 0; i < A.m; i++) {
    T s = 0;
    for (int j = 0; j < A.n; j++) s += A(i, j) * v[j];
    out[i] = s;
  }
}

// op(A) for the product body below: entry (i, k) of A, or of A^T when Trans
// is set (BLAS's op argument, resolved at compile time).
template <bool Trans, class T>
inline T op_entry(MatView<const T> A, int i, int k) {
  if constexpr (Trans) return A(k, i);
  else                 return A(i, k);
}

// @anchor matmul-ijk-order
// out := op(A) * op(B), packed row-major, op = transpose where the flag is
// set.  ONE body for the products: the i-j-k loop with a k-ascending scalar
// accumulator, a transposition being an index swap in the operand read and
// nothing else.  LOAD-BEARING: the Alexandrov solver's byte gates are
// validated against this association, matrix<double>::operator* delegates
// to the plain instance, and dense-linalg-test's frozen product oracle pins
// it.  Do not reassociate, block, or vectorize this loop.
// The two inputs may be the same matrix (R^T R); out is disjoint from both.
// @pre op(A) is m x kk and op(B) is kk x n; out.size() >= m * n
template <bool TransA, bool TransB, class T>
inline void matmul_op(MatView<const T> A, std::type_identity_t<MatView<const T>> B,
                      out_<T> out) {
  const int m  = TransA ? A.n : A.m;
  const int kk = TransA ? A.m : A.n;
  const int n  = TransB ? B.m : B.n;
  const MatView<T> C{out, m, n, n};
  for (int i = 0; i < m; i++)
    for (int j = 0; j < n; j++) {
      T x = 0;
      for (int k = 0; k < kk; k++)
        x += op_entry<TransA>(A, i, k) * op_entry<TransB>(B, k, j);
      C(i, j) = x;
    }
}

// out := A * B  (the historical matmul).  @pre A.n == B.m, out.size() >= A.m * B.n
template <class T = double>
inline void matmul(MatView<const T> A, std::type_identity_t<MatView<const T>> B,
                   out_<T> out) {
  matmul_op<false, false, T>(A, B, out);
}
// out := A * B^T.  @pre A.n == B.n, out.size() >= A.m * B.m
template <class T = double>
inline void matmul_nt(MatView<const T> A, std::type_identity_t<MatView<const T>> B,
                      out_<T> out) {
  matmul_op<false, true, T>(A, B, out);
}
// out := A^T * B.  @pre A.m == B.m, out.size() >= A.n * B.n
template <class T = double>
inline void matmul_tn(MatView<const T> A, std::type_identity_t<MatView<const T>> B,
                      out_<T> out) {
  matmul_op<true, false, T>(A, B, out);
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

// --- The elimination, step by step.  row_reduce is the sequential
//     composition of the three words below; a lane-parallel lowering (each
//     row of a step in its own lane -- parallel-primitives' par::linsolve)
//     composes the SAME words, so its per-element arithmetic -- the
//     operands, the roundings, the order in which an element receives its
//     updates -- is that of the sequential body, and its results are
//     byte-equal to it by construction. ---

// pivot_row: the partial-pivot choice for column c -- the first row q >= c
// whose |M(q, c)| is largest.  A NaN entry is never chosen (no comparison
// with it succeeds); when M(c, c) itself is NaN nothing displaces row c.
// The historical scan, verbatim.
template <class T>
inline int pivot_row(MatView<const T> M, int c) {
  int p = c;
  for (int q = c + 1; q < M.n; q++)
    if (std::fabs(M(q, c)) > std::fabs(M(p, c))) p = q;
  return p;
}

// pivot_key: the SAME choice as a total order on the rows q >= c, so that a
// collective evaluation -- the largest key, then the smallest row carrying
// it -- selects pivot_row's row whatever order the rows are visited in.
// The key of a row is the bit pattern of |M(q, c)|: for non-negative finite
// values and +inf the unsigned pattern is ordered as the value is, and -0.0
// maps to +0.0.  NaN is stated explicitly: at q > c the minimum (never
// chosen; a zero entry ties it, and the smaller row -- never a NaN row
// while row c has a smaller index -- wins), at q == c the maximum (row c
// stays, as in the scan).
// @post argmax over q >= c of (pivot_key(M(q, c), q == c), then the smaller q)
//           == pivot_row(M, c)   (dense-linalg-test pins this on the NaN, zero,
//           infinity and tie columns no random matrix produces)
template <class T>
inline std::uint64_t pivot_key(T entry, bool is_row_c) {
  if (std::isnan(entry)) return is_row_c ? ~std::uint64_t{0} : std::uint64_t{0};
  const T a = std::fabs(entry);
  if constexpr (sizeof(T) == 4) return std::bit_cast<std::uint32_t>(a);
  else                          return std::bit_cast<std::uint64_t>(a);
}

// swap_entry: column j of a row swap -- the element statement swap_rows
// loops, and the one a column-parallel swap runs one per lane.
template <class T>
inline void swap_entry(MatView<T> M, int c, int p, int j) { std::swap(M(c, j), M(p, j)); }

// swap_rows: rows c and p from column c on, and the right-hand side entries
// (the columns before c are dead after step c: U lives on and above the
// diagonal).
template <class T>
inline void swap_rows(MatView<T> M, int c, int p, out_<T> b) {
  for (int j = c; j < M.n; j++) swap_entry(M, c, p, j);
  if (!b.empty()) std::swap(b[c], b[p]);
}

// The two arithmetic statements of an elimination step, on VALUES, stated
// once so that any loop order composes them -- by row (eliminate_row below,
// the sequential body) or by column (a lane per column on a GPU, which
// supplies the pivot-row value from a register): the multiplier of a row,
// and the update of one entry -- of the matrix, or of the right-hand side.
// A row whose multiplier is exactly zero is skipped by every composition
// (the update is not the identity when the pivot row holds an infinity).
template <class T>
inline T eliminate_multiplier(T entry_qc, T pivot) { return entry_qc / pivot; }
template <class T>
inline void eliminate_entry(T& entry_qj, T mult, T entry_cj) { entry_qj -= mult * entry_cj; }

// eliminate_row: row q of step c -- the multiplier mult = M(q, c) / M(c, c)
// is STORED in the entry it eliminates (the strict lower triangle thus holds
// the multipliers, as in-place LU does), then mult times the pivot row is
// subtracted from row q on the columns after c and from b.  Reads rows c and
// q, writes row q and b[q] only, so the rows of one step are independent of
// each other.
template <class T>
inline void eliminate_row(MatView<T> M, int c, int q, out_<T> b) {
  const T mult = eliminate_multiplier(M(q, c), M(c, c));
  M(q, c) = mult;
  if (mult == 0) return;
  for (int j = c + 1; j < M.n; j++) eliminate_entry(M(q, j), mult, M(c, j));
  if (!b.empty()) eliminate_entry(b[q], mult, b[c]);
}

// Partial-pivot row reduction of a COPY of A staged in the packed scratch
// M (capacity n*n): M holds U on and above the diagonal and each step's
// multipliers below it (in the entries they eliminated); the returned sign
// is sign(det A).  When b is non-empty the same swaps and forward
// elimination are applied to it, leaving the triangular system U x = b
// ready for back_substitute.  The permutation is not recorded and a swap
// does not permute the multipliers of earlier columns, so M is NOT an LU
// factorisation of P A -- only U, sign det, and the eliminated RHS are
// outputs (hence not "lu_decompose": a second RHS cannot be solved from
// the result).  On Singular, M, sign, and b are only partially reduced.
// @pre A.m == A.n (square); M.size() >= n*n; b empty or b.size() >= n
template <class T = double>
inline LuReduction row_reduce(MatView<const T> A, out_<T> Mbuf, out_<T> b) {
  const int n = A.n;
  const MatView<T> M{Mbuf, n, n, n};
  for (int i = 0; i < n; i++)
    for (int j = 0; j < n; j++) M(i, j) = A(i, j);
  LuReduction lu;
  for (int c = 0; c < n; c++) {
    const int p = pivot_row(MatView<const T>{M}, c);
    if (M(p, c) == 0) { lu.status = LuStatus::Singular; return lu; }
    if (p != c) { swap_rows(M, c, p, b); lu.sign = -lu.sign; }   // row swap parity
    if (M(c, c) < 0) lu.sign = -lu.sign;                          // sign of diag(U)
    for (int q = c + 1; q < n; q++) eliminate_row(M, c, q, b);
  }
  return lu;
}

// Back-substitution on the reduced system (Mbuf = packed U from
// row_reduce, x the forward-eliminated RHS).  Each x[c] is a fold in
// ascending j whose first term needs x[c+1], the value finished last, so
// the chain of subtractions is sequential by data flow: no lowering can
// spread it over lanes without changing the association.
template <class T = double>
inline void back_substitute(in_<T> Mbuf, int n, out_<T> x) {
  const MatView<const T> M{Mbuf, n, n, n};
  for (int c = n - 1; c >= 0; c--) {
    T s = x[c];
    for (int j = c + 1; j < n; j++) s -= M(c, j) * x[j];
    x[c] = s / M(c, c);
  }
}

// Solve A x = b over caller scratch M.  On Singular, x is the ZERO vector
// -- the owner solve()'s failure value, so callers that ignore the status
// still compute with the reference result (a zero step: guaranteed
// reject, radius shrink).
// @pre A.m == A.n; M.size() >= n*n; b.size() >= n; x.size() >= n
template <class T = double>
inline LuStatus solve(MatView<const T> A, in_<T> b, out_<T> M, out_<T> x) {
  const int n = A.n;
  for (int i = 0; i < n; i++) x[i] = b[i];
  const LuReduction lu = row_reduce<T>(A, M, x.first(n));
  if (lu.status != LuStatus::Ok) {
    for (int i = 0; i < n; i++) x[i] = T(0);
    return lu.status;
  }
  back_substitute<T>(M, n, x);
  return LuStatus::Ok;
}

// Solve (A + lambda I) x = b: the shifted copy staged in M2, then solve()
// on it with LU scratch M.  Copy-then-add-diagonal, the owner's exact
// arithmetic: adding 0.0 to off-diagonals instead would flip the sign of
// negative zeros.
// @pre A.m == A.n; M2.size() >= n*n; M.size() >= n*n
template <class T = double>
inline LuStatus solve_shifted(MatView<const T> A, in_<T> b, T lambda,
                              out_<T> M2, out_<T> M, out_<T> x) {
  const int n = A.n;
  const MatView<T> S{M2, n, n, n};
  for (int i = 0; i < n; i++)
    for (int j = 0; j < n; j++) S(i, j) = A(i, j);
  for (int i = 0; i < n; i++) S(i, i) += lambda;
  return solve<T>(MatView<const T>{S}, b, M, x);
}

// --- Symmetric eigendecomposition: the cyclic Jacobi family, step by step.
//     jacobi_eig (dense_linalg.cc) is the sequential composition of these
//     words in its historical order.  A lane-parallel lowering that rotates
//     DISJOINT pairs concurrently (the tournament-scheduled sweep of a
//     batched eigensolver) composes the SAME value statements, so its
//     per-element arithmetic is that of the sequential body; but rotations
//     on disjoint pairs, which commute as matrices, meet in the entries at
//     the intersection of one rotation's rows and the other's columns, and
//     the two orders of applying them round differently -- hence agreement
//     to a tolerance, not to the bit.  Conventions: the accumulator V holds
//     eigenvector m as ROW m; the rotation J of a pair (p, q) acts as
//     A <- J^T A J and V <- J^T V. ---

// jacobi_negligible: |x| <= tol, the one predicate of the family -- it
// decides which pairs a sweep leaves alone and when the off-diagonal is
// done.  @pre tol >= 0 (so an exact zero is always negligible)
template <class T>
inline bool jacobi_negligible(T x, T tol) { return std::fabs(x) <= tol; }

template <class T>
struct Rotation { T c, s; };

// jacobi_rotation: the stable plane rotation annihilating a_pq, as (c, s):
// theta = (a_qq - a_pp) / (2 a_pq), t = sign(theta) / (|theta| + sqrt(theta^2 + 1)),
// c = 1 / sqrt(t^2 + 1), s = t c.  @pre a_pq != 0
template <class T>
inline Rotation<T> jacobi_rotation(T a_pp, T a_qq, T a_pq) {
  const T theta = (a_qq - a_pp) / (2 * a_pq);
  const T t = (theta >= 0 ? T(1) : T(-1)) /
              (std::fabs(theta) + std::sqrt(theta * theta + 1));
  const T c = T(1) / std::sqrt(t * t + 1);
  return {c, t * c};
}

// rotate_pair: the one arithmetic statement of a plane rotation, on the
// coordinate pair (x, y) at positions (p, q):  x' = c x - s y,  y' = s x + c y.
template <class T>
inline void rotate_pair(T& x, T& y, T c, T s) {
  const T x0 = x, y0 = y;
  x = c * x0 - s * y0;
  y = s * x0 + c * y0;
}

// rotate_columns: A <- A J, the rotation applied to columns p and q of every
// row.  rotate_rows: A <- J^T A, applied to rows p and q over every column;
// the eigenvector accumulator's update is this word on V.
// @pre 0 <= p, q < the rotated dimension and p != q (p == q would write one
//      entry twice)
template <class T>
inline void rotate_columns(MatView<T> A, int p, int q, T c, T s) {
  for (int i = 0; i < A.m; i++) rotate_pair(A(i, p), A(i, q), c, s);
}
template <class T>
inline void rotate_rows(MatView<T> A, int p, int q, T c, T s) {
  for (int i = 0; i < A.n; i++) rotate_pair(A(p, i), A(q, i), c, s);
}

// off_diagonal_max: max |a_pq| over the strict upper triangle, the
// historical convergence measure.  Unlike max_abs, a NaN entry is IGNORED
// here (std::max keeps the running value), so a NaN matrix reports
// convergence; byte-identity to the frozen jacobi_eig oracle fixes this.
// @pre A.m == A.n, A symmetric
template <class T>
inline T off_diagonal_max(MatView<const T> A) {
  T off = 0;
  for (int p = 0; p < A.m; p++)
    for (int q = p + 1; q < A.n; q++) off = std::max(off, std::fabs(A(p, q)));
  return off;
}

// jacobi_sweep: one cyclic sweep, every pair (p, q) with p < q in row-major
// order.  A pair whose entry is negligible is skipped and WRITTEN NOT AT
// ALL -- the historical continue; applying the identity rotation instead
// would not be byte-safe, since 0 * x + y turns a negative zero positive.
// (With tol >= 0 an exact zero a_pq always skips, so jacobi_rotation never
// divides by zero.)  The rest are rotated: columns, then rows, then the
// accumulator when V is non-empty (V.m == 0 means no vectors are wanted).
// @pre A.m == A.n; V.m == 0 or (V.m == A.n and V.n == A.n); tol >= 0
template <class T>
inline void jacobi_sweep(MatView<T> A, MatView<T> V, T tol) {
  const int n = A.n;
  for (int p = 0; p < n; p++)
    for (int q = p + 1; q < n; q++) {
      const T a_pq = A(p, q);
      if (jacobi_negligible(a_pq, tol)) continue;
      const Rotation<T> r = jacobi_rotation(A(p, p), A(q, q), a_pq);
      rotate_columns(A, p, q, r.c, r.s);
      rotate_rows(A, p, q, r.c, r.s);
      if (V.m) rotate_rows(V, p, q, r.c, r.s);
    }
}

// jacobi_diagonalize: cyclic Jacobi, the historical sweep loop.  Before
// each sweep the off-diagonal maximum is tested; converged when it is
// negligible, with the diagonal of A the eigenvalues (UNSORTED: entry m
// stays paired with row m of V) and V the accumulated rotations applied to
// whatever the caller put in it -- the identity for a fresh
// eigendecomposition, a previous eigenbasis for a warm start.  tol is the
// caller's negligibility scale (the owner jacobi_eig uses
// 1e-15 * max(max|A|, 1e-300)).  The guard is the historical one: the
// sweep that would be the max_sweeps-th is refused (at most max_sweeps - 1
// sweeps run), and a trip means a non-symmetric or non-finite input.
// sweeps counts the sweeps performed.  Two conventions to know: a NaN
// off-diagonal is not seen by off_diagonal_max, so a NaN matrix reports
// convergence with NaN eigenvalues (the historical behaviour, unlike
// max_abs); and max_sweeps <= 0 reports not converged, where the
// historical loop fell through as converged (unreachable from the owner,
// whose guard is 60).
// @pre as jacobi_sweep
struct JacobiResult { bool converged = false; int sweeps = 0; };
template <class T>
inline JacobiResult jacobi_diagonalize(MatView<T> A, MatView<T> V, T tol, int max_sweeps) {
  for (int sweep = 0; sweep < max_sweeps; sweep++) {
    if (jacobi_negligible(off_diagonal_max(MatView<const T>{A}), tol)) return {true, sweep};
    if (sweep == max_sweeps - 1) return {false, sweep};
    jacobi_sweep(A, V, tol);
  }
  return {false, 0};   // max_sweeps <= 0: nothing was tested
}

}  // namespace LinAlg
