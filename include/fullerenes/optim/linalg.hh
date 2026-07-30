#pragma once

// =====================================================================
// optim::la -- the CPU LinAlg backend of the framework.
//
// Two layers:
//
//   1. Flat-span vector kernels (dot, inf_norm, axpy, scal, copy) --
//      written in EXACTLY the left-to-right accumulation order of
//      minimize::detail, because the wu re-expression gate requires
//      bit-identical trajectories to the parent minimize::lbfgs.
//      Deliberately NOT the parallel-primitives D1 fixed-shape monoid
//      reduction; that shape arrives with the step-6 lowering, gated by
//      its own parity run.
//
//   2. Dense small-n operations, re-exported from the parent library's
//      BLAS-free LinAlg (fullerenes/dense_linalg.hh): solve /
//      solve_shifted / matvec / energy / norm / max_abs / is_usable_step and
//      the SymEigen truncated pseudoinverse.  The dense matrix type is
//      the library's matrix<double>.  The BLAS-free rationale (silently
//      wrong deployed OpenBLAS) is that header's hard rule and applies
//      here unchanged.
//
// The sparse and KKT/saddle shapes are migration step 5 (deferred by
// sign-off); nothing here precludes them -- they will be additional
// overload sets behind the same names.
// =====================================================================

#include "fullerenes/dense_linalg.hh"
#include "fullerenes/matrix.hh"

#include <cmath>
#include <span>

namespace optim {

using Mat = matrix<double>;

namespace la {

// --- Flat-span kernels (lib-order floating point; see header note) ---

inline double dot(std::span<const double> a, std::span<const double> b) {
  double s = 0;
  for (std::size_t i = 0; i < a.size(); ++i) s += a[i] * b[i];
  return s;
}

inline double inf_norm(std::span<const double> a) {
  double m = 0;
  for (double v : a) m = std::max(m, std::fabs(v));
  return m;
}

inline double nrm2(std::span<const double> a) { return std::sqrt(dot(a, a)); }

inline void axpy(std::span<double> y, double a, std::span<const double> x) {
  for (std::size_t i = 0; i < y.size(); ++i) y[i] += a * x[i];
}

inline void scal(std::span<double> v, double a) {
  for (double& x : v) x *= a;
}

inline void copy(std::span<const double> src, std::span<double> dst) {
  for (std::size_t i = 0; i < src.size(); ++i) dst[i] = src[i];
}

// --- Dense small-n backend: the parent library's BLAS-free module ---

using ::LinAlg::V;
using ::LinAlg::solve;
using ::LinAlg::solve_with_sign;
using ::LinAlg::solve_shifted;
using ::LinAlg::matvec;
using ::LinAlg::energy;
using ::LinAlg::norm;
using ::LinAlg::max_abs;
using ::LinAlg::is_usable_step;
using ::LinAlg::jacobi_eig;
using ::LinAlg::sym_eigvals;
namespace SymEigen = ::LinAlg::SymEigen;

}  // namespace la
}  // namespace optim
