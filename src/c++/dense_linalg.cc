#include "fullerenes/dense_linalg.hh"

#include <algorithm>
#include <cmath>
#include <numeric>

using namespace std;

namespace LinAlg {

// (Vector reductions: view-level only -- vectors convert to spans, so the
// view functions serve every owner caller with no wrapper.)

// --- Symmetric eigendecomposition (cyclic Jacobi) ---

bool jacobi_eig(std::vector<double> A, int n, std::vector<double>& lam,
                std::vector<double>* V_out)
{
  constexpr int MAX_SWEEPS = 60;

  // Eigenvector accumulator (row m = eigvec m); only when requested.
  // Named Vacc, NOT V: LinAlg::V is the vector type alias in this scope.
  std::vector<double> Vacc;
  if (V_out) {
    Vacc.assign(size_t(n) * n, 0.0);
    for (int i = 0; i < n; i++) Vacc[size_t(i) * n + i] = 1.0;
  }

  // The negligibility scale: 1e-15 of the largest entry, floored at 1e-300.
  // A hand loop rather than max_abs, whose NaN policy (+inf) would declare
  // instant convergence on NaN input and change this body's bytes there.
  double anorm = 0;
  for (size_t x = 0; x < A.size(); x++) anorm = std::max(anorm, std::fabs(A[x]));
  const double tol = 1e-15 * std::max(anorm, 1e-300);

  // The sweep loop is the view-level jacobi_diagonalize
  // (dense_linalg_view.hh), composed in the historical order -- columns,
  // rows, accumulator -- so this body is byte-identical to its
  // pre-promotion self (the frozen oracle in dense-linalg-test).  Vv is
  // EMPTY (m = 0) when no vectors are requested; the accumulator is then
  // left alone.
  const MatView<double> Av{A, n, n, n};
  const MatView<double> Vv{Vacc, V_out ? n : 0, n, n};
  if (!jacobi_diagonalize(Av, Vv, tol, MAX_SWEEPS).converged)
    return false;   // guard trip = bug (non-symmetric input)

  // Sort ascending, permuting eigenvector rows along.
  std::vector<int> order(n);
  std::iota(order.begin(), order.end(), 0);
  std::sort(order.begin(), order.end(),
            [&](int x, int y) { return Av(x, x) < Av(y, y); });
  lam.assign(n, 0);
  for (int m = 0; m < n; m++) lam[m] = Av(order[m], order[m]);
  if (V_out) {
    V_out->assign(size_t(n) * n, 0.0);
    for (int m = 0; m < n; m++)
      for (int i = 0; i < n; i++)
        (*V_out)[size_t(m) * n + i] = Vacc[size_t(order[m]) * n + i];
  }
  return true;
}

// Symmetric part of A as a flat row-major buffer: ½(A + Aᵀ).  The Jacobi
// eigensolver assumes a symmetric input, so both eigen entry points funnel
// their matrix through here first.
static std::vector<double> symmetrize_flat(const matrix<double>& A)
{
  int n = A.m;
  std::vector<double> Af(size_t(n)*n);
  for (int i = 0; i < n; i++)
    for (int j = 0; j < n; j++)
      Af[size_t(i)*n + j] = 0.5 * (A(i,j) + A(j,i));
  return Af;
}

std::vector<double> sym_eigvals(const matrix<double>& A)
{
  int n = A.m;
  std::vector<double> w;
  if (!jacobi_eig(symmetrize_flat(A), n, w, nullptr)) return {};
  return w;
}

// --- LU with partial pivoting (bodies: dense_linalg_view.hh's shared
// row_reduce core; these wrappers allocate the scratch and the result) ---

std::expected<LuSolved, LuFail>
solve_with_sign(const matrix<double>& A, const V& b)
{
  const int n = A.m;
  vector<double> M(size_t(n) * n);
  V x(b);
  const LuReduction lu = row_reduce(view_of(A), M, x);
  if (lu.status != LuStatus::Ok) return std::unexpected(LuFail::Singular);
  back_substitute(M, n, x);
  return LuSolved{std::move(x), lu.sign};
}

double det(const matrix<double>& A)
{
  const int n = A.m;
  vector<double> M(size_t(n) * n);
  // An exact zero pivot means A is singular — det = 0 is the true value here,
  // not an error, so we return 0.0 rather than propagating LuFail::Singular.
  const LuReduction lu = row_reduce(view_of(A), M, {});
  if (lu.status != LuStatus::Ok) return 0.0;
  // |det A| = ∏|U_ii|; sign already carries (−1)^{#swaps} × ∏ sign(U_ii).
  double mag = 1.0;
  for (int i = 0; i < n; i++) mag *= fabs(M[size_t(i)*n + i]);
  return lu.sign * mag;
}

V solve(const matrix<double>& A, const V& b)
{
  const int n = A.m;
  vector<double> M(size_t(n) * n);
  V x(n, 0.0);
  // The view solve zero-fills x on a singular A — the documented value.
  solve(view_of(A), b, M, x);
  return x;
}

V solve_shifted(const matrix<double>& A, const V& b, double lambda)
{
  const int n = A.m;
  vector<double> M2(size_t(n) * n), M(size_t(n) * n);
  V x(n, 0.0);
  solve_shifted(view_of(A), b, lambda, M2, M, x);
  return x;
}

V matvec(const matrix<double>& A, const V& v)
{
  V r(A.m, 0.0);
  matvec(view_of(A), v, r);
  return r;
}


// --- Symmetric truncated pseudoinverse ---

namespace SymEigen {

Decomp decompose(const matrix<double>& A)
{
  int n = A.m;
  V lambda;
  std::vector<double> Vrows;                        // row m = eigvec m
  if (!jacobi_eig(symmetrize_flat(A), n, lambda, &Vrows))
    return Decomp{matrix<double>(0, 0, 0.0), {}};
  matrix<double> Q(n, n, 0.0);
  for (int i = 0; i < n; i++)
    for (int j = 0; j < n; j++)
      Q(j, i) = Vrows[size_t(i) * n + j];           // eigvec i is column i
  return {std::move(Q), std::move(lambda)};
}

Signature signature(const V& lambda, double tol)
{
  Signature s{0, 0, 0};
  for (double l : lambda) {
    if      (l >  tol) s.pos++;
    else if (l < -tol) s.neg++;
    else               s.zero++;
  }
  return s;
}

ApplyResult apply(const Decomp& d, const V& b, double cutoff)
{
  int n = (int)d.lambda.size();
  V   x(n, 0.0);
  int rank = 0;
  double lam_min = 0;
  for (int i = 0; i < n; i++) {
    double l = d.lambda[i];
    if (std::fabs(l) <= cutoff) continue;
    rank++;
    if (lam_min == 0 || std::fabs(l) < lam_min) lam_min = std::fabs(l);
    double bq = 0;
    for (int j = 0; j < n; j++) bq += d.Q(j, i) * b[j];
    double coeff = bq / l;
    for (int j = 0; j < n; j++) x[j] += coeff * d.Q(j, i);
  }
  return {std::move(x), rank, lam_min};
}

Solution solve(const matrix<double>& A, const V& b, double rcond)
{
  auto d = decompose(A);
  if (d.lambda.empty()) return {};
  double lam_max = 0;
  for (double l : d.lambda) lam_max = std::max(lam_max, std::fabs(l));
  double cutoff = rcond * lam_max;
  auto   r      = apply(d, b, cutoff);
  return {std::move(r.x), r.rank, r.lambda_min_kept, lam_max,
          signature(d.lambda, cutoff)};
}

}  // namespace SymEigen

}  // namespace LinAlg
