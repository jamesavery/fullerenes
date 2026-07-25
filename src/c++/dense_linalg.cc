#include "fullerenes/dense_linalg.hh"

#include <algorithm>
#include <cmath>
#include <numeric>

using namespace std;

namespace LinAlg {

// --- Vector reductions ---

double dot(const V& a, const V& b) { double s=0; for (size_t i=0;i<a.size();i++) s+=a[i]*b[i]; return s; }
double norm(const V& v)            { return sqrt(dot(v, v)); }
double max_abs(const V& v)         { double m=0; for (double x:v) { if (std::isnan(x)) return HUGE_VAL; m=max(m,fabs(x)); } return m; }
double sum_sq(const V& v)          { return dot(v, v); }
bool   is_valid(const V& v)        { double n=sum_sq(v); return isfinite(n) && n > 1e-30; }
double energy(const V& v)          { return 0.5 * sum_sq(v); }

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

  double anorm = 0;
  for (size_t x = 0; x < A.size(); x++) anorm = std::max(anorm, std::fabs(A[x]));
  const double tol = 1e-15 * std::max(anorm, 1e-300);

  auto at = [&](int i, int j) -> double& { return A[size_t(i) * n + j]; };

  for (int sweep = 0; sweep < MAX_SWEEPS; sweep++) {
    double off = 0;
    for (int p = 0; p < n; p++)
      for (int q = p + 1; q < n; q++) off = std::max(off, std::fabs(at(p, q)));
    if (off <= tol) break;
    if (sweep == MAX_SWEEPS - 1) return false;   // guard trip = bug

    for (int p = 0; p < n; p++)
      for (int q = p + 1; q < n; q++) {
        double apq = at(p, q);
        if (std::fabs(apq) <= tol) continue;
        double theta = (at(q, q) - at(p, p)) / (2 * apq);
        double t = (theta >= 0 ? 1.0 : -1.0) /
                   (std::fabs(theta) + std::sqrt(theta * theta + 1));
        double c = 1.0 / std::sqrt(t * t + 1), s = t * c;

        for (int i = 0; i < n; i++) {      // rotate rows/cols p,q of A
          double aip = at(i, p), aiq = at(i, q);
          at(i, p) = c * aip - s * aiq;
          at(i, q) = s * aip + c * aiq;
        }
        for (int i = 0; i < n; i++) {
          double api = at(p, i), aqi = at(q, i);
          at(p, i) = c * api - s * aqi;
          at(q, i) = s * api + c * aqi;
        }
        if (V_out)
          for (int i = 0; i < n; i++) {    // accumulate: Vacc row = eigvecᵀ
            double vpi = Vacc[size_t(p) * n + i], vqi = Vacc[size_t(q) * n + i];
            Vacc[size_t(p) * n + i] = c * vpi - s * vqi;
            Vacc[size_t(q) * n + i] = s * vpi + c * vqi;
          }
      }
  }

  // Sort ascending, permuting eigenvector rows along.
  std::vector<int> order(n);
  std::iota(order.begin(), order.end(), 0);
  std::sort(order.begin(), order.end(),
            [&](int x, int y) { return at(x, x) < at(y, y); });
  lam.assign(n, 0);
  for (int m = 0; m < n; m++) lam[m] = at(order[m], order[m]);
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

// --- LU with partial pivoting ---

// Partial-pivot LU factorization, shared by solve_with_sign and det.  Copies A
// into the row-major working buffer M (n×n) and reduces it so M holds U on and
// above the diagonal; sign := (−1)^{#row swaps} × ∏ sign(diag U).  When b is
// non-null the same swaps and forward elimination are applied to it, leaving the
// triangular system U·x = b ready for back-substitution.  Returns false on an
// exact zero pivot (A singular); M and sign are then only partially reduced.
static bool lu_decompose(const matrix<double>& A, vector<double>& M,
                         int& sign, V* b)
{
  int n = A.m;
  M.assign(size_t(n)*n, 0.0);                  // row-major working copy
  for (int i = 0; i < n; i++)
    for (int j = 0; j < n; j++)
      M[size_t(i)*n + j] = A(i, j);
  sign = 1;
  for (int c = 0; c < n; c++) {
    int p = c;
    for (int q = c+1; q < n; q++)
      if (fabs(M[size_t(q)*n + c]) > fabs(M[size_t(p)*n + c])) p = q;
    if (M[size_t(p)*n + c] == 0) return false;
    if (p != c) {
      for (int j = c; j < n; j++) swap(M[size_t(c)*n + j], M[size_t(p)*n + j]);
      if (b) swap((*b)[c], (*b)[p]);
      sign = -sign;                            // row swap parity
    }
    if (M[size_t(c)*n + c] < 0) sign = -sign;  // sign of diag(U)
    for (int q = c+1; q < n; q++) {
      double m = M[size_t(q)*n + c] / M[size_t(c)*n + c];
      if (m == 0) continue;
      for (int j = c+1; j < n; j++) M[size_t(q)*n + j] -= m * M[size_t(c)*n + j];
      if (b) (*b)[q] -= m * (*b)[c];
    }
  }
  return true;
}

std::expected<LuSolved, LuFail>
solve_with_sign(const matrix<double>& A, const V& b)
{
  int n = A.m;
  vector<double> M;
  V   x(b);
  int sign;
  if (!lu_decompose(A, M, sign, &x)) return std::unexpected(LuFail::Singular);
  for (int c = n-1; c >= 0; c--) {             // back substitution
    double s = x[c];
    for (int j = c+1; j < n; j++) s -= M[size_t(c)*n + j] * x[j];
    x[c] = s / M[size_t(c)*n + c];
  }
  return LuSolved{std::move(x), sign};
}

double det(const matrix<double>& A)
{
  int n = A.m;
  vector<double> M;
  int sign;
  // An exact zero pivot means A is singular — det = 0 is the true value here,
  // not an error, so we return 0.0 rather than propagating LuFail::Singular.
  if (!lu_decompose(A, M, sign, nullptr)) return 0.0;
  // |det A| = ∏|U_ii|; sign already carries (−1)^{#swaps} × ∏ sign(U_ii).
  double mag = 1.0;
  for (int i = 0; i < n; i++) mag *= fabs(M[size_t(i)*n + i]);
  return sign * mag;
}

V solve(const matrix<double>& A, const V& b)
{
  auto r = solve_with_sign(A, b);
  return r ? std::move(r->x) : V(A.m, 0.0);
}

V solve_shifted(const matrix<double>& A, const V& b, double lambda)
{
  matrix<double> Al = A;
  for (int i = 0; i < A.m; i++) Al(i,i) += lambda;
  return solve(Al, b);
}

V matvec(const matrix<double>& A, const V& v)
{
  int n = A.m;
  V r(n, 0.0);
  for (int i = 0; i < n; i++)
    for (int j = 0; j < n; j++)
      r[i] += A(i,j) * v[j];
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
