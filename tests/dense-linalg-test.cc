// dense_linalg bit-exactness tests.  This module exists because a deployed
// OpenBLAS silently returned wrong dgesv solutions from n ~ 60 (info = 0),
// so the in-house LU and Jacobi eigensolver carry their own known-solution
// tests across the sizes where the corruption was observed.

#include "fullerenes/dense_linalg.hh"

#include <gtest/gtest.h>
#include <algorithm>
#include <cmath>
#include <cstring>
#include <limits>
#include <numeric>
#include <random>
#include <vector>

using namespace std;
using LinAlg::V;

// ============================================================================
// frozen: the pre-promotion owner bodies VERBATIM (dense_linalg.cc and
// matrix.hh operator* @ 262fb63e) -- the independent historical oracle for
// the view bodies.  FROZEN 2026-07-27: never refactor, never modernize; its
// entire value is being a second, unchanging implementation.  A frozen BODY
// compiled in this TU (rather than stored hex) cancels FP-contraction
// differences across compilers and architectures.
// FROZEN 2026-09-24 in addition: jacobi_eig @ e10921f3, the pre-promotion
// cyclic Jacobi body verbatim -- the oracle for the view-level rotation
// family (dense_linalg_view.hh) that the owner jacobi_eig now composes.
// (What the in-TU body cancels is cross-compiler contraction; the library
// body compiles in its own translation unit, so its agreement with this
// one also rests on the two sharing the one global CMAKE_CXX_FLAGS.)
// ============================================================================
namespace frozen {

static bool lu_decompose(const matrix<double>& A, vector<double>& M,
                         int& sign, V* b)
{
  int n = A.m;
  M.assign(size_t(n)*n, 0.0);
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
      sign = -sign;
    }
    if (M[size_t(c)*n + c] < 0) sign = -sign;
    for (int q = c+1; q < n; q++) {
      double m = M[size_t(q)*n + c] / M[size_t(c)*n + c];
      if (m == 0) continue;
      for (int j = c+1; j < n; j++) M[size_t(q)*n + j] -= m * M[size_t(c)*n + j];
      if (b) (*b)[q] -= m * (*b)[c];
    }
  }
  return true;
}

static std::expected<LinAlg::LuSolved, LinAlg::LuFail>
solve_with_sign(const matrix<double>& A, const V& b)
{
  int n = A.m;
  vector<double> M;
  V   x(b);
  int sign;
  if (!lu_decompose(A, M, sign, &x)) return std::unexpected(LinAlg::LuFail::Singular);
  for (int c = n-1; c >= 0; c--) {
    double s = x[c];
    for (int j = c+1; j < n; j++) s -= M[size_t(c)*n + j] * x[j];
    x[c] = s / M[size_t(c)*n + c];
  }
  return LinAlg::LuSolved{std::move(x), sign};
}

static double det(const matrix<double>& A)
{
  int n = A.m;
  vector<double> M;
  int sign;
  if (!lu_decompose(A, M, sign, nullptr)) return 0.0;
  double mag = 1.0;
  for (int i = 0; i < n; i++) mag *= fabs(M[size_t(i)*n + i]);
  return sign * mag;
}

static V solve(const matrix<double>& A, const V& b)
{
  auto r = solve_with_sign(A, b);
  return r ? std::move(r->x) : V(A.m, 0.0);
}

static V solve_shifted(const matrix<double>& A, const V& b, double lambda)
{
  matrix<double> Al = A;
  for (int i = 0; i < A.m; i++) Al(i,i) += lambda;
  return solve(Al, b);
}

static V matvec(const matrix<double>& A, const V& v)
{
  int n = A.m;
  V r(n, 0.0);
  for (int i = 0; i < n; i++)
    for (int j = 0; j < n; j++)
      r[i] += A(i,j) * v[j];
  return r;
}

static matrix<double> matmul(const matrix<double>& A, const matrix<double>& B)
{
  // the historical matrix.hh operator* loop, verbatim
  matrix<double> C(A.m, B.n);
  for (int i = 0; i < A.m; i++)
    for (int j = 0; j < B.n; j++) {
      double x = 0;
      for (int k = 0; k < A.n; k++) x += A[size_t(i)*A.n + k] * B[size_t(k)*B.n + j];
      C(i, j) = x;
    }
  return C;
}

// The pre-promotion jacobi_eig, verbatim (dense_linalg.cc @ e10921f3).
static bool jacobi_eig(std::vector<double> A, int n, std::vector<double>& lam,
                       std::vector<double>* V_out)
{
  constexpr int MAX_SWEEPS = 60;

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

}  // namespace frozen

namespace {

// Random m x n entries in [-1, 1), row-major draw order.  (Its own name:
// an overload of random_matrix would be ambiguous against the (n, seed,
// symmetric) form for an int seed and a bool literal.)
matrix<double> random_rectangular(int m, int n, unsigned seed)
{
  mt19937 rng(seed);
  uniform_real_distribution<double> u(-1.0, 1.0);
  matrix<double> A(m, n, 0.0);
  for (int i = 0; i < m; i++)
    for (int j = 0; j < n; j++) A(i, j) = u(rng);
  return A;
}

// Square; symmetric draws the upper triangle only and mirrors it (a
// different draw sequence from the full square, kept as the historical
// inputs of the tests below).
matrix<double> random_matrix(int n, unsigned seed, bool symmetric)
{
  if (!symmetric) return random_rectangular(n, n, seed);
  mt19937 rng(seed);
  uniform_real_distribution<double> u(-1.0, 1.0);
  matrix<double> A(n, n, 0.0);
  for (int i = 0; i < n; i++)
    for (int j = i; j < n; j++) {
      A(i, j) = u(rng);
      A(j, i) = A(i, j);
    }
  return A;
}

matrix<double> transpose_of(const matrix<double>& A)
{
  matrix<double> At(A.n, A.m, 0.0);
  for (int i = 0; i < A.m; i++) for (int j = 0; j < A.n; j++) At(j, i) = A(i, j);
  return At;
}

// The eigenpair residual of eigenvectors held as ROWS of V:
// max over m, j of | (row m of V) A - lambda_m (row m of V) |_j, in double.
template <class TV, class TA>
double row_residual(LinAlg::MatView<TV> V, const std::vector<double>& lambda,
                    LinAlg::MatView<TA> A)
{
  const int n = V.m;
  double res = 0;
  for (int m = 0; m < n; m++)
    for (int j = 0; j < n; j++) {
      double acc = 0;
      for (int i = 0; i < n; i++) acc += (double)V(m, i) * (double)A(i, j);
      res = std::max(res, std::fabs(acc - lambda[m] * (double)V(m, j)));
    }
  return res;
}
// The bar every double eigenpair residual in this file is held to.
constexpr double kEigenpairResidualBar = 1e-12;

V random_vector(int n, unsigned seed)
{
  mt19937 rng(seed);
  uniform_real_distribution<double> u(-1.0, 1.0);
  V x(n);
  for (double& v : x) v = u(rng);
  return x;
}

// A = Q·diag(D)·Qᵀ with Q orthogonal (a product of random Givens rotations).
// Then det(A) = det(Q)²·∏D = ∏D exactly, independent of whether Q is a rotation
// or reflection, and A is symmetric with eigenvalues D — a dense, well-
// conditioned matrix of known determinant when the |D_k| are bounded away from
// zero and each other.
matrix<double> qdqt(const V& D, unsigned seed)
{
  int n = (int)D.size();
  matrix<double> Q(n, n, 0.0);
  for (int i = 0; i < n; i++) Q(i, i) = 1.0;
  mt19937 rng(seed);
  uniform_real_distribution<double> ang(-M_PI, M_PI);
  uniform_int_distribution<int> pick(0, n - 1);
  for (int r = 0; r < 4 * n; r++) {                 // rotate columns p,q of Q
    int p = pick(rng), q = pick(rng);
    if (p == q) continue;
    double th = ang(rng), c = cos(th), s = sin(th);
    for (int k = 0; k < n; k++) {
      double qp = Q(k, p), qq = Q(k, q);
      Q(k, p) = c * qp - s * qq;
      Q(k, q) = s * qp + c * qq;
    }
  }
  matrix<double> A(n, n, 0.0);
  for (int i = 0; i < n; i++)
    for (int j = 0; j < n; j++) {
      double s = 0;
      for (int k = 0; k < n; k++) s += Q(i, k) * D[k] * Q(j, k);
      A(i, j) = s;
    }
  return A;
}

}  // namespace

TEST(DenseLinalg, LuSolvesKnownSolution)
{
  for (int n : {12, 30, 60, 90, 128, 200}) {
    matrix<double> A = random_matrix(n, 42 + n, /*symmetric=*/false);
    V x_true = random_vector(n, 7 + n);
    V b = LinAlg::matvec(A, x_true);
    V x = LinAlg::solve(A, b);
    double err = 0, scale = 0;
    for (int i = 0; i < n; i++) {
      err = max(err, fabs(x[i] - x_true[i]));
      scale = max(scale, fabs(x_true[i]));
    }
    EXPECT_LT(err / scale, 1e-9) << "n = " << n;
  }
}

TEST(DenseLinalg, LuDetSignMatchesEigenvalueSigns)
{
  // For a symmetric matrix, sign(det) = (-1)^{#negative eigenvalues}.
  for (int n : {12, 60}) {
    matrix<double> A = random_matrix(n, 1000 + n, /*symmetric=*/true);
    auto lu = LinAlg::solve_with_sign(A, random_vector(n, 5));
    ASSERT_TRUE(lu.has_value());
    auto lam = LinAlg::sym_eigvals(A);
    ASSERT_EQ((int)lam.size(), n);
    int n_neg = 0;
    for (double l : lam) n_neg += (l < 0);
    EXPECT_EQ(lu->det_sign, (n_neg % 2 == 0) ? 1 : -1) << "n = " << n;
  }
}

TEST(DenseLinalg, DeterminantKnownValues)
{
  // Hand-computable small cases, including a negative determinant.
  EXPECT_NEAR(LinAlg::det(matrix<double>(2, 2, V{1, 2, 3, 4})),     -2.0, 1e-12);
  EXPECT_NEAR(LinAlg::det(matrix<double>(2, 2, V{2, 0, 0, 3})),      6.0, 1e-12);
  // 3x3, det = -306 (classic worked example).
  EXPECT_NEAR(LinAlg::det(matrix<double>(3, 3, V{6, 1, 1,
                                                 4, -2, 5,
                                                 2, 8, 7})),    -306.0, 1e-11);
  // Upper-triangular 4x4: det = product of the diagonal, no elimination needed.
  EXPECT_NEAR(LinAlg::det(matrix<double>(4, 4, V{1, 2, 3, 4,
                                                 0, 5, 6, 7,
                                                 0, 0, 8, 9,
                                                 0, 0, 0, 10})),  400.0, 1e-9);
}

TEST(DenseLinalg, DeterminantSingularIsExactZero)
{
  // A zero column forces an exact zero pivot -> det == 0.0 exactly (a value,
  // not a failure).
  EXPECT_EQ(LinAlg::det(matrix<double>(2, 2, V{0, 1, 0, 2})), 0.0);
  EXPECT_EQ(LinAlg::det(matrix<double>(3, 3, V{0, 1, 2,
                                               0, 3, 4,
                                               0, 5, 6})), 0.0);
  // Two identical rows: the second reduces to all zeros -> exact zero pivot.
  EXPECT_EQ(LinAlg::det(matrix<double>(3, 3, V{1, 2, 3,
                                               1, 2, 3,
                                               4, 5, 6})), 0.0);
}

TEST(DenseLinalg, DeterminantPermutationSign)
{
  // Permutation matrices have det = sign of the permutation, exactly +-1.
  EXPECT_DOUBLE_EQ(LinAlg::det(matrix<double>(2, 2, V{0, 1,
                                                      1, 0})), -1.0);  // 1 swap
  EXPECT_DOUBLE_EQ(LinAlg::det(matrix<double>(3, 3, V{0, 0, 1,
                                                      1, 0, 0,
                                                      0, 1, 0})), 1.0); // 3-cycle
  // 4x4 single transposition (swap rows 0,1): odd -> det = -1.
  EXPECT_DOUBLE_EQ(LinAlg::det(matrix<double>(4, 4, V{0, 1, 0, 0,
                                                      1, 0, 0, 0,
                                                      0, 0, 1, 0,
                                                      0, 0, 0, 1})), -1.0);
}

TEST(DenseLinalg, DeterminantMatchesKnownConstruction)
{
  // det(Q diag(D) Qᵀ) = prod(D), with D of mixed sign and bounded magnitude so
  // the matrix stays well-conditioned across the sizes where OpenBLAS corrupts.
  for (int n : {12, 30, 60}) {
    V D(n);
    for (int k = 0; k < n; k++)
      D[k] = ((k % 2) ? -1.0 : 1.0) * (0.5 + 0.25 * (k % 4));  // in [0.5, 1.25]
    double prodD = 1.0;
    for (double d : D) prodD *= d;
    matrix<double> A = qdqt(D, 1234 + n);
    double det = LinAlg::det(A);
    EXPECT_NEAR(det, prodD, 1e-9 * fabs(prodD)) << "n = " << n;
  }
}

TEST(DenseLinalg, JacobiEigenpairsSatisfyDefinition)
{
  for (int n : {12, 60, 200}) {
    matrix<double> A = random_matrix(n, 4242 + n, /*symmetric=*/true);
    auto d = LinAlg::SymEigen::decompose(A);
    ASSERT_EQ((int)d.lambda.size(), n);
    // Eigenvalues ascending; ||A q_i - lambda_i q_i||_inf small.  SymEigen's
    // Q holds the eigenvectors as COLUMNS, row_residual takes rows.
    for (int i = 0; i + 1 < n; i++) EXPECT_LE(d.lambda[i], d.lambda[i + 1]);
    const double res = row_residual(LinAlg::view_of(transpose_of(d.Q)), d.lambda,
                                    LinAlg::view_of(A));
    EXPECT_LT(res, kEigenpairResidualBar) << "n = " << n;
  }
}

TEST(DenseLinalg, PseudoinverseSolvesFullRankSystem)
{
  const int n = 60;
  matrix<double> A = random_matrix(n, 99, /*symmetric=*/true);
  V x_true = random_vector(n, 3);
  V b = LinAlg::matvec(A, x_true);
  auto sol = LinAlg::SymEigen::solve(A, b, /*rcond=*/1e-12);
  ASSERT_EQ((int)sol.x.size(), n);
  EXPECT_EQ(sol.rank, n);
  double err = 0;
  for (int i = 0; i < n; i++) err = max(err, fabs(sol.x[i] - x_true[i]));
  EXPECT_LT(err, 1e-8);
}

// ============================================================================
// View-level conformance (dense_linalg_view.hh).  The owner API and
// matrix::operator* both DELEGATE to the view bodies, so owner-vs-view is no
// pin; the independent instruments are:
//   - the FROZEN HISTORICAL ORACLE above (the pre-promotion bodies,
//     byte-compared -- kills reassociation and pivot-rule drift the
//     value-level tests are measured-blind to),
//   - the strided live-block path vs the packed one (the batch-port shape,
//     exercised by no owner call), and
//   - the singular zero-fill contract, incl. a LATE zero pivot (partially
//     forward-eliminated RHS).
// ============================================================================

namespace {

// Hilbert matrix: dense, ill-conditioned, exactly representable inputs.
matrix<double> hilbert(int n)
{
  matrix<double> A(n, n, 0.0);
  for (int i = 0; i < n; i++)
    for (int j = 0; j < n; j++) A(i, j) = 1.0 / (i + j + 1);
  return A;
}

// Exact-magnitude pivot ties in column 0 -- the discriminator for the
// pivot selection rule (> vs >=), which no random matrix produces.
matrix<double> tie3()
{
  return matrix<double>(3, 3, V{1, 2, 0, -1, 5, 1, 1, 0, 7});
}

// Near-singular: row 3 = row 0 + row 1 over random doubles.  NOT exactly
// singular in floating point (the eliminated pivot is ~1e-16 rounding
// residue, and the solve "succeeds" with huge components) -- a conformance
// input for the frozen oracle, not a zero-fill case.
matrix<double> near_singular(unsigned seed)
{
  matrix<double> A = random_matrix(5, seed, false);
  for (int j = 0; j < 5; j++) A(3, j) = A(0, j) + A(1, j);
  return A;
}

// EXACTLY singular at the LAST pivot: integer entries with row 2 a bitwise
// duplicate of row 0.  No swap at c = 0 (row 0 holds the strict max), the
// duplicate eliminates with multiplier exactly 1.0 to an all-zero row, and
// the zero pivot appears at c = 2 -- after the RHS has been forward-
// eliminated through two columns.
matrix<double> singular_late_exact()
{
  return matrix<double>(3, 3, V{4, 1, 7,  2, 3, 5,  4, 1, 7});
}

void expect_bits_equal(const V& a, const V& b, const char* tag)
{
  ASSERT_EQ(a.size(), b.size()) << tag;
  if (a.empty()) return;                     // equal, and memcmp wants real pointers
  EXPECT_EQ(memcmp(a.data(), b.data(), a.size() * sizeof(double)), 0) << tag;
}

// The full solve family, new vs frozen, bitwise.
void expect_matches_frozen(const matrix<double>& A, const V& b, const char* tag)
{
  expect_bits_equal(LinAlg::solve(A, b), frozen::solve(A, b), tag);
  for (double lam : {0.37, 1e-14, 0.0, -0.5})
    expect_bits_equal(LinAlg::solve_shifted(A, b, lam),
                      frozen::solve_shifted(A, b, lam), tag);
  expect_bits_equal(LinAlg::matvec(A, b), frozen::matvec(A, b), tag);

  const double dn = LinAlg::det(A), df = frozen::det(A);
  EXPECT_EQ(memcmp(&dn, &df, sizeof(double)), 0) << tag << " det";

  auto sn = LinAlg::solve_with_sign(A, b);
  auto sf = frozen::solve_with_sign(A, b);
  ASSERT_EQ(sn.has_value(), sf.has_value()) << tag;
  if (sn) {
    EXPECT_EQ(sn->det_sign, sf->det_sign) << tag;
    expect_bits_equal(sn->x, sf->x, tag);
  }

  const double dotn = LinAlg::dot(b, b);
  double dotf = 0;
  for (size_t i = 0; i < b.size(); i++) dotf += b[i] * b[i];   // frozen dot
  EXPECT_EQ(memcmp(&dotn, &dotf, sizeof(double)), 0) << tag << " dot";
}

}  // namespace

TEST(DenseLinalgFrozenOracle, SolveFamilyBitIdenticalToHistory)
{
  for (int n : {5, 12, 17, 60, 128})
    expect_matches_frozen(random_matrix(n, 300 + n, /*symmetric=*/false),
                          random_vector(n, 400 + n), "random");
  for (int n : {5, 12, 17})
    expect_matches_frozen(hilbert(n), random_vector(n, 500 + n), "hilbert");
  expect_matches_frozen(tie3(), random_vector(3, 42), "pivot-tie");
  expect_matches_frozen(near_singular(77), random_vector(5, 78), "near-singular");
  expect_matches_frozen(singular_late_exact(), random_vector(3, 79),
                        "singular-late-exact");
}

TEST(DenseLinalgFrozenOracle, MatmulBitIdenticalToHistory)
{
  // matrix::operator* now DELEGATES to LinAlg::matmul, so operator*-vs-view
  // is wiring, not a pin; the frozen historical product is the oracle.
  for (int n : {3, 7, 17, 60, 128}) {
    matrix<double> A = random_matrix(n, 101 + n, /*symmetric=*/false);
    matrix<double> B = random_matrix(n, 202 + n, /*symmetric=*/false);
    matrix<double> Cf = frozen::matmul(A, B);
    matrix<double> Cd = A * B;                       // the delegation wiring
    std::vector<double> out((size_t)n * n);
    LinAlg::matmul(LinAlg::view_of(A), LinAlg::view_of(B), out);
    ASSERT_EQ(memcmp(&Cf[0], out.data(), (size_t)n * n * sizeof(double)), 0)
        << "view matmul must be bit-identical to the frozen product, n = " << n;
    ASSERT_EQ(memcmp(&Cf[0], &Cd[0], (size_t)n * n * sizeof(double)), 0)
        << "operator* must delegate to the same body, n = " << n;
  }
  // Rectangular operands pin the m x kk x n bookkeeping of the one product
  // body against the frozen loop (the square cases above cannot tell m
  // from kk from n).
  struct Shape { int m, k, n; };
  for (Shape s : {Shape{3, 4, 5}, Shape{17, 9, 13}, Shape{1, 7, 1}}) {
    matrix<double> A = random_rectangular(s.m, s.k, 301 + s.m);
    matrix<double> B = random_rectangular(s.k, s.n, 302 + s.n);
    matrix<double> Cf = frozen::matmul(A, B);
    std::vector<double> out((size_t)s.m * s.n);
    LinAlg::matmul(LinAlg::view_of(A), LinAlg::view_of(B), out);
    ASSERT_EQ(memcmp(&Cf[0], out.data(), out.size() * sizeof(double)), 0)
        << "view matmul must be bit-identical to the frozen product, m k n = "
        << s.m << ' ' << s.k << ' ' << s.n;
  }
}

TEST(DenseLinalgView, StridedLiveBlockMatchesPacked)
{
  // The batch-port shape: the live n x n block embedded in an lda-strided
  // buffer must compute bit-identically to the packed owner path.  Scratch
  // is POISONED, so a slot the copy loops fail to overwrite is detected.
  const int n   = 17;
  const int lda = 23;
  matrix<double> A = random_matrix(n, 7, /*symmetric=*/false);
  matrix<double> B = random_matrix(n, 9, /*symmetric=*/false);
  V b = random_vector(n, 8);
  auto widen = [&](const matrix<double>& S) {
    std::vector<double> wide((size_t)n * lda, -999.0);
    for (int i = 0; i < n; i++)
      for (int j = 0; j < n; j++) wide[(size_t)i * lda + j] = S(i, j);
    return wide;
  };
  std::vector<double> Awide = widen(A), Bwide = widen(B);
  const LinAlg::MatConstView Av{Awide, n, n, lda};
  const LinAlg::MatConstView Bv{Bwide, n, n, lda};

  std::vector<double> M((size_t)n * n, -999.0), M2((size_t)n * n, -999.0);
  std::vector<double> x(n), xs(n), mv(n), mm((size_t)n * n);
  ASSERT_EQ(LinAlg::solve(Av, b, M, x), LinAlg::LuStatus::Ok);
  expect_bits_equal(V(x.begin(), x.end()), LinAlg::solve(A, b), "strided solve");

  ASSERT_EQ(LinAlg::solve_shifted(Av, b, 0.37, M2, M, xs), LinAlg::LuStatus::Ok);
  expect_bits_equal(V(xs.begin(), xs.end()), LinAlg::solve_shifted(A, b, 0.37),
                    "strided solve_shifted");

  LinAlg::matvec(Av, b, mv);
  expect_bits_equal(V(mv.begin(), mv.end()), LinAlg::matvec(A, b), "strided matvec");

  // Strided matmul -- BOTH operands strided: exactly the Alexandrov batch
  // shape (matmul(Jview, Jview, JtJ) with lda = nv_cap != n).
  LinAlg::matmul(Av, Bv, mm);
  std::vector<double> mm_packed((size_t)n * n);
  LinAlg::matmul(LinAlg::view_of(A), LinAlg::view_of(B), mm_packed);
  EXPECT_EQ(memcmp(mm.data(), mm_packed.data(), (size_t)n * n * sizeof(double)), 0)
      << "strided matmul must match packed";
}

TEST(DenseLinalgView, SingularSolveZeroFillsThroughBothAPIs)
{
  {   // pivot-0 singular: x is the raw b copy when the guard fires
    const int n = 5;
    matrix<double> A(n, n, 0.0);
    V b = random_vector(n, 9);
    std::vector<double> M((size_t)n * n), x(n, 7.0);
    EXPECT_EQ(LinAlg::solve(LinAlg::view_of(A), b, M, x), LinAlg::LuStatus::Singular);
    for (int i = 0; i < n; i++) EXPECT_EQ(x[i], 0.0);
    V xo = LinAlg::solve(A, b);
    for (int i = 0; i < n; i++) EXPECT_EQ(xo[i], 0.0);
  }
  {   // LATE zero pivot: the RHS is partially forward-eliminated before the
      // trip -- the zero-fill must still win, through every entry point.
    matrix<double> A = singular_late_exact();
    V b = random_vector(3, 32);
    std::vector<double> M(9), x(3, 7.0);
    EXPECT_EQ(LinAlg::solve(LinAlg::view_of(A), b, M, x), LinAlg::LuStatus::Singular);
    for (int i = 0; i < 3; i++) EXPECT_EQ(x[i], 0.0);
    V xo = LinAlg::solve(A, b);
    for (int i = 0; i < 3; i++) EXPECT_EQ(xo[i], 0.0);
    EXPECT_FALSE(LinAlg::solve_with_sign(A, b).has_value());
    EXPECT_EQ(LinAlg::det(A), 0.0);
  }
}

TEST(DenseLinalgView, VectorAssignmentsBitIdenticalToOwnerExpressions)
{
  // Each span word restates an owner V expression elementwise with one
  // rounding per element -- so the results must be BIT-identical, not just
  // close.  These are the expressions the Alexandrov solver family writes
  // (r = r - dr, r_trial = r + delta, F -= t1*kappa1, result += r_j*basis).
  const int n = 17;
  const V a = random_vector(n, 3), b = random_vector(n, 5);
  const double s = 0.37;

  V x(a);
  LinAlg::copy_into(x, b);
  expect_bits_equal(x, b, "copy_into");

  LinAlg::neg_into(x, a);
  expect_bits_equal(x, -a, "neg_into");

  x = a; LinAlg::add_into(x, b);
  expect_bits_equal(x, a + b, "add_into");

  x = a; LinAlg::sub_into(x, b);
  expect_bits_equal(x, a - b, "sub_into");

  x = a; LinAlg::add_scaled(x, s, b);
  expect_bits_equal(x, a + b * s, "add_scaled");

  x = a; LinAlg::sub_scaled(x, s, b);
  expect_bits_equal(x, a - b * s, "sub_scaled");

  const double lhs = LinAlg::sum(a);
  const double rhs = LinAlg::dot(a, V(a.size(), 1.0));
  EXPECT_EQ(memcmp(&lhs, &rhs, sizeof lhs), 0) << "sum must be dot(v, ones) bitwise";
}

// ============================================================================
// The pivot as a total order (dense_linalg_view.hh pivot_key): a collective
// evaluation -- the largest key over the rows q >= c, then the smallest row
// carrying it -- must select the row the sequential scan (pivot_row) selects,
// on exactly the columns no random matrix produces: NaN at the diagonal, NaN
// below it, signed zeros, infinity, equal magnitudes of opposite sign, ties,
// subnormals.  This is what lets a lane-parallel elimination reproduce the
// sequential one byte for byte.
// ============================================================================

namespace {

// The collective's selection rule, evaluated sequentially: max key, ties to
// the smaller row (the order is what is pinned; any evaluation of it agrees).
template <class T>
int pivot_by_key(LinAlg::MatView<const T> M, int c)
{
  int p = c;
  uint64_t best = LinAlg::pivot_key<T>(M(c, c), true);
  for (int q = c + 1; q < M.n; q++) {
    const uint64_t k = LinAlg::pivot_key<T>(M(q, c), false);
    if (k > best) { best = k; p = q; }
  }
  return p;
}

template <class T>
void expect_pivot_key_matches_scan(std::vector<T> column, const char* tag)
{
  const LinAlg::MatView<const T> M{std::span<const T>(column), (int)column.size(), 1, 1};
  EXPECT_EQ(pivot_by_key<T>(M, 0), LinAlg::pivot_row<T>(M, 0)) << tag;
}

template <class T>
void pivot_key_columns()
{
  const T nan = std::numeric_limits<T>::quiet_NaN();
  const T inf = std::numeric_limits<T>::infinity();
  const T sub = std::numeric_limits<T>::denorm_min();
  expect_pivot_key_matches_scan<T>({nan, 1, 2}, "NaN at the diagonal keeps row c");
  expect_pivot_key_matches_scan<T>({1, nan, 3}, "NaN below the diagonal is skipped");
  expect_pivot_key_matches_scan<T>({0, nan}, "zero pivot with a NaN below stays at c");
  expect_pivot_key_matches_scan<T>({nan, nan}, "all NaN keeps row c");
  expect_pivot_key_matches_scan<T>({T(0), T(-0.0), T(0)}, "signed zeros tie to the first");
  expect_pivot_key_matches_scan<T>({T(-0.0), T(0)}, "negative zero at the diagonal");
  expect_pivot_key_matches_scan<T>({1, inf, 2}, "infinity wins");
  expect_pivot_key_matches_scan<T>({2, -2, 2}, "equal magnitudes tie to the first");
  expect_pivot_key_matches_scan<T>({1, -3, 3}, "opposite signs, first of the tie");
  expect_pivot_key_matches_scan<T>({sub, 0, T(-2) * sub}, "subnormals order by magnitude");
  expect_pivot_key_matches_scan<T>({T(0.5), T(-0.5), T(0.25)}, "exact halves tie");
}

}  // namespace

TEST(DenseLinalgView, PivotKeyOrderMatchesScan)
{
  pivot_key_columns<double>();
  pivot_key_columns<float>();
  // Every column of every step of random reductions, and one with a NaN
  // planted on the diagonal midway: the order agrees at c > 0 as at c == 0.
  for (int n = 2; n <= 9; n++) {
    matrix<double> A = random_matrix(n, 700 + n, /*symmetric=*/false);
    if (n == 6) A(3, 3) = std::numeric_limits<double>::quiet_NaN();
    const LinAlg::MatConstView M = LinAlg::view_of(A);
    for (int c = 0; c < n; c++)
      EXPECT_EQ(pivot_by_key<double>(M, c), LinAlg::pivot_row<double>(M, c))
          << "random n = " << n << ", column " << c;
  }
}

// ============================================================================
// The Jacobi rotation family and the transposed products (promoted
// 2026-09-24): jacobi_eig against its frozen self, bitwise; A B^T and A^T B
// against the plain product on materialised transposes, bitwise; the float
// instantiation against the eigenpair definition to a counted bound; the
// accumulator's warm start.
// ============================================================================
namespace {

// (matrix<double> is a std::vector<double>, so it binds to jacobi_eig's
// by-value flat argument as a copy of its row-major storage.)
void expect_jacobi_matches_frozen(const matrix<double>& A, bool vectors, const char* tag)
{
  const int n = A.m;
  std::vector<double> lam_n, lam_f, V_n, V_f;
  const bool ok_n = LinAlg::jacobi_eig(A, n, lam_n, vectors ? &V_n : nullptr);
  const bool ok_f = frozen::jacobi_eig(A, n, lam_f, vectors ? &V_f : nullptr);
  EXPECT_EQ(ok_n, ok_f) << tag;
  expect_bits_equal(lam_n, lam_f, tag);
  if (vectors) expect_bits_equal(V_n, V_f, tag);
}

// Two random symmetric blocks on the diagonal, exact zeros elsewhere: a
// rotation inside a block leaves the cross-block entries exactly zero, so
// every cross-block pair takes the skip in every sweep.
matrix<double> block_diagonal(int n1, int n2, unsigned seed)
{
  matrix<double> A(n1 + n2, n1 + n2, 0.0);
  const matrix<double> B1 = random_matrix(n1, seed, true);
  const matrix<double> B2 = random_matrix(n2, seed + 1, true);
  for (int i = 0; i < n1; i++) for (int j = 0; j < n1; j++) A(i, j) = B1(i, j);
  for (int i = 0; i < n2; i++) for (int j = 0; j < n2; j++) A(n1 + i, n1 + j) = B2(i, j);
  return A;
}

}  // namespace

TEST(DenseLinalgFrozenOracle, JacobiEigBitIdenticalToHistory)
{
  for (int n : {3, 12, 60, 128, 240}) {
    expect_jacobi_matches_frozen(random_matrix(n, 700 + n, true), /*vectors=*/true,  "random+V");
    expect_jacobi_matches_frozen(random_matrix(n, 800 + n, true), /*vectors=*/false, "random");
  }
  for (int n : {8, 30}) expect_jacobi_matches_frozen(hilbert(n), true, "hilbert");
  {  // exactly diagonal: converged at the first test, no sweep
    matrix<double> D(20, 20, 0.0);
    for (int i = 0; i < 20; i++) D(i, i) = 20 - i;
    expect_jacobi_matches_frozen(D, true, "diagonal");
  }
  expect_jacobi_matches_frozen(block_diagonal(20, 20, 900), true, "block-diagonal");
  {
    V D(50);
    for (int k = 0; k < 50; k++) D[k] = (k % 2 ? 1.0 : -1.0) * (1.0 + k);
    expect_jacobi_matches_frozen(qdqt(D, 950), true, "qdqt");
  }
  // Non-symmetric inputs trip the sweep guard: both bodies return false and
  // neither writes lam or V (compared as the empty vectors they stay).
  for (int n : {2, 4, 6})
    expect_jacobi_matches_frozen(random_matrix(n, 990 + n, /*symmetric=*/false), true, "guard-trip");
  // No pair to rotate: converged at the first test.
  expect_jacobi_matches_frozen(matrix<double>(0, 0, 0.0), true, "n=0");
  expect_jacobi_matches_frozen(matrix<double>(1, 1, V{2.5}), true, "n=1");
}

TEST(DenseLinalgView, TransposedProductsBitIdenticalToMatmul)
{
  // A B^T against A (B^T) and A^T B against (A^T) B: the same products in
  // the same k order, so the claim is equal bits.  The plain product's own
  // pin is MatmulBitIdenticalToHistory.
  struct Shape { int m, k, n; };
  for (Shape s : {Shape{3, 4, 5}, Shape{17, 9, 13}, Shape{60, 60, 60}, Shape{1, 7, 1}}) {
    std::vector<double> got((size_t)s.m * s.n), ref((size_t)s.m * s.n);

    const matrix<double> A  = random_rectangular(s.m, s.k, 1000 + s.m);   // m x k
    const matrix<double> Bt = random_rectangular(s.n, s.k, 2000 + s.n);   // n x k, B = Bt^T
    LinAlg::matmul_nt(LinAlg::view_of(A), LinAlg::view_of(Bt), got);
    LinAlg::matmul(LinAlg::view_of(A), LinAlg::view_of(transpose_of(Bt)), ref);
    EXPECT_EQ(memcmp(got.data(), ref.data(), got.size() * sizeof(double)), 0)
        << "A B^T, m k n = " << s.m << ' ' << s.k << ' ' << s.n;

    const matrix<double> At = random_rectangular(s.k, s.m, 3000 + s.m);   // k x m, A = At^T
    const matrix<double> B  = random_rectangular(s.k, s.n, 4000 + s.n);   // k x n
    LinAlg::matmul_tn(LinAlg::view_of(At), LinAlg::view_of(B), got);
    LinAlg::matmul(LinAlg::view_of(transpose_of(At)), LinAlg::view_of(B), ref);
    EXPECT_EQ(memcmp(got.data(), ref.data(), got.size() * sizeof(double)), 0)
        << "A^T B, m k n = " << s.m << ' ' << s.k << ' ' << s.n;
  }
}

TEST(DenseLinalgView, JacobiFamilyFloatSatisfiesDefinition)
{
  // The float instantiation has no history to freeze; the claim is the
  // eigenpair definition.  A rigorous rounding count gives no usable bound
  // at these sizes (every entry receives 2(n - 1) rotation writes per
  // sweep, and errors migrate between entries under later rotations), so
  // the bar has the first-order FORM n (tol + K S eps max|A|), S the sweeps
  // performed, with K = 10 an EMPIRICAL constant: the residual measured at
  // n = 12, 60, 128 and 240 lies 18 to 31 times below it.  The bar's job
  // is to catch a broken instantiation, whose residual is of order
  // n max|A|, not to certify precision.
  using LinAlg::MatView;
  const float eps = std::numeric_limits<float>::epsilon();
  for (int n : {12, 60, 240}) {
    const matrix<double> Ad = random_matrix(n, 5000 + n, true);
    std::vector<float> A((size_t)n * n), Vf((size_t)n * n, 0.0f);
    float anorm = 0;
    for (int i = 0; i < n; i++)
      for (int j = 0; j < n; j++) {
        A[(size_t)i * n + j] = (float)Ad(i, j);
        anorm = std::max(anorm, std::fabs(A[(size_t)i * n + j]));
      }
    const std::vector<float> A0 = A;                  // the matrix decomposed
    for (int i = 0; i < n; i++) Vf[(size_t)i * n + i] = 1.0f;
    const float tol = 4.5f * eps * anorm;             // jacobi_eig's 1e-15 is 4.5 eps64
    const MatView<float> Av{A, n, n, n}, Vv{Vf, n, n, n};
    const LinAlg::JacobiResult run = LinAlg::jacobi_diagonalize(Av, Vv, tol, 60);
    ASSERT_TRUE(run.converged) << "n = " << n;
    EXPECT_LE(LinAlg::off_diagonal_max<float>(MatView<const float>{Av}), tol);
    std::vector<double> lambda(n);
    for (int m = 0; m < n; m++) lambda[m] = Av(m, m);
    const MatView<const float> A0v{A0, n, n, n};
    const double res   = row_residual(Vv, lambda, A0v);
    const double bound = n * ((double)tol + 10.0 * run.sweeps * eps * anorm);
    EXPECT_LE(res, bound) << "n = " << n << ", sweeps = " << run.sweeps;
  }
}

TEST(DenseLinalgView, JacobiAccumulatorWarmStart)
{
  // jacobi_diagonalize accumulates into whatever V holds.  From an
  // eigenbasis V of A (rows), M = V A V^T is diagonal to rounding, and a
  // second run on (M, V) must leave V an eigenbasis of A within a sweep or
  // two -- the warm start a self-consistent-field iteration performs,
  // re-diagonalising each cycle's matrix in the previous cycle's
  // eigenbasis.  M is formed through the plain and the transposed product.
  using LinAlg::MatView;
  const int n = 60;
  const matrix<double> Am = random_matrix(n, 6000, true);
  std::vector<double> A = Am, Vd((size_t)n * n, 0.0);
  for (int i = 0; i < n; i++) Vd[(size_t)i * n + i] = 1.0;
  double anorm = 0;
  for (double x : A) anorm = std::max(anorm, std::fabs(x));
  const double tol = 1e-15 * anorm;
  const MatView<double> Av{A, n, n, n}, Vv{Vd, n, n, n};
  ASSERT_TRUE(LinAlg::jacobi_diagonalize(Av, Vv, tol, 60).converged);

  std::vector<double> T1((size_t)n * n), M((size_t)n * n);
  LinAlg::matmul(MatView<const double>{Vv}, LinAlg::view_of(Am), T1);
  LinAlg::matmul_nt(MatView<const double>{T1, n, n, n}, MatView<const double>{Vv}, M);
  const MatView<double> Mv{M, n, n, n};
  const LinAlg::JacobiResult again = LinAlg::jacobi_diagonalize(Mv, Vv, tol, 60);
  ASSERT_TRUE(again.converged);
  EXPECT_LE(again.sweeps, 3) << "a warm start from an eigenbasis should need at most a few sweeps";

  std::vector<double> lambda(n);
  for (int m = 0; m < n; m++) lambda[m] = Mv(m, m);
  EXPECT_LT(row_residual(Vv, lambda, LinAlg::view_of(Am)), kEigenpairResidualBar);
}
