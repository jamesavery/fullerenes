// dense_linalg bit-exactness tests.  This module exists because a deployed
// OpenBLAS silently returned wrong dgesv solutions from n ~ 60 (info = 0),
// so the in-house LU and Jacobi eigensolver carry their own known-solution
// tests across the sizes where the corruption was observed.

#include "fullerenes/dense_linalg.hh"

#include <gtest/gtest.h>
#include <cmath>
#include <cstring>
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

}  // namespace frozen

namespace {

matrix<double> random_matrix(int n, unsigned seed, bool symmetric)
{
  mt19937 rng(seed);
  uniform_real_distribution<double> u(-1.0, 1.0);
  matrix<double> A(n, n, 0.0);
  for (int i = 0; i < n; i++)
    for (int j = symmetric ? i : 0; j < n; j++) {
      A(i, j) = u(rng);
      if (symmetric) A(j, i) = A(i, j);
    }
  return A;
}

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
    // ||A q_i - lambda_i q_i||_inf small, eigenvalues ascending.
    double res = 0;
    for (int i = 0; i < n; i++) {
      if (i + 1 < n) EXPECT_LE(d.lambda[i], d.lambda[i + 1]);
      for (int r = 0; r < n; r++) {
        double acc = 0;
        for (int c = 0; c < n; c++) acc += A(r, c) * d.Q(c, i);
        res = max(res, fabs(acc - d.lambda[i] * d.Q(r, i)));
      }
    }
    EXPECT_LT(res, 1e-12) << "n = " << n;
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
