// dense_linalg bit-exactness tests.  This module exists because a deployed
// OpenBLAS silently returned wrong dgesv solutions from n ~ 60 (info = 0),
// so the in-house LU and Jacobi eigensolver carry their own known-solution
// tests across the sizes where the corruption was observed.

#include "fullerenes/dense_linalg.hh"

#include <gtest/gtest.h>
#include <cmath>
#include <random>
#include <vector>

using namespace std;
using LinAlg::V;

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
