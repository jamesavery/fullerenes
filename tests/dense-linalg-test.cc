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
