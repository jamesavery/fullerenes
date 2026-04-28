// Analytical validator for the icosahedral series of Goldberg-Coxeter
// fullerenes — those with full I_h symmetry: GC(m,n) where m·n = 0 OR
// m = n.  For these, the 12 cone points are arranged with full
// icosahedral symmetry, so the Alexandrov polytope is a *regular
// icosahedron* (up to overall scale):
//
//   - all 30 iDT edges equal length
//   - all 20 face triangles equilateral
//   - vol_norm = V / ⟨ℓ⟩³ = 5(3 + √5) / 12 ≈ 2.181694990624912
//
// The intrinsic edge length d_GC(m,n) (= surface geodesic between any
// two adjacent cone points on the all-edge-1 piecewise-flat metric)
// has a closed form: d² = m² + m·n + n² in Eisenstein integers, so
// d = √(m² + m·n + n²).  But since vol_norm and "all edges equal" are
// scale-invariant, we don't need to predict d; the regular-icosahedron
// shape is the assertion.
//
// This test is the strongest check we have for solver correctness:
// for the entire Ih GC family, the analytical answer is known and
// must match.  Failure on any size means the solver is buggy.

#include <gtest/gtest.h>
#include <cmath>
#include <vector>

#include "fullerenes/triangulation.hh"
#include "fullerenes/delaunay.hh"
#include "fullerenes/delaunay_alexandrov12.hh"

using std::vector;

// Regular icosahedron normalized volume: V / a³ for unit edge length a.
//   V_icosa = (5/12)(3 + √5) a³
//   ⟨ℓ⟩ = a, ⟨ℓ⟩³ = a³
// so vol_norm = (5/12)(3 + √5).
static constexpr double VOL_NORM_ICOSA = 5.0 * (3.0 + 1.4142135623730951 + 0.7320508075688772) / 12.0;

// More precise: 5*(3+sqrt(5))/12 — compute at runtime.
static double vol_norm_icosa() {
  return 5.0 * (3.0 + std::sqrt(5.0)) / 12.0;
}

namespace {

// Build the GC(k,l) icosahedral fullerene's dual triangulation.
Triangulation make_gc_dual(int k, int l) {
  // The icosahedron itself is the dual of C20 (12 deg-5 vertices).
  Triangulation icosa(vector<int>(12, 5));
  return icosa.GCtransform(k, l);
}

struct IcoCheckResult {
  AlexandrovSolver::ValidationStatus status;
  double vol_norm;
  double edge_min, edge_max, edge_mean;
  double edge_relative_spread;     // (max - min) / mean — robust to cancellation
};

IcoCheckResult run_solver_and_measure(int k, int l) {
  Triangulation T = make_gc_dual(k, l);
  auto D = DelaunayTriangulation::compute(T);
  AlexandrovSolver solver;
  solver.D = D;
  solver.solve();

  IcoCheckResult r;
  r.status = solver.stats_status;
  r.vol_norm = solver.stats_volume_norm;

  // Edge-length stats over alive non-bigon iDT edges.  Use (max-min)/mean
  // as the spread metric — std-via-sum-of-squares loses ~half the digits
  // to catastrophic cancellation when values are nearly equal.
  double sum = 0;
  int n = 0;
  r.edge_min = std::numeric_limits<double>::infinity();
  r.edge_max = -std::numeric_limits<double>::infinity();
  for (int h = 0; h < solver.D.nh; h += 2) {
    if (!solver.D.alive(h)) continue;
    double L = solver.D.he_length[h];
    sum += L; n++;
    if (L < r.edge_min) r.edge_min = L;
    if (L > r.edge_max) r.edge_max = L;
  }
  r.edge_mean = (n > 0) ? sum / n : 0;
  r.edge_relative_spread = (r.edge_mean > 0)
                            ? (r.edge_max - r.edge_min) / r.edge_mean
                            : 0;
  return r;
}

void expect_regular_icosahedron(const char* label, int k, int l,
                                  double vol_tol = 1e-6,
                                  double edge_spread_tol = 1e-12) {
  SCOPED_TRACE(label);
  auto r = run_solver_and_measure(k, l);
  EXPECT_EQ(r.status, AlexandrovSolver::ValidationStatus::OK)
    << "solver failed; status=" << AlexandrovSolver::status_str(r.status);
  EXPECT_NEAR(r.vol_norm, vol_norm_icosa(), vol_tol)
    << "vol_norm mismatch: got " << r.vol_norm
    << ", expected " << vol_norm_icosa();
  EXPECT_LT(r.edge_relative_spread, edge_spread_tol)
    << "edge lengths not all equal: (max−min)/mean="
    << r.edge_relative_spread
    << " (min=" << r.edge_min << ", max=" << r.edge_max
    << ", mean=" << r.edge_mean << ")";
}

}  // namespace

TEST(AlexandrovIcosahedral, GC_1_0_C20) {
  expect_regular_icosahedron("GC(1,0) = C20-Ih", 1, 0);
}

TEST(AlexandrovIcosahedral, GC_1_1_C60) {
  expect_regular_icosahedron("GC(1,1) = C60-Ih (Buckminsterfullerene)", 1, 1);
}

TEST(AlexandrovIcosahedral, GC_2_0_C80) {
  expect_regular_icosahedron("GC(2,0) = C80-Ih", 2, 0);
}

TEST(AlexandrovIcosahedral, GC_2_2_C240) {
  expect_regular_icosahedron("GC(2,2) = C240-Ih", 2, 2);
}

TEST(AlexandrovIcosahedral, GC_3_0_C180) {
  expect_regular_icosahedron("GC(3,0) = C180-Ih", 3, 0);
}

TEST(AlexandrovIcosahedral, GC_3_3_C540) {
  expect_regular_icosahedron("GC(3,3) = C540-Ih", 3, 3);
}

TEST(AlexandrovIcosahedral, GC_4_0_C320) {
  expect_regular_icosahedron("GC(4,0) = C320-Ih", 4, 0);
}

TEST(AlexandrovIcosahedral, GC_5_0_C500) {
  expect_regular_icosahedron("GC(5,0) = C500-Ih", 5, 0);
}
