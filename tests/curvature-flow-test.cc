// curvature-flow-test -- confine_curvature (fullerenes/curvature-flow.hh) and the DCEL
// primitives it graduated onto DelaunayTriangulation (split_face, geodesic_disks,
// curvature, incident). Asserts the intrinsic invariants the curvature-flow sub-project
// validates on real tomogram meshes, here on equilateral fullerene-dual metrics:
//   Gauss-Bonnet      Sigma K* = 4*pi
//   metric preserved  every original vertex keeps its angle defect; cones are flat
//   DCEL consistent   well-formed + consistent after the 1->3 splits

#include <gtest/gtest.h>
#include "fullerenes/curvature-flow.hh"
#include "fullerenes/delaunay.hh"
#include "fullerenes/buckygen-wrapper.hh"

#include <cmath>
#include <set>
#include <vector>

// The first fullerene dual triangulation of size N (buckygen; needs no database).
static Triangulation first_dual(int N) {
  BuckyGen::buckygen_queue Q = BuckyGen::start(N, false, false);
  Triangulation T;
  BuckyGen::next_fullerene(Q, T);
  BuckyGen::stop(Q);
  return T;
}

// Greedily pick up to `want` faces that pairwise share no vertex, as cone sites at the
// face centroid -- so the 1->3 splits never touch a face an earlier split disturbed.
static std::vector<ConeSite> spread_cone_sites(const DelaunayTriangulation& D, int want) {
  std::vector<ConeSite> sites;
  std::set<int> used;
  for (int f = 0; f < D.nf && (int)sites.size() < want; f++) {
    int h = D.f_he[f];
    if (h < 0) continue;                                  // dead face slot
    int a = D.he_origin[h], b = D.he_origin[D.he_next[h]], c = D.he_origin[D.he_next[D.he_next[h]]];
    if (used.count(a) || used.count(b) || used.count(c)) continue;
    sites.push_back({tri_t(a, b, c), {1.0/3, 1.0/3, 1.0/3}});
    used.insert(a); used.insert(b); used.insert(c);
  }
  return sites;
}

TEST(CurvatureFlow, ConfinesAllCurvatureToDisks) {
  DelaunayTriangulation D = DelaunayTriangulation::from_triangulation(first_dual(60));
  const std::vector<double> K_before = D.curvature();
  const int nv_before = D.nv;

  std::vector<ConeSite> sites = spread_cone_sites(D, 8);
  ASSERT_GE(sites.size(), 1u);

  ConePrescription pre = confine_curvature(D, sites, 2.0, CurvatureMode::Uniform, DiskMetric::Edge);

  // Gauss-Bonnet: all 4*pi confined, regardless of cone count (the per-disk pi/3 is
  // rescaled to a total of 4*pi).
  double sum = 0; for (double k : pre.kstar) sum += k;
  EXPECT_NEAR(sum, 4.0 * M_PI, 1e-9);

  // DCEL stays well-formed and consistent after the 1->3 splits.
  EXPECT_TRUE(pre.surface.is_well_formed());
  EXPECT_TRUE(pre.surface.check_consistency());

  // The inserted cones are flat: split_face preserves the metric exactly.
  EXPECT_EQ(pre.cones.size(), sites.size());
  for (int cone : pre.cones)
    EXPECT_NEAR(pre.surface.vertex_angle_sum(cone), 2.0 * M_PI, 1e-9);

  // Metric preserved: every original vertex keeps its angle defect.
  const std::vector<double> K_after = pre.surface.curvature();
  double max_dK = 0;
  for (int v = 0; v < nv_before; v++)
    max_dK = std::max(max_dK, std::fabs(K_after[v] - K_before[v]));
  EXPECT_LT(max_dK, 1e-9);
}

// Both disk metrics partition disjointly and confine 4*pi; Unfold's fast-marching must
// not break Gauss-Bonnet or DCEL consistency.
TEST(CurvatureFlow, BothDiskMetricsConfine4Pi) {
  DelaunayTriangulation D = DelaunayTriangulation::from_triangulation(first_dual(80));
  std::vector<ConeSite> sites = spread_cone_sites(D, 6);
  ASSERT_GE(sites.size(), 1u);

  for (DiskMetric m : {DiskMetric::Edge, DiskMetric::Unfold}) {
    ConePrescription pre = confine_curvature(D, sites, 2.0, CurvatureMode::Gaussian, m);
    double sum = 0; for (double k : pre.kstar) sum += k;
    EXPECT_NEAR(sum, 4.0 * M_PI, 1e-9);
    EXPECT_TRUE(pre.surface.check_consistency());
  }
}

// split_face is metric-preserving and degree-correct in isolation: one split of a
// single face yields a degree-3 flat vertex and three live faces, total +2.
TEST(CurvatureFlow, SplitFaceIsMetricPreserving) {
  DelaunayTriangulation D = DelaunayTriangulation::from_triangulation(first_dual(60));
  const std::vector<double> K_before = D.curvature();
  const int nv_before = D.nv;

  int h0 = D.f_he[0];
  // Spoke lengths to the centroid of an equilateral (unit-edge) face.
  std::array<double,3> spoke = { 1.0/std::sqrt(3.0), 1.0/std::sqrt(3.0), 1.0/std::sqrt(3.0) };
  int P = D.split_face(h0, spoke);

  EXPECT_EQ(P, nv_before);
  EXPECT_EQ(D.vertex_degree(P), 3);
  EXPECT_NEAR(D.vertex_angle_sum(P), 2.0 * M_PI, 1e-9);   // flat
  EXPECT_TRUE(D.is_well_formed());
  EXPECT_TRUE(D.check_consistency());

  const std::vector<double> K_after = D.curvature();
  for (int v = 0; v < nv_before; v++)
    EXPECT_NEAR(K_after[v], K_before[v], 1e-9);
}
