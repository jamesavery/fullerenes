// AlexandrovIDTCubic closed-form regression: the cubic polyhedral metric's
// unique convex realization is known exactly for the two maximally
// symmetric fullerenes,
//
//   C20     -> regular dodecahedron   (20 cones, kappa = 3*pi/15 = pi/5)
//   C60-Ih  -> truncated icosahedron  (60 cones, kappa = pi/15)
//
// and both must be reproduced at machine precision: cone census,
// circumradius, enclosed volume, unit polytope edges, and Tbar(0) equal to
// the cubic graph's face lattice.

#include "fullerenes/delaunay_alexandrov.hh"
#include "fullerenes/buckygen-wrapper.hh"

#include <gtest/gtest.h>
#include <cmath>
#include <vector>

using namespace std;

namespace {

Triangulation nth_dual(int N, bool IPR, int idx)
{
  auto Q = BuckyGen::start(N, IPR, false);
  Graph G;
  Triangulation T;
  bool found = false;
  for (int i = 0; BuckyGen::next_fullerene(Q, G); i++)
    if (i == idx) { T = Triangulation(G); found = true; break; }
  BuckyGen::stop(Q);
  EXPECT_TRUE(found) << "isomer C" << N << " #" << idx << " not found";
  return T;
}

// Signed volume of the triangulated closed surface (CCW-outward), centered.
double volume(const DelaunayTriangulation& D, const vector<coord3d>& pos)
{
  coord3d c(0, 0, 0);
  for (const auto& p : pos) c += p;
  c *= 1.0 / pos.size();
  double V = 0;
  for (int f = 0; f < D.nf; f++) {
    if (D.f_he[f] < 0) continue;
    int h0 = D.f_he[f], h1 = D.he_next[h0], h2 = D.he_next[h1];
    coord3d a = pos[D.he_origin[h0]] - c;
    coord3d b = pos[D.he_origin[h1]] - c;
    coord3d d = pos[D.he_origin[h2]] - c;
    V += a.dot(b.cross(d)) / 6.0;
  }
  return V;
}

void expect_closed_form(int N, bool IPR, int n_cones, int npent,
                        double R_exact, double vol_exact)
{
  Triangulation T = nth_dual(N, IPR, 0);

  AlexandrovIDTCubic AC;
  auto P = AC.solve_polytope(T);

  ASSERT_TRUE(AC.solver.valid())
      << AlexandrovSolver::status_str(AC.solver.stats_status);
  EXPECT_EQ(AC.solver.D.nv, n_cones);
  for (int k : AC.cone_npent) EXPECT_EQ(k, npent);
  EXPECT_LT(AC.solver.stats_final_kappa, 1e-10);

  // Circumradius: all cones equidistant from the centroid.
  coord3d c(0, 0, 0);
  for (const auto& p : P.positions) c += p;
  c *= 1.0 / P.positions.size();
  for (const auto& p : P.positions)
    EXPECT_NEAR((p - c).norm(), R_exact, 1e-6);

  // Enclosed volume.
  EXPECT_NEAR(volume(AC.solver.D, P.positions), vol_exact, 1e-6 * vol_exact);

  // Tbar(0) is exactly the cubic graph's face lattice, with unit edges.
  EXPECT_EQ(P.tesselation.n_cells(), T.N);
  auto census = AC.flat_face_census(T, P.tesselation);
  EXPECT_TRUE(census.face_lattice());
  EXPECT_EQ(census.pent_flat, 12);
  EXPECT_EQ(census.n_hex, T.N - 12);
  for (const auto& cell : P.tesselation.cells) {
    const int n = cell.size();
    for (int i = 0; i < n; i++)
      EXPECT_NEAR((P.positions[cell[i].first] -
                   P.positions[cell[(i + 1) % n].first]).norm(), 1.0, 1e-6);
  }
}

}  // namespace

TEST(AlexandrovCubic, C20IsRegularDodecahedron)
{
  const double S5 = sqrt(5.0);
  expect_closed_form(20, false, 20, 3,
                     sqrt(3.0) * (1 + S5) / 4,   // circumradius, unit edge
                     (15 + 7 * S5) / 4);         // volume, unit edge
}

TEST(AlexandrovCubic, C60IhIsTruncatedIcosahedron)
{
  const double S5 = sqrt(5.0);
  expect_closed_form(60, true, 60, 1,
                     sqrt(58 + 18 * S5) / 4,
                     (125 + 43 * S5) / 4);
}

// The Gauss-Newton trust-region fallback -- the solver path through the
// shared JtJ product (dense_linalg_view.hh matmul) -- fires only when the
// pure Newton step is rejected, which C20 and C60-Ih never do: measured
// coverage puts the smallest provoking isomer at C40 buckygen idx 18
// (36/1,625 solves across C20-C50 reach the branch).  This leg gives the
// JtJ integration path a parent-side gate.
TEST(AlexandrovCubic, C40Idx18ReachesGaussNewtonAndConverges)
{
  Triangulation T = nth_dual(40, false, 18);
  {   // dual-metric path
    AlexandrovSolver s;
    s.D = DelaunayTriangulation::compute(T);
    const vector<coord3d> pos = s.solve();
    EXPECT_TRUE(s.valid());
    EXPECT_LT(s.stats_final_kappa, 1e-10);
    EXPECT_FALSE(pos.empty());
  }
  {   // cubic/kis path
    AlexandrovIDTCubic c;
    const vector<coord3d> pos = c.solve(T);
    EXPECT_TRUE(c.solver.valid());
    EXPECT_LT(c.solver.stats_final_kappa, 1e-10);
    EXPECT_FALSE(pos.empty());
  }
}
