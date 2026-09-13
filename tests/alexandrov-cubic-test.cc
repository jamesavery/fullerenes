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

enum class Regime { Banded, Exact };

// Build through one regime (class banner, step 2), then solve.
AlexandrovSolver::AlexandrovPolytope build_and_solve(AlexandrovIDTCubic& AC,
                                                     const Triangulation& T,
                                                     Regime regime)
{
  if (regime == Regime::Exact) AC.build(T);
  else                         AC.build_banded(T);
  return AC.solver.solve_polytope();
}

vector<int> sorted(vector<int> v) { sort(v.begin(), v.end()); return v; }

void expect_closed_form(int N, bool IPR, int n_cones, int npent,
                        double R_exact, double vol_exact,
                        Regime regime = Regime::Exact)
{
  Triangulation T = nth_dual(N, IPR, 0);

  AlexandrovIDTCubic AC;
  auto P = build_and_solve(AC, T, regime);

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
                     (15 + 7 * S5) / 4,          // volume, unit edge
                     Regime::Banded);
}

TEST(AlexandrovCubic, C60IhIsTruncatedIcosahedron)
{
  const double S5 = sqrt(5.0);
  expect_closed_form(60, true, 60, 1,
                     sqrt(58 + 18 * S5) / 4,
                     (125 + 43 * S5) / 4, Regime::Banded);
}

// The same closed forms through the exact regime (build, the default:
// exact algebraic signs plus the canonical completion).  C20's kis surface
// is twelve cocircular pentagon cells, so this case also exercises the
// completion: every pentagon must be fanned (flips > 0) and the polytope
// unchanged.
TEST(AlexandrovCubicExact, C20IsRegularDodecahedron)
{
  const double S5 = sqrt(5.0);
  expect_closed_form(20, false, 20, 3,
                     sqrt(3.0) * (1 + S5) / 4,
                     (15 + 7 * S5) / 4);
}

TEST(AlexandrovCubicExact, C60IhIsTruncatedIcosahedron)
{
  const double S5 = sqrt(5.0);
  expect_closed_form(60, true, 60, 1,
                     sqrt(58 + 18 * S5) / 4,
                     (125 + 43 * S5) / 4);
}

// Both regimes keep the same cones (the reduced triangulations are
// compared corpus-wide by tools/bench_cubic_regimes); the exact build
// reports its completion, which refuses no cell and, these kis surfaces
// carrying cocircular pentagon cells, always fans.
TEST(AlexandrovCubicExact, RegimesAgreeOnConesAndCompletionFans)
{
  for (int N : {20, 36, 40}) {
    Triangulation T = nth_dual(N, false, 0);
    AlexandrovIDTCubic banded, exact;
    banded.build_banded(T);
    const auto completion = exact.build(T);
    EXPECT_EQ(banded.solver.D.nv, exact.solver.D.nv) << "C" << N;
    EXPECT_EQ(sorted(banded.cone_kis_vertex), sorted(exact.cone_kis_vertex)) << "C" << N;
    EXPECT_EQ(completion.ambiguous, 0) << "C" << N;
    EXPECT_EQ(completion.nondisk, 0) << "C" << N;
    EXPECT_GT(completion.flips, 0) << "C" << N;
  }
}

// Point tracking through the exact build: every removed kis vertex (face
// centres and hexagon-only cubic vertices) rides the removal and the
// completion's flips as a tracked point, and ends in a live face with
// normalized, non-negative barycentric coordinates, its label intact.
TEST(AlexandrovCubicExact, TrackedRemovalSurvivesCompletion)
{
  for (int N : {36, 60}) {
    Triangulation T = nth_dual(N, N == 60, 0);
    AlexandrovIDTCubic exact;
    exact.track_removed = true;
    const auto completion = exact.build(T);
    EXPECT_GT(completion.flips, 0) << "C" << N;
    const DelaunayTriangulation& D = exact.solver.D;
    ASSERT_TRUE(D.tracker.active) << "C" << N;
    const int kis_nv = T.N + 2 * T.N - 4;
    EXPECT_EQ(D.tracker.view.n, kis_nv - D.nv) << "C" << N;
    vector<int> labels, expected;
    for (int i = 0; i < D.tracker.view.n; i++) {
      labels.push_back(D.tracker.view.label[i]);
      const int f = D.tracker.view.face[i];
      EXPECT_TRUE(f >= 0 && f < D.nf && D.f_he[f] >= 0) << "C" << N << " point " << i;
      const double* b = D.tracker.view.b_of(i);
      EXPECT_NEAR(b[0] + b[1] + b[2], 1.0, 1e-9) << "C" << N << " point " << i;
      EXPECT_TRUE(b[0] >= 0 && b[1] >= 0 && b[2] >= 0) << "C" << N << " point " << i;
    }
    const vector<int> cones = sorted(exact.cone_kis_vertex);
    for (int v = 0; v < kis_nv; v++)
      if (!binary_search(cones.begin(), cones.end(), v)) expected.push_back(v);
    EXPECT_EQ(sorted(labels), expected) << "C" << N;
  }
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
  {   // cubic/kis path (exact, the default)
    AlexandrovIDTCubic c;
    const vector<coord3d> pos = c.solve(T);
    EXPECT_TRUE(c.solver.valid());
    EXPECT_LT(c.solver.stats_final_kappa, 1e-10);
    EXPECT_FALSE(pos.empty());
  }
  {   // cubic/kis path, tolerance-based regime
    AlexandrovIDTCubic c;
    c.build_banded(T);
    const vector<coord3d> pos = c.solver.solve();
    EXPECT_TRUE(c.solver.valid());
    EXPECT_LT(c.solver.stats_final_kappa, 1e-10);
    EXPECT_FALSE(pos.empty());
  }
}
