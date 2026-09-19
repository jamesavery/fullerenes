// dwu-forcefield-test -- the derived-geometry Wu force field (fullerenes/dwu_forcefield.hh).
//
// The claim this field rests on is that a fullerene's rest geometry is DERIVED,
// not tabulated: two bond lengths fix the corner angles and the vertex dihedrals
// by exact construction.  So the tests are mostly consistency statements about
// that construction rather than regressions against stored numbers:
//
//   cyclic polygon   equal sides give the regular polygon exactly; general sides
//                    give a polygon that actually closes on a single circle
//   derived values   equal bonds reproduce the published extwu table (and expose
//                    its one typo)
//   flatness         zero on a planar face, (m-1)d^2/m for one atom d off-plane
//   gradient         central differences, all three term groups together
//   relaxation       a real cage from the library's own paint start lands on a
//                    geometry with fullerene bond lengths and a small gradient

#include <gtest/gtest.h>

#include "fullerenes/dwu_forcefield.hh"
#include "fullerenes/buckygen-wrapper.hh"
#include "fullerenes/triangulation.hh"

#include <cmath>
#include <numeric>
#include <vector>

namespace {

constexpr double DEG = 180.0 / M_PI;

// The first fullerene of size N (buckygen; needs no isomer database).
FullereneGraph first_cage(int N) {
  BuckyGen::buckygen_queue Q = BuckyGen::start(N, false, false);
  Triangulation T;
  BuckyGen::next_fullerene(Q, T);
  BuckyGen::stop(Q);
  return FullereneGraph(T.dual_graph());
}

// Every bond class one length, so every derived quantity must collapse to the
// regular-polygon value the published table holds.
dwu::Parameters uniform_bonds(double r) {
  dwu::Parameters P = dwu::Parameters::fitted();
  P.R55 = P.R56 = P.R66 = dwu::Field::constant(r);
  return P;
}

}  // namespace

// --- the cyclic polygon ------------------------------------------------------

TEST(DwuCyclicPolygon, EqualSidesGiveTheRegularPolygon) {
  for (int m : {5, 6}) {
    const std::vector<double> s(m, 1.43);
    const std::vector<double> a = dwu::cyclic_polygon_angles(s);
    ASSERT_EQ(int(a.size()), m);
    for (double x : a) EXPECT_NEAR(x, 180.0 * (m - 2) / m, 1e-9) << "m = " << m;
  }
  // Scale invariance: the angles are shape, not size.
  EXPECT_NEAR(dwu::cyclic_polygon_angles({1.0, 1.0, 1.0, 1.0, 1.0})[0],
              dwu::cyclic_polygon_angles({7.0, 7.0, 7.0, 7.0, 7.0})[0], 1e-9);
}

// The construction's real content: walking the returned angles must trace a
// closed polygon whose vertices lie on ONE circle.  This checks the answer
// against its own definition rather than against a stored number.
TEST(DwuCyclicPolygon, GeneralSidesClosOnASingleCircle) {
  const std::vector<std::vector<double>> cases = {
      {1.39, 1.45, 1.39, 1.45, 1.39, 1.45},   // a hexagon between pentagons
      {1.47, 1.43, 1.43, 1.43, 1.43},         // a pentagon with one long edge
      {1.30, 1.55, 1.40, 1.48, 1.35, 1.44},   // no symmetry at all
  };
  for (const std::vector<double>& s : cases) {
    const int m = int(s.size());
    const std::vector<double> a = dwu::cyclic_polygon_angles(s);

    EXPECT_NEAR(std::accumulate(a.begin(), a.end(), 0.0), 180.0 * (m - 2), 1e-7);

    // Walk side i, then turn by the exterior angle at the vertex it ends on.
    double x = 0, y = 0, heading = 0;
    std::vector<double> px(m), py(m);
    for (int i = 0; i < m; ++i) {
      px[i] = x; py[i] = y;
      x += s[i] * std::cos(heading);
      y += s[i] * std::sin(heading);
      heading += M_PI - a[(i + 1) % m] / DEG;      // exterior angle at the far end
    }
    EXPECT_NEAR(std::hypot(x - px[0], y - py[0]), 0.0, 1e-7) << "polygon does not close";

    // Circumcentre from the first three vertices; every vertex must share it.
    const double ax = px[0], ay = py[0], bx = px[1], by = py[1], cx = px[2], cy = py[2];
    const double d = 2 * (ax * (by - cy) + bx * (cy - ay) + cx * (ay - by));
    ASSERT_GT(std::fabs(d), 1e-12);
    const double ux = ((ax*ax + ay*ay) * (by - cy) + (bx*bx + by*by) * (cy - ay)
                       + (cx*cx + cy*cy) * (ay - by)) / d;
    const double uy = ((ax*ax + ay*ay) * (cx - bx) + (bx*bx + by*by) * (ax - cx)
                       + (cx*cx + cy*cy) * (bx - ax)) / d;
    const double R = std::hypot(ax - ux, ay - uy);
    for (int i = 0; i < m; ++i)
      EXPECT_NEAR(std::hypot(px[i] - ux, py[i] - uy), R, 1e-7) << "vertex " << i << " off the circle";
  }
}

TEST(DwuCyclicPolygon, RejectsSidesThatFormNoPolygon) {
  // One side longer than the rest combined: no polygon, cyclic or otherwise.
  EXPECT_THROW(dwu::cyclic_polygon_angles({10.0, 1.0, 1.0, 1.0, 1.0}), std::domain_error);
}

// --- derived rest values -----------------------------------------------------

// The published extwu dihedral table is not an independent measurement: it is
// this construction, evaluated at extwu's own bond lengths.  Three of its four
// entries come back to the precision they are printed at; the fourth does not,
// and the way it misses says what happened to it.
TEST(DwuDerivedRestValues, ReproducesThePublishedTableAndExposesItsTypo) {
  const wu::Parameters e = wu::Parameters::extwu();
  const double R55 = e.R55, R56 = e.R56, R66 = e.R66;

  // Each vertex class takes the edges it actually has: a 555 vertex is bounded by
  // three 55 bonds, a 556 vertex by 56, 55, 56, a 566 vertex by 56, 66, 56.
  EXPECT_NEAR(coord3d::ideal_dihedral(5, 5, 5, R55, R55, R55) * DEG, e.D555, 5e-3);
  EXPECT_NEAR(coord3d::ideal_dihedral(5, 5, 6, R56, R55, R56) * DEG, e.D556, 5e-3);
  EXPECT_EQ  (coord3d::ideal_dihedral(6, 6, 6, R66, R66, R66),       e.D666);

  // 23.9430 against a printed 23.49: the published constant has the digits of its
  // own construction transposed.
  const double D566 = coord3d::ideal_dihedral(6, 6, 5, R56, R66, R56) * DEG;
  EXPECT_NEAR(D566, 23.9430, 1e-4);
  EXPECT_GT(std::fabs(D566 - e.D566), 0.4);
}

TEST(DwuDerivedRestValues, EqualBondsGiveRegularPolygonCorners) {
  const FullereneGraph G = first_cage(60);
  const dwu::ForceField FF = dwu::forcefield(G, uniform_bonds(1.43));

  int n5 = 0, n6 = 0;
  for (const wu::Corner& c : FF.wu.corners) {
    const double deg = c.q0 * DEG;
    if (std::fabs(deg - 108.0) < 1e-9) ++n5;
    else if (std::fabs(deg - 120.0) < 1e-9) ++n6;
    else ADD_FAILURE() << "corner rest angle " << deg << " is neither 108 nor 120";
  }
  EXPECT_EQ(n5, 60);                 // 12 pentagons x 5 corners, on any C60 isomer
  EXPECT_EQ(n6, 120);                // 20 hexagons  x 6 corners
  EXPECT_EQ(n5 + n6, 3 * G.N);
}

// Classification: every vertex must get the dihedral of the faces it actually
// carries.  Equal bond lengths make the answer depend on the face multiset alone,
// so this isolates the classification from the slot convention tested below.
TEST(DwuDerivedRestValues, EveryVertexGetsItsOwnClass) {
  const double r = 1.43;
  const double D[4] = {0.0,                                          // 666
                       coord3d::ideal_dihedral(6, 6, 5, r, r, r),    // 566
                       coord3d::ideal_dihedral(5, 5, 6, r, r, r),    // 556
                       coord3d::ideal_dihedral(5, 5, 5, r, r, r)};   // 555
  for (int N : {60, 80, 180}) {
    const FullereneGraph G = first_cage(N);
    const dwu::ForceField FF = dwu::forcefield(G, uniform_bonds(r));
    ASSERT_EQ(int(FF.wu.dihedrals.size()), G.N);
    int seen[4] = {0, 0, 0, 0};
    for (const wu::Dihedral& d : FF.wu.dihedrals) {
      const node_t u = d.atoms[0];
      int npent = 0;
      for (node_t v : G.nbrs(u)) npent += (G.face_size(u, v) == 5);
      ASSERT_GE(npent, 0); ASSERT_LE(npent, 3);
      EXPECT_NEAR(d.q0, D[npent], 1e-12) << "N = " << N << ", vertex " << u
                                         << " has " << npent << " pentagons";
      ++seen[npent];
    }
    // Every pentagon corner is a vertex with at least one pentagon: 12 x 5 = 60.
    EXPECT_EQ(seen[1] + 2 * seen[2] + 3 * seen[3], 60) << "N = " << N;
  }
}

// The slot convention.  ideal_dihedral's face A sits between edges ur and us, and
// the quadruple it evaluates is not symmetric under every rotation of the three:
// posing a 566 vertex as (5,6,6) gives 19.77 deg where the correct (6,6,5) gives
// 23.94.  (Permuting only the edge LENGTHS, faces fixed, moves it by 0.04 deg --
// so this is the posing of the faces, not the order of the lengths.)  Equal bonds
// cannot see the mistake; distinct ones can.
TEST(DwuDerivedRestValues, OddFaceGoesIntoSlotC) {
  dwu::Parameters P = dwu::Parameters::fitted();
  const double a = 1.44, b = 1.38, c = 1.49;      // R56, R66, R55 -- deliberately apart
  P.R56 = dwu::Field::constant(a);
  P.R66 = dwu::Field::constant(b);
  P.R55 = dwu::Field::constant(c);

  // The test has teeth only if the wrong posing gives a different answer.
  ASSERT_GT(std::fabs(coord3d::ideal_dihedral(6, 6, 5, a, b, a)
                    - coord3d::ideal_dihedral(5, 6, 6, a, a, b)) * DEG, 1.0);

  const double D[4] = {0.0,
                       coord3d::ideal_dihedral(6, 6, 5, a, b, a),    // 566: 56, 66, 56
                       coord3d::ideal_dihedral(5, 5, 6, a, c, a),    // 556: 56, 55, 56
                       coord3d::ideal_dihedral(5, 5, 5, c, c, c)};   // 555: 55, 55, 55
  for (int N : {60, 180}) {
    const FullereneGraph G = first_cage(N);
    const dwu::ForceField FF = dwu::forcefield(G, P);
    for (const wu::Dihedral& d : FF.wu.dihedrals) {
      const node_t u = d.atoms[0];
      int npent = 0;
      for (node_t v : G.nbrs(u)) npent += (G.face_size(u, v) == 5);
      EXPECT_NEAR(d.q0, D[npent], 1e-12) << "N = " << N << ", vertex " << u;
    }
  }
}

TEST(DwuDerivedRestValues, BondLengthsFollowTheirClassAndShell) {
  const FullereneGraph G = first_cage(180);
  dwu::Parameters P = dwu::Parameters::fitted();
  const dwu::ForceField FF = dwu::forcefield(G, P);

  const std::vector<face_t> faces = G.compute_faces_oriented(6);
  const std::vector<int> pd = dwu::pentagon_distance(G, faces);
  ASSERT_EQ(int(pd.size()), G.N);
  EXPECT_EQ(*std::min_element(pd.begin(), pd.end()), 0);
  EXPECT_GT(*std::max_element(pd.begin(), pd.end()), 0) << "C180 has vertices off the pentagons";

  for (const wu::Bond& b : FF.wu.bonds) {
    const node_t u = b.atoms[0], v = b.atoms[1];
    const int np = (G.face_size(u, v) == 5) + (G.face_size(v, u) == 5);
    const dwu::Field& f = (np == 2 ? P.R55 : np == 1 ? P.R56 : P.R66);
    EXPECT_DOUBLE_EQ(b.q0, f.at(pd[u] + pd[v]));
    EXPECT_GT(b.q0, 1.30);
    EXPECT_LT(b.q0, 1.55);
  }
}

// The disclination field's shape, asserted against the field itself rather than
// against transcribed numbers: the constants are provisional and will move when the
// full-corpus fit lands, but the shape is the model and must not.
TEST(DwuField, FreeHeadThenMonotoneDecayToTheAsymptote) {
  const dwu::Field R66 = dwu::Parameters::fitted().R66;
  ASSERT_EQ(R66.shell.size(), R66.value.size());
  ASSERT_FALSE(R66.shell.empty()) << "R66 is the field with free head values";

  // The listed shells return their own value, exactly and independently of the tail.
  for (std::size_t i = 0; i < R66.shell.size(); ++i)
    EXPECT_DOUBLE_EQ(R66.at(R66.shell[i]), R66.value[i]);

  // Past the head the tail takes over and decays monotonically to q_inf.
  for (int p = R66.p0; p < 60; ++p)
    EXPECT_GT(R66.at(p), R66.at(p + 1)) << "the tail must decay, p = " << p;
  EXPECT_NEAR(R66.at(400), R66.q_inf, 1e-12);

  // The physics the field exists to carry: the 66 bond is SHORTEST against a
  // pentagon and relaxes outward, so the head must sit below the bulk value.
  EXPECT_LT(R66.at(0), R66.at(R66.p0)) << "no disclination strain to describe";
  EXPECT_GT(R66.delta, 0.0) << "the tail must approach q_inf from above";

  // A constant field ignores p entirely.
  EXPECT_DOUBLE_EQ(dwu::Field::constant(1.42).at(0), 1.42);
  EXPECT_DOUBLE_EQ(dwu::Field::constant(1.42).at(17), 1.42);
}

// --- the flatness term -------------------------------------------------------

TEST(DwuFlatness, VanishesOnAPlanarFaceAndGrowsAsTheSquareOffIt) {
  for (int m : {4, 5, 6}) {
    dwu::FlatnessField F;
    face_t f;
    for (int i = 0; i < m; ++i) f.push_back(i);
    F.faces = {f};
    F.k = {1.0};

    std::vector<coord3d> x(m), g(m);
    for (int i = 0; i < m; ++i) {
      const double t = 2 * M_PI * i / m;
      x[i] = coord3d(std::cos(t), std::sin(t), 0.0);
    }
    EXPECT_NEAR(F.energy_gradient(x, g), 0.0, 1e-18) << "planar face, m = " << m;
    for (const coord3d& v : g) EXPECT_NEAR(v.norm(), 0.0, 1e-9);

    // One atom lifted by d.  The least-squares plane does not stay horizontal: it
    // tilts toward the lifted atom, and the tilt is first order in d, so it changes
    // the leading coefficient.  Fitting z = p + q x + r y over a regular m-gon of
    // circumradius R (sum x = sum y = 0, sum x^2 = sum y^2 = m R^2 / 2) gives
    // p = d/m and q = 2d/(mR), whence
    //     lambda_min = d^2 - d^2/m - 2 d^2/m = d^2 (m - 3)/m,
    // independent of R -- and exactly zero at m = 3, as three points must give.
    const double d = 1e-4;
    x[0] = coord3d(x[0][0], x[0][1], d);
    for (coord3d& v : g) v = coord3d(0, 0, 0);
    const double E = F.energy_gradient(x, g);
    const double want = 0.5 * d * d * (m - 3) / m;
    EXPECT_NEAR(E, want, 1e-3 * want) << "m = " << m;
  }
}

// --- the gradient ------------------------------------------------------------

TEST(DwuGradient, MatchesCentralDifferences) {
  for (int N : {60, 80, 180}) {
    const FullereneGraph G = first_cage(N);
    const dwu::ForceField FF = dwu::forcefield(G, dwu::Parameters::fitted());

    // A geometry off the minimum, so every term is active: the paint start with
    // a deterministic, reproducible perturbation.
    std::vector<coord3d> x = G.zero_order_geometry(1.42);
    ASSERT_EQ(int(x.size()), G.N);
    for (int i = 0; i < G.N; ++i)
      for (int j = 0; j < 3; ++j)
        x[i][j] += 0.02 * std::sin(3.1 * i + 1.7 * j);

    std::vector<coord3d> g(G.N);
    FF.energy_gradient(x, {g.data(), g.size()});

    // Every 37th coordinate, so the sample spreads over the whole cage.
    const double h = 1e-6;
    int checked = 0;
    for (int c = 0; c < 3 * G.N; c += 37) {
      const int i = c / 3, j = c % 3;
      const double x0 = x[i][j];
      x[i][j] = x0 + h; const double Ep = FF.energy(x);
      x[i][j] = x0 - h; const double Em = FF.energy(x);
      x[i][j] = x0;
      const double fd = (Ep - Em) / (2 * h);
      EXPECT_NEAR(g[i][j], fd, 1e-4 * std::max(1.0, std::fabs(fd)))
          << "N = " << N << ", atom " << i << ", axis " << j;
      ++checked;
    }
    EXPECT_GT(checked, 4);
  }
}

// --- end to end --------------------------------------------------------------

TEST(DwuOptimize, RelaxesThePaintStartToAFullereneGeometry) {
  for (int N : {60, 80, 180}) {
    const FullereneGraph G = first_cage(N);
    const dwu::ForceField FF = dwu::forcefield(G, dwu::Parameters::fitted());

    std::vector<coord3d> x = G.zero_order_geometry(1.42);
    const double E0 = FF.energy(x);

    const minimize::Outcome out = dwu::optimize(FF, {x.data(), x.size()}, 1e-8);
    EXPECT_TRUE(out.converged) << "N = " << N << ": hit the iteration safeguard";
    EXPECT_LT(out.f, E0) << "N = " << N << ": relaxation did not lower the energy";

    // Bond lengths land in the range real fullerene bonds occupy.  This is the
    // statement that matters: the derived rest values are self-consistent, so
    // the minimiser is not fighting contradictory targets.
    double lo = 1e9, hi = 0;
    for (const wu::Bond& b : FF.wu.bonds) {
      const double r = (x[b.atoms[0]] - x[b.atoms[1]]).norm();
      lo = std::min(lo, r); hi = std::max(hi, r);
    }
    EXPECT_GT(lo, 1.30) << "N = " << N;
    EXPECT_LT(hi, 1.55) << "N = " << N;

    // The cage must stay a cage: every atom on a shell, none collapsed to the centre.
    coord3d c(0, 0, 0);
    for (const coord3d& p : x) c += p;
    c /= double(G.N);
    double rmin = 1e9;
    for (const coord3d& p : x) rmin = std::min(rmin, (p - c).norm());
    EXPECT_GT(rmin, 1.0) << "N = " << N << ": the cage collapsed";
  }
}
