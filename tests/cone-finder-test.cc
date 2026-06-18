// DeltahedronView<double>::find_cones (deltahedron-cones.cc): locate the curvature
// cones of a closed triangulation.
//
// Fixture: a GC(k,0)-subdivided icosahedron. With flat-face (harmonic) subdivision
// the Gaussian curvature concentrates exactly at the 12 original degree-5 vertices
// (angle defect pi/3 each, summing to 4*pi by Gauss-Bonnet); the degree-6 vertices
// are flat (K=0). Taubin smoothing spreads each cone into a smooth bump, giving the
// detector 12 unique peaks. find_cones must recover all 12, each integrating ~pi/3.

#include <gtest/gtest.h>
#include <cmath>
#include <limits>
#include <vector>
#include "fullerenes/deltahedron.hh"
#include "fullerenes/polyhedron.hh"

using namespace std;

// Icosahedron = dual of the C20 dodecahedron (12 vertices, 20 triangular faces).
static Deltahedron make_icosahedron() {
  Polyhedron P20 = Polyhedron::C20();
  return Deltahedron(P20.dual());
}

// The degree-5 (cone) vertices, which a GC(k,0) icosahedron has exactly 12 of.
static vector<int> cone_vertices(const Deltahedron& D) {
  vector<int> deg5;
  for (int v = 0; v < D.N; v++) if (D.degree(v) == 5) deg5.push_back(v);
  return deg5;
}

// Minimum spacing between a given vertex set -- the scale that bounds the disk
// radius so each cone's R-disk holds exactly one cone.
static double min_spacing(const Deltahedron& D, const vector<int>& vs) {
  double m = numeric_limits<double>::infinity();
  for (size_t i = 0; i < vs.size(); i++)
    for (size_t j = i+1; j < vs.size(); j++)
      m = std::min(m, (D.points[vs[i]] - D.points[vs[j]]).norm());
  return m;
}

TEST(ConeFinder, GCIcosahedronTwelveCones) {
  Deltahedron geo = make_icosahedron().GCtransform(5, 0);   // 252 vertices, 12 deg-5 cones
  const vector<int> cones5 = cone_vertices(geo);
  ASSERT_EQ(cones5.size(), 12u) << "a GC(5,0) icosahedron must have exactly 12 degree-5 vertices";
  const double L = min_spacing(geo, cones5);

  ConeFinderParams params(0.5 * L);            // disks disjoint (adjacent cones L apart)

  std::vector<SurfaceCone> cones = geo.find_cones(params);

  ASSERT_EQ(cones.size(), 12u) << "12 pentagon cones expected";

  double sum = 0;
  for (const auto& c : cones) {
    sum += c.integrated_K;
    EXPECT_NEAR(c.integrated_K, M_PI/3.0, 0.10) << "each cone integrates ~pi/3";
    EXPECT_GE(c.face, 0);
    EXPECT_LT(c.face, geo.n_triangles());
    EXPECT_NEAR(c.bary[0] + c.bary[1] + c.bary[2], 1.0, 1e-9);
  }
  // Gauss-Bonnet: the 12 cone integrals recover the total curvature 4*pi.
  EXPECT_NEAR(sum, 4.0 * M_PI, 0.15);
}

// max_centres bounds the result (and exercises the early-break / under-segmentation path).
TEST(ConeFinder, MaxCentresBound) {
  Deltahedron geo = make_icosahedron().GCtransform(5, 0);
  const double L = min_spacing(geo, cone_vertices(geo));
  ConeFinderParams params(0.5 * L);
  params.max_centres = 5;
  std::vector<SurfaceCone> cones = geo.find_cones(params);
  EXPECT_LE(cones.size(), 5u);
  EXPECT_EQ(cones.size(), 5u) << "12 cones exist, so the cap of 5 should be reached";
}

// The optional smoothed_points out-param receives the full Taubin-smoothed copy.
TEST(ConeFinder, SmoothedPointsOut) {
  Deltahedron geo = make_icosahedron().GCtransform(5, 0);
  const double L = min_spacing(geo, cone_vertices(geo));
  std::vector<coord3d> smoothed;
  geo.find_cones(ConeFinderParams(0.5 * L), &smoothed);
  ASSERT_EQ((int)smoothed.size(), geo.N);
  // Smoothing must have moved at least one vertex off the original geometry.
  double max_move = 0;
  for (int v = 0; v < geo.N; v++) max_move = std::max(max_move, (smoothed[v] - geo.points[v]).norm());
  EXPECT_GT(max_move, 0.0);
}

// Precondition violations throw std::invalid_argument: bad radius / fraction, and
// a non-finite vertex coordinate.
TEST(ConeFinder, RejectsBadParams) {
  Deltahedron geo = make_icosahedron().GCtransform(3, 0);
  EXPECT_THROW(geo.find_cones(ConeFinderParams(0.0)), std::invalid_argument);
  EXPECT_THROW(geo.find_cones(ConeFinderParams(-1.0)), std::invalid_argument);
  EXPECT_THROW(geo.find_cones(ConeFinderParams(std::numeric_limits<double>::infinity())),
               std::invalid_argument);
  ConeFinderParams bad(1.0); bad.threshold_frac = 1.5;
  EXPECT_THROW(geo.find_cones(bad), std::invalid_argument);

  // A NaN vertex coordinate is rejected, not silently propagated.
  std::vector<coord3d> pts(geo.points.begin(), geo.points.end());
  pts[0] = coord3d(std::numeric_limits<double>::quiet_NaN(), 0, 0);
  geo.set_points(pts);
  EXPECT_THROW(geo.find_cones(ConeFinderParams(1.0)), std::invalid_argument);
}
