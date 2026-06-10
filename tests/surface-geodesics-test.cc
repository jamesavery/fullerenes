// Tests for the surface_distances geodesics extension.
//
// Verifies, on a fullerene dual:
//  - Every simple_geodesic (U, V) is reproducible: end_of_the_line of
//    the recorded (axis, a, b) starting at U reaches V.
//  - Every (a, b) is in the first quadrant (a, b >= 0) per the spec.
//  - The sum of segment lengths matches sqrt(surface_distance_squared)
//    returned by surface_distances.
//  - Each segment's endpoint chains correctly into the next segment's
//    start.

#include <gtest/gtest.h>
#include <cmath>

#include "fullerenes/buckygen-wrapper.hh"
#include "fullerenes/triangulation.hh"

namespace {

Triangulation load_first(int N, bool IPR=false) {
  BuckyGen::buckygen_queue Q = BuckyGen::start(N, IPR, false);
  Graph G; bool ok = BuckyGen::next_fullerene(Q, G); BuckyGen::stop(Q);
  EXPECT_TRUE(ok);
  return Triangulation(G);
}

void verify_geodesics(const Triangulation& T) {
  const int n = T.N;

  matrix<Triangulation::geodesic> Gs(0, 0, Triangulation::geodesic());
  matrix<double> Dsq = T.surface_distances({}, false, &Gs);

  for (int U = 0; U < n; U++) {
    for (int V = 0; V < n; V++) {
      const auto& segs = Gs(U, V).segments;

      int cur = U;
      double len = 0;
      for (size_t k = 0; k < segs.size(); k++) {
        const auto& s = segs[k];
        EXPECT_GE(s.axis, 0)            << "U=" << U << " V=" << V << " k=" << k;
        EXPECT_LT(s.axis, T.degree(cur)) << "U=" << U << " V=" << V << " k=" << k;
        EXPECT_GE(s.g.first,  0)        << "U=" << U << " V=" << V << " k=" << k;
        EXPECT_GE(s.g.second, 0)        << "U=" << U << " V=" << V << " k=" << k;

        int v_end = T.end_of_the_line(cur, s.axis, s.g.first, s.g.second);
        int expected = (k + 1 == segs.size()) ? V : -1;
        if (expected != -1) {
          EXPECT_EQ(v_end, expected)
            << "U=" << U << " V=" << V << " k=" << k
            << " axis=" << s.axis << " (a,b)=(" << s.g.first << "," << s.g.second << ")";
        }
        len += std::sqrt((double)s.g.norm2());
        cur = v_end;
      }

      double got_dsq = len * len;
      double tol = 1e-9 * std::max(1.0, Dsq(U, V));
      EXPECT_NEAR(Dsq(U, V), got_dsq, tol)
        << "U=" << U << " V=" << V << " Dsq=" << Dsq(U, V)
        << " segments=" << segs.size();
    }
  }
}

}  // namespace

TEST(SurfaceGeodesics, C20) {
  verify_geodesics(load_first(20));
}

TEST(SurfaceGeodesics, C28) {
  verify_geodesics(load_first(28));
}

TEST(SurfaceGeodesics, C60) {
  verify_geodesics(load_first(60));
}
