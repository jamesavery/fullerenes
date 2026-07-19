// Robustness tests for convex_hull_tris (geometry.hh, impl in deltahedron.cc).
//
// Focus: exactly-coplanar / on-edge inputs (the midpoint-subdivided-mesh
// degeneracy class that used to send the incremental hull into a non-
// terminating blow-up). All data is deterministic -- a small LCG for the
// pseudo-random cases, closed-form for the structured ones -- so a failure
// always reproduces. The function's own guards abort() on impossible
// intermediate states, so "the test simply completes" is itself the
// termination check (a regression would hang or abort, not silently pass).

#include <gtest/gtest.h>
#include "fullerenes/geometry.hh"

#include <array>
#include <cmath>
#include <map>
#include <set>
#include <span>
#include <vector>

using namespace std;

namespace {

// ---- deterministic PRNG (splitmix64-style LCG) ------------------------------
struct LCG {
  uint64_t s;
  explicit LCG(uint64_t seed) : s(seed) {}
  uint64_t next(){
    s = s * 6364136223846793005ULL + 1442695040888963407ULL;
    uint64_t z = s;
    z = (z ^ (z >> 30)) * 0xbf58476d1ce4e5b9ULL;
    z = (z ^ (z >> 27)) * 0x94d049bb133111ebULL;
    return z ^ (z >> 31);
  }
  double uniform(){ return (next() >> 11) * (1.0 / 9007199254740992.0); }  // [0,1)
};

// A uniform-ish point on a sphere of given radius (deterministic given LCG).
coord3d sphere_point(LCG& g, double radius){
  double z   = 2.0 * g.uniform() - 1.0;
  double phi = 2.0 * M_PI * g.uniform();
  double r   = sqrt(max(0.0, 1.0 - z*z));
  return coord3d(r * cos(phi) * radius, r * sin(phi) * radius, z * radius);
}

// De-duplicating point accumulator: shared subdivision midpoints collapse to a
// single index so the hull vertex set stays crisp.
struct PointSet {
  vector<coord3d> P;
  int add(const coord3d& q){
    for(int i = 0; i < (int)P.size(); i++)
      if((P[i] - q).norm() < 1e-6) return i;   // points are >=~1 apart in these tests
    P.push_back(q);
    return (int)P.size() - 1;
  }
};

double bbox_diag(span<const coord3d> pts){
  coord3d lo = pts[0], hi = pts[0];
  for(size_t i = 1; i < pts.size(); i++)
    for(int k = 0; k < 3; k++){ lo[k] = min(lo[k], pts[i][k]); hi[k] = max(hi[k], pts[i][k]); }
  return (hi - lo).norm();
}

set<int> hull_vertices(const vector<array<int,3>>& tris){
  set<int> v;
  for(const auto& t : tris){ v.insert(t[0]); v.insert(t[1]); v.insert(t[2]); }
  return v;
}

// One round of 1->4 midpoint subdivision. Consumes a triangle list indexing
// into ps, emits the refined list, and registers every edge midpoint in ps.
vector<array<int,3>> subdivide_once(const vector<array<int,3>>& tris, PointSet& ps){
  vector<array<int,3>> out;
  out.reserve(tris.size() * 4);
  for(const auto& t : tris){
    coord3d A = ps.P[t[0]], B = ps.P[t[1]], C = ps.P[t[2]];
    int ab = ps.add((A + B) / 2.0);
    int bc = ps.add((B + C) / 2.0);
    int ca = ps.add((C + A) / 2.0);
    out.push_back({t[0], ab, ca});
    out.push_back({ab, t[1], bc});
    out.push_back({ca, bc, t[2]});
    out.push_back({ab, bc, ca});
  }
  return out;
}

// 12 icosahedron vertices at the requested circumradius (0,+/-1,+/-phi) cyclic.
vector<coord3d> icosahedron(double circumradius){
  const double phi = (1.0 + sqrt(5.0)) / 2.0;
  const double s   = circumradius / sqrt(1.0 + phi*phi);
  vector<coord3d> v;
  for(double a : {-1.0, 1.0})
    for(double b : {-phi, phi}){
      v.push_back(coord3d(0, a*s, b*s));
      v.push_back(coord3d(a*s, b*s, 0));
      v.push_back(coord3d(b*s, 0, a*s));
    }
  return v;   // 12 vertices
}

// ---- shared invariant checks ------------------------------------------------

// Watertight + oriented: every directed edge appears exactly once and its twin
// (reverse) appears exactly once.
void expect_watertight(const vector<array<int,3>>& tris){
  ASSERT_FALSE(tris.empty());
  map<pair<int,int>,int> cnt;
  for(const auto& t : tris)
    for(int j = 0; j < 3; j++)
      cnt[{t[j], t[(j+1)%3]}]++;
  for(const auto& [e, c] : cnt){
    EXPECT_EQ(c, 1) << "directed edge (" << e.first << "," << e.second
                    << ") used " << c << " times";
    auto it = cnt.find({e.second, e.first});
    ASSERT_NE(it, cnt.end()) << "directed edge (" << e.first << "," << e.second
                             << ") has no twin";
    EXPECT_EQ(it->second, 1);
  }
}

// Convexity: every input point lies at most 10*eps outside every (outward-
// oriented) face plane.
void expect_convex(span<const coord3d> pts, const vector<array<int,3>>& tris){
  const double eps = 1e-10 * bbox_diag(pts);
  double worst = -1;
  for(const auto& t : tris){
    coord3d fn = (pts[t[1]] - pts[t[0]]).cross(pts[t[2]] - pts[t[0]]);
    double fnl = fn.norm();
    ASSERT_GT(fnl, 0.0) << "degenerate (zero-area) output face";
    coord3d nhat = fn / fnl;
    for(const auto& q : pts){
      double d = nhat.dot(q - pts[t[0]]);
      worst = max(worst, d);
    }
  }
  EXPECT_LE(worst, 10.0 * eps) << "a point sits " << worst
                               << " outside a face (eps=" << eps << ")";
}

}  // namespace

// -----------------------------------------------------------------------------
// The formerly failing degeneracy class: a mesh whose new points lie exactly in
// parent face/edge planes. The clean icosahedron, refined twice, must hull back
// to exactly its 12 vertices / 20 faces.
TEST(ConvexHullTris, SubdividedIcosahedron){
  vector<coord3d> ico = icosahedron(100.0);
  ASSERT_EQ(ico.size(), 12u);

  vector<array<int,3>> base = convex_hull_tris(ico);
  ASSERT_EQ(base.size(), 20u) << "clean icosahedron hull is not 20 triangles";

  PointSet ps;
  for(const auto& v : ico) ps.add(v);          // indices 0..11 are the corners
  vector<array<int,3>> tri = base;
  tri = subdivide_once(tri, ps);
  tri = subdivide_once(tri, ps);               // exactly-coplanar interior points

  ASSERT_GT(ps.P.size(), 100u);                // ~162 points, most on faces/edges

  vector<array<int,3>> hull = convex_hull_tris(ps.P);

  EXPECT_EQ(hull.size(), 20u) << "subdivided icosahedron hull is not 20 faces";
  set<int> verts = hull_vertices(hull);
  set<int> corners = {0,1,2,3,4,5,6,7,8,9,10,11};
  EXPECT_EQ(verts, corners) << "hull vertex set is not exactly the 12 corners";
  expect_watertight(hull);
  expect_convex(ps.P, hull);
}

// A dense cube lattice: coplanar faces, collinear edges, interior points. The 8
// corners must be on the hull; every hull vertex must lie on the cube surface.
TEST(ConvexHullTris, CubeGrid){
  const int G = 11;                 // 11x11x11
  const double h = 10.0;            // spacing -> side 100
  PointSet ps;
  set<int> corner_idx;
  for(int i = 0; i < G; i++)
    for(int j = 0; j < G; j++)
      for(int k = 0; k < G; k++){
        coord3d q(i*h, j*h, k*h);
        int id = ps.add(q);
        bool xc = (i == 0 || i == G-1), yc = (j == 0 || j == G-1), zc = (k == 0 || k == G-1);
        if(xc && yc && zc) corner_idx.insert(id);
      }
  ASSERT_EQ(corner_idx.size(), 8u);

  vector<array<int,3>> hull = convex_hull_tris(ps.P);
  ASSERT_FALSE(hull.empty());
  expect_watertight(hull);
  expect_convex(ps.P, hull);

  set<int> verts = hull_vertices(hull);
  for(int c : corner_idx)
    EXPECT_TRUE(verts.count(c)) << "corner " << c << " missing from hull";
  // Every hull vertex must lie on the cube surface (min or max in some coord).
  const double lo = 0.0, hi = (G-1)*h, tol = 1e-6;
  for(int v : verts){
    const coord3d& q = ps.P[v];
    bool on_surface =
      fabs(q[0]-lo) < tol || fabs(q[0]-hi) < tol ||
      fabs(q[1]-lo) < tol || fabs(q[1]-hi) < tol ||
      fabs(q[2]-lo) < tol || fabs(q[2]-hi) < tol;
    EXPECT_TRUE(on_surface) << "hull vertex " << v << " is interior to the cube";
  }
}

// A random polytope, refined twice: the subdivision points lie exactly on the
// polytope's planar faces, so the hull must return to (a subset of) the
// original 30 generators.
TEST(ConvexHullTris, SubdividedRandomPolytope){
  LCG g(0x1234abcdULL);
  const int NORIG = 30;
  vector<coord3d> orig;
  for(int i = 0; i < NORIG; i++) orig.push_back(sphere_point(g, 100.0));

  vector<array<int,3>> base = convex_hull_tris(orig);
  ASSERT_FALSE(base.empty());

  PointSet ps;
  for(const auto& v : orig) ps.add(v);         // indices 0..NORIG-1
  vector<array<int,3>> tri = base;
  tri = subdivide_once(tri, ps);
  tri = subdivide_once(tri, ps);

  vector<array<int,3>> hull = convex_hull_tris(ps.P);
  expect_watertight(hull);
  expect_convex(ps.P, hull);

  set<int> verts = hull_vertices(hull);
  for(int v : verts)
    EXPECT_LT(v, NORIG) << "hull vertex " << v << " is a subdivision point";
}

// A generic well-separated point cloud: all 500 sphere points are extreme, so
// the hull is a full triangulated 2-sphere with F = 2V-4.
TEST(ConvexHullTris, GenericSpherePoints){
  LCG g(0xdeadbeefULL);
  const int NP = 500;
  vector<coord3d> pts;
  for(int i = 0; i < NP; i++) pts.push_back(sphere_point(g, 100.0));

  vector<array<int,3>> hull = convex_hull_tris(pts);
  expect_watertight(hull);
  expect_convex(pts, hull);

  set<int> verts = hull_vertices(hull);
  EXPECT_EQ((int)verts.size(), NP) << "not every sphere point is on the hull";
  EXPECT_EQ((int)hull.size(), 2*NP - 4) << "F != 2V-4 for a triangulated sphere";
}

// Complexity regression: the conflict-graph construction does expected
// O(n log n) point-face conflict tests and creates O(n) faces over the
// fixed-seed shuffled insertion order. The all-points-extreme sphere cloud is
// the adversarial input class for the quadratic full-scan implementations this
// replaced (each of which needed ~n^2/2 = 2e8 tests here, an order of
// magnitude over the bound below). The bound is a machine-independent COUNT
// (HullStats), so this test cannot flake on wall-clock and fails loudly if the
// construction ever regresses to a per-insertion full scan.
TEST(ConvexHullTris, ConflictGraphComplexity){
  LCG g(0xfeedfaceULL);
  const int NP = 20000;
  vector<coord3d> pts;
  pts.reserve(NP);
  for(int i = 0; i < NP; i++) pts.push_back(sphere_point(g, 100.0));

  HullStats st;
  vector<array<int,3>> hull = convex_hull_tris(pts, &st);
  EXPECT_EQ((int)hull_vertices(hull).size(), NP);
  EXPECT_EQ((int)hull.size(), 2*NP - 4);
  expect_watertight(hull);

  const double nlogn = double(NP) * (log(double(NP)) / log(2.0));
  EXPECT_LT((double)st.conflict_tests, 60.0 * nlogn)
      << "conflict tests " << st.conflict_tests << " exceed 60*n*log2(n) = "
      << 60.0 * nlogn << " -- quadratic regression?";
  EXPECT_LT((double)st.faces_created, 30.0 * NP)
      << "faces created " << st.faces_created << " exceed 30n";
  fprintf(stderr, "[stats] n=%d conflict_tests=%zu (%.1f per n*log2 n) "
          "faces_created=%zu (%.2f per n) inserted=%zu\n",
          NP, st.conflict_tests, st.conflict_tests / nlogn,
          st.faces_created, st.faces_created / double(NP), st.points_inserted);
}
