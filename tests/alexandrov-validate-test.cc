// The polytope acceptance gate (delaunay_polytope.hh), which the solver's
// validate_polytope and a GPU batch validator both run.
//
// The gate was promoted out of the solver so a device kernel can apply it:
// its 2-skeleton census counts cells and tests simplicity WITHOUT building
// them, since a kernel has neither heap nor a marker array over the
// half-edges.  That rewrite is what this test gates: over a real corpus, on
// both polyhedral metrics, the census must agree cell for cell with the
// tesselation it replaces, and the assembled ladder must agree with the
// solver's own verdict.
//
//   - census.n_cells == polytope_tesselation(...).n_cells(), and
//     census.simple == is_simple_polygonal(that tesselation), on every
//     isomer of C20..C50 and the icosahedral C60, both metrics.  The corpus
//     deliberately includes the multi-edge delta-complexes and the flat
//     hexagonal caps (C24-D6d, C36-D6h), where cells are polygons rather
//     than triangles and a repeated corner is a real possibility.
//   - the gate's verdict agrees with the solver's reported status.
//   - a realized polytope passes, and a deliberately broken one fails with
//     the property that was broken: a pushed-in vertex is not convex, a
//     collapsed one is volume-degenerate.

#include "fullerenes/buckygen-wrapper.hh"
#include "fullerenes/delaunay_alexandrov.hh"
#include "fullerenes/delaunay_polytope.hh"
#include "fullerenes/spiral.hh"

#include <gtest/gtest.h>
#include <numeric>
#include <string>
#include <vector>

using namespace std;

namespace {

// The interior-edge predicate of the 2-skeleton, in the host's arithmetic:
// the solver's inessential mask, read as a predicate.
struct Tight {
  const vector<bool>& mask;
  bool operator()(int h) const { return mask[h]; }
};

// Every isomer of one size, as dual triangulations.
vector<Triangulation> duals_of(int N) {
  vector<Triangulation> out;
  auto Q = BuckyGen::start(N, false, false);
  Graph G;
  while (BuckyGen::next_fullerene(Q, G)) out.push_back(Triangulation(G));
  BuckyGen::stop(Q);
  return out;
}

struct Solved {
  DelaunayTriangulation D;
  vector<double> r;
  vector<coord3d> pos;
  AlexandrovSolver::ValidationStatus status;
};

Solved solve_dual(const Triangulation& T) {
  AlexandrovSolver S;
  S.D = DelaunayTriangulation::compute(T);
  vector<coord3d> pos = S.solve();
  return {S.D, S.r, pos, S.stats_status};
}

Solved solve_cubic(const Triangulation& T) {
  AlexandrovIDTCubic AC;
  const auto P = AC.solve_polytope(T);
  return {AC.solver.D, AC.solver.r, P.positions, P.status};
}

// The census against the tesselation it replaces, plus the gate against the
// solver's own verdict, for one solved isomer.
void check_agreement(const Solved& s, const string& what) {
  SCOPED_TRACE(what);
  if (s.pos.empty()) return;                     // reconstruction refused: nothing realized
  const vector<bool> mask = AlexandrovSolver::inessential_edges(s.D, s.r);
  const polytope::CellCensus c = polytope::cell_census(s.D, Tight{mask});
  ASSERT_TRUE(c.closed) << "a cell-boundary walk did not close";

  vector<int> labels(s.D.nv);
  iota(labels.begin(), labels.end(), 0);
  const CanonicalTesselation tbar =
      AlexandrovSolver::polytope_tesselation(s.D, s.r, labels);
  EXPECT_EQ(c.n_cells, tbar.n_cells()) << "the census counts the tesselation's cells";
  EXPECT_EQ(c.simple, AlexandrovSolver::is_simple_polygonal(tbar))
      << "the census decides simplicity as the tesselation does";

  polytope::Record rec;
  const polytope::Verdict v = polytope::validate(s.D, s.pos, Tight{mask}, &rec);
  const bool gate_ok = v == polytope::Verdict::Ok;
  EXPECT_EQ(gate_ok, s.status == AlexandrovSolver::ValidationStatus::OK)
      << "gate says " << polytope::verdict_str(v) << ", solver says "
      << AlexandrovSolver::status_str(s.status);
}

}  // namespace

TEST(PolytopeGate, CensusAgreesWithTesselationOnBothMetrics) {
  for (int N : {20, 24, 26, 28, 30, 32, 34, 36, 38, 40, 42, 44, 46, 48, 50}) {
    const vector<Triangulation> duals = duals_of(N);
    ASSERT_FALSE(duals.empty()) << "C" << N << " enumerated no isomer";
    for (size_t i = 0; i < duals.size(); i++) {
      check_agreement(solve_dual(duals[i]), "C" + to_string(N) + " #" + to_string(i) + " dual");
      check_agreement(solve_cubic(duals[i]), "C" + to_string(N) + " #" + to_string(i) + " cubic");
    }
  }
}

TEST(PolytopeGate, IcosahedralC60PassesOnBothMetrics) {
  const string name = "[1,7,9,11,13,15,18,20,22,24,26,32]-fullerene";
  const Triangulation T{spiral_nomenclature(name)};
  for (const Solved& s : {solve_dual(T), solve_cubic(T)}) {
    ASSERT_FALSE(s.pos.empty());
    const vector<bool> mask = AlexandrovSolver::inessential_edges(s.D, s.r);
    polytope::Record rec;
    const polytope::Verdict v = polytope::validate(s.D, s.pos, Tight{mask}, &rec);
    EXPECT_EQ(v, polytope::Verdict::Ok) << polytope::verdict_str(v);
    EXPECT_GE(rec.n_cells, 3);
    EXPECT_GT(rec.volume_norm, polytope::kVolumeFloor);
    EXPECT_TRUE(rec.convex);
    EXPECT_TRUE(rec.no_self_cross);
  }
}

TEST(PolytopeGate, BrokenGeometryFailsWithItsOwnProperty) {
  const string name = "[1,7,9,11,13,15,18,20,22,24,26,32]-fullerene";
  const Triangulation T{spiral_nomenclature(name)};
  const Solved s = solve_dual(T);
  ASSERT_FALSE(s.pos.empty());
  const vector<bool> mask = AlexandrovSolver::inessential_edges(s.D, s.r);

  // A vertex pushed through the body: no longer convex.
  {
    vector<coord3d> pos = s.pos;
    pos[0] = pos[0] * -0.5;
    const polytope::Verdict v = polytope::validate(s.D, pos, Tight{mask});
    EXPECT_EQ(v, polytope::Verdict::NotConvex) << polytope::verdict_str(v);
  }
  // The body collapsed onto a plane: degenerate volume.
  {
    vector<coord3d> pos = s.pos;
    for (coord3d& p : pos) p[2] = 0;
    const polytope::Verdict v = polytope::validate(s.D, pos, Tight{mask});
    EXPECT_EQ(v, polytope::Verdict::VolumeDegenerate) << polytope::verdict_str(v);
  }
}
