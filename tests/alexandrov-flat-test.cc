// The doubling witness (AlexandrovSolver::doubled_polygon): the exact
// certificate that a twelve-cone dual metric is a convex polygon doubled
// along its boundary, read off the solver's own floor state.
//
//   - The five known flat fullerenes -- C96 and IPR C384 of delaunay-
//     geometry's counterexample note, and C120, C170, C180 found by the
//     whole-space sweeps -- refuse the dual solve, and the twelve edges of
//     smallest dihedral angle at the floor state verify as the boundary of
//     a doubled polygon of N/2 unit triangles (the dual has N faces, two
//     sheets of N/2 each) with the twelve cones as its corners.
//   - Ordinary isomers (the dodecahedral C20, the icosahedral C60) deliver,
//     and their proposal fails the verification.
//   - The verdict is the proposal's, not the state's: a verified proposal
//     with one edge exchanged for a sheet diagonal fails.

#include "fullerenes/delaunay_alexandrov.hh"
#include "fullerenes/spiral.hh"

#include <gtest/gtest.h>
#include <algorithm>
#include <array>
#include <string>
#include <vector>

using namespace std;

namespace {

// The dual triangulation of a canonically named isomer ("C<N>-[...]-fullerene").
Triangulation dual_of(const string& name) {
  const string spiral = name.substr(name.find('-') + 1);
  return Triangulation(spiral_nomenclature(spiral));
}

// The solver on the isomer's twelve-cone dual metric, left at its final
// state (D flipped, r the final radii).
AlexandrovSolver solve_dual(const Triangulation& T) {
  AlexandrovSolver S;
  S.D = DelaunayTriangulation::compute(T);
  S.solve();
  return S;
}

struct Flat { const char* name; int carbons; };
const Flat kFlat[] = {
    {"C96-[1,2,12,17,23,29,35,40,41,45,49,50]-fullerene", 96},
    {"C120-[1,2,12,17,23,35,41,52,53,57,61,62]-fullerene", 120},
    {"C170-[1,2,12,17,23,54,61,67,78,82,86,87]-fullerene", 170},
    {"C180-[1,11,16,30,38,54,62,78,80,85,90,92]-fullerene", 180},
    {"C384-[1,16,40,71,83,119,131,164,166,174,192,194]-fullerene", 384},
};
const char* const kOrdinary[] = {
    "C20-[1,2,3,4,5,6,7,8,9,10,11,12]-fullerene",
    "C60-[1,7,9,11,13,15,18,20,22,24,26,32]-fullerene",
};

}  // namespace

// The SEARCH decides flatness on the input complex: no solve, no radii, no
// coordinates.  The fold is already an edge set of the intrinsic Delaunay
// complex on every fullerene known to be flat, so the solver's floor state
// is not needed to find it.
TEST(AlexandrovFlat, TheInputComplexDecidesFlatness) {
  for (const Flat& f : kFlat) {
    SCOPED_TRACE(f.name);
    const DelaunayTriangulation D = DelaunayTriangulation::compute(dual_of(f.name));
    const AlexandrovSolver::DoubledPolygon W = AlexandrovSolver::doubled_polygon(D);
    ASSERT_TRUE(W.ok()) << AlexandrovSolver::doubling_verdict_str(W.verdict);
    EXPECT_EQ(W.area, f.carbons / 2);
    array<int, 12> corners = W.corner;
    sort(corners.begin(), corners.end());
    for (int k = 0; k < 12; k++) EXPECT_EQ(corners[k], k) << "the corners are the twelve cones";
  }
}

TEST(AlexandrovFlat, OrdinaryInputComplexesHoldNoFold) {
  for (const char* name : kOrdinary) {
    SCOPED_TRACE(name);
    const DelaunayTriangulation D = DelaunayTriangulation::compute(dual_of(name));
    const AlexandrovSolver::DoubledPolygon W = AlexandrovSolver::doubled_polygon(D);
    EXPECT_FALSE(W.ok()) << "a realizable metric has no fold";
    EXPECT_EQ(W.verdict, AlexandrovSolver::DoublingVerdict::NoFoldFound);
  }
}

TEST(AlexandrovFlat, KnownFlatFullerenesAreCertified) {
  for (const Flat& f : kFlat) {
    SCOPED_TRACE(f.name);
    const AlexandrovSolver S = solve_dual(dual_of(f.name));
    EXPECT_NE(S.stats_status, AlexandrovSolver::ValidationStatus::OK)
        << "a flat metric has no three-dimensional realization to deliver";
    const AlexandrovSolver::DoubledPolygon W = AlexandrovSolver::doubled_polygon(S.D, S.r);
    ASSERT_TRUE(W.ok()) << AlexandrovSolver::doubling_verdict_str(W.verdict);
    EXPECT_EQ(W.area, f.carbons / 2);
    array<int, 12> corners = W.corner;
    sort(corners.begin(), corners.end());
    for (int k = 0; k < 12; k++) EXPECT_EQ(corners[k], k) << "the corners are the twelve cones";
    for (int k = 0; k < 12; k++) {
      EXPECT_EQ(S.D.he_origin[W.boundary[k]], W.corner[k]);
      EXPECT_EQ(S.D.dest(W.boundary[k]), W.corner[(k + 1) % 12]) << "the fold is one cycle";
    }
  }
}

TEST(AlexandrovFlat, OrdinaryIsomersAreNot) {
  for (const char* name : kOrdinary) {
    SCOPED_TRACE(name);
    const AlexandrovSolver S = solve_dual(dual_of(name));
    EXPECT_EQ(S.stats_status, AlexandrovSolver::ValidationStatus::OK);
    const AlexandrovSolver::DoubledPolygon W = AlexandrovSolver::doubled_polygon(S.D, S.r);
    EXPECT_FALSE(W.ok()) << "a delivered polytope has no doubling witness";
  }
}

TEST(AlexandrovFlat, TheVerdictIsTheProposals) {
  const AlexandrovSolver S = solve_dual(dual_of(kFlat[1].name));   // C120
  const AlexandrovSolver::DoubledPolygon W = AlexandrovSolver::doubled_polygon(S.D, S.r);
  ASSERT_TRUE(W.ok());
  array<int, 12> edges{};
  for (int k = 0; k < 12; k++) edges[k] = S.D.edge(W.boundary[k]);
  // The same twelve edges in another order verify to the same polygon.
  reverse(edges.begin(), edges.end());
  const AlexandrovSolver::DoubledPolygon W2 = AlexandrovSolver::verify_doubled_polygon(S.D, edges);
  ASSERT_TRUE(W2.ok());
  EXPECT_EQ(W2.area, W.area);
  // One fold edge exchanged for a diagonal of the first sheet fails.
  const auto hs = S.D.face_halfedges(W.sheet_a[0]);
  int diagonal = -1;
  for (int h : hs)
    if (find(edges.begin(), edges.end(), S.D.edge(h)) == edges.end()) diagonal = S.D.edge(h);
  ASSERT_GE(diagonal, 0);
  edges[0] = diagonal;
  EXPECT_FALSE(AlexandrovSolver::verify_doubled_polygon(S.D, edges).ok());
}
