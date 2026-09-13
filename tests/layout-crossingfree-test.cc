// Unit coverage for layout2d::layout_is_crossingfree, one case per defect of
// the slope-intercept form it used to be written in (refactor-debt entry
// 2026-08-09-layout-is-crossingfree-slope-form).
//
// The predicate is now the orientation determinant, so both cases below are
// decided without a division.  Each fails under the old form:
//
//   1. A VERTICAL edge made the slope (ay-by)/(ax-bx) divide by zero, and the
//      inf/NaN comparisons that followed decided nothing, so the pair fell
//      through to "intersect".  This is the defect that actually reddened
//      orientation-test: the Tutte layout of C60-GC(1,1) puts both ends of
//      edge [8,9] at x = 0.36562987420987297.
//   2. The straddle tests were STRICT, so a point exactly ON the other line
//      satisfied neither branch and was likewise reported as an intersection.
//
// The graphs here are hand-built quadrilaterals rather than fullerenes: the
// point is to pin the predicate on a named configuration, which a 60-vertex
// barycentric solution cannot do.

#include <gtest/gtest.h>
#include "fullerenes/planargraph.hh"
#include "fullerenes/layout2d.hh"

using namespace std;

// A 4-cycle 0-1-2-3 plus the two diagonals' worth of vertices is more than we
// need; what the predicate reads is the edge list and the coordinates, so a
// graph carrying exactly the two edges under test is the honest fixture.
static PlanarGraph two_edge_graph()
{
  Graph g(4, 3);
  g.push_back(0, 1); g.push_back(1, 0);   // edge A
  g.push_back(2, 3); g.push_back(3, 2);   // edge B
  return PlanarGraph(g);
}

// Defect 1: edge B is exactly vertical and the two edges do not cross.
// Old form: dx == 0 -> infinite slope -> reported as intersecting.
TEST(LayoutCrossingFree, VerticalEdgeThatDoesNotCross)
{
  PlanarGraph G = two_edge_graph();
  vector<coord2d> layout(4);
  layout[0] = coord2d(0.0, 0.0);   // edge A: a short horizontal segment
  layout[1] = coord2d(1.0, 0.0);
  layout[2] = coord2d(2.0, -1.0);  // edge B: exactly vertical, x == 2.0,
  layout[3] = coord2d(2.0,  1.0);  //         and well clear of A
  EXPECT_TRUE(layout2d::layout_is_crossingfree(G, layout));
}

// The same vertical edge, now genuinely crossed: the verdict must still be no.
// Guards against "fixing" defect 1 by making the predicate permissive.
TEST(LayoutCrossingFree, VerticalEdgeThatDoesCross)
{
  PlanarGraph G = two_edge_graph();
  vector<coord2d> layout(4);
  layout[0] = coord2d(1.0, 0.0);
  layout[1] = coord2d(3.0, 0.0);
  layout[2] = coord2d(2.0, -1.0);  // vertical, straddled by A
  layout[3] = coord2d(2.0,  1.0);
  EXPECT_FALSE(layout2d::layout_is_crossingfree(G, layout));
}

// Defect 2: two collinear edges lying on a common line but NOT overlapping.
// Old form: the strict straddle tests rejected the exact zeros and called it an
// intersection.  Disjoint segments on one line do not intersect.
TEST(LayoutCrossingFree, CollinearButDisjoint)
{
  PlanarGraph G = two_edge_graph();
  vector<coord2d> layout(4);
  layout[0] = coord2d(0.0, 0.0);   // [0,1] spans x in [0,1] on the x-axis
  layout[1] = coord2d(1.0, 0.0);
  layout[2] = coord2d(2.0, 0.0);   // [2,3] spans x in [2,3] on the same line
  layout[3] = coord2d(3.0, 0.0);
  EXPECT_TRUE(layout2d::layout_is_crossingfree(G, layout));
}

// Collinear AND overlapping is a degenerate drawing, reported as such rather
// than as a crossing -- but still not a crossing-free embedding.
TEST(LayoutCrossingFree, CollinearOverlapping)
{
  PlanarGraph G = two_edge_graph();
  vector<coord2d> layout(4);
  layout[0] = coord2d(0.0, 0.0);
  layout[1] = coord2d(2.0, 0.0);
  layout[2] = coord2d(1.0, 0.0);   // sits inside [0,1]'s span
  layout[3] = coord2d(3.0, 0.0);
  EXPECT_FALSE(layout2d::layout_is_crossingfree(G, layout));
}

// An endpoint landing exactly on the other edge's interior: a T-junction
// between edges that share no vertex.  Also degenerate, also not crossing-free.
TEST(LayoutCrossingFree, EndpointOnOtherEdgeInterior)
{
  PlanarGraph G = two_edge_graph();
  vector<coord2d> layout(4);
  layout[0] = coord2d(0.0, 0.0);
  layout[1] = coord2d(2.0, 0.0);
  layout[2] = coord2d(1.0, 0.0);   // exactly on the interior of [0,1]
  layout[3] = coord2d(1.0, 1.0);
  EXPECT_FALSE(layout2d::layout_is_crossingfree(G, layout));
}
