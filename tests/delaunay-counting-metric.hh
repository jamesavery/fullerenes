// Counting instrumentation for the exact-metric Delaunay driver: a
// forwarding wrapper around ExactIntegerMetric whose only addition is a
// counter on first_tie_side -- the observable that the cocircular
// tie-break was actually CONSULTED (a self-loop resisted both direct
// flips, i.e. a genuine exact tie).
//
// Single source shared by the parent falsifier test
// (tests/delaunay-test.cc, DCELExact.ConstructedTieExercisesTieBreak) and
// the research seed-search tool
// (claude-projects/delaunay/tools/diag_tie_exact.cc): the search's
// verdict transfers to the test only if both count the SAME event
// through the SAME driver, so neither may carry its own copy.
#pragma once

#include "fullerenes/delaunay.hh"

struct CountingExactMetric {
  ExactIntegerMetric inner;
  int* n_tie;
  // Last-consultation observation: the exact rule's side, and the side the
  // banded float rule (theta[0] <= theta[1]) would have taken there.  At an
  // exact tie both thetas are pi in reals, so the float side is rounding
  // noise -- report it, never assert on it; the exact side is the
  // deterministic one.
  int* last_side;
  int* last_float_side;
  bool is_flat(const DelaunayView& V, int v) const { return inner.is_flat(V, v); }
  bool delaunay(const DelaunayView& V, int h) const { return inner.delaunay(V, h); }
  bool convex(const DelaunayView& V, int h) const { return inner.convex(V, h); }
  bool cocircular(const DelaunayView& V, int h) const { return inner.cocircular(V, h); }
  std::optional<Length> flipped(DelaunayView& V, int h) const { return inner.flipped(V, h); }
  void set_edge_length(DelaunayView& V, int h, Length l) const { inner.set_edge_length(V, h, l); }
  void prepare_star(DelaunayView& V, DelaunayWorkspace& ws, int v) const {
    inner.prepare_star(V, ws, v);
  }
  Length ear(DelaunayView& V, const FanPolygon& fan, int pp, int pi, int pn) const {
    return inner.ear(V, fan, pp, pi, pn);
  }
  int first_tie_side(DelaunayView& V, DelaunayWorkspace& ws, int h_loop,
                     const double* theta) const {
    ++*n_tie;
    const int side = inner.first_tie_side(V, ws, h_loop, theta);
    *last_side = side;
    *last_float_side = theta[0] <= theta[1] ? 0 : 1;
    return side;
  }
};
