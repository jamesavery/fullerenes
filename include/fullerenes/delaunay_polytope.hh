#pragma once
// The acceptance gate of a realized convex polytope, as DEVICE-LEGAL bodies:
// given a metric complex (DelaunayView), an interior-edge mask and the
// vertex positions, decide whether the result is a well-formed convex
// polytope.  One body per property, callable from a host loop and from a GPU
// kernel over one isomer's workspace alike -- the parent's
// AlexandrovSolver::validate_polytope is the host wrapper around the ladder
// below, and a batched device validator runs the same functions on the same
// data (promotion pattern: DelaunayView::compact_vertices,
// eisenstein_paint_tables::cell_developments_into).
//
// DEVICE-LEGAL means: no allocation, no std::vector/std::set, no exceptions,
// no I/O, and no transcendental functions.  The last one is why the
// INTERIOR-EDGE MASK is a caller-supplied predicate rather than computed
// here: deciding whether the dihedral at an edge is pi needs an arctangent
// per pyramid, which each caller already has in the arithmetic it runs in
// (the parent solver's GCP::theta on the host; the port's field-surface
// Geometry::theta in a kernel, gated equal to it).  Everything the gate
// decides from that mask and the positions is rational arithmetic plus
// square roots, and lives here.
//
// THE GATE, in the order a failure is reported (Bobenko-Izmestiev's
// conditions for the deformation's endpoint to be the polytope):
//
//   SIMPLICITY       the polytope's 2-skeleton -- the cells of the complex
//                    after collapsing every interior edge -- is a simple
//                    polygonal tesselation (every cell at least a triangle,
//                    no cell visiting a vertex twice) with at least three
//                    cells.  A non-simple cell means the "polytope" folds
//                    onto itself; fewer than three cells is a degenerate
//                    two-sided or one-sided body.
//   WELL-FORMEDNESS  the enclosed volume, in units of the mean edge length
//                    cubed, is not degenerate, AND no two non-adjacent
//                    triangles cross in space.  Convexity does imply the
//                    second, but it is part of the definition of a closed
//                    embedded surface and is enforced on its own so that a
//                    failure is reported with its true cause.
//   CONVEXITY        every vertex lies on the inner side of every face's
//                    plane, with the outward normal taken from the CCW
//                    half-edge convention.
//
// No spherical or radial assumption enters: normals come from the face
// orientation, never from a centre.

#include "fullerenes/delaunay_view.hh"
#include "fullerenes/geometry.hh"

#include <cmath>
#include <span>

namespace polytope {

// The verdict of the gate below.  A solver's own status vocabulary adds its
// convergence outcomes to these (AlexandrovSolver::ValidationStatus); this
// enum is exactly the properties of a realized polytope, which is what the
// device decides.
enum class Verdict {
  Ok,                 // a well-formed convex polytope
  NotSimple,          // the 2-skeleton is not a simple tesselation, or has < 3 cells
  WalkUnclosed,       // a cell boundary did not close: the complex is corrupt
  VolumeDegenerate,   // flat or near-flat: volume / <edge>^3 below the floor
  SelfIntersecting,   // two non-adjacent triangles cross in space
  NotConvex,          // a vertex lies outside a face plane
};

inline const char* verdict_str(Verdict v) {
  switch (v) {
    case Verdict::Ok:               return "OK";
    case Verdict::NotSimple:        return "NOT_SIMPLE";
    case Verdict::WalkUnclosed:     return "CELL_WALK_UNCLOSED";
    case Verdict::VolumeDegenerate: return "VOLUME_DEGENERATE";
    case Verdict::SelfIntersecting: return "SELF_INTERSECTING";
    case Verdict::NotConvex:        return "NOT_CONVEX";
  }
  return "?";
}

// The gate's constants, named once so host and device apply the same ones.
//   kVolumeFloor  volume / <edge>^3 below this is degenerate.  A 1.03M-isomer
//                 scan separates healthy polytopes (>= 0.12, median 1.01)
//                 from degenerate ones (<= 1e-6) by five orders of
//                 magnitude; the floor sits in the empty gap.
//   kConvexTol    slack of the vertex-inside-plane test, relative to the
//                 mean edge length.
//   kCrossTol     slack of the signed-distance tests of the triangle
//                 crossing predicate.
// The two tolerances are the values the parent's AlexandrovSolver::is_convex
// and ::has_self_intersection have always defaulted to; they are named here
// so the host gate and a device batch apply the same numbers.
inline constexpr double kVolumeFloor = 0.01;
inline constexpr double kConvexTol   = 1e-3;
inline constexpr double kCrossTol    = 1e-6;

// ── The 2-skeleton's census ────────────────────────────────────────────────
// The cells of D after collapsing every edge the mask calls interior, counted
// and tested for simplicity WITHOUT materializing them: each cell is walked
// from every one of its boundary half-edges, and counted once, at the walk's
// minimal slot (DelaunayView::CellWalk::hmin, invariant under interior
// flips).  That is what makes the census allocation-free -- a "seen" array
// over the half-edges is exactly what a kernel cannot have.
struct CellCensus {
  int  n_cells = 0;
  bool simple  = false;   // every cell has >= 3 corners, all distinct
  bool closed  = true;    // every cell-boundary walk closed
  int  witness = -1;      // the half-edge of the first non-closing walk, or
                          // of the first non-simple cell
};

// @pre  tight(h) == tight(twin(h)) for every live h (an edge is interior
//       from both sides)
// @post result.n_cells is the number of cells; result.simple iff every cell
//       has at least three corners and no repeated corner
template <class Tight>
inline CellCensus cell_census(const DelaunayView& D, Tight&& tight) {
  CellCensus c;
  c.simple = true;
  for (int h0 = 0; h0 < D.nh; h0++) {
    if (!D.alive(h0) || tight(h0)) continue;
    const typename DelaunayView::CellWalk C = D.walk_cell(h0, tight);
    if (C.closure != DelaunayView::CellWalk::Closure::Closed) {
      c.closed = false;
      c.witness = C.witness;
      return c;
    }
    if (C.hmin != h0) continue;             // this cell is counted at its minimal slot
    c.n_cells++;
    if (C.d < 3) { c.simple = false; if (c.witness < 0) c.witness = h0; continue; }
    // No corner twice.  The cell is walked once per corner and compared
    // against the corners before it -- quadratic in the cell's degree (a
    // handful), and free of the marker array a kernel cannot allocate.
    int i = 0;
    bool distinct = true;
    D.visit_cell(h0, tight, [&](int h) {
      int j = 0;
      D.visit_cell(h0, tight, [&](int g) {
        if (j < i && D.he_origin[g] == D.he_origin[h]) distinct = false;
        j++;
      });
      i++;
    });
    if (!distinct) { c.simple = false; if (c.witness < 0) c.witness = h0; }
  }
  return c;
}

// ── Extrinsic geometry of a realized complex ───────────────────────────────

// The mean length of D's live edges (the scale every relative tolerance of
// the gate is measured in); 1 for an edgeless complex, so a ratio is never
// divided by zero.
inline double mean_edge_length(const DelaunayView& D) {
  double sum = 0;
  int n = 0;
  for (int h = 0; h < D.nh; h += 2) {
    if (!D.alive(h)) continue;
    sum += D.he_length[h];
    n++;
  }
  return n > 0 ? sum / (double)n : 1.0;
}

// Six times the signed volume enclosed by pos under D's face structure:
// sum over live faces of a . (b x c), the signed tetrahedra from the origin.
// Positive exactly when the CCW half-edge convention gives outward normals
// (the divergence theorem on a closed oriented surface).
inline double signed_volume6(const DelaunayView& D, std::span<const coord3d> pos) {
  double vol6 = 0;
  for (int f = 0; f < D.nf; f++) {
    if (D.f_he[f] < 0) continue;
    const auto v = D.face_vertices(f);
    vol6 += pos[v[0]].dot(pos[v[1]].cross(pos[v[2]]));
  }
  return vol6;
}

// Do the triangles (a0,b0,c0) and (a1,b1,c1) cross with non-empty common
// interior?  Moller-style: classify each triangle's corners by signed
// distance to the other's plane, and when both planes interleave, intersect
// the two crossing intervals along the planes' common line.  `tol` is the
// slack of the signed-distance tests.  Triangles sharing a corner are the
// caller's to filter (they are reported as not crossing).
//
// COPLANAR PAIRS ARE NOT CROSSINGS.  Two coplanar non-adjacent triangles on a
// valid polytope are sub-triangulations of one flat 2-face (the complex
// triangulates a flat polygon into several triangles), which share interior
// without the surface crossing itself.  The other case, two distinct 2-faces
// in one plane, is a flat body, which the volume floor already rejects.
inline bool tri_tri_intersect(const coord3d& a0, const coord3d& b0, const coord3d& c0,
                              const coord3d& a1, const coord3d& b1, const coord3d& c1,
                              double tol = kCrossTol) {
  auto plane_dist = [](const coord3d& a, const coord3d& b, const coord3d& c,
                       const coord3d& p) {
    const coord3d n = (b - a).cross(c - a);
    return n.dot(p - a);            // unscaled signed distance, times twice the area
  };
  const double da0 = plane_dist(a1, b1, c1, a0);
  const double db0 = plane_dist(a1, b1, c1, b0);
  const double dc0 = plane_dist(a1, b1, c1, c0);
  if ((da0 > tol && db0 > tol && dc0 > tol) ||
      (da0 < -tol && db0 < -tol && dc0 < -tol)) return false;
  const double da1 = plane_dist(a0, b0, c0, a1);
  const double db1 = plane_dist(a0, b0, c0, b1);
  const double dc1 = plane_dist(a0, b0, c0, c1);
  if ((da1 > tol && db1 > tol && dc1 > tol) ||
      (da1 < -tol && db1 < -tol && dc1 < -tol)) return false;

  const coord3d n0 = (b0 - a0).cross(c0 - a0);
  const coord3d n1 = (b1 - a1).cross(c1 - a1);
  const coord3d L  = n0.cross(n1);
  const double L2 = L.dot(L), n0_2 = n0.dot(n0), n1_2 = n1.dot(n1);
  // |L|^2 = |n0|^2 |n1|^2 sin^2(angle): coplanar within a microradian.
  if (L2 < 1e-12 * n0_2 * n1_2) return false;

  struct Interval { double lo, hi; };
  auto interval = [&](const coord3d& a, const coord3d& b, const coord3d& c,
                      double da, double db, double dc) {
    auto edge_t = [&](const coord3d& p, const coord3d& q, double dp, double dq) {
      const double s = dp / (dp - dq);
      return (p + (q - p) * s).dot(L);
    };
    Interval r{1e300, -1e300};
    auto upd = [&](double t) { if (t < r.lo) r.lo = t; if (t > r.hi) r.hi = t; };
    if ((da > 0) != (db > 0)) upd(edge_t(a, b, da, db));
    if ((db > 0) != (dc > 0)) upd(edge_t(b, c, db, dc));
    if ((dc > 0) != (da > 0)) upd(edge_t(c, a, dc, da));
    if (std::fabs(da) <= tol) upd(a.dot(L));   // a corner on the plane counts
    if (std::fabs(db) <= tol) upd(b.dot(L));
    if (std::fabs(dc) <= tol) upd(c.dot(L));
    return r;
  };
  const Interval i0 = interval(a0, b0, c0, da0, db0, dc0);
  const Interval i1 = interval(a1, b1, c1, da1, db1, dc1);
  if (i0.lo > i0.hi || i1.lo > i1.hi) return false;   // one triangle misses the other's plane
  return !(i0.hi < i1.lo - tol || i1.hi < i0.lo - tol);
}

// Does the realized surface cross itself?  Every pair of live faces that
// share no corner is tested; adjacent faces meet along an edge by
// construction and are skipped.
inline bool has_self_intersection(const DelaunayView& D, std::span<const coord3d> pos,
                                  double tol = kCrossTol) {
  for (int f = 0; f < D.nf; f++) {
    if (D.f_he[f] < 0) continue;
    const auto A = D.face_vertices(f);
    for (int g = f + 1; g < D.nf; g++) {
      if (D.f_he[g] < 0) continue;
      const auto B = D.face_vertices(g);
      bool share = false;
      for (int i = 0; i < 3; i++)
        for (int j = 0; j < 3; j++)
          if (A[i] == B[j]) share = true;
      if (share) continue;
      if (tri_tri_intersect(pos[A[0]], pos[A[1]], pos[A[2]],
                            pos[B[0]], pos[B[1]], pos[B[2]], tol)) return true;
    }
  }
  return false;
}

// Is every vertex on the inner side of every face's plane?  The outward
// normal of a face is (b - a) x (c - a) over its corners in CCW half-edge
// order, which points outward exactly when the enclosed signed volume is
// positive -- checked first, so a globally inverted placement is rejected
// here rather than passing vacuously.  A flat body fails the same check.
// `tol` is relative to the mean edge length.
inline bool is_convex(const DelaunayView& D, std::span<const coord3d> pos,
                      double tol = kConvexTol) {
  const double mean_l = mean_edge_length(D);
  const double abs_tol = tol * mean_l;
  const double vol6 = signed_volume6(D, pos);
  const double vol_threshold6 = 1e-6 * mean_l * mean_l * mean_l;
  if (!std::isfinite(vol6) || vol6 < vol_threshold6) return false;
  for (int f = 0; f < D.nf; f++) {
    if (D.f_he[f] < 0) continue;
    const auto [va, vb, vc] = D.face_vertices(f);
    const coord3d a = pos[va], b = pos[vb], c = pos[vc];
    const coord3d nf_raw = (b - a).cross(c - a);
    const double nlen = std::sqrt(nf_raw.dot(nf_raw));
    if (nlen < 1e-15) continue;                       // degenerate triangle
    const coord3d nf = nf_raw * (1.0 / nlen);
    for (int v = 0; v < D.nv; v++) {
      if (v == va || v == vb || v == vc) continue;
      if ((pos[v] - a).dot(nf) > abs_tol) return false;
    }
  }
  return true;
}

// ── The gate ───────────────────────────────────────────────────────────────

// What each rung of the gate found.  A rung a failure skipped keeps its
// default, so the record says which check produced the verdict.
struct Record {
  bool   simple        = false;   // the 2-skeleton is a simple tesselation
  int    n_cells       = 0;       // its cells
  bool   closed        = true;    // every cell-boundary walk closed
  double volume_norm   = 0;       // volume / <edge>^3
  bool   no_self_cross = false;
  bool   convex        = false;
  int    witness       = -1;      // the offending half-edge, when there is one
};

// The gate, in the banner's order.  `tight(h)` says whether edge h is
// interior to a cell of the 2-skeleton (the caller's dihedral test).
// @pre  pos.size() >= D.nv, and the complex is the one pos realizes
// @post the verdict names the first property that failed; `out` holds every
//       rung that ran
template <class Tight>
inline Verdict validate(const DelaunayView& D, std::span<const coord3d> pos,
                        Tight&& tight, Record* out = nullptr) {
  Record v;
  auto done = [&](Verdict s) { if (out) *out = v; return s; };

  const CellCensus c = cell_census(D, tight);
  v.simple = c.simple && c.n_cells >= 3;
  v.n_cells = c.n_cells;
  v.closed = c.closed;
  v.witness = c.witness;
  if (!c.closed) return done(Verdict::WalkUnclosed);
  if (!v.simple) return done(Verdict::NotSimple);

  const double mean_l = mean_edge_length(D);
  const double volume = std::fabs(signed_volume6(D, pos)) / 6.0;
  v.volume_norm = volume / (mean_l * mean_l * mean_l);
  if (!(std::isfinite(v.volume_norm) && v.volume_norm >= kVolumeFloor))
    return done(Verdict::VolumeDegenerate);

  v.no_self_cross = !has_self_intersection(D, pos);
  if (!v.no_self_cross) return done(Verdict::SelfIntersecting);

  v.convex = is_convex(D, pos);
  if (!v.convex) return done(Verdict::NotConvex);

  return done(Verdict::Ok);
}

}  // namespace polytope
