#pragma once

// ============================================================================
// The DCEL point tracker's FLIP-TAPE TRANSPORT: the tracker state as a span
// view, and the Transport policy that carries it through the two
// topology-changing operations.
//
// DelaunayView's surgery bodies (flip_edge, remove_flat_vertex,
// delaunay_view.hh) call a Transport policy's five hooks; TrackerTransport
// below is THE policy that transports tracked points, and it is written over
// caller-owned spans -- so the host owner (DelaunayTriangulation, whose
// DelaunayPointTracker owns the vectors behind the view) and a device kernel
// (whose arena carves them) run ONE body.
//
// The mathematics: for a surgery S re-triangulating a flat region R, the
// transport is the composition
//     T_S  =  express_post^-1  o  express_pre
// where express_pre maps (face, barycentric) -> R^2 through the PRE-surgery
// triangles' isometric development of R, and express_post^-1 locates the
// planar point in the POST-surgery triangles of the SAME development.  Both
// chart sets develop R isometrically, so T_S fixes the intrinsic surface
// point (the transport lemma).  Both halves are computed before the first
// write -- the post-charts are determined by the planned surgery.
// @post (transport) each point's intrinsic surface position is unchanged:
//       the developments are isometries, re-expression only changes
//       coordinates.
//
// Flip development: u = (0,0), v = (e,0), B above, D below; positions keyed
// per HALF-EDGE (slot anchoring), so repeated-corner faces transport the
// same.  The post-flip charts are (B, D, v) -> fh and (D, B, u) -> ft.
// Star development: the FanPolygon IS the isometric development (apex at the
// origin); the removed vertex seeds a tracked point with label = v.
//
// FAILURE CHANNEL.  A kernel cannot throw, so every transport failure trips
// the DCEL's own Status latch (the run-path convention of delaunay_view.hh),
// which is TERMINAL for the complex; the owner converts it to its documented
// throw at its boundary (throw_on_status).  The latch keeps the removal's
// commit-or-nothing property -- plan_charts trips => splice_fan is a no-op =>
// commit_removal is skipped -- while a flip whose plan hook trips still
// rewires (the view's flip_edge has no status check between the hook and the
// first write): the commit hook then refuses and the complex stays poisoned,
// per the file's deep-invariant convention.
//
// BOUNDEDNESS.  No allocation, no exceptions, no std::function on the run
// path; the two capacity formulas below size the tracked points and the
// plan-to-commit carry.
// ============================================================================

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <limits>
#include <span>

#include "fullerenes/delaunay_geometry.hh"   // FanPolygon, FanTriangulation
#include "fullerenes/delaunay_view.hh"       // DelaunayView, NoTransport

// ============================================================================
// Capacities.  Removal seeds exactly one tracked point per removed flat
// vertex and nothing else creates them, so the pre-removal vertex count
// bounds the points; the carry additionally holds the star's seeded apex.
// ============================================================================
constexpr long tracked_points_cap(int nv0) { return nv0; }
constexpr long transport_carry_cap(long points_cap) { return points_cap + 1; }

// One tracked point as a value record: the caller's label, the containing
// live face, and the barycentric triple over that face's half-edge SLOTS
// (slot i of face f is the origin of the i-th half-edge of the cycle
// starting at f_he[f]).  The tracker stores this SoA (see the view below);
// this is the projection consumers read.
struct DelaunayTrackedPoint {
  int    label = -1;
  int    face  = -1;
  double b[3]  = {0, 0, 0};       // b >= 0, sum == 1
};

// ============================================================================
// DelaunayPointTrackerView: the tracker's state as SoA spans over caller
// storage -- the points (label / face / barycentric) plus the by-face
// partition as intrusive singly-linked lists (head/tail per face, next per
// point), appended at the tail so bucket order is insertion order.
//
// @inv partition: the bucket of face f is exactly { i : face[i] == f } --
//      every live point indexed exactly once, under its own live face.
//      Holds after every hooked operation.
// Iteration idiom:  for (int32_t i = head[f]; i >= 0; i = next[i]) ...
// ============================================================================
struct DelaunayPointTrackerView {
  std::span<int32_t> label;    // [P]
  std::span<int32_t> face;     // [P]
  std::span<double>  bary;     // [3P] barycentric per slot
  std::span<int32_t> head;     // [NF] bucket first point (-1 = empty)
  std::span<int32_t> tail;     // [NF] bucket last point
  std::span<int32_t> next;     // [P]  successor within the bucket (-1 = end)
  int32_t n         = 0;       // live point count
  int32_t n_clamped = 0;       // clamp accounting (see clamp_barycentric)
  double  max_clamp = 0;

  // FP policy: after transport, a barycentric in [-CLAMP_TOL, 0) is clamped
  // to 0 and the triple renormalized (accounted above); below -CLAMP_TOL the
  // transport refuses -- that is a wrong-side transport bug, not roundoff.
  // Never widen the band to make a case pass.
  static constexpr double CLAMP_TOL = 1e-9;

  int32_t cap() const { return (int32_t)label.size(); }

  // The empty tracker (no points, every bucket empty).  EVERY bucket slot is
  // cleared, not just the live prefix: a surgery may take a face slot the
  // caller's counts did not cover, and a stale head there would silently
  // splice a dead list onto a live bucket.
  void reset() {
    n = 0; n_clamped = 0; max_clamp = 0;
    for (std::size_t f = 0; f < head.size(); f++) { head[f] = -1; tail[f] = -1; }
  }

  // @pre 0 <= f < head.size() -- callers own the range check (the transport
  //      trips its Status latch; see TrackerTransport::bucket_face_ok).
  // @pre bucket_push: i is not already in bucket f.  A repeated push would
  //      set next[i] = i (a self-cycle) and every later walk of the bucket
  //      would not terminate; the transport refuses the repeat upstream, at
  //      the one place it can arise (a plan naming one face slot twice).
  void bucket_clear(int f) { head[f] = -1; tail[f] = -1; }
  void bucket_push(int f, int32_t i) {
    next[i] = -1;
    if (head[f] < 0) head[f] = i; else next[tail[f]] = i;
    tail[f] = i;
  }
  const double* b_of(int32_t i) const { return &bary[(std::size_t)3 * i]; }
  double*       b_of(int32_t i)       { return &bary[(std::size_t)3 * i]; }

  DelaunayTrackedPoint point(int32_t i) const {
    return { label[i], face[i], { b_of(i)[0], b_of(i)[1], b_of(i)[2] } };
  }

  // *this := src (points, buckets, counters) -- the tracker half of "copies
  // of a tracking complex snapshot the tracker" (the owner's copy does this
  // by copying its vectors; a device trial copy calls this).  Buckets are
  // copied verbatim, so the snapshot's iteration order is the source's.
  void copy_from(const DelaunayPointTrackerView& src) {
    n = src.n; n_clamped = src.n_clamped; max_clamp = src.max_clamp;
    for (int32_t i = 0; i < src.n; i++) {
      label[i] = src.label[i];
      face[i]  = src.face[i];
      next[i]  = src.next[i];
      for (int k = 0; k < 3; k++)
        bary[(std::size_t)3 * i + k] = src.bary[(std::size_t)3 * i + k];
    }
    for (std::size_t f = 0; f < src.head.size(); f++) {
      head[f] = src.head[f]; tail[f] = src.tail[f];
    }
  }
};

// The clamp policy (see CLAMP_TOL): a barycentric in [-CLAMP_TOL, 0) clamps
// to 0 with accounting and the triple renormalizes; below -CLAMP_TOL, or any
// non-finite coordinate (a degenerate development), the triple is REFUSED
// (false) -- the caller supplies the failure channel.
inline bool clamp_barycentric(DelaunayPointTrackerView& tk, double b[3]) {
  if (!std::isfinite(b[0]) || !std::isfinite(b[1]) || !std::isfinite(b[2]))
    return false;
  double neg = 0;
  for (int i = 0; i < 3; i++) if (b[i] < neg) neg = b[i];
  if (neg < -DelaunayPointTrackerView::CLAMP_TOL) return false;
  if (neg < 0) {
    tk.n_clamped++;
    if (-neg > tk.max_clamp) tk.max_clamp = -neg;
  }
  double s = 0;
  for (int i = 0; i < 3; i++) { if (b[i] < 0) b[i] = 0; s += b[i]; }
  for (int i = 0; i < 3; i++) b[i] /= s;
  return true;
}

// ============================================================================
// TransportCarry: the plan hooks' scratch, drained by the commit hooks.
// cidx == -1 marks the point being SEEDED (the removed flat vertex at the
// development's origin).
// ============================================================================
struct TransportCarry {
  std::span<int32_t> cidx;    // [C] tracked-point index, -1 = seeded
  std::span<double>  cx, cy;  // [C] development coordinates
  std::span<int32_t> rcand;   // [C] winning candidate chart
  std::span<double>  rb;      // [3C] barycentric in that chart
  int32_t cap() const { return (int32_t)cidx.size(); }
};

namespace transport_detail {

// One post-operation chart: three corner positions in the operation's
// isometric planar development, in the slot order the face's half-edge cycle
// will have.  The face id is resolved at commit time by the caller.
struct DevelopedTriangle { double x0, y0, x1, y1, x2, y2; };

// The flip's two post-charts: (B, D, v) -> fh and (D, B, u) -> ft.
struct FlipCharts {
  DevelopedTriangle c[2];
  int               size() const { return 2; }
  DevelopedTriangle get(int i) const { return c[i]; }
};

// The star's post-charts: the ear triangulation read straight off
// (fan, tri) -- the developed corners of each new triangle.
struct StarCharts {
  const FanPolygon*       fan;
  const FanTriangulation* tri;
  int                     size() const { return tri->n_triangles; }
  DevelopedTriangle get(int i) const {
    const FanTriangulation::Triangle& t = tri->triangles[i];
    return {fan->x(t.v0), fan->y(t.v0), fan->x(t.v1), fan->y(t.v1),
            fan->x(t.v2), fan->y(t.v2)};
  }
};

// Barycentric coordinates of p in the planar triangle (A, B, C), CCW.
// Returns the minimum coordinate (negative iff p lies outside).
inline double planar_barycentric(double px, double py,
                                 double ax, double ay,
                                 double bx, double by,
                                 double cx, double cy,
                                 double out[3]) {
  const double d = (bx - ax) * (cy - ay) - (cx - ax) * (by - ay);  // 2*area > 0 (CCW)
  out[0] = ((bx - px) * (cy - py) - (cx - px) * (by - py)) / d;
  out[1] = ((cx - px) * (ay - py) - (ax - px) * (cy - py)) / d;
  out[2] = 1.0 - out[0] - out[1];
  return std::min({out[0], out[1], out[2]});
}

}  // namespace transport_detail

// ============================================================================
// TrackerTransport -- THE Transport policy carrying the point tracker
// through the surgery bodies (see the banner).  ONE object per surgery
// call: the plan hook fills the carry, the commit hook drains it.
// ============================================================================
struct TrackerTransport {
  DelaunayView&             T;
  DelaunayPointTrackerView& tk;
  TransportCarry            carry;
  int32_t                   n_carry = 0;
  // The star protocol (plan_star .. commit_removal) holds the carry across
  // the surgery; a flip planned INSIDE that window would clobber it.  The
  // surgery bodies never flip (asserted here).
  bool                      star_in_flight = false;

  // The owner instantiates this policy only for tracked runs, so the trait
  // is compile-time true (the single gate is the owner's selection).
  static constexpr bool tracking() { return true; }

  bool trip(const char* site, int witness) {
    return T.trip(DelaunayView::Status::InvariantViolated, site, witness);
  }

  bool clamp(double b[3], const char* who, int witness) {
    return clamp_barycentric(tk, b) ? true : trip(who, witness);
  }

  // Every bucket WRITE is range-checked against the tracker's own bucket
  // count.  The READ path (develop_bucket) guards already; an unchecked
  // write past the span would corrupt whatever follows it -- on device, the
  // neighbouring isomer's tracker slab -- instead of refusing, which is what
  // the boundedness contract promises.
  bool bucket_face_ok(int f, const char* who) {
    if (f >= 0 && f < (int)tk.head.size()) return true;
    trip(who, f);
    return false;
  }

  bool push_carry(int32_t idx, double x, double y) {
    if (n_carry >= carry.cap())
      return T.trip(DelaunayView::Status::CapacityExceeded,
                    "TrackerTransport: transport carry", n_carry);
    carry.cidx[n_carry] = idx; carry.cx[n_carry] = x; carry.cy[n_carry] = y;
    n_carry++;
    return true;
  }

  // Express every tracked point of face f in the operation's development:
  // the barycentric combination of f's slot corners in the CURRENT cycle
  // order, where `corner` maps each half-edge of the cycle to its origin's
  // development position.  Appends to the carry; reads the tracker only.
  template <class CornerFn>
  void develop_bucket(int f, CornerFn corner) {
    if (f < 0 || f >= (int)tk.head.size()) return;
    if (tk.head[f] < 0) return;
    const std::array<int, 3> h = T.face_halfedges(f);
    double x0, y0, x1, y1, x2, y2;
    corner(h[0], x0, y0); corner(h[1], x1, y1); corner(h[2], x2, y2);
    for (int32_t i = tk.head[f]; i >= 0; i = tk.next[i]) {
      const double* b = tk.b_of(i);
      if (!push_carry(i, b[0] * x0 + b[1] * x1 + b[2] * x2,
                         b[0] * y0 + b[1] * y1 + b[2] * y2)) return;
    }
  }

  // Re-express every carried point in the post-operation charts: choose the
  // candidate triangle where the point's minimum barycentric is largest,
  // then clamp.  Deterministic on shared edges -- the winning chart depends
  // only on the development coordinates, never on which source face supplied
  // the point -- and boundary points cannot fall through.  Pure computation
  // over the development: a refusal here fires BEFORE the caller has mutated
  // anything.  Every candidate scoring NaN leaves no winner and refuses.
  template <class Charts>
  void relocate_all(const Charts& ch, const char* who) {
    const int nc = ch.size();
    for (int32_t i = 0; i < n_carry; i++) {
      const double px = carry.cx[i], py = carry.cy[i];
      int    cand   = -1;
      double best_m = -std::numeric_limits<double>::infinity();
      double bb[3]  = {0, 0, 0};
      for (int c = 0; c < nc; c++) {
        const transport_detail::DevelopedTriangle t = ch.get(c);
        double b[3];
        const double m = transport_detail::planar_barycentric(
            px, py, t.x0, t.y0, t.x1, t.y1, t.x2, t.y2, b);
        if (m > best_m) {
          best_m = m; cand = c;
          bb[0] = b[0]; bb[1] = b[1]; bb[2] = b[2];
        }
      }
      if (cand < 0) { trip(who, (int)i); return; }
      if (!clamp(bb, who, (int)i)) return;
      carry.rcand[i] = cand;
      for (int k = 0; k < 3; k++) carry.rb[(std::size_t)3 * i + k] = bb[k];
    }
  }

  // Commit the precomputed relocations: point idx >= 0 moves to its new
  // chart; idx == -1 seeds a new tracked point with label seed_label.
  // face_of maps each candidate chart to its post-operation face id.  Cannot
  // fail except on tracked-point capacity: every relocation was verified
  // before any mutation.
  template <class FaceOf>
  void commit(FaceOf face_of, int seed_label) {
    for (int32_t i = 0; i < n_carry; i++) {
      const int f = face_of(carry.rcand[i]);
      int32_t idx = carry.cidx[i];
      if (idx < 0) {
        if (tk.n >= tk.cap()) {
          T.trip(DelaunayView::Status::CapacityExceeded,
                 "TrackerTransport: tracked points", tk.n);
          return;
        }
        idx = tk.n++;
        tk.label[idx] = seed_label;
      }
      if (!bucket_face_ok(f, "TrackerTransport: commit face out of range"))
        return;
      tk.face[idx] = f;
      for (int k = 0; k < 3; k++) tk.b_of(idx)[k] = carry.rb[(std::size_t)3 * i + k];
      tk.bucket_push(f, idx);
    }
  }

  // ── The five hooks. ──

  void plan_flip(int h, int fh, int ft) {
    if (star_in_flight) {
      trip("TrackerTransport: flip planned inside a star protocol", h);
      return;
    }
    n_carry = 0;
    // The two diamond faces must be distinct slots: one slot named twice
    // would carry its points twice and commit would push each index into
    // that bucket twice (a self-cycling list).  Refuse rather than guess
    // which of the two developments the shared face should take.
    if (fh == ft) {
      trip("TrackerTransport: flip diamond names one face slot twice", h);
      return;
    }
    const int t  = h ^ 1;
    const int h1 = T.he_next[h], h2 = T.he_next[h1];
    const int h4 = T.he_next[t], h5 = T.he_next[h4];
    const double e_dev = T.he_length[h];
    const double lvB = T.he_length[h1], lBu = T.he_length[h2];
    const double luD = T.he_length[h4], lDv = T.he_length[h5];
    const double xB = (e_dev * e_dev + lBu * lBu - lvB * lvB) / (2 * e_dev);
    const double yB =  std::sqrt(std::max(0.0, lBu * lBu - xB * xB));
    const double xD = (e_dev * e_dev + luD * luD - lDv * lDv) / (2 * e_dev);
    const double yD = -std::sqrt(std::max(0.0, luD * luD - xD * xD));
    auto corner = [&](int he, double& x, double& y) {
      if      (he == h)  { x = 0;     y = 0;  }
      else if (he == h1) { x = e_dev; y = 0;  }
      else if (he == h2) { x = xB;    y = yB; }
      else if (he == t)  { x = e_dev; y = 0;  }
      else if (he == h4) { x = 0;     y = 0;  }
      else               { x = xD;    y = yD; }   // h5
    };
    develop_bucket(fh, corner);
    develop_bucket(ft, corner);
    if (n_carry == 0) return;
    const transport_detail::FlipCharts ch{{{xB, yB, xD, yD, e_dev, 0.0},
                                           {xD, yD, xB, yB, 0.0,   0.0}}};
    relocate_all(ch, "TrackerTransport: flip_edge transport");
    // Clear the source buckets only once every re-expression is verified: on
    // a refusal the tracker is left untouched.
    if (T.status != DelaunayView::Status::Ok) return;
    if (!bucket_face_ok(fh, "TrackerTransport: flip face out of range") ||
        !bucket_face_ok(ft, "TrackerTransport: flip face out of range")) return;
    tk.bucket_clear(fh);
    tk.bucket_clear(ft);
  }

  void commit_flip(int fh, int ft) {
    // f_he[fh] = h and f_he[ft] = t were set by the rewire, so the chart
    // slot orders (B, D, v) / (D, B, u) match the anchoring convention.
    // A non-empty carry IS the planned-flip fact.
    if (n_carry == 0) return;
    if (T.status != DelaunayView::Status::Ok) { n_carry = 0; return; }
    const int fc[2] = {fh, ft};
    commit([&](int c) { return fc[c]; }, /*seed_label*/ -1);
    n_carry = 0;
  }

  void plan_star(int /*v*/, const FanPolygon& fan) {
    star_in_flight = true;
    n_carry = 0;
    if (!push_carry(-1, 0.0, 0.0)) return;              // the seeded apex
    for (int i = 0; i < fan.k; i++) {
      // A delta-complex fan can traverse ONE face slot twice (a face with the
      // apex at two corners).  Its points would then carry twice, at two
      // different planar placements, and commit would self-cycle the bucket.
      // Refuse: the development is genuinely ambiguous, so there is nothing
      // to pick.  (The pre-span tracker duplicated the entry instead.)
      const int fi = T.he_face[fan.spoke_he[i]];
      for (int p = 0; p < i; p++)
        if (T.he_face[fan.spoke_he[p]] == fi) {
          trip("TrackerTransport: star names one face slot twice", fi);
          return;
        }
      const int j = (i + 1) % fan.k;
      auto corner = [&](int he, double& x, double& y) {
        if      (he == fan.spoke_he[i])  { x = 0;        y = 0; }
        else if (he == fan.inner_rim[i]) { x = fan.x(i); y = fan.y(i); }
        else                             { x = fan.x(j); y = fan.y(j); }
      };
      develop_bucket(fi, corner);
    }
  }

  void plan_charts(const FanPolygon& fan, const FanTriangulation& tri) {
    relocate_all(transport_detail::StarCharts{&fan, &tri},
                 "TrackerTransport: remove_flat_vertex transport");
    // All re-expressions verified: clear the source buckets before the
    // surgery recycles their face slots.
    if (T.status != DelaunayView::Status::Ok) return;
    for (int i = 0; i < fan.k; i++) {
      const int fi = T.he_face[fan.spoke_he[i]];
      if (!bucket_face_ok(fi, "TrackerTransport: star face out of range")) return;
      tk.bucket_clear(fi);
    }
  }

  void commit_removal(std::span<const int> faces, int v) {
    if (T.status != DelaunayView::Status::Ok) {
      n_carry = 0; star_in_flight = false; return;
    }
    commit([&](int c) { return faces[(std::size_t)c]; }, /*seed_label*/ v);
    n_carry = 0;
    star_in_flight = false;
  }
};
