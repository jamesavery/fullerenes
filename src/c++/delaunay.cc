#include "fullerenes/delaunay.hh"

#include <cmath>
#include <algorithm>
#include <map>
#include <climits>
#include <limits>
#include <queue>
#include <stdexcept>
#include <string>
#include <sstream>
#include <cstdio>

// The intrinsic geometry primitives and the Diamond predicates are
// header-inline in delaunay_view.hh (delaunay_detail::).

// Old adjacency-list FulleroidDelaunay + IDTAudit implementation retired to
// src/c++/attic/delaunay_old.cc.attic (superseded by the DCEL
// DelaunayTriangulation below; zero consumers as of 2026-07-09).

// ============================================================================
// DelaunayTriangulation — DCEL-based iDT (delta-complex)
// (The Rule of 5 is header-inline: one-liners over the DelaunayStorage base.)
// ============================================================================

// --- Allocation ---
// (The capacity-guarded allocators live on DelaunayView; the owner's
// growth-shadowing versions are header-inline in delaunay.hh.)

void DelaunayTriangulation::throw_on_status(const char* op) const
{
  if (status == Status::Ok) return;
  const char* kind =
      status == Status::BudgetExceeded   ? "budget exhausted (fail-loud backstop)"
    : status == Status::CapacityExceeded ? "workspace capacity exceeded"
                                         : "DCEL invariant violated";
  throw std::runtime_error(std::string(op) + ": " + kind + " at " +
                           (status_site ? status_site : "?") +
                           " (witness=" + std::to_string(status_witness) + ")");
}

// ============================================================================
// Weighted graph algorithms on the 1-skeleton (he_length as edge weight)
// ============================================================================

std::vector<double>
DelaunayTriangulation::single_source_shortest_paths(int src) const
{
  const double kInf = std::numeric_limits<double>::infinity();
  std::vector<double> dist(nv, kInf);
  if (src < 0 || src >= nv || v_out[src] < 0) return dist;
  dist[src] = 0.0;

  using Item = std::pair<double, int>;     // (dist, vertex)
  std::priority_queue<Item, std::vector<Item>, std::greater<Item>> pq;
  pq.push({0.0, src});

  while (!pq.empty()) {
    auto [d, u] = pq.top(); pq.pop();
    if (d > dist[u]) continue;             // stale entry
    // incident(u) is empty when v_out[u] < 0, subsuming the old h0<0 guard.
    for (int h : incident(u)) {
      int    v  = dest(h);
      double nd = d + he_length[h];
      if (nd < dist[v]) { dist[v] = nd; pq.push({nd, v}); }
    }
  }
  return dist;
}

double DelaunayTriangulation::diameter_upper_bound() const
{
  if (nv == 0) return 0.0;
  // Pick the first live vertex as the seed.
  int seed = -1;
  for (int v = 0; v < nv; v++) if (v_out[v] >= 0) { seed = v; break; }
  if (seed < 0) return 0.0;

  auto farthest = [&](int src) {
    std::vector<double> d = single_source_shortest_paths(src);
    int    far = src;
    double m   = 0.0;
    for (int v = 0; v < nv; v++)
      if (std::isfinite(d[v]) && d[v] > m) { m = d[v]; far = v; }
    return std::pair<int, double>{far, m};
  };

  auto [u1, _ ] = farthest(seed);
  auto [u2, m1] = farthest(u1);
  // Double-sweep returns m1 with m1 >= D/2 (lower bound on true diameter D);
  // 2 * m1 is therefore a safe upper bound on D.
  return 2.0 * m1;
}

// --- Construction ---

// Owner storage for the view-level build (delaunay_view.hh carries the
// algorithm), sized at dcel_build_capacities(degree_sum(T)) -- the input's
// own sizes, byte-identical to the historical grow-as-discovered build for
// every closed oriented input, connected genus-0 or not.  ensure_* pre-fills
// the dead representation (its growth contract) and the build's reset
// re-establishes it: the double pass is deliberate, keeping the two
// contracts independent.
static DelaunayTriangulation owner_sized_for_build(const Triangulation& T)
{
  const DcelCapacities c = dcel_build_capacities(degree_sum(T));
  DelaunayTriangulation D;
  D.ensure_vertices(T.N);
  D.ensure_halfedges((int)c.nh_cap);
  D.ensure_faces((int)c.nf_cap);
  return D;
}

DelaunayTriangulation DelaunayTriangulation::from_triangulation(const Triangulation& T)
{
  // Equilateral metric: every edge length 1, every angle pi/3, cone angle
  // deg*pi/3. (The special case length == 1 of from_intrinsic_metric, kept
  // explicit because the constants are exact and this is the hot path.)
  DelaunayTriangulation D = owner_sized_for_build(T);
  std::vector<int> he_of_arc(arc_space(T));
  DelaunayView::from_triangulation(D, T, he_of_arc);
  D.throw_on_status("from_triangulation");
  return D;
}

DelaunayTriangulation
DelaunayTriangulation::from_intrinsic_metric(const Triangulation& T,
                                             const EdgeLengthFn& length)
{
  // Prescribed intrinsic metric: edge lengths from `length`, angles and cone
  // angles derived. The topology is identical to the equilateral case; only
  // the metric differs.
  DelaunayTriangulation D = owner_sized_for_build(T);
  std::vector<int> he_of_arc(arc_space(T));
  D.build_topology(T, he_of_arc);
  D.throw_on_status("from_intrinsic_metric");
  for (int h = 0; h < D.nh; h += 2) {
    double L = length(D.he_origin[h], D.he_origin[D.twin(h)]);
    D.he_length[h] = D.he_length[D.twin(h)] = L;
  }
  D.recompute_all_angles();   // fills he_angle AND refreshes the v_cone_angle cache
  return D;
}

// (The intrinsic-geometry methods -- diamond, recompute_face_angles /
// recompute_all_angles / recompute_all_face_angles / recompute_cone_angles,
// vertex_angle_sum, is_flat -- are header-inline on DelaunayView.)

// ============================================================================
// Point transport (flip-tape) helpers
// See claude-projects/delaunay-fillin/DESIGN-cubic-exact-paint.md and the
// tracker banner in delaunay.hh.  Transport happens inside the two
// topology-changing operations (flip_edge, remove_flat_vertex); everything
// here is planar geometry over the operations' own isometric developments.
// ============================================================================

namespace {

// Barycentric coordinates of p in the planar triangle (A, B, C), CCW.
// Returns the minimum coordinate (negative iff p lies outside).
double planar_barycentric(double px, double py,
                          double ax, double ay,
                          double bx, double by,
                          double cx, double cy,
                          double out[3])
{
  const double d = (bx - ax) * (cy - ay) - (cx - ax) * (by - ay);  // 2*area > 0 (CCW)
  out[0] = ((bx - px) * (cy - py) - (cx - px) * (by - py)) / d;
  out[1] = ((cx - px) * (ay - py) - (ax - px) * (cy - py)) / d;
  out[2] = 1.0 - out[0] - out[1];
  return std::min({out[0], out[1], out[2]});
}

// Apply the clamp policy (tracker banner): [-CLAMP_TOL, 0) clamps to 0 with
// accounting; below -CLAMP_TOL, or any non-finite coordinate (a degenerate
// development), throws (wrong-side transport bug, not noise).
void clamp_barycentric(DelaunayTriangulation::PointTracker& tk,
                       double b[3], const char* who)
{
  if (!std::isfinite(b[0]) || !std::isfinite(b[1]) || !std::isfinite(b[2]))
    throw std::runtime_error(std::string(who) +
        ": non-finite tracked-point barycentric (degenerate development)");
  double neg = 0;
  for (int i = 0; i < 3; i++) if (b[i] < neg) neg = b[i];
  if (neg < -DelaunayTriangulation::PointTracker::CLAMP_TOL)
    throw std::runtime_error(std::string(who) +
        ": tracked-point barycentric " + std::to_string(neg) +
        " below -CLAMP_TOL (transport bug, not roundoff)");
  if (neg < 0) {
    tk.n_clamped++;
    if (-neg > tk.max_clamp) tk.max_clamp = -neg;
  }
  double s = 0;
  for (int i = 0; i < 3; i++) { if (b[i] < 0) b[i] = 0; s += b[i]; }
  for (int i = 0; i < 3; i++) b[i] /= s;
}

// One tracked point EXPRESSED IN THE DEVELOPMENT: its index in
// tracker.points and its coordinates in the operation's isometric planar
// development.  idx == -1 marks a point being seeded (the removed flat
// vertex itself, at the development's origin).
struct DevelopedPoint { int idx; double x, y; };

// The development of one POST-operation face: its three corner positions,
// in the slot order the face's half-edge cycle will have.  The face id is
// resolved at commit time by the caller.
struct DevelopedTriangle { double x0, y0, x1, y1, x2, y2; };

// The re-expression of a developed point in the post-operation charts:
// which candidate triangle hosts it, and its barycentric coordinates there.
struct Relocation { int cand; double b[3]; };

// Re-express one developed point: choose the candidate triangle where the
// point's minimum barycentric is largest, then clamp.  Deterministic on
// shared edges -- the winning chart depends only on the development
// coordinates, never on which source face supplied the point -- and
// boundary points cannot fall through.  Pure computation over the
// development: a clamp throw here fires BEFORE the caller has mutated
// anything (the commit-or-nothing property of both transport hooks).
Relocation relocate(DelaunayTriangulation::PointTracker& tk,
                    double px, double py,
                    const std::vector<DevelopedTriangle>& cands, const char* who)
{
  Relocation out{-1, {0, 0, 0}};
  double best_m = -std::numeric_limits<double>::infinity();
  for (size_t i = 0; i < cands.size(); i++) {
    const DevelopedTriangle& t = cands[i];
    double b[3];
    const double m = planar_barycentric(px, py, t.x0, t.y0, t.x1, t.y1,
                                        t.x2, t.y2, b);
    if (m > best_m) {
      best_m = m; out.cand = (int)i;
      out.b[0] = b[0]; out.b[1] = b[1]; out.b[2] = b[2];
    }
  }
  clamp_barycentric(tk, out.b, who);
  return out;
}

// Express every tracked point of face f in the operation's development:
// the barycentric combination of f's slot corners (current cycle order),
// where `corner` maps each half-edge of the cycle to its origin's
// development position.  Appends to `out`; reads the tracker only.
template <class CornerFn>
void develop_bucket(const DelaunayTriangulation& D, int f, CornerFn corner,
                    std::vector<DevelopedPoint>& out)
{
  if ((int)D.tracker.by_face.size() <= f) return;
  const std::vector<int>& bucket = D.tracker.by_face[f];
  if (bucket.empty()) return;
  const auto h = D.face_halfedges(f);
  double x0, y0, x1, y1, x2, y2;
  corner(h[0], x0, y0); corner(h[1], x1, y1); corner(h[2], x2, y2);
  for (int idx : bucket) {
    const DelaunayTriangulation::TrackedPoint& p = D.tracker.points[idx];
    out.push_back({idx, p.b[0]*x0 + p.b[1]*x1 + p.b[2]*x2,
                        p.b[0]*y0 + p.b[1]*y1 + p.b[2]*y2});
  }
}

// Commit precomputed relocations: point idx >= 0 moves to its new chart;
// idx == -1 seeds a new tracked point with label seed_label.  face_of_cand
// maps each candidate triangle to its (post-operation) face id.  Cannot
// fail: every Relocation was verified by relocate() before any mutation.
void commit_relocations(DelaunayTriangulation& D,
                        const std::vector<DevelopedPoint>& pts,
                        const std::vector<Relocation>& rel,
                        const std::vector<int>& face_of_cand,
                        int seed_label)
{
  for (size_t i = 0; i < pts.size(); i++) {
    const int f = face_of_cand[rel[i].cand];
    if (pts[i].idx < 0) {
      const int idx = (int)D.tracker.points.size();
      D.tracker.points.push_back({seed_label, f,
          {rel[i].b[0], rel[i].b[1], rel[i].b[2]}});
      D.tracker.bucket(f).push_back(idx);
    } else {
      DelaunayTriangulation::TrackedPoint& p = D.tracker.points[pts[i].idx];
      p.face = f;
      p.b[0] = rel[i].b[0]; p.b[1] = rel[i].b[1]; p.b[2] = rel[i].b[2];
      D.tracker.bucket(f).push_back(pts[i].idx);
    }
  }
}

}  // namespace

std::vector<int>& DelaunayPointTracker::bucket(int f)
{
  if ((int)by_face.size() <= f) by_face.resize(f + 1);
  return by_face[f];
}

void DelaunayTriangulation::enable_point_tracking()
{
  if (!tracker.points.empty())
    throw std::runtime_error(
        "enable_point_tracking: tracker already carries points "
        "(one tracking session per complex; stale state would be transported)");
  tracker.active = true;
  if ((int)tracker.by_face.size() < nf) tracker.by_face.resize(nf);
}

int DelaunayTriangulation::track_point(int label, int face,
                                       double b0, double b1, double b2)
{
  if (!tracker.active)
    throw std::runtime_error("track_point: tracker not active (enable_point_tracking first)");
  if (face < 0 || face >= nf || f_he[face] < 0)
    throw std::runtime_error("track_point: face " + std::to_string(face) + " is not live");
  if (!std::isfinite(b0) || !std::isfinite(b1) || !std::isfinite(b2) ||
      std::fabs(b0 + b1 + b2 - 1.0) > 1e-9)
    throw std::runtime_error("track_point: barycentrics must be finite and sum to 1");
  double b[3] = {b0, b1, b2};
  clamp_barycentric(tracker, b, "track_point");
  const int idx = (int)tracker.points.size();
  tracker.points.push_back({label, face, {b[0], b[1], b[2]}});
  tracker.bucket(face).push_back(idx);
  return idx;
}

// ============================================================================
// TrackerTransport: the host-side transport policy.  Reproduces the
// flip-tape transport (DESIGN-cubic-exact-paint.md 4.1/4.2) through the
// view's surgery bodies: plan hooks develop and re-express every tracked
// point BEFORE any mutation (a clamp throw fires with complex and tracker
// untouched -- commit-or-nothing), commit hooks store the precomputed
// relocations after the surgery.
//
// The mathematics: for a surgery S re-triangulating a flat region R, the
// transport is the composition
//     T_S  =  express_post^-1  o  express_pre
// where express_pre maps (face, barycentric) -> R^2 through the
// PRE-surgery triangles' isometric development of R, and express_post^-1
// locates the planar point in the POST-surgery triangles of the SAME
// development.  Both chart sets develop R isometrically, so T_S fixes the
// intrinsic surface point (the transport lemma).  Both halves are computed
// before the first write -- the post-charts are determined by the planned
// surgery -- which is what makes the failure contract commit-or-nothing.
// (The exact-transport rework will type this composition over lattice
// charts, unifying it with the atlas's point-location vocabulary.)
//
// Flip development: u = (0,0), v = (e,0), B above, D below; positions keyed
// per HALF-EDGE (slot anchoring), so repeated-corner faces transport the
// same.  The post-flip charts are (B, D, v) -> fh and (D, B, u) -> ft.
// Star development: the FanPolygon IS the isometric development (apex at
// the origin); the removed vertex seeds a tracked point with label = v.
// @post (transport) each point's intrinsic surface position is unchanged:
//       the developments are isometries, re-expression only changes
//       coordinates.
// ============================================================================
namespace {

struct TrackerTransport {
  DelaunayTriangulation& T;
  std::vector<DevelopedPoint>    carry;
  std::vector<Relocation>        rel;
  std::vector<DevelopedTriangle> charts;
  std::vector<std::array<int,3>> new_tris;
  // The star protocol (plan_star .. commit_removal) holds carry/rel across
  // the surgery; a flip planned INSIDE that window would clobber them.  The
  // surgery bodies never flip (asserted here), so the two protocols can
  // share the buffers.
  bool star_in_flight = false;

  // The owner instantiates this policy only when the tracker is active, so
  // the trait is compile-time true (the single gate is the owner's
  // selection).
  static constexpr bool tracking() { return true; }

  void plan_flip(int h, int fh, int ft) {
    if (star_in_flight)
      throw std::runtime_error("TrackerTransport: flip planned inside a star "
                               "protocol (surgery bodies must not flip)");
    carry.clear(); rel.clear();
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
    develop_bucket(T, fh, corner, carry);
    develop_bucket(T, ft, corner, carry);
    if (carry.empty()) return;
    const std::vector<DevelopedTriangle> fcharts = {
      {xB, yB, xD, yD, e_dev, 0.0},   // cand 0 -> face fh (B, D, v)
      {xD, yD, xB, yB, 0.0,   0.0}};  // cand 1 -> face ft (D, B, u)
    rel.reserve(carry.size());
    for (const DevelopedPoint& c : carry)
      rel.push_back(relocate(T.tracker, c.x, c.y, fcharts, "flip_edge transport"));
    T.tracker.bucket(fh).clear();
    T.tracker.bucket(ft).clear();
  }

  void commit_flip(int fh, int ft) {
    // f_he[fh] = h and f_he[ft] = t were set by the rewire, so the chart
    // slot orders (B, D, v) / (D, B, u) match the anchoring convention.
    // A non-empty carry IS the planned-flip fact (HEAD's own condition).
    if (carry.empty()) return;
    commit_relocations(T, carry, rel, {fh, ft}, /*seed_label*/ -1);
    carry.clear(); rel.clear();
  }

  void plan_star(int /*v*/, const FanPolygon& fan) {
    star_in_flight = true;
    carry.clear(); rel.clear(); new_tris.clear(); charts.clear();
    carry.push_back({-1, 0.0, 0.0});                    // the seeded apex
    for (int i = 0; i < fan.k; i++) {
      const int j = (i + 1) % fan.k;
      auto corner = [&](int he, double& x, double& y) {
        if      (he == fan.spoke_he[i])  { x = 0;        y = 0; }
        else if (he == fan.inner_rim[i]) { x = fan.x(i); y = fan.y(i); }
        else                             { x = fan.x(j); y = fan.y(j); }
      };
      develop_bucket(T, T.he_face[fan.spoke_he[i]], corner, carry);
    }
  }

  void plan_charts(const FanPolygon& fan, const FanTriangulation& tri) {
    for (int ti = 0; ti < tri.n_triangles; ti++)
      new_tris.push_back({tri.triangles[ti].v0, tri.triangles[ti].v1,
                          tri.triangles[ti].v2});
    charts.reserve(new_tris.size());
    for (const auto& tv : new_tris)
      charts.push_back({fan.x(tv[0]), fan.y(tv[0]),
                        fan.x(tv[1]), fan.y(tv[1]),
                        fan.x(tv[2]), fan.y(tv[2])});
    rel.reserve(carry.size());
    for (const DevelopedPoint& c : carry)
      rel.push_back(relocate(T.tracker, c.x, c.y, charts,
                             "remove_flat_vertex transport"));
    // All re-expressions verified: clear the source buckets before the
    // surgery recycles their face slots.
    for (int i = 0; i < fan.k; i++)
      T.tracker.bucket(T.he_face[fan.spoke_he[i]]).clear();
  }

  void commit_removal(std::span<const int> faces, int v) {
    std::vector<int> f(faces.begin(), faces.end());
    commit_relocations(T, carry, rel, f, /*seed_label*/ v);
    carry.clear(); rel.clear();
    star_in_flight = false;
  }
};

}  // namespace

// --- Delaunay operations: owner wrappers over the view machinery ---

bool DelaunayTriangulation::flip_edge(int h)
{
  // Owner wrappers serve the general (float-metric) API surface -- weighted
  // flips, post-reduction surgery -- so they instantiate the banded policy.
  bool r;
  if (tracker.active) {
    TrackerTransport tr{*this};
    r = DelaunayView::flip_edge(h, BandedFloatMetric{}, tr);
  } else {
    r = DelaunayView::flip_edge(h);
  }
  throw_on_status("flip_edge");
  return r;
}

int DelaunayTriangulation::lawson_sweep()
{
  // Sweep-only workspace: k_max = 0 (no fan machinery), O(nh) bytes.
  HostDelaunayWorkspace ws({.nv0 = nv, .k_max = 0, .nh_explicit = nh});
  ws.sweep_in_stack.clear_all();
  int flips;
  if (tracker.active) {
    TrackerTransport tr{*this};
    flips = DelaunayView::flip_to_delaunay(ws.lawson_stack, ws.sweep_in_stack,
                                           BandedFloatMetric{}, tr);
  } else {
    flips = DelaunayView::flip_to_delaunay(ws.lawson_stack, ws.sweep_in_stack);
  }
  throw_on_status("lawson_sweep");
  return flips;
}

int DelaunayTriangulation::flip_to_delaunay()
{
  return lawson_sweep();
}

// ============================================================================
// Vertex removal
// (FanPolygon / FanTriangulation and the whole removal machinery --
// extract_fan, ear clipping, splice_fan, the cocircular tie-break,
// flip_away_self_loops, and the remove_flat_vertices driver -- are
// header-inline on DelaunayView over DelaunayWorkspace; the owner
// wrappers below keep the historical signatures.)
// ============================================================================

// (Removal helper functions: header-inline on DelaunayView.)

// ============================================================================
// Vertex removal: main entry point
// ============================================================================

void DelaunayTriangulation::remove_flat_vertex(int v)
{
  // k_max = nh: the provable star-degree bound on a delta-complex (deg(v)
  // is bounded by the outgoing-arc count, not by nv -- multi-edges and
  // self-loops can push it past the vertex count).  O(nh) bytes.
  HostDelaunayWorkspace ws({.nv0 = nv, .k_max = nh, .nh_explicit = nh});
  if (tracker.active) {
    TrackerTransport tr{*this};
    DelaunayView::remove_flat_vertex(v, ws, BandedFloatMetric{}, tr);
  } else {
    DelaunayView::remove_flat_vertex(v, ws);
  }
  throw_on_status("remove_flat_vertex");
}

// (Self-loop cleanup + cocircular tie-break: header-inline on DelaunayView.)

namespace {

// The shared owner glue of the removal driver, for BOTH metric regimes: the
// work-list driver itself is header-inline on DelaunayView (see its doc for
// the algorithm and fixed-point structure).  This materializes the
// workspace, selects the transport policy by tracker.active, adapts the
// on_pop / verbose_removal diagnostics to the observer policy, and converts
// a non-Ok status to the documented throws.
template <class Metric>
void run_flat_removal(DelaunayTriangulation& D, Metric&& m,
                      const std::function<void(int)>& on_pop)
{
  auto count_live = [&]() {
    int n = 0;
    for (int v = 0; v < D.nv; v++) if (D.v_out[v] >= 0) n++;
    return n;
  };
  if (D.verbose_removal) {
    std::fprintf(stderr, "[remove_flat] start: %d live vertices\n", count_live());
    std::fflush(stderr);
  }

  struct OwnerObserver {
    DelaunayTriangulation& T;
    const std::function<void(int)>& on_pop_fn;
    long long removed = 0;
    void on_pop(int v) { if (on_pop_fn) on_pop_fn(v); }
    void on_removed(int) {
      if (T.verbose_removal && (++removed % T.verbose_removal == 0)) {
        int live = 0;
        for (int v = 0; v < T.nv; v++) if (T.v_out[v] >= 0) live++;
        std::fprintf(stderr, "[remove_flat] removed %lld, %d live remain\n",
                     removed, live);
        std::fflush(stderr);
      }
    }
  } obs{D, on_pop};

  HostDelaunayWorkspace ws({.nv0 = D.nv, .k_max = D.nh, .nh_explicit = D.nh});
  if (D.tracker.active) {
    TrackerTransport tr{D};
    D.DelaunayView::remove_flat_vertices(ws, m, tr, obs);
  } else {
    D.DelaunayView::remove_flat_vertices(ws, m, NoTransport{}, obs);
  }
  D.throw_on_status("remove_flat_vertices");

  if (D.verbose_removal) {
    std::fprintf(stderr, "[remove_flat] done: removed %lld, %d live remain\n",
                 obs.removed, count_live());
    std::fflush(stderr);
  }
}

}  // namespace

void DelaunayTriangulation::remove_flat_vertices(double flat_tol,
                                                 const std::function<void(int)>& on_pop)
{
  run_flat_removal(*this, BandedFloatMetric{flat_tol}, on_pop);
}

void DelaunayTriangulation::remove_flat_vertices_exact(const std::function<void(int)>& on_pop)
{
  // The exact-regime driver (paper sec:exactness).  The Lsq carry is
  // DERIVED from the current metric and VERIFIED -- entering the exact
  // regime on a metric it does not describe would be a silent wrong answer,
  // so both preconditions are loud:
  //   - every live he_length must square to a positive in-envelope integer
  //     within the integrality band (a from_intrinsic_metric complex, or
  //     the half-lengths bisect_multi_edges introduces, fail here);
  //   - the float curvature must agree with the integer cone excess at
  //     every live vertex (a complex with split_face / bisect-inserted
  //     vertices fails here: their v_orig_degree is not the equilateral
  //     labeling the exact flatness predicate reads).
  // The carry is transient (sized nh_cap: the machinery may allocate up to
  // capacity): the surviving he_length values are square roots of integers,
  // which the post-hoc utilities (canonical tesselation) recover exactly
  // through the same integrality boundary.
  using delaunay_detail::lsq_integrality_band;
  const double pi_3 = std::numbers::pi_v<double> / 3.0;
  std::vector<long long> Lsq((std::size_t)nh_cap, 0);
  for (int h = 0; h < nh; h++) {
    if (!alive(h)) continue;
    const double sq = he_length[h] * he_length[h];
    const long long n = std::llround(sq);
    if (n < 1 || n > exact_lsq_max ||
        std::abs(sq - (double)n) > lsq_integrality_band * std::max(1.0, sq))
      throw std::runtime_error(
          "remove_flat_vertices_exact: he_length[" + std::to_string(h) +
          "] does not square to a positive in-envelope integer (not an exact metric)");
    Lsq[h] = n;
  }
  for (int v = 0; v < nv; v++) {
    if (v_out[v] < 0) continue;
    // DelaunayView:: qualification: the owner's vector-returning curvature()
    // name-hides the view's per-vertex form.
    if (std::abs(DelaunayView::curvature(v) - cone_excess(v) * pi_3) > 1e-6)
      throw std::runtime_error(
          "remove_flat_vertices_exact: curvature at v=" + std::to_string(v) +
          " disagrees with its integer cone excess (not the equilateral labeling)");
  }
  run_flat_removal(*this, ExactIntegerMetric{std::span<long long>(Lsq)}, on_pop);
}

std::vector<int> DelaunayTriangulation::compact_vertices()
{
  // Live vertices (v_out >= 0) get fresh contiguous indices in their current
  // order; dead vertices are dropped. Half-edge origins and the per-vertex
  // arrays are rewritten; nv shrinks to the live count.
  std::vector<int> new_of_old(nv, -1), new_to_old;
  for (int v = 0; v < nv; v++)
    if (v_out[v] >= 0) { new_of_old[v] = (int)new_to_old.size(); new_to_old.push_back(v); }
  int nlive = (int)new_to_old.size();

  for (int h = 0; h < nh; h++)
    if (he_origin[h] >= 0) he_origin[h] = new_of_old[he_origin[h]];

  std::vector<int>    nout(nlive);
  std::vector<double> ncone(nlive);
  std::vector<int>    ndeg(nlive);
  for (int w = 0; w < nlive; w++) {
    int o = new_to_old[w];
    nout[w]  = v_out[o];
    ncone[w] = v_cone_angle[o];
    ndeg[w]  = v_orig_degree[o];
  }
  owned_v_out = std::move(nout);
  owned_v_cone_angle = std::move(ncone);
  owned_v_orig_degree = std::move(ndeg);
  nv = nlive;
  repoint();
  return new_to_old;
}

// --- Full algorithm ---

// (min_live_degree: header-inline on DelaunayView.)

bool DelaunayTriangulation::is_simplicial() const
{
  // Simplicial iff the arc map  h |-> (origin(h), dest(h))  is injective
  // on live half-edges.  Self-loops fail because both twins encode (v,v);
  // multi-edges fail because two non-twin half-edges encode the same
  // arc (u,v).  Both pathologies show up as a duplicate insertion.
  std::set<std::pair<int,int>> arcs;
  for (int h = 0; h < nh; h++) {
    if (!alive(h)) continue;
    if (!arcs.insert({he_origin[h], dest(h)}).second) return false;
  }
  return true;
}

bool DelaunayTriangulation::is_well_formed() const
{
  // Delegates to the view's walk over a temporary bit workspace (single
  // source of truth; the view form is the device-facing one).
  std::vector<std::uint64_t> words(Spanify::BitSpan::words_for(nh));
  Spanify::BitSpan visited{std::span<std::uint64_t>(words)};
  return DelaunayView::is_well_formed(visited);
}

// The shared prelude of both equilateral pipelines: flat vertices sorted
// last, DCEL built with the unit metric.
static DelaunayTriangulation equilateral_seed(const Triangulation& T)
{
  Triangulation sorted = T;
  sorted.apply_permutation(T.sort_flat_last());
  return DelaunayTriangulation::from_triangulation(sorted);
}

DelaunayTriangulation DelaunayTriangulation::compute(const Triangulation& T)
{
  // Sort flat vertices last, then build DCEL and run the algorithm.
  // Returns the unique intrinsic Delaunay triangulation of the input
  // surface (Bobenko-Springborn 2007).  The output is generally a
  // delta-complex: it may contain multi-edges, self-loops at cones,
  // and (rarely) bigons around cones (deg-2 cone vertices) -- all
  // legitimate features of the iDT object.  Callers that need a
  // strictly simplicial output should query min_live_degree() and
  // post-process (e.g. via bisect_multi_edges()) when needed.
  //
  // Equilateral input selects the EXACT integer metric: every decision
  // through the reduction -- flips, ears, flatness, the tie-break side
  // order included -- is a function of integer arithmetic alone, and the
  // output he_length values are square roots of integers.
  DelaunayTriangulation D = equilateral_seed(T);
  D.remove_flat_vertices_exact();
  return D;
}

DelaunayTriangulation DelaunayTriangulation::compute_banded(const Triangulation& T)
{
  // The general-metric (banded float) predicates applied to equilateral
  // input: the identical prelude, the float regime throughout.  The frozen
  // cross-reference for the parallel-primitives mirrors (which pin the
  // banded machinery byte-for-byte) and for the exact-vs-banded corpus
  // verdict; production callers use compute(T).
  DelaunayTriangulation D = equilateral_seed(T);
  D.remove_flat_vertices();
  return D;
}

DelaunayTriangulation DelaunayTriangulation::compute(const Triangulation& T,
                                                     const EdgeLengthFn& length,
                                                     double flat_tol,
                                                     std::vector<int>* new_to_old,
                                                     bool track_removed)
{
  // Prescribed-metric iDT. Unlike the equilateral compute(T), we do NOT
  // sort_flat_last() (that classifies flatness by degree, which is wrong for
  // a general metric, and would also permute the vertices the caller's
  // `length` oracle is keyed on). remove_flat_vertices() identifies flat
  // vertices by cone angle (is_flat) and its full-scan loop is order-robust,
  // so no pre-sort is needed.
  DelaunayTriangulation D = from_intrinsic_metric(T, length);
  // With track_removed, every removed flat vertex is seeded as a tracked
  // point (label = its T-label: removals never renumber) and transported
  // through all removal/flip surgery -- see the tracker banner (delaunay.hh).
  if (track_removed) D.enable_point_tracking();
  D.remove_flat_vertices(flat_tol);
  // Cones are scattered (no flat-last sort); compact to nv = cones.  The
  // surviving-cone labels in T's numbering are reported through new_to_old
  // when requested (callers that annotate cones, e.g. AlexandrovIDTCubic).
  // compact_vertices renumbers only vertices: tracked points (face +
  // barycentric slots) are untouched.
  std::vector<int> n2o = D.compact_vertices();
  if (new_to_old) *new_to_old = std::move(n2o);
  return D;
}

// ============================================================================
// Canonical Delaunay tesselation
// See CANONICAL-TESSELATION.md for the algorithm derivation and the
// Bobenko-Springborn uniqueness rationale.
// ============================================================================

bool DelaunayTriangulation::is_cocircular_edge(int h) const
{
  if (!alive(h)) return false;
  return diamond(h).is_cocircular_exact();
}

// Per-half-edge mask from a per-edge predicate, symmetric on twins.
template <class Pred>
static vector<bool> edge_mask(const DelaunayTriangulation& D, Pred pred)
{
  vector<bool> mask(D.nh, false);
  for (int h : D.edges())
    mask[h] = mask[h ^ 1] = pred(h);
  return mask;
}

vector<bool> DelaunayTriangulation::cocircular_edges() const
{
  return edge_mask(*this, [&](int h) { return diamond(h).is_cocircular_exact(); });
}

vector<bool> DelaunayTriangulation::cocircular_edges(double tol) const
{
  return edge_mask(*this, [&](int h) { return diamond(h).is_cocircular(tol); });
}

// Lex-min cyclic rotation of a polygon boundary (oriented surface, no reverse).
static CanonicalTesselation::Polygon
min_rotation(const CanonicalTesselation::Polygon& p)
{
  int d = (int)p.size();
  if (d <= 1) return p;
  CanonicalTesselation::Polygon best = p;
  CanonicalTesselation::Polygon rot(d);
  for (int r = 1; r < d; r++) {
    for (int i = 0; i < d; i++) rot[i] = p[(r + i) % d];
    if (rot < best) best = rot;
  }
  return best;
}

CanonicalTesselation
DelaunayTriangulation::canonical_tesselation(const vector<int>& vertex_labels) const
{
  return canonical_tesselation(vertex_labels, cocircular_edges());
}

CanonicalTesselation
DelaunayTriangulation::canonical_tesselation(const vector<int>& vertex_labels,
                                             const vector<bool>& tight) const
{
  // Walk cell boundaries.  Each non-tight half-edge h sits on exactly one
  // cell (the one to its left in the DCEL CCW orientation).  Within a cell,
  // tight edges are interior; we step across them with `next(twin(.))`.
  vector<bool> visited(nh, false);
  CanonicalTesselation T;
  for (int h_start = 0; h_start < nh; h_start++) {
    if (!alive(h_start) || visited[h_start] || tight[h_start]) continue;
    CanonicalTesselation::Polygon poly;
    int h = h_start;
    do {
      visited[h] = true;
      int u = he_origin[h];
      long long L = (long long)std::llround(he_length[h] * he_length[h]);
      poly.push_back({(u >= 0 && u < (int)vertex_labels.size()) ? vertex_labels[u] : u, L});
      // Advance: walk the cell boundary CCW.  Within the cell, tight edges
      // are interior, so step into the adjacent triangle until the next
      // boundary edge.
      int h_next = he_next[h];
      int safety = 0;
      while (tight[h_next]) {
        h_next = he_next[h_next ^ 1];
        if (++safety > nh)
          // Deep invariant failure: two silently-empty results would
          // compare equal, so fail loud instead of returning a sentinel.
          throw std::runtime_error(
              "canonical_tesselation: cell-boundary walk from half-edge " +
              std::to_string(h_start) + " failed to close after " +
              std::to_string(nh) + " interior-edge (tight) steps, stuck at "
              "half-edge " + std::to_string(h_next) + "; the tight mask "
              "encloses a cell -- a well-formed iDT tesselation always closes. "
              "(Thrown rather than returning an empty tesselation, which a "
              "legitimately empty iDT also yields and so could not signal this.)");
      }
      h = h_next;
    } while (h != h_start);
    T.cells.push_back(min_rotation(poly));
  }
  std::sort(T.cells.begin(), T.cells.end());
  return T;
}

size_t CanonicalTesselation::fingerprint() const
{
  size_t acc = 0xcbf29ce484222325ULL;          // FNV-1a basis, used as seed
  std::hash<long long> hh;
  auto mix = [&](long long x){
    acc ^= hh(x) + 0x9e3779b97f4a7c15ULL + (acc << 6) + (acc >> 2);
  };
  for (auto& cell : cells) {
    mix((long long)cell.size());
    for (auto& [v, L] : cell) { mix((long long)v); mix(L); }
    mix(0x5bd1e9955LL);                        // cell separator
  }
  return acc;
}

string CanonicalTesselation::to_string() const
{
  ostringstream os;
  os << "CanonicalTesselation{n_cells=" << cells.size() << ", cells=[\n";
  for (size_t i = 0; i < cells.size(); i++) {
    os << "  [n=" << cells[i].size() << "]";
    for (auto& [v, L] : cells[i]) os << " (" << v << "," << L << ")";
    os << "\n";
  }
  os << "]}";
  return os.str();
}

// Bisect a multi-edge by inserting a midpoint vertex.
// h: one half-edge of the multi-edge to bisect (u -> v with length L).
// The midpoint w is at distance L/2 from both u and v.
// Each of the two adjacent faces (u,v,B) and (v,u,D) is split into two.
// Returns the index of the new vertex w.
//
// Before:
//      B                      B
//     / \                    /|\
//    a   b       After:     a | b'
//   /     \                /  |  \
//  u---e---v     =>    u--e/2-w-e/2-v
//   \     /                \  |  /
//    c   d                  c | d'
//     \ /                    \|/
//      D                      D
//
// New edges: u-w (L/2), w-v (L/2), w-B (computed), w-D (computed).
// New faces: (u,w,B), (w,v,B), (v,w,D), (w,u,D).
static int bisect_edge(DelaunayTriangulation& D, int h) {
  int t = D.twin(h);
  int u = D.he_origin[h], v = D.dest(h);
  double L = D.he_length[h];

  // Face of h: u -> v -> B -> u
  int h_vB = D.he_next[h];        // v -> B
  int h_Bu = D.he_next[h_vB];     // B -> u
  int B = D.dest(h_vB);
  double a = D.he_length[h_Bu];   // |B-u|
  double b = D.he_length[h_vB];   // |v-B|

  // Face of t: v -> u -> D -> v
  int h_uD = D.he_next[t];        // u -> D
  int h_Dv = D.he_next[h_uD];     // D -> v
  int Dv = D.dest(h_uD);
  double c = D.he_length[h_uD];   // |u-D|
  double d = D.he_length[h_Dv];   // |D-v|

  double half = L / 2.0;

  // Stewart's theorem for a cevian to the midpoint:  m² = (a² + b²)/2 - L²/4.
  double wB2 = (a*a + b*b) / 2.0 - L*L / 4.0;
  double wB = wB2 > 0 ? sqrt(wB2) : 0;

  // Similarly for |w-D|.
  double wD2 = (c*c + d*d) / 2.0 - L*L / 4.0;
  double wD = wD2 > 0 ? sqrt(wD2) : 0;

  // Allocate new vertex w (flat, original degree 4).
  int w = D.alloc_vertex(2.0 * M_PI, 4);

  // Delete original edge and its two faces.
  D.dealloc_face(D.he_face[h]);
  D.dealloc_face(D.he_face[t]);
  D.dealloc_edge(h);

  // Allocate four new directed edges incident to w.
  int uw    = D.alloc_directed_edge(u, w,  half);
  int wv    = D.alloc_directed_edge(w, v,  half);
  int wB_he = D.alloc_directed_edge(w, B,  wB);
  int wD_he = D.alloc_directed_edge(w, Dv, wD);

  D.v_out[w] = D.twin(uw);  // w -> u

  // Wire four new CCW faces around w.
  D.wire_triangle(uw,          wB_he,   h_Bu);          // (u, w, B)
  D.wire_triangle(wv,          h_vB,    D.twin(wB_he)); // (w, v, B)
  D.wire_triangle(D.twin(wv),  wD_he,   h_Dv);          // (v, w, D)
  D.wire_triangle(D.twin(uw),  h_uD,    D.twin(wD_he)); // (w, u, D)

  // Fix v_out for u, v if they pointed at the deleted edge.
  if (D.v_out[u] == h) D.v_out[u] = uw;
  if (D.v_out[v] == t) D.v_out[v] = D.twin(wv);

  return w;
}

int DelaunayTriangulation::bisect_multi_edges() {
  // Not transport-hooked: bisect_edge deallocs/rewires faces without moving
  // tracked points, and a recycled face slot would silently inherit stale
  // bucket entries.  Fail loud instead of corrupting (hook it if a tracking
  // caller ever needs it).
  if (tracker.active)
    throw std::runtime_error("bisect_multi_edges: point tracking is active; "
                             "this operation is not transport-hooked");
  // Find multi-edges: vertex pairs with >1 edge.
  map<pair<int,int>, vector<int>> pair_to_hes;
  for (int h : edges()) {
    int u = he_origin[h], v = dest(h);
    pair_to_hes[{min(u,v), max(u,v)}].push_back(h);
  }

  int n_inserted = 0;
  for (auto& [vp, hes] : pair_to_hes) {
    if (hes.size() <= 1) continue;
    // Bisect each edge in this multi-edge group.
    for (int h : hes) {
      if (!alive(h)) continue;
      bisect_edge(*this, h);
      n_inserted++;
    }
  }
  throw_on_status("bisect_multi_edges");
  return n_inserted;
}

// --- Vertex / face surgery ---

int DelaunayTriangulation::alloc_vertex(double cone_angle, int orig_degree) {
  int v = nv++;
  ensure_vertices(nv);
  v_out[v] = -1;
  v_cone_angle[v] = cone_angle;
  v_orig_degree[v] = orig_degree;
  return v;
}

int DelaunayTriangulation::split_face(int h0, std::array<double,3> spokes) {
  // Not transport-hooked (same rationale as bisect_multi_edges).
  if (tracker.active)
    throw std::runtime_error("split_face: point tracking is active; "
                             "this operation is not transport-hooked");
  int h1 = he_next[h0], h2 = he_next[h1];
  int a = he_origin[h0], b = he_origin[h1], c = he_origin[h2];
  int f = he_face[h0];

  int P = alloc_vertex(2.0 * M_PI, 3);
  dealloc_face(f);                                  // only the face; the 3 edges stay
  int sa = alloc_directed_edge(P, a, spokes[0]);
  int sb = alloc_directed_edge(P, b, spokes[1]);
  int sc = alloc_directed_edge(P, c, spokes[2]);
  wire_triangle(sa, h0, twin(sb));                  // (P, a, b)
  wire_triangle(sb, h1, twin(sc));                  // (P, b, c)
  wire_triangle(sc, h2, twin(sa));                  // (P, c, a)

  v_out[P] = sa;
  v_out[a] = twin(sa); v_out[b] = twin(sb); v_out[c] = twin(sc);
  // wire_triangle set the new faces' angles; refresh the cached cone angle at the four
  // vertices whose incident faces changed.
  for (int v : {P, a, b, c}) recompute_cone_angle(v);
  throw_on_status("split_face");
  return P;
}

// --- Curvature ---

std::vector<double> DelaunayTriangulation::curvature() const {
  std::vector<double> K(nv);
  for (int v = 0; v < nv; v++) K[v] = DelaunayView::curvature(v);
  return K;
}

// --- Geodesic disks (bounded multi-source Voronoi) ---

// Fast-marching update: the wavefront reaches the two ends of an edge (length `edge`)
// at times d0, d1; return its arrival time at the triangle's opposite vertex, whose
// edges to those two ends are e0, e1, travelling straight across the unfolded
// interior. +inf when the straight crossing is invalid (obtuse triangle, or d0,d1
// inconsistent) so the caller's edge relaxation is the fallback.
static double fast_march(double d0, double d1, double edge, double e0, double e1) {
  const double kInf = std::numeric_limits<double>::infinity();
  if (edge <= 0.0) return kInf;
  double xt  = (edge*edge + e0*e0 - e1*e1) / (2.0*edge);   // target, unfolded above
  double yt2 = e0*e0 - xt*xt;
  if (yt2 <= 0.0) return kInf;                             // degenerate face
  double yt  = sqrt(yt2);
  double xs  = (d0*d0 - d1*d1 + edge*edge) / (2.0*edge);   // virtual source, below
  double ys2 = d0*d0 - xs*xs;
  if (ys2 < 0.0) return kInf;                              // |d0-d1| > edge
  double ys  = -sqrt(ys2);
  double cross = xs + (xt - xs) * (-ys) / (yt - ys);       // line source->target meets edge
  if (cross < 0.0 || cross > edge) return kInf;            // wave rounds a vertex
  double d = hypot(xt - xs, yt - ys);
  return d >= std::max(d0, d1) ? d : kInf;                 // causality (monotone front)
}

std::vector<GeodesicDisk>
DelaunayTriangulation::geodesic_disks(const std::vector<int>& sources, double R,
                                      DiskMetric metric) const {
  const double kInf = std::numeric_limits<double>::infinity();
  std::vector<GeodesicDisk> disks;
  for (int s : sources) disks.push_back({s, {}});   // ORDER CONTRACT: disks[i].source == sources[i]
  std::vector<double> dist(nv, kInf);
  std::vector<int>    owner(nv, -1);    // index into sources/disks of the claiming source
  std::vector<char>   frozen(nv, 0);
  std::priority_queue<std::pair<double,int>, std::vector<std::pair<double,int>>,
                      std::greater<>> frontier;
  for (int i = 0; i < (int)sources.size(); i++) {
    dist[sources[i]] = 0.0; owner[sources[i]] = i; frontier.push({0.0, sources[i]});
  }

  auto relax = [&](int v, double cand, int own) {
    if (!frozen[v] && cand <= R && cand < dist[v]) { dist[v] = cand; owner[v] = own; frontier.push({cand, v}); }
  };
  while (!frontier.empty()) {
    auto [d, u] = frontier.top(); frontier.pop();
    if (frozen[u]) continue;
    frozen[u] = 1;
    int o = owner[u];
    disks[o].members.push_back({u, d});
    for (int h : incident(u)) {                     // each incident face (u, x, y)
      int    h1 = he_next[h], h2 = he_next[h1];
      int    x = dest(h), y = he_origin[h2];
      double L_ux = he_length[h], L_xy = he_length[h1], L_uy = he_length[h2];
      relax(x, d + L_ux, o);                         // edge u -> x
      if (metric == DiskMetric::Unfold) {            // wavefront across the face, within the cell
        if (frozen[x] && owner[x] == o) relax(y, fast_march(d, dist[x], L_ux, L_uy, L_xy), o);
        if (frozen[y] && owner[y] == o) relax(x, fast_march(d, dist[y], L_uy, L_ux, L_xy), o);
      }
    }
  }
  // disks[i].source == sources[i] holds throughout: disks were built in `sources` order and
  // never reordered; owner[] indexes into them. Callers rely on this 1:1 ordering.
  return disks;
}

// --- Validation ---
// (check_consistency is header-inline on DelaunayView: the nine numbered
// class invariants, allocation-free.)

// ============================================================================
// Serialization -- the ".idt" text format
// ============================================================================
//
// Format (text, versioned):
//   iDT-DCEL 1
//   <nv> <nf> <ne>
//   <v_orig_degree[v]>              (nv lines)
//   <o0> <o1> <n0> <n1> <length>    (ne lines)
//
// ALIVE elements only, reindexed dense (nh/nf keep dead slots after remove_flat_vertices +
// compact_vertices, which shrinks only nv). Twin convention: half-edges 2k and 2k+1 are twins.
// An edge line gives the two half-edges' origin vertices (o0 = origin(2k), o1 = origin(2k+1) =
// dest(2k)), their he_next successors (n0, n1 as reindexed half-edge ids), and the shared length.
// Faces (he_face / f_he), v_out and angles are derived on read, so only the topology + metric is
// stored, at full double precision (%.17g = DBL_DECIMAL_DIG). The intrinsic geometry round-trips
// exactly. The format is faithful to a non-simplicial delta-complex (multi-edges / self-loops).

bool DelaunayTriangulation::to_ascii(const DelaunayTriangulation& D, FILE* file)
{
  for (int v = 0; v < D.nv; v++)
    if (D.v_out[v] < 0)
      throw std::runtime_error("DelaunayTriangulation::to_ascii: vertex " + std::to_string(v)
                               + " is dead (non-compacted DCEL; compact_vertices first)");

  // old edge id -> dense new edge index, assigned in ascending old order; -1 for dead.
  std::vector<int> new_edge(D.nh / 2, -1);
  int ne = 0;
  for (int h : D.edges()) new_edge[D.edge(h)] = ne++;
  auto new_he = [&](int h) { return 2 * new_edge[D.edge(h)] + (h & 1); };

  int nf_alive = 0;
  for (int f = 0; f < D.nf; f++)
    if (D.f_he[f] >= 0) nf_alive++;

  std::fprintf(file, "iDT-DCEL 1\n%d %d %d\n", D.nv, nf_alive, ne);
  for (int v = 0; v < D.nv; v++)
    std::fprintf(file, "%d\n", D.v_orig_degree[v]);
  for (int h : D.edges())
    std::fprintf(file, "%d %d %d %d %.17g\n", D.he_origin[h], D.he_origin[h + 1],
                 new_he(D.he_next[h]), new_he(D.he_next[h + 1]), D.he_length[h]);
  return std::ferror(file) == 0;
}

DelaunayTriangulation DelaunayTriangulation::from_ascii(FILE* file)
{
  char magic[16];
  int version = 0;
  if (std::fscanf(file, "%15s %d", magic, &version) != 2 || std::string(magic) != "iDT-DCEL")
    throw std::runtime_error("DelaunayTriangulation::from_ascii: bad header (expected 'iDT-DCEL <version>')");
  if (version != 1)
    throw std::runtime_error("DelaunayTriangulation::from_ascii: unsupported version " + std::to_string(version));

  int nv = 0, nf = 0, ne = 0;
  // ne is bounded so nh = 2*ne cannot overflow int (it is then used as a vector size).
  if (std::fscanf(file, "%d %d %d", &nv, &nf, &ne) != 3 || nv < 0 || nf < 0 || ne < 0 || ne > INT_MAX / 2)
    throw std::runtime_error("DelaunayTriangulation::from_ascii: bad counts line");

  DelaunayTriangulation D;
  D.nv = nv;
  D.ensure_vertices(nv);      // v_out -1, cone 0, degree 0 (read below)
  for (int v = 0; v < nv; v++)
    if (std::fscanf(file, "%d", &D.v_orig_degree[v]) != 1 || D.v_orig_degree[v] < 0)
      throw std::runtime_error("DelaunayTriangulation::from_ascii: bad vertex table (truncated or "
                               "negative degree) at vertex " + std::to_string(v));

  const int nh = 2 * ne;
  D.nh = nh;
  D.ensure_halfedges(nh);     // dead-slot fill == the historical assigns
  auto in_verts = [&](int x) { return x >= 0 && x < nv; };
  auto in_hes   = [&](int x) { return x >= 0 && x < nh; };
  for (int k = 0; k < ne; k++) {
    int o0, o1, n0, n1;
    double L;
    if (std::fscanf(file, "%d %d %d %d %lf", &o0, &o1, &n0, &n1, &L) != 5)
      throw std::runtime_error("DelaunayTriangulation::from_ascii: truncated edge table at edge " + std::to_string(k));
    if (!in_verts(o0) || !in_verts(o1) || !in_hes(n0) || !in_hes(n1) || !(L > 0))
      throw std::runtime_error("DelaunayTriangulation::from_ascii: out-of-range field at edge " + std::to_string(k));
    D.he_origin[2 * k] = o0;  D.he_origin[2 * k + 1] = o1;
    D.he_next[2 * k]   = n0;  D.he_next[2 * k + 1]   = n1;
    D.he_length[2 * k] = D.he_length[2 * k + 1] = L;
  }

  // Faces: each he_next cycle is one face. Assign he_face + a representative f_he per cycle.
  D.owned_f_he.clear();
  for (int h = 0; h < nh; h++) {
    if (D.he_face[h] >= 0) continue;
    const int f = (int)D.owned_f_he.size();
    D.owned_f_he.push_back(h);
    int g = h, steps = 0;
    do {
      D.he_face[g] = f;
      g = D.he_next[g];
      if (++steps > nh)
        throw std::runtime_error("DelaunayTriangulation::from_ascii: he_next cycle from half-edge "
                                 + std::to_string(h) + " did not close");
    } while (g != h);
  }
  D.nf = (int)D.owned_f_he.size();
  D.owned_free_faces.resize(D.owned_f_he.size());
  D.repoint();
  if (D.nf != nf)
    throw std::runtime_error("DelaunayTriangulation::from_ascii: derived " + std::to_string(D.nf)
                             + " faces but header declared " + std::to_string(nf));

  // v_out: any outgoing half-edge per vertex (the cw ring from it covers the fan).
  for (int h = 0; h < nh; h++)
    if (D.v_out[D.he_origin[h]] < 0) D.v_out[D.he_origin[h]] = h;

  // An isolated declared vertex (no incident edge) would keep v_out == -1, so recompute_cone_angles
  // leaves its v_cone_angle at 0 and curvature()/is_flat() silently misreport it. Fail loud.
  for (int v = 0; v < nv; v++)
    if (D.v_out[v] < 0)
      throw std::runtime_error("DelaunayTriangulation::from_ascii: vertex " + std::to_string(v)
                               + " has no incident edge (isolated / malformed input)");

  D.recompute_all_angles();   // he_angle from the law of cosines, THEN refresh v_cone_angle

  // Full structural + metric validation: triangular faces (he_next^3 == h), origin chaining
  // (dest(h) == origin(next(h))), positive twin-consistent lengths, and the triangle inequality per
  // face (which also rejects a +inf length that slips past the per-field L>0 guard). Without this a
  // malformed stream -- a non-triangular he_next permutation whose cycle count happens to match nf,
  // a +inf length, a broken chain -- would load as a silently-wrong triangulation.
  if (!D.check_consistency())
    throw std::runtime_error("DelaunayTriangulation::from_ascii: reconstructed DCEL failed "
                             "check_consistency (non-triangular face, broken chaining, or a length / "
                             "triangle-inequality violation)");
  return D;
}

