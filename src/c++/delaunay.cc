#include "fullerenes/delaunay.hh"

#include <stack>
#include <cmath>
#include <cassert>
#include <algorithm>
#include <map>
#include <climits>
#include <limits>
#include <queue>
#include <stdexcept>
#include <string>
#include <sstream>
#include <cstdio>
#include <limits>
#include <queue>

// ============================================================================
// Intrinsic geometry primitives
// ============================================================================

// Heron product: H(a,b,c) = (a+b+c)(-a+b+c)(a-b+c)(a+b-c) = 16*Area^2.
// Returns 0 if triangle inequality is violated.
static double heron(double a, double b, double c)
{
  double s1 = -a + b + c;
  double s2 =  a - b + c;
  double s3 =  a + b - c;
  if (s1 < 0 || s2 < 0 || s3 < 0) return 0;
  return (a + b + c) * s1 * s2 * s3;
}

// Cotangent of angle opposite side `opp` in triangle with sides (opp, b, c).
// cot(alpha) = (b^2 + c^2 - opp^2) / sqrt(H).
static double cot_opposite(double opp, double b, double c)
{
  double H = heron(opp, b, c);
  double num = b*b + c*c - opp*opp;
  if (H <= 0) return (num >= 0) ? 1e15 : -1e15;
  return num / sqrt(H);
}

// Angle of the corner adjacent to sides `adj1`, `adj2` in a triangle
// whose opposite side is `opp`.  Law of cosines, clamped for floating-point
// safety at triangle-inequality boundaries.
static double triangle_angle(double adj1, double adj2, double opp)
{
  double c = (adj1*adj1 + adj2*adj2 - opp*opp) / (2 * adj1 * adj2);
  return acos(std::clamp(c, -1.0, 1.0));
}

// ============================================================================
// Diamond geometry
// ============================================================================

bool Diamond::is_delaunay() const
{
  return cot_opposite(e, a, b) + cot_opposite(e, c, d) >= -1e-10;
}

bool Diamond::is_convex() const
{
  // sin(angle_at_u) proportional to sqrt(Ha)*Q + P*sqrt(Hd), must be > 0.
  // sin(angle_at_v) proportional to sqrt(Ha)*Qv + Pv*sqrt(Hd), must be > 0.
  double e2 = e*e;
  double Ha = heron(e, a, b), Hd = heron(e, c, d);
  double sHa = (Ha > 0) ? sqrt(Ha) : 0;
  double sHd = (Hd > 0) ? sqrt(Hd) : 0;

  double Pu = e2 + a*a - b*b, Qu = e2 + c*c - d*d;
  if (sHa * Qu + Pu * sHd <= 1e-12) return false;

  double Pv = e2 + b*b - a*a, Qv = e2 + d*d - c*c;
  return sHa * Qv + Pv * sHd > 1e-12;
}

double Diamond::flipped_length() const
{
  // f^2 = a^2 + c^2 - (PQ - sqrt(Ha*Hd)) / (2e^2)
  double e2 = e*e, a2 = a*a, b2 = b*b, c2 = c*c, d2 = d*d;
  double P = e2 + a2 - b2;
  double Q = e2 + c2 - d2;
  double Ha = heron(e, a, b), Hd = heron(e, c, d);
  double sqrtHH = (Ha > 0 && Hd > 0) ? sqrt(Ha * Hd) : 0;
  double f2 = a2 + c2 - (P * Q - sqrtHH) / (2 * e2);
  return (f2 > 0) ? sqrt(f2) : 0;
}

bool Diamond::is_cocircular() const
{
  // Tight Delaunay: cot(angle_B) + cot(angle_D) == 0 exactly.  Equivalent
  // to s1 * area_2 + s2 * area_1 = 0 where s1 = a^2+b^2-e^2, s2 = c^2+d^2-e^2.
  // Squaring after sign-check: tight iff sign(s1) != sign(s2) and
  // s1^2 * H2 == s2^2 * H1 (with H = 16*area^2).  Or: both s1, s2 == 0.
  // Done in integer length-squared arithmetic so the predicate is exact for
  // equilateral triangulations and any sequence of flips.
  long long Le = (long long)std::llround(e * e);
  long long La = (long long)std::llround(a * a);
  long long Lb = (long long)std::llround(b * b);
  long long Lc = (long long)std::llround(c * c);
  long long Ld = (long long)std::llround(d * d);
  long long s1 = La + Lb - Le;
  long long s2 = Lc + Ld - Le;
  if (s1 == 0 && s2 == 0) return true;
  if (s1 == 0 || s2 == 0) return false;
  if ((s1 > 0) == (s2 > 0)) return false;          // same sign: not tight
  auto H = [](long long x, long long y, long long z) {
    return 2*(x*y + y*z + x*z) - (x*x + y*y + z*z); // 16 * area^2
  };
  return s1 * s1 * H(Le, Lc, Ld) == s2 * s2 * H(Le, La, Lb);
}

bool Diamond::is_cocircular(double tol) const
{
  // Float tight-Delaunay test for non-equilateral metrics: cot(angle_B) +
  // cot(angle_D) == 0, where angle_B is opposite e in triangle (e,a,b) and
  // angle_D opposite e in (e,c,d). cot_opposite returns +/-1e15 on a
  // degenerate triangle, which never lands within a sane tol, so degenerate
  // diamonds are correctly reported non-tight.
  double cotB = cot_opposite(e, a, b);
  double cotD = cot_opposite(e, c, d);
  return std::abs(cotB + cotD) < tol;
}

// Old FulleroidDelaunay + IDTAudit implementation moved to delaunay_old.cc.

// ============================================================================
// DelaunayTriangulation — DCEL-based iDT (delta-complex)
// ============================================================================

// --- Allocation ---

int DelaunayTriangulation::alloc_edge()
{
  int eid;
  if (!free_edges.empty()) {
    eid = free_edges.back();
    free_edges.pop_back();
  } else {
    eid = nh / 2;
    nh += 2;
    he_next.resize(nh, -1);
    he_origin.resize(nh, -1);
    he_face.resize(nh, -1);
    he_length.resize(nh, 0);
    he_angle.resize(nh, 0);
  }
  he_origin[2*eid] = -1;
  he_origin[2*eid+1] = -1;
  return 2 * eid;
}

int DelaunayTriangulation::alloc_face()
{
  int fid;
  if (!free_faces.empty()) {
    fid = free_faces.back();
    free_faces.pop_back();
  } else {
    fid = nf++;
    f_he.push_back(-1);
  }
  f_he[fid] = -1;
  return fid;
}

void DelaunayTriangulation::dealloc_edge(int h)
{
  int eid = h >> 1;
  he_origin[2*eid] = -1;
  he_origin[2*eid+1] = -1;
  he_next[2*eid] = -1;
  he_next[2*eid+1] = -1;
  he_face[2*eid] = -1;
  he_face[2*eid+1] = -1;
  he_length[2*eid] = 0;
  he_length[2*eid+1] = 0;
  free_edges.push_back(eid);
}

void DelaunayTriangulation::dealloc_face(int f)
{
  f_he[f] = -1;
  free_faces.push_back(f);
}

int DelaunayTriangulation::alloc_directed_edge(int u, int v, double L)
{
  int h = alloc_edge();
  he_origin[h]     = u;
  he_origin[h ^ 1] = v;
  he_length[h] = he_length[h ^ 1] = L;
  return h;
}

int DelaunayTriangulation::wire_triangle(int h0, int h1, int h2)
{
  int fid = alloc_face();
  he_next[h0] = h1; he_next[h1] = h2; he_next[h2] = h0;
  he_face[h0] = he_face[h1] = he_face[h2] = fid;
  f_he[fid] = h0;
  recompute_face_angles(fid);
  return fid;
}

void DelaunayTriangulation::ensure_v_out(int v)
{
  int h = v_out[v];
  if (h >= 0 && alive(h) && he_origin[h] == v) return;
  // Scan all half-edges for a live outgoing one from v.  O(nh) fallback;
  // only used to recover after structural edits that may have killed
  // the previously-recorded outgoing pointer.
  for (int hh = 0; hh < nh; hh++)
    if (alive(hh) && he_origin[hh] == v) { v_out[v] = hh; return; }
  v_out[v] = -1;
}

int DelaunayTriangulation::vertex_degree(int v) const
{
  if (v_out[v] < 0) return 0;
  int h0 = v_out[v], h = h0, deg = 0;
  do { deg++; h = cw(h); } while (h != h0);
  return deg;
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
    const int h0 = v_out[u];
    if (h0 < 0) continue;
    int h = h0;
    do {
      if (alive(h)) {
        int v = he_origin[h ^ 1];
        if (v >= 0 && v < nv) {
          double nd = d + he_length[h];
          if (nd < dist[v]) { dist[v] = nd; pq.push({nd, v}); }
        }
      }
      h = cw(h);
    } while (h != h0);
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

// Build the metric-independent half-edge topology (connectivity + faces +
// v_orig_degree) from a combinatorial triangulation. The metric arrays
// (he_length, he_angle, v_cone_angle) are sized but left for the caller to
// fill -- this is the shared core of from_triangulation (equilateral) and
// from_intrinsic_metric (prescribed lengths).
static DelaunayTriangulation build_topology(const Triangulation& T)
{
  DelaunayTriangulation D;
  D.nv = T.N;

  // Phase 1: Assign half-edge IDs to directed arcs.
  // For each undirected edge {u,v} with u<v: half-edges 2k (u->v) and 2k+1 (v->u).
  map<pair<int,int>, int> arc_to_he;
  int eid = 0;
  for (node_t u = 0; u < T.N; u++)
    for (node_t v : T[u])
      if (u < v) {
        arc_to_he[{u,v}] = 2 * eid;
        arc_to_he[{v,u}] = 2 * eid + 1;
        eid++;
      }

  D.nh = 2 * eid;
  D.he_next.resize(D.nh);
  D.he_origin.resize(D.nh);
  D.he_face.resize(D.nh, -1);
  D.he_length.resize(D.nh);   // metric: caller fills
  D.he_angle.resize(D.nh);    // metric: caller fills (recompute_all_angles)

  // Phase 2: Set origins.
  for (auto& [arc, hid] : arc_to_he)
    D.he_origin[hid] = arc.first;

  // Phase 3: Set next pointers and face assignments.
  // For vertex u with CCW neighbors [..., v, w, ...]:
  //   arc u->v is in face (u, v, w) where w = next neighbor after v.
  //   next(u->v) = v->w.
  D.nf = 0;
  D.v_out.assign(D.nv, -1);

  for (node_t u = 0; u < T.N; u++) {
    auto row = T[u];
    int deg = row.size();
    for (int j = 0; j < deg; j++) {
      node_t v = row[j], w = row[(j+1) % deg];
      int h_uv = arc_to_he.at({u, v});
      int h_vw = arc_to_he.at({v, w});
      D.he_next[h_uv] = h_vw;

      // Assign face when u is the smallest vertex (canonical representative).
      if (u < v && u < w) {
        int fid = D.nf++;
        int h_wu = arc_to_he.at({w, u});
        D.he_face[h_uv] = fid;
        D.he_face[h_vw] = fid;
        D.he_face[h_wu] = fid;
        D.f_he.push_back(h_uv);
      }
    }
    if (deg > 0) D.v_out[u] = arc_to_he.at({u, row[0]});
  }

  // Per-vertex original degree (topological; metric-independent).
  D.v_orig_degree.resize(D.nv);
  D.v_cone_angle.resize(D.nv);  // caller fills
  for (node_t v = 0; v < T.N; v++) D.v_orig_degree[v] = T.degree(v);

  return D;
}

DelaunayTriangulation DelaunayTriangulation::from_triangulation(const Triangulation& T)
{
  // Equilateral metric: every edge length 1, every angle pi/3, cone angle
  // deg*pi/3. (The special case length == 1 of from_intrinsic_metric, kept
  // explicit because the constants are exact and this is the hot path.)
  DelaunayTriangulation D = build_topology(T);
  D.he_length.assign(D.nh, 1.0);
  D.he_angle.assign(D.nh, M_PI / 3.0);
  for (node_t v = 0; v < D.nv; v++)
    D.v_cone_angle[v] = D.v_orig_degree[v] * M_PI / 3.0;
  return D;
}

DelaunayTriangulation
DelaunayTriangulation::from_intrinsic_metric(const Triangulation& T,
                                             const EdgeLengthFn& length)
{
  // Prescribed intrinsic metric: edge lengths from `length`, angles and cone
  // angles derived. The topology is identical to the equilateral case; only
  // the metric differs.
  DelaunayTriangulation D = build_topology(T);
  for (int h = 0; h < D.nh; h += 2) {
    double L = length(D.he_origin[h], D.he_origin[h ^ 1]);
    D.he_length[h] = D.he_length[h ^ 1] = L;
  }
  D.recompute_all_angles();   // fills he_angle AND refreshes the v_cone_angle cache
  return D;
}

// --- Geometry ---

Diamond DelaunayTriangulation::diamond(int h) const
{
  // h: u->v.  Face left of h has third vertex B = dest(next(h)).
  // Twin face has third vertex D = dest(next(twin(h))).
  int t = h ^ 1;
  int u = he_origin[h], v = he_origin[t];
  int B = dest(he_next[h]);
  int D = dest(he_next[t]);

  double e_len = he_length[h];
  // a = edge from u to B, b = edge from v to B (in face of h)
  int h_vB = he_next[h];          // v->B
  int h_Bu = he_next[h_vB];       // B->u
  double a = he_length[h_Bu];     // side adjacent to u (B-u)
  double b = he_length[h_vB];     // side adjacent to v (v-B)

  // c = edge from u to D, d = edge from v to D (in face of twin)
  int h_uD = he_next[t];          // u->D
  int h_Dv = he_next[h_uD];       // D->v
  double c = he_length[h_uD];     // side adjacent to u (u-D)
  double d = he_length[h_Dv];     // side adjacent to v (D-v)

  return {e_len, a, b, c, d};
}

void DelaunayTriangulation::recompute_face_angles(int f)
{
  if (f_he[f] < 0) return;
  // h_i: u_i -> u_{i+1} with length L_i.  Angle at origin(h_i) is the
  // corner between sides L_i (outgoing) and L_{i-1} (incoming), opposite
  // to L_{i+1}.
  int h[3] = { f_he[f], he_next[f_he[f]], he_next[he_next[f_he[f]]] };
  double L[3] = { he_length[h[0]], he_length[h[1]], he_length[h[2]] };
  for (int i = 0; i < 3; i++)
    he_angle[h[i]] = triangle_angle(L[i], L[(i + 2) % 3], L[(i + 1) % 3]);
}

// Recompute every corner angle (he_angle) from the current edge lengths, then
// refresh the per-vertex cone-angle cache (v_cone_angle) that curvature() /
// is_flat() / remove_flat_vertices() read. Both are derived from he_length, so
// this is THE entry point that re-establishes angle coherence after any change
// to the metric (e.g. a conformal rescale of the lengths). Delaunay flips do NOT
// need it: the cone angle is flip-invariant (the quad's interior angle at each
// diamond vertex is independent of which diagonal is present), so flip_edge keeps
// v_cone_angle correct on its own.
void DelaunayTriangulation::recompute_all_angles()
{
  recompute_all_face_angles();
  recompute_cone_angles();
}

// Recompute he_angle for every face, WITHOUT refreshing the v_cone_angle cache.
// For a hot caller (a line search) that reads only he_angle per trial; pair with
// recompute_cone_angles once the kept metric needs the cone cache.
void DelaunayTriangulation::recompute_all_face_angles()
{
  for (int f = 0; f < nf; f++)
    recompute_face_angles(f);
}

// Refresh the cone-angle cache v_cone_angle[v] = vertex_angle_sum(v) for every
// live vertex. Requires he_angle current (call after recompute_face_angles /
// recompute_all_angles). O(sum of degrees) = O(nh).
void DelaunayTriangulation::recompute_cone_angles()
{
  for (int v = 0; v < nv; v++)
    if (v_out[v] >= 0)
      v_cone_angle[v] = vertex_angle_sum(v);
}

double DelaunayTriangulation::vertex_angle_sum(int v) const
{
  // Sum the corner angle at v over all incident faces. he_angle[h] is the
  // angle at origin(h); cw() walks the outgoing half-edges (one per incident
  // face) around v.
  int h0 = v_out[v];
  if (h0 < 0) return 0.0;
  double sum = 0.0;
  int h = h0;
  do { sum += he_angle[h]; h = cw(h); } while (h != h0);
  return sum;
}

bool DelaunayTriangulation::is_flat(int v, double flat_tol) const
{
  // Flat == zero cone curvature == cone angle 2*pi. flat_tol must separate the
  // metric's noise floor from the smallest real cone curvature (pi/3 ~ 1.047);
  // for the equilateral metric any tol in (0, pi/3) agrees exactly with the
  // degree-6 test (deg-6 cone angle = 6*pi/3 = 2*pi).
  return std::abs(v_cone_angle[v] - 2 * M_PI) < flat_tol;
}

// --- Delaunay operations ---

bool DelaunayTriangulation::is_delaunay_edge(int h) const
{
  return diamond(h).is_delaunay();
}

// Flip the diagonal of the diamond around edge h.  Returns true on success.
// Accepts any non-Delaunay edge with a convex diamond, including the
// B == D case (which produces a self-loop edge at B, strictly Delaunay
// by Lemma 1 of CORRECTNESS-PROOF.md).
bool DelaunayTriangulation::flip_edge(int h)
{
  int t = h ^ 1;
  int h1 = he_next[h], h2 = he_next[h1];
  int h4 = he_next[t], h5 = he_next[h4];
  int u = he_origin[h],  v = he_origin[t];
  int B = he_origin[h2], D = he_origin[h5];

  Diamond dm = diamond(h);
  if (!dm.is_convex()) return false;
  double f_len = dm.flipped_length();
  if (!std::isfinite(f_len) || f_len <= 0) return false;

  // Rewire the diagonal: h becomes B->D, t becomes D->B.  Reuse the two
  // face slots fh, ft by rewiring in place (avoids dealloc/realloc).
  int fh = he_face[h], ft = he_face[t];
  he_origin[h] = B;
  he_origin[t] = D;
  he_length[h] = he_length[t] = f_len;

  // Face left of h: (B, D, u) via half-edges h -> h5 -> h1.
  he_next[h]  = h5;  he_next[h5] = h1;  he_next[h1] = h;
  he_face[h]  = he_face[h5] = he_face[h1] = fh;
  f_he[fh] = h;

  // Face left of t: (D, B, v) via half-edges t -> h2 -> h4.
  he_next[t]  = h2;  he_next[h2] = h4;  he_next[h4] = t;
  he_face[t]  = he_face[h2] = he_face[h4] = ft;
  f_he[ft] = t;

  recompute_face_angles(fh);
  recompute_face_angles(ft);

  // u and v lost their incident diagonal; find a new outgoing half-edge.
  if (v_out[u] == h) v_out[u] = h4;
  if (v_out[v] == t) v_out[v] = h1;
  // B == D case: h and t are now self-loops at B; ensure_v_out anchors
  // B's outgoing pointer at a live half-edge.
  if (B == D) ensure_v_out(B);
  return true;
}

bool DelaunayTriangulation::is_delaunay() const
{
  for (int h = 0; h < nh; h += 2)
    if (alive(h) && !is_delaunay_edge(h))
      return false;
  return true;
}

int DelaunayTriangulation::count_non_delaunay() const
{
  int count = 0;
  for (int h = 0; h < nh; h += 2)
    if (alive(h) && !is_delaunay_edge(h))
      count++;
  return count;
}

int DelaunayTriangulation::lawson_sweep()
{
  int flips = 0;

  // Mark all live edges for checking.
  vector<bool> in_stack(nh / 2, false);
  stack<int> S;
  for (int h = 0; h < nh; h += 2)
    if (alive(h)) {
      S.push(h);
      in_stack[h >> 1] = true;
    }

  // The Bobenko-Springborn discrete Dirichlet energy strictly decreases on
  // every flip (Theorem 1, CORRECTNESS-PROOF.md), so the sweep always
  // terminates in finitely many flips on a metrized delta-complex.  The
  // 200*nv budget is a defensive guard, not part of the algorithm: if it
  // fires we have a counterexample to Theorem 1 and must abort loudly
  // rather than return a non-Delaunay output.
  int budget = 200 * nv;
  while (!S.empty()) {
    if (budget <= 0)
      throw std::runtime_error(
          "lawson_sweep: budget exhausted (200*nv = " +
          std::to_string(200 * nv) + " flips); Theorem 1 should preclude "
          "this -- see CORRECTNESS-PROOF.md");

    int h = S.top(); S.pop();
    in_stack[h >> 1] = false;

    if (!alive(h) || is_delaunay_edge(h)) continue;

    // Record rim edges before flipping (they'll be checked next).
    int h1 = he_next[h], h2 = he_next[h1];
    int h4 = he_next[h ^ 1], h5 = he_next[h4];

    if (!flip_edge(h))
      throw std::runtime_error(
          "lawson_sweep: flip_edge rejected a non-Delaunay edge h=" +
          std::to_string(h) + "; Lemma 2 (convex-flippable diamond) "
          "violated -- see CORRECTNESS-PROOF.md");
    flips++; budget--;

    // Push the 4 rim edges.
    for (int rim : {h1, h2, h4, h5}) {
      int eid = rim >> 1;
      if (!in_stack[eid]) { S.push(rim & ~1); in_stack[eid] = true; }
    }
  }
  return flips;
}

int DelaunayTriangulation::flip_to_delaunay()
{
  return lawson_sweep();
}

// ============================================================================
// Vertex removal: data structures
// ============================================================================

// Fan polygon: isometric 2D embedding of a flat vertex's star.
//
// A flat vertex (cone angle = 2pi) has a star that unfolds without overlap
// into a planar disk.  The cumulative angle parameterization gives polar
// coordinates (spokes[i], cum[i]) for each boundary vertex.
//
//         nb[1]
//        / | \
//       /  |  \     spokes[i] = |v - nb[i]|
//      /   |   \    rims[i]   = |nb[i] - nb[(i+1)%k]|
//     v----+----    cum[i]    = sum of fan angles 0..i-1
//      \   |   /
//       \  |  /
//        \ | /
//         nb[0]
//
struct FanPolygon {
  int k;                    // number of fan vertices (= star degree)
  vector<int> nb;           // neighbor vertex IDs in CCW order
  vector<int> spoke_he;     // spoke half-edges: v -> nb[i]
  vector<int> inner_rim;    // inner rim half-edges: nb[i] -> nb[(i+1)%k]
  vector<double> spokes;    // spoke lengths
  vector<double> rims;      // rim edge lengths
  vector<double> cum;       // cumulative fan angles [0, ..., 2pi]

  // 2D fan coordinates of boundary vertex i.
  double x(int i) const { return spokes[i] * cos(cum[i]); }
  double y(int i) const { return spokes[i] * sin(cum[i]); }

  // Diagonal length between fan boundary vertices i and j,
  // computed as Euclidean distance in the isometric 2D embedding.
  double diag_length(int from, int to) const {
    double angle = (to > from) ? cum[to] - cum[from]
                                : (cum[k] - cum[from]) + cum[to];
    double sf = spokes[from], st = spokes[to];
    double len2 = sf*sf + st*st - 2*sf*st*cos(angle);
    return (len2 > 0) ? sqrt(len2) : 0;
  }

  // Signed area of triangle (pp, pi, pn) in fan coordinates.
  // Positive means CCW orientation (valid ear).
  double ear_area(int pp, int pi, int pn) const {
    double rp = spokes[pp], ri = spokes[pi], rn = spokes[pn];
    double tp = cum[pp], ti = cum[pi], tn = cum[pn];
    return rp*ri*sin(ti - tp) + ri*rn*sin(tn - ti) + rn*rp*sin(tp - tn);
  }
};

// Fan triangulation: the result of ear-clipping a fan polygon.
struct FanTriangulation {
  struct Diagonal { int from, ear, to; double length; };
  struct Triangle { int v0, v1, v2; };

  vector<Diagonal> diagonals;   // k-3 ear diagonals
  vector<Triangle> triangles;   // k-2 ear triangles
};

// ============================================================================
// Vertex removal: helper functions
// ============================================================================

// Extract the fan polygon around vertex v.
static FanPolygon extract_fan(const DelaunayTriangulation& D, int v) {
  FanPolygon fan;
  fan.k = D.vertex_degree(v);

  fan.spoke_he.resize(fan.k);
  fan.spoke_he[0] = D.v_out[v];
  for (int i = 1; i < fan.k; i++)
    fan.spoke_he[i] = D.ccw(fan.spoke_he[i-1]);

  fan.nb.resize(fan.k);
  fan.spokes.resize(fan.k);
  fan.rims.resize(fan.k);
  fan.inner_rim.resize(fan.k);
  for (int i = 0; i < fan.k; i++) {
    fan.nb[i] = D.dest(fan.spoke_he[i]);
    fan.spokes[i] = D.he_length[fan.spoke_he[i]];
    fan.inner_rim[i] = D.he_next[fan.spoke_he[i]];
    fan.rims[i] = D.he_length[fan.inner_rim[i]];
  }

  fan.cum.resize(fan.k + 1, 0);
  for (int i = 0; i < fan.k; i++)
    fan.cum[i+1] = fan.cum[i]
                 + triangle_angle(fan.spokes[i], fan.spokes[(i+1)%fan.k], fan.rims[i]);

  return fan;
}

// Ear acceptance test for a candidate ear (pp, pi, pn) in the fan polygon.
// Purely geometric: Meisters-style area/convexity/length.  Self-loop and
// multi-edge diagonals are legal delta-complex edges; splice_fan wires
// them independently via their polygon-position keys.
// Returns the diagonal length if acceptable, else 0.
static double ear_length_if_acceptable(
    const FanPolygon& fan, int pp, int pi, int pn)
{
  if (fan.ear_area(pp, pi, pn) <= 1e-10) return 0;
  double sub = (pn > pp) ? fan.cum[pn] - fan.cum[pp]
                         : (fan.cum[fan.k] - fan.cum[pp]) + fan.cum[pn];
  if (sub > M_PI + 1e-10) return 0;
  double len = fan.diag_length(pp, pn);
  return (len > 1e-15) ? len : 0;
}

// Ear-clip the fan polygon into triangles.  By Meisters' theorem
// (Lemma 4 of the proof), a simple polygon with k >= 4 always has an
// ear; the scan below therefore terminates at k = 3.  If a pass finds no
// ear we throw rather than return silently: that would be a Meisters
// counterexample, not a recoverable condition.
static FanTriangulation ear_clip_fan(const FanPolygon& fan) {
  int k = fan.k;
  FanTriangulation tri;

  vector<int> poly(k);
  for (int i = 0; i < k; i++) poly[i] = i;

  while ((int)poly.size() > 3) {
    int n = poly.size();
    bool clipped = false;

    for (int j = 0; j < n; j++) {
      int pp = poly[(j - 1 + n) % n], pi = poly[j], pn = poly[(j + 1) % n];
      double len = ear_length_if_acceptable(fan, pp, pi, pn);
      if (len <= 0) continue;
      tri.diagonals.push_back({pp, pi, pn, len});
      poly.erase(poly.begin() + j);
      clipped = true;
      break;
    }

    if (!clipped)
      throw std::runtime_error(
          "ear_clip_fan: no ear found at remaining polygon size " +
          std::to_string(poly.size()) + " (Meisters / Lemma 4 violated) "
          "-- see CORRECTNESS-PROOF.md");
  }

  // Compose diagonals into triangle list, appending the final base triangle.
  vector<int> rpoly(k);
  for (int i = 0; i < k; i++) rpoly[i] = i;
  for (auto& ear : tri.diagonals) {
    tri.triangles.push_back({ear.from, ear.ear, ear.to});
    rpoly.erase(std::find(rpoly.begin(), rpoly.end(), ear.ear));
  }
  if (rpoly.size() != 3)
    throw std::runtime_error("ear_clip_fan: residual polygon size != 3");
  tri.triangles.push_back({rpoly[0], rpoly[1], rpoly[2]});

  return tri;
}

static void splice_fan(DelaunayTriangulation& D, int v,
                        const FanPolygon& fan, const FanTriangulation& tri) {
  int k = fan.k;

  // --- Fix v_out for neighbors before deallocation ---
  for (int i = 0; i < k; i++) {
    int spoke_twin = fan.spoke_he[i] ^ 1;
    if (D.v_out[fan.nb[i]] == spoke_twin)
      D.v_out[fan.nb[i]] = fan.inner_rim[(i - 1 + k) % k] ^ 1;
  }

  // --- Deallocate fan faces and spoke edges ---
  for (int i = 0; i < k; i++) {
    D.dealloc_face(D.he_face[fan.spoke_he[i]]);
    D.dealloc_edge(fan.spoke_he[i]);
  }
  D.v_out[v] = -1;

  // --- Build arc-to-halfedge map (inner rim + new diagonals) ---
  map<pair<int,int>, int> local_arc;
  for (int i = 0; i < k; i++)
    local_arc[{i, (i + 1) % k}] = fan.inner_rim[i];
  for (auto& d : tri.diagonals) {
    int h_d = D.alloc_directed_edge(fan.nb[d.from], fan.nb[d.to], d.length);
    local_arc[{d.from, d.to}] = h_d;
    local_arc[{d.to, d.from}] = h_d ^ 1;
  }

  // --- Wire each ear triangle ---
  for (auto& t : tri.triangles)
    D.wire_triangle(
      local_arc.at({t.v0, t.v1}),
      local_arc.at({t.v1, t.v2}),
      local_arc.at({t.v2, t.v0}));

  for (int i = 0; i < k; i++) D.ensure_v_out(fan.nb[i]);
}

// Remove a degree-3 vertex: three fan faces merge into one triangle.
static void remove_degree3(DelaunayTriangulation& D, int v) {
  int h0 = D.v_out[v], h1 = D.ccw(h0), h2 = D.ccw(h1);
  int f0 = D.he_face[h0], f1 = D.he_face[h1], f2 = D.he_face[h2];
  int inner0 = D.he_next[h0], inner1 = D.he_next[h1], inner2 = D.he_next[h2];

  // Snapshot the neighbour vertex ids BEFORE deallocation (dest reads he_origin
  // of the twin half-edge, which dealloc_edge clears to -1).
  int nb0 = D.dest(h0), nb1 = D.dest(h1), nb2 = D.dest(h2);

  D.dealloc_face(f0); D.dealloc_face(f1); D.dealloc_face(f2);
  D.dealloc_edge(h0); D.dealloc_edge(h1); D.dealloc_edge(h2);
  D.v_out[v] = -1;

  // Repair v_out for each neighbour whose outgoing pointer was a killed spoke.
  for (int nbr : {nb0, nb1, nb2}) D.ensure_v_out(nbr);

  D.wire_triangle(inner0, inner1, inner2);
}

// ============================================================================
// Vertex removal: main entry point
// ============================================================================

void DelaunayTriangulation::remove_flat_vertex(int v)
{
  int deg = vertex_degree(v);
  if (deg == 0) return;

  if (deg >= 4) {
    // Phase 2: Ear clipping + DCEL surgery.
    //   extract_fan:  read star geometry -> FanPolygon (isometric 2D embedding)
    //   ear_clip_fan: triangulate the fan polygon (Meisters; throws on failure)
    //   splice_fan:   replace the star with the triangulation (DCEL surgery)
    FanPolygon fan = extract_fan(*this, v);
    FanTriangulation tri = ear_clip_fan(fan);
    splice_fan(*this, v, fan, tri);
  } else if (deg == 3) {
    // Phase 3: Direct removal (three faces merge into one triangle).
    remove_degree3(*this, v);
  }
}

// Flip away all self-loops at vertex v.
// Self-loops at a flat vertex arise from ear diagonals in previous
// removals; they must be cleared before remove_flat_vertex, otherwise
// extract_fan sees v in its own polygon and splice_fan would wire a
// live edge to the about-to-be-dead vertex.  Correctness obligation
// (CORRECTNESS-PROOF.md, Theorem 3): every self-loop at a flat v is
// flippable, i.e. its diamond is convex at v.  Empirically true on
// 1.94B+ fullerene isomers; in the adversarial case of a strictly
// Delaunay self-loop at a flat vertex the assert below catches it.
static void flip_away_self_loops(DelaunayTriangulation& D, int v) {
  if (D.v_out[v] < 0) return;
  bool flipped_any = true;
  while (flipped_any) {
    flipped_any = false;
    int h0 = D.v_out[v];
    if (h0 < 0) break;
    int h = h0;
    do {
      if (D.dest(h) == v) {
        if (D.flip_edge(h))     { flipped_any = true; break; }
        if (D.flip_edge(h ^ 1)) { flipped_any = true; break; }
      }
      h = D.cw(h);
    } while (h != h0);
  }
  // Invariant check (runtime, not assert: asserts compile out with -DNDEBUG):
  // no self-loop survives at a flat vertex.  If one does, splice_fan would
  // double-deallocate the self-loop edge and corrupt the DCEL silently.
  // Theorem in CORRECTNESS-PROOF.md proves this never fires on
  // non-Delaunay or freshly-created self-loops; the residual is the
  // rim-flip-evolution case, which is empirically clean across all
  // tested inputs.
  int h0 = D.v_out[v];
  if (h0 >= 0) {
    int h = h0;
    do {
      if (D.dest(h) == v)
        throw std::runtime_error(
            "flip_away_self_loops: un-flippable self-loop at flat v=" +
            std::to_string(v) + " (Obligation 1 violated; "
            "see CORRECTNESS-PROOF.md)");
      h = D.cw(h);
    } while (h != h0);
  }
}

void DelaunayTriangulation::remove_flat_vertices(double flat_tol)
{
  // Scan all live flat vertices top-down; any single removable vertex keeps
  // the algorithm moving.  Repeat until a pass produces no removals, then
  // apply a Delaunay restructuring sweep and try once more.  Only after
  // two consecutive fruitless passes do we declare the graph truly stuck.
  //
  // Rationale: removing in strict descending-index order aborts as soon as
  // the highest-index flat vertex resists, even when other flat vertices
  // are still removable — which strands labeling-dependent multi-edge
  // clusters and leaves the iDT with > 12 vertices.  A full-scan pass is
  // O(N) per iteration and has no worse asymptotic cost than the old loop.
  // Progress reporting (when verbose_removal > 0). Counts live vertices and
  // reports every `verbose_removal` removals. fprintf(stderr) is unbuffered;
  // we fflush each line so it shows on a slow run (don't pipe — piping
  // re-introduces block buffering).
  auto count_live = [&]() {
    int n = 0;
    for (int v = 0; v < nv; v++) if (v_out[v] >= 0) n++;
    return n;
  };
  long long removed = 0;
  if (verbose_removal) {
    std::fprintf(stderr, "[remove_flat] start: %d live vertices\n", count_live());
    std::fflush(stderr);
  }

  auto remove_any_flat = [&]() {
    bool progress = false;
    for (int v = nv - 1; v >= 0; v--) {
      if (v_out[v] < 0 || !is_flat(v, flat_tol)) continue;
      flip_away_self_loops(*this, v);
      remove_flat_vertex(v);
      if (v_out[v] < 0) {
        progress = true;
        while (nv > 0 && v_out[nv - 1] < 0) nv--;
        lawson_sweep();
        if (verbose_removal && (++removed % verbose_removal == 0)) {
          std::fprintf(stderr, "[remove_flat] removed %lld, %d live remain\n",
                       removed, count_live());
          std::fflush(stderr);
        }
      }
    }
    return progress;
  };

  while (true) {
    if (remove_any_flat()) continue;
    // Nothing removable this pass; restructure via full Delaunay sweep and retry.
    flip_to_delaunay();
    if (!remove_any_flat()) break;
  }

  if (verbose_removal) {
    std::fprintf(stderr, "[remove_flat] done: removed %lld, %d live remain\n",
                 removed, count_live());
    std::fflush(stderr);
  }

  // Final Lawson.  With guards dropped, this converges unconditionally by
  // Theorem 1 (Bobenko-Springborn energy strictly decreases on every flip,
  // including B == D self-loop-creating flips; the new self-loop is
  // strictly Delaunay by Lemma 1).
  flip_to_delaunay();
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
  v_out = std::move(nout);
  v_cone_angle = std::move(ncone);
  v_orig_degree = std::move(ndeg);
  nv = nlive;
  return new_to_old;
}

// --- Full algorithm ---

int DelaunayTriangulation::min_live_degree() const
{
  int m = INT_MAX;
  for (int v = 0; v < nv; v++) {
    if (v_out[v] < 0) continue;
    int d = vertex_degree(v);
    if (d < m) m = d;
  }
  return m;
}

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
  // Every live half-edge belongs to exactly one face cycle (walked via
  // he_next), and every such cycle has length 3.  DCEL counterpart of
  // Graph::is_consistently_oriented: there an arc u->v is consistent
  // if it appears in exactly one face walked via next(v, u); here a
  // half-edge h is consistent if it appears in exactly one face
  // walked via he_next.
  std::vector<bool> visited(nh, false);
  for (int h0 = 0; h0 < nh; h0++) {
    if (!alive(h0) || visited[h0]) continue;
    int h = h0;
    for (int step = 0; step < 3; step++) {
      if (!alive(h) || visited[h]) return false;
      visited[h] = true;
      h = he_next[h];
    }
    if (h != h0) return false;             // cycle longer than 3
  }
  for (int h = 0; h < nh; h++)
    if (alive(h) && !visited[h]) return false;
  return true;
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
  Triangulation sorted = T;
  sorted.apply_permutation(T.sort_flat_last());
  DelaunayTriangulation D = from_triangulation(sorted);
  D.remove_flat_vertices();
  return D;
}

DelaunayTriangulation DelaunayTriangulation::compute(const Triangulation& T,
                                                     const EdgeLengthFn& length,
                                                     double flat_tol)
{
  // Prescribed-metric iDT. Unlike the equilateral compute(T), we do NOT
  // sort_flat_last() (that classifies flatness by degree, which is wrong for
  // a general metric, and would also permute the vertices the caller's
  // `length` oracle is keyed on). remove_flat_vertices() identifies flat
  // vertices by cone angle (is_flat) and its full-scan loop is order-robust,
  // so no pre-sort is needed.
  DelaunayTriangulation D = from_intrinsic_metric(T, length);
  D.remove_flat_vertices(flat_tol);
  D.compact_vertices();  // cones are scattered (no flat-last sort); compact to nv=cones
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
  return diamond(h).is_cocircular();
}

vector<bool> DelaunayTriangulation::cocircular_edges() const
{
  vector<bool> tight(nh, false);
  for (int h = 0; h < nh; h += 2) {
    if (!alive(h)) continue;
    bool t = diamond(h).is_cocircular();
    tight[h]     = t;
    tight[h ^ 1] = t;
  }
  return tight;
}

vector<bool> DelaunayTriangulation::cocircular_edges(double tol) const
{
  vector<bool> tight(nh, false);
  for (int h = 0; h < nh; h += 2) {
    if (!alive(h)) continue;
    bool t = diamond(h).is_cocircular(tol);
    tight[h]     = t;
    tight[h ^ 1] = t;
  }
  return tight;
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
        if (++safety > nh) {
          // Topology gone wrong -- bail out with an empty cell so the
          // caller can detect and report the problem.
          T.cells.clear();
          return T;
        }
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
  int t = h ^ 1;
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

  D.v_out[w] = uw ^ 1;  // w -> u

  // Wire four new CCW faces around w.
  D.wire_triangle(uw,     wB_he,   h_Bu);     // (u, w, B)
  D.wire_triangle(wv,     h_vB,    wB_he^1);  // (w, v, B)
  D.wire_triangle(wv^1,   wD_he,   h_Dv);     // (v, w, D)
  D.wire_triangle(uw^1,   h_uD,    wD_he^1);  // (w, u, D)

  // Fix v_out for u, v if they pointed at the deleted edge.
  if (D.v_out[u] == h) D.v_out[u] = uw;
  if (D.v_out[v] == t) D.v_out[v] = wv ^ 1;

  return w;
}

int DelaunayTriangulation::bisect_multi_edges() {
  // Find multi-edges: vertex pairs with >1 edge.
  map<pair<int,int>, vector<int>> pair_to_hes;
  for (int h = 0; h < nh; h += 2) {
    if (!alive(h)) continue;
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
  return n_inserted;
}

// --- Vertex / face surgery ---

int DelaunayTriangulation::alloc_vertex(double cone_angle, int orig_degree) {
  int v = nv++;
  if (v >= (int)v_out.size()) {
    v_out.push_back(-1);
    v_cone_angle.push_back(cone_angle);
    v_orig_degree.push_back(orig_degree);
  } else {
    v_out[v] = -1;
    v_cone_angle[v] = cone_angle;
    v_orig_degree[v] = orig_degree;
  }
  return v;
}

int DelaunayTriangulation::split_face(int h0, std::array<double,3> spoke) {
  int h1 = he_next[h0], h2 = he_next[h1];
  int a = he_origin[h0], b = he_origin[h1], c = he_origin[h2];
  int f = he_face[h0];

  int P = alloc_vertex(2.0 * M_PI, 3);
  dealloc_face(f);                                  // only the face; the 3 edges stay
  int sa = alloc_directed_edge(P, a, spoke[0]);
  int sb = alloc_directed_edge(P, b, spoke[1]);
  int sc = alloc_directed_edge(P, c, spoke[2]);
  wire_triangle(sa, h0, twin(sb));                  // (P, a, b)
  wire_triangle(sb, h1, twin(sc));                  // (P, b, c)
  wire_triangle(sc, h2, twin(sa));                  // (P, c, a)

  v_out[P] = sa;
  v_out[a] = twin(sa); v_out[b] = twin(sb); v_out[c] = twin(sc);
  // wire_triangle set the new faces' angles; refresh the cached cone angle at the four
  // vertices whose incident faces changed.
  for (int v : {P, a, b, c}) v_cone_angle[v] = vertex_angle_sum(v);
  return P;
}

// --- Curvature ---

std::vector<double> DelaunayTriangulation::curvature() const {
  std::vector<double> K(nv);
  for (int v = 0; v < nv; v++) K[v] = 2.0 * M_PI - v_cone_angle[v];
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
  for (int s : sources) disks.push_back({s, {}});
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
  return disks;
}

// --- Validation ---

bool DelaunayTriangulation::check_consistency() const
{
  // 1. Twin pairs: twin(twin(h)) == h.
  for (int h = 0; h < nh; h++)
    if (alive(h) && (h ^ 1) >= nh) return false;

  // 2. Next-cycle closure: following next 3 times returns to start (triangulation).
  for (int h = 0; h < nh; h++)
    if (alive(h) && he_next[he_next[he_next[h]]] != h) return false;

  // 3. dest(h) == origin(next(h)).
  for (int h = 0; h < nh; h++)
    if (alive(h) && dest(h) != he_origin[he_next[h]]) return false;

  // 4. Face consistency: all half-edges in a face have the same face ID.
  for (int h = 0; h < nh; h++)
    if (alive(h) && he_face[he_next[h]] != he_face[h]) return false;

  // 5. v_out: for each live vertex, v_out points to a live outgoing half-edge.
  for (int v = 0; v < nv; v++)
    if (v_out[v] >= 0 && (!alive(v_out[v]) || he_origin[v_out[v]] != v))
      return false;

  // 6. f_he: for each live face, f_he points to a half-edge in that face.
  for (int f = 0; f < nf; f++)
    if (f_he[f] >= 0 && (!alive(f_he[f]) || he_face[f_he[f]] != f))
      return false;

  // 7. Positive edge lengths for live half-edges.
  for (int h = 0; h < nh; h++)
    if (alive(h) && he_length[h] <= 0) return false;

  // 8. Twin length consistency.
  for (int h = 0; h < nh; h += 2)
    if (alive(h) && he_length[h] != he_length[h ^ 1]) return false;

  // 9. Triangle inequality on every live face.
  for (int h = 0; h < nh; h++) {
    if (!alive(h)) continue;
    int h1 = he_next[h], h2 = he_next[h1];
    double L0 = he_length[h], L1 = he_length[h1], L2 = he_length[h2];
    double eps = 1e-9 * (L0 + L1 + L2);
    if (L0 > L1 + L2 + eps) return false;
    if (L1 > L0 + L2 + eps) return false;
    if (L2 > L0 + L1 + eps) return false;
  }

  return true;
}

