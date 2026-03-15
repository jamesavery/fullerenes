#include "fullerenes/delaunay.hh"
#include "fullerenes/eisenstein.hh"

#include <stack>
#include <queue>
#include <cmath>
#include <cassert>
#include <algorithm>
#include <numeric>
#include <map>
#include <array>
#include <unordered_set>
#include <unordered_map>

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

// ============================================================================
// Constructor
// ============================================================================

FulleroidDelaunay::FulleroidDelaunay(const Triangulation& T)
  : Triangulation(T.sort_flat_last()), edge_lengths(N, N, 0)
{
  // Edge flips during vertex removal can temporarily push vertex degrees
  // well above 6 (the max for fullerene duals). Restride to give headroom.
  if (dmax < 20) {
    auto restrided = restride(20);
    owned_values = std::move(restrided.owned_values);
    owned_deg = std::move(restrided.owned_deg);
    dmax = 20;
    repoint();
  }
  for (node_t u = 0; u < N; u++)
    for (node_t v : (*this)[u])
      edge_lengths(u, v) = 1.0;
}

// ============================================================================
// Diamond extraction and Delaunay operations
// ============================================================================

Diamond FulleroidDelaunay::diamond(node_t u, node_t v) const
{
  node_t B = next(v, u), D = next(u, v);
  Diamond d = {get_length(u,v), get_length(u,B), get_length(v,B),
               get_length(u,D), get_length(v,D)};
  assert(d.e > 0 && d.a > 0 && d.b > 0 && d.c > 0 && d.d > 0);
  return d;
}

bool FulleroidDelaunay::is_delaunay_edge(node_t u, node_t v) const
{
  return diamond(u, v).is_delaunay();
}

bool FulleroidDelaunay::flip_edge(node_t u, node_t v, bool verbose)
{
  node_t B = next(v, u), D = next(u, v);

  // Topological guards
  if (B == D) { if (verbose) fprintf(stderr, "  flip(%d,%d): B==D==%d\n", u, v, B); return false; }
  if (edge_exists(edge_t(B, D))) { if (verbose) fprintf(stderr, "  flip(%d,%d): multi-edge (%d,%d)\n", u, v, B, D); return false; }

  // Geometric guards
  auto d = diamond(u, v);
  if (!d.is_convex()) { if (verbose) fprintf(stderr, "  flip(%d,%d): not convex\n", u, v); return false; }
  double f = d.flipped_length();
  if (!std::isfinite(f) || f <= 0) { if (verbose) fprintf(stderr, "  flip(%d,%d): bad length f=%g\n", u, v, f); return false; }

  // Execute flip: remove diagonal u-v, insert diagonal B-D.
  Graph::remove_edge(edge_t(u, v));
  set_length(u, v, 0);
  Graph::insert_edge(arc_t(B, D), u, v);
  set_length(B, D, f);

  if (audit) audit->after_flip(*this, u, v);
  return true;
}

int FulleroidDelaunay::flip_to_delaunay()
{
  int flips = lawson_sweep();
  if (!is_delaunay())
    flips += delaunay_resolve();
  return flips;
}

bool FulleroidDelaunay::is_delaunay() const
{
  for (node_t u = 0; u < N; u++)
    for (node_t v : (*this)[u])
      if (u < v && !is_delaunay_edge(u, v))
        return false;
  return true;
}

int FulleroidDelaunay::count_non_delaunay() const
{
  int count = 0;
  for (node_t u = 0; u < N; u++)
    for (node_t v : (*this)[u])
      if (u < v && !is_delaunay_edge(u, v))
        count++;
  return count;
}

int FulleroidDelaunay::delaunay_resolve()
{
  // Search for a sequence of flips that reduces the non-Delaunay edge count.
  //
  // On some surfaces, the iDT is non-simplicial (has multi-edges between
  // the same pair of cone points).  In that case, no simplicial triangulation
  // is fully Delaunay, and this function cannot reach 0 non-Delaunay edges.
  // It still improves the result when possible (e.g., reducing 3 → 1).

  int total_flips = 0;

  for (int round = 0; round < 20; round++) {
    int count0 = count_non_delaunay();
    if (count0 == 0) break;

    FulleroidDelaunay saved = *this;
    bool improved = false;

    // Collect all edges from saved state.
    vector<edge_t> edges;
    for (node_t u = 0; u < N; u++)
      for (node_t v : saved[u])
        if (u < v) edges.push_back(edge_t(u, v));

    // Depth 1: try every single edge flip.
    for (size_t i = 0; i < edges.size() && !improved; i++) {
      *this = saved;
      if (!flip_edge(edges[i].first, edges[i].second)) continue;
      int f = 1 + lawson_sweep();
      if (count_non_delaunay() < count0) {
        total_flips += f;
        improved = true;
      }
    }

    // Depth 2: try all pairs of flips.
    for (size_t i = 0; i < edges.size() && !improved; i++) {
      for (size_t j = i + 1; j < edges.size() && !improved; j++) {
        *this = saved;
        if (!flip_edge(edges[i].first, edges[i].second)) continue;
        if (!edge_exists(edges[j]) ||
            !flip_edge(edges[j].first, edges[j].second)) continue;
        int f = 2 + lawson_sweep();
        if (count_non_delaunay() < count0) {
          total_flips += f;
          improved = true;
        }
      }
    }

    if (!improved) {
      *this = saved;  // restore
      break;
    }
  }

  return total_flips;
}

// ============================================================================
// Vertex removal
// ============================================================================

void FulleroidDelaunay::remove_flat_vertex(node_t v)
{
  // Phase 1: Reduce v's degree via incident edge flips.
  // Stop at degree 4 — the ear clipping in Phase 2 handles degree >= 4.
  while (degree(v) > 4) {
    bool progress = false;
    vector<node_t> vnbrs((*this)[v].begin(), (*this)[v].end());
    for (node_t u : vnbrs)
      if (flip_edge(v, u)) { progress = true; break; }
    if (!progress) break;
  }

  int deg = degree(v);

  // Phase 2: Star retriangulation for degree >= 4 via intrinsic ear clipping.
  // Unfolds the fan of triangles around v into a planar polygon using
  // cumulative angles and spoke lengths, then clips convex ears.
  if (deg >= 4) {
    int k = deg;
    auto row = (*this)[v];
    vector<node_t> nb(row.begin(), row.end());

    // Gather spoke and rim lengths.
    vector<double> spokes(k), rims(k);
    for (int i = 0; i < k; i++) {
      spokes[i] = get_length(v, nb[i]);
      rims[i]   = get_length(nb[i], nb[(i + 1) % k]);
    }

    // Cumulative angle at v: cum[i] = sum of fan-triangle angles from 0 to i-1.
    vector<double> cum(k + 1, 0);
    for (int i = 0; i < k; i++) {
      double si = spokes[i], sn = spokes[(i+1)%k], r = rims[i];
      double cos_a = std::clamp((si*si + sn*sn - r*r) / (2*si*sn), -1.0, 1.0);
      cum[i+1] = cum[i] + acos(cos_a);
    }

    // Diagonal length from nb[from] to nb[to] through the CCW sub-fan.
    auto diag_len = [&](int from, int to) -> double {
      double angle = (to > from) ? cum[to] - cum[from]
                                  : (cum[k] - cum[from]) + cum[to];
      double sf = spokes[from], st = spokes[to];
      double len2 = sf*sf + st*st - 2*sf*st*cos(angle);
      return (len2 > 0) ? sqrt(len2) : 0;
    };

    // Signed area of triangle (P_pp, P_pi, P_pn) in the 2D fan unfolding,
    // where P_j = (spokes[j]*cos(cum[j]), spokes[j]*sin(cum[j])).
    // Positive means the polygon vertex at pi is convex (valid ear).
    auto ear_signed_area = [&](int pp, int pi, int pn) -> double {
      double rp = spokes[pp], ri = spokes[pi], rn = spokes[pn];
      double tp = cum[pp], ti = cum[pi], tn = cum[pn];
      return rp*ri*sin(ti - tp) + ri*rn*sin(tn - ti) + rn*rp*sin(tp - tn);
    };

    // Ear clipping with intrinsic geometry.
    struct EarDiag { int prev_idx, ear_idx, next_idx; double len; };
    vector<EarDiag> diagonals;

    // Polygon as ring indices + edge lengths.
    vector<int> poly(k);
    vector<double> poly_edge(k);
    for (int i = 0; i < k; i++) { poly[i] = i; poly_edge[i] = rims[i]; }

    while ((int)poly.size() > 3) {
      int n = poly.size();
      bool found = false;
      for (int j = 0; j < n; j++) {
        int jm = (j - 1 + n) % n, jp = (j + 1) % n;
        int pp = poly[jm], pi = poly[j], pn = poly[jp];

        // Topological checks.
        if (nb[pp] == nb[pn]) continue;
        if (edge_exists(edge_t(nb[pp], nb[pn]))) continue;
        bool dup = false;
        for (auto& d : diagonals)
          if (edge_t(nb[d.prev_idx], nb[d.next_idx]) == edge_t(nb[pp], nb[pn]))
            { dup = true; break; }
        if (dup) continue;

        // Convexity check: polygon vertex must be convex (positive signed area).
        if (ear_signed_area(pp, pi, pn) <= 1e-10) continue;

        // Sub-fan angle check: the ear's fan arc from pp to pn (through pi)
        // must be <= pi.  This ensures we always clip a "small" ear.
        double sub_angle = (pn > pp) ? cum[pn] - cum[pp]
                                      : (cum[k] - cum[pp]) + cum[pn];
        if (sub_angle > M_PI + 1e-10) continue;

        // Diagonal length via sub-fan angle + law of cosines.
        double len = diag_len(pp, pn);
        if (len <= 1e-15) continue;

        diagonals.push_back({pp, pi, pn, len});
        poly.erase(poly.begin() + j);
        poly_edge.erase(poly_edge.begin() + j);
        poly_edge[(j > 0) ? j - 1 : (int)poly.size() - 1] = len;
        found = true;
        break;
      }
      if (!found) {
        // Blocker-flip: try flipping edges that block ear diagonals.
        bool flipped_blocker = false;
        int n2 = poly.size();
        for (int j = 0; j < n2 && !flipped_blocker; j++) {
          int jm = (j - 1 + n2) % n2, jp = (j + 1) % n2;
          int pp = poly[jm], pn = poly[jp];
          if (nb[pp] != nb[pn] && edge_exists(edge_t(nb[pp], nb[pn]))) {
            if (flip_edge(nb[pp], nb[pn]))
              flipped_blocker = true;
          }
        }
        if (!flipped_blocker) return;  // truly stuck
      }
    }

    // Remove all spoke edges.
    for (int i = 0; i < k; i++) set_length(v, nb[i], 0);
    for (int i = k - 1; i >= 0; i--) Graph::remove_edge(edge_t(v, nb[i]));

    // Insert ear diagonals.
    for (size_t di = 0; di < diagonals.size(); di++) {
      auto& d = diagonals[di];
      node_t p_prev = nb[d.prev_idx];
      node_t p_ear  = nb[d.ear_idx];
      node_t p_next = nb[d.next_idx];
      // Insert diagonal p_prev--p_next so that the ear face
      // (p_prev, p_ear, p_next) is created with correct CCW orientation.
      node_t suc_uv = next(p_prev, p_ear);
      node_t suc_vu = p_ear;
      assert(find(p_next, p_ear) >= 0);
      Graph::insert_edge(arc_t(p_prev, p_next), suc_uv, suc_vu);
      set_length(p_prev, p_next, d.len);
    }

    assert(v == N - 1 && degree(v) == 0);
    pop_back();
    if (audit) audit->after_removal(*this, v);
    return;
  }

  // Phase 3: Degree 3 — the three surrounding triangles collapse into one.
  assert(deg == 3);
  auto vrow = (*this)[v];
  node_t a = vrow[0], b = vrow[1], c = vrow[2];
  set_length(v, a, 0); set_length(v, b, 0); set_length(v, c, 0);
  Graph::remove_edge(edge_t(v, a));
  Graph::remove_edge(edge_t(v, b));
  Graph::remove_edge(edge_t(v, c));

  assert(v == N - 1);
  pop_back();
  if (audit) audit->after_removal(*this, v);
}

// Standard Lawson sweep: flip all flippable non-Delaunay edges.
// Terminates because the Delaunay energy E_D strictly decreases with each flip
// (Bobenko-Springborn 2005) and there are finitely many triangulations.
int FulleroidDelaunay::lawson_sweep()
{
  int flips = 0;
  map<edge_t, bool> in_stack;
  stack<edge_t> S;
  for (node_t u = 0; u < N; u++)
    for (node_t v : (*this)[u])
      if (u < v) {
        edge_t e(u, v);
        S.push(e); in_stack[e] = true;
      }

  int budget = 200 * N;
  while (!S.empty() && budget > 0) {
    edge_t e = S.top(); S.pop();
    in_stack[e] = false;
    node_t u = e.first, v = e.second;
    if (!edge_exists(e) || is_delaunay_edge(u, v)) continue;

    node_t B = next(v, u), D = next(u, v);
    if (!flip_edge(u, v)) continue;

    flips++; budget--;
    for (edge_t ec : {edge_t(u,B), edge_t(B,v), edge_t(v,D), edge_t(D,u)})
      if (!in_stack[ec]) { S.push(ec); in_stack[ec] = true; }
  }
  if (audit) audit->after_sweep(*this, flips);
  return flips;
}

void FulleroidDelaunay::remove_flat_vertices()
{
  vector<int> original_degrees(N);
  for (node_t v = 0; v < N; v++)
    original_degrees[v] = degree(v);

  while (N > 0 && original_degrees[N - 1] == 6) {
    int old_N = N;
    remove_flat_vertex(N - 1);

    if (N == old_N) {
      // Vertex removal stuck — restructure via Delaunay flips and retry.
      bool removed = false;
      for (int retry = 0; retry < 5; retry++) {
        flip_to_delaunay();
        remove_flat_vertex(N - 1);
        if (N < old_N) { removed = true; break; }
      }
      if (!removed) break;  // truly stuck
    }

    lawson_sweep();  // maintain near-Delaunay between removals
  }

  // Final: lawson_sweep + delaunay_resolve on the cone-points-only graph.
  flip_to_delaunay();
}

// ============================================================================
// Validation
// ============================================================================

bool FulleroidDelaunay::edge_lengths_are_symmetric() const
{
  for (int i = 0; i < N; i++)
    for (int j = i; j < N; j++)
      if (edge_lengths(i,j) != edge_lengths(j,i))
        return false;
  return true;
}

// ============================================================================
// IDTAudit — postcondition checking for iDT operations
// ============================================================================

IDTAudit::IDTAudit(const FulleroidDelaunay& D)
  : original_degrees(D.N), original_faces(0), expected_area(0)
{
  for (node_t v = 0; v < D.N; v++)
    original_degrees[v] = D.degree(v);

  // Count faces: each triangle (u,v,w) with u<v<w is counted once.
  for (node_t u = 0; u < D.N; u++)
    for (size_t j = 0; j < D[u].size(); j++) {
      node_t v = D[u][j], w = D[u][(j+1) % D[u].size()];
      if (u < v && u < w) original_faces++;
    }

  expected_area = original_faces * sqrt(3.0) / 4.0;
}

void IDTAudit::after_flip(const FulleroidDelaunay& D, node_t u, node_t v) {
  char buf[64];
  snprintf(buf, sizeof(buf), "flip(%d,%d)", u, v);
  check_all(D, buf);
}

void IDTAudit::after_removal(const FulleroidDelaunay& D, node_t removed) {
  char buf[64];
  snprintf(buf, sizeof(buf), "remove(%d)", removed);
  check_all(D, buf);
}

void IDTAudit::after_sweep(const FulleroidDelaunay& D, int n_flips) {
  char buf[64];
  snprintf(buf, sizeof(buf), "sweep(%d flips)", n_flips);
  check_all(D, buf);
}

void IDTAudit::check_all(const FulleroidDelaunay& D, const char* ctx)
{
  verify_euler(D, ctx);
  verify_orientation(D, ctx);
  verify_edge_consistency(D, ctx);
  verify_triangle_inequality(D, ctx);
  verify_positive_area(D, ctx);
  verify_cone_angles(D, ctx);
  verify_total_area(D, ctx);
  verify_loeschian(D, ctx);
  verify_no_multi_edges(D, ctx);
}

void IDTAudit::fail(const char* invariant, const char* ctx, const string& detail) {
  n_failures++;
  fprintf(stderr, "IDT INVARIANT FAILURE: %s\n  Context: %s\n  Detail: %s\n",
          invariant, ctx, detail.c_str());
  if (stop_on_failure) abort();
}

// --- Individual invariant checks ---

bool IDTAudit::verify_euler(const FulleroidDelaunay& D, const char* ctx) {
  n_checks++;
  int V = D.N;
  int E = 0;
  for (node_t u = 0; u < D.N; u++) E += D.degree(u);
  E /= 2;
  int F = 0;
  for (node_t u = 0; u < D.N; u++)
    for (size_t j = 0; j < D[u].size(); j++) {
      node_t v = D[u][j], w = D[u][(j+1) % D[u].size()];
      if (u < v && u < w) F++;
    }
  if (V - E + F != 2) {
    char buf[128];
    snprintf(buf, sizeof(buf), "V=%d E=%d F=%d, V-E+F=%d (expected 2)", V, E, F, V-E+F);
    fail("euler", ctx, buf);
    return false;
  }
  return true;
}

bool IDTAudit::verify_orientation(const FulleroidDelaunay& D, const char* ctx) {
  n_checks++;
  if (!D.is_consistently_oriented()) {
    fail("orientation", ctx, "is_consistently_oriented() returned false");
    return false;
  }
  return true;
}

bool IDTAudit::verify_edge_consistency(const FulleroidDelaunay& D, const char* ctx) {
  n_checks++;
  // Symmetry
  for (node_t u = 0; u < D.N; u++)
    for (node_t v = u+1; v < D.N; v++)
      if (D.edge_lengths(u,v) != D.edge_lengths(v,u)) {
        char buf[128];
        snprintf(buf, sizeof(buf), "edge_lengths(%d,%d)=%g != edge_lengths(%d,%d)=%g",
                 u, v, D.edge_lengths(u,v), v, u, D.edge_lengths(v,u));
        fail("edge_symmetry", ctx, buf);
        return false;
      }
  // Adjacency match
  for (node_t u = 0; u < D.N; u++)
    for (node_t v : D[u])
      if (D.get_length(u,v) <= 0) {
        char buf[128];
        snprintf(buf, sizeof(buf), "edge(%d,%d) exists but length=%g", u, v, D.get_length(u,v));
        fail("edge_positive", ctx, buf);
        return false;
      }
  for (node_t u = 0; u < D.N; u++)
    for (node_t v = 0; v < D.N; v++)
      if (u != v && D.edge_lengths(u,v) > 0 &&
          std::find(D[u].begin(), D[u].end(), v) == D[u].end()) {
        char buf[128];
        snprintf(buf, sizeof(buf), "edge_lengths(%d,%d)=%g but not in neighbour list",
                 u, v, D.edge_lengths(u,v));
        fail("edge_phantom", ctx, buf);
        return false;
      }
  return true;
}

bool IDTAudit::verify_triangle_inequality(const FulleroidDelaunay& D, const char* ctx) {
  n_checks++;
  for (node_t u = 0; u < D.N; u++)
    for (size_t j = 0; j < D[u].size(); j++) {
      node_t v = D[u][j], w = D[u][(j+1) % D[u].size()];
      if (u < v && u < w) {
        double a = D.get_length(u,v), b = D.get_length(v,w), c = D.get_length(w,u);
        if (a + b <= c || b + c <= a || c + a <= b) {
          char buf[128];
          snprintf(buf, sizeof(buf), "triangle(%d,%d,%d) sides %g,%g,%g", u, v, w, a, b, c);
          fail("triangle_inequality", ctx, buf);
          return false;
        }
      }
    }
  return true;
}

bool IDTAudit::verify_positive_area(const FulleroidDelaunay& D, const char* ctx) {
  n_checks++;
  for (node_t u = 0; u < D.N; u++)
    for (size_t j = 0; j < D[u].size(); j++) {
      node_t v = D[u][j], w = D[u][(j+1) % D[u].size()];
      if (u < v && u < w) {
        double a = D.get_length(u,v), b = D.get_length(v,w), c = D.get_length(w,u);
        double s = (a+b+c)/2;
        double H = s*(s-a)*(s-b)*(s-c);
        if (H <= 0) {
          char buf[128];
          snprintf(buf, sizeof(buf), "triangle(%d,%d,%d) heron=%g (sides %g,%g,%g)",
                   u, v, w, H, a, b, c);
          fail("positive_area", ctx, buf);
          return false;
        }
      }
    }
  return true;
}

bool IDTAudit::verify_cone_angles(const FulleroidDelaunay& D, const char* ctx) {
  n_checks++;
  for (node_t u = 0; u < D.N; u++) {
    double total = 0;
    auto nbrs = D[u];
    int k = nbrs.size();
    for (int j = 0; j < k; j++) {
      node_t v = nbrs[j], w = nbrs[(j+1)%k];
      double a = D.get_length(v,w), b = D.get_length(u,v), c = D.get_length(u,w);
      double cos_u = (b*b + c*c - a*a) / (2.0*b*c);
      cos_u = std::max(-1.0, std::min(1.0, cos_u));
      total += acos(cos_u);
    }
    double expected = original_degrees[u] * M_PI / 3.0;
    double err = fabs(total - expected);
    if (err > 1e-8) {
      char buf[128];
      snprintf(buf, sizeof(buf), "vertex %d: cone_angle=%g expected=%g (orig_deg=%d) err=%g",
               u, total, expected, original_degrees[u], err);
      fail("cone_angle", ctx, buf);
      return false;
    }
  }
  return true;
}

bool IDTAudit::verify_total_area(const FulleroidDelaunay& D, const char* ctx) {
  n_checks++;
  double area = 0;
  for (node_t u = 0; u < D.N; u++)
    for (size_t j = 0; j < D[u].size(); j++) {
      node_t v = D[u][j], w = D[u][(j+1) % D[u].size()];
      if (u < v && u < w) {
        double a = D.get_length(u,v), b = D.get_length(v,w), c = D.get_length(w,u);
        double s = (a+b+c)/2;
        double H = s*(s-a)*(s-b)*(s-c);
        area += (H > 0) ? sqrt(H) : 0;
      }
    }
  double err = fabs(area - expected_area);
  if (err > 1e-8) {
    char buf[128];
    snprintf(buf, sizeof(buf), "total_area=%g expected=%g err=%g", area, expected_area, err);
    fail("total_area", ctx, buf);
    return false;
  }
  return true;
}

// Check if n is a Loeschian number (a^2 + ab + b^2 for non-negative integers a,b).
static bool is_loeschian(int n) {
  if (n < 0) return false;
  if (n == 0) return true;
  for (int a = 0; 3*a*a <= 4*n; a++) {
    int disc = 4*n - 3*a*a;
    int s = (int)round(sqrt((double)disc));
    if (s*s != disc) continue;
    if ((s - a) % 2 != 0) continue;
    int b = (s - a) / 2;
    if (b >= 0) return true;
  }
  return false;
}

bool IDTAudit::verify_loeschian(const FulleroidDelaunay& D, const char* ctx) {
  n_checks++;
  for (node_t u = 0; u < D.N; u++)
    for (node_t v : D[u])
      if (u < v) {
        double L = D.get_length(u,v);
        double L2 = L*L;
        int L2i = (int)round(L2);
        if (fabs(L2 - L2i) > 1e-6) {
          char buf[128];
          snprintf(buf, sizeof(buf), "edge(%d,%d) L=%g L^2=%g not integer (nearest=%d, err=%g)",
                   u, v, L, L2, L2i, L2-L2i);
          fail("loeschian_integer", ctx, buf);
          return false;
        }
        if (!is_loeschian(L2i)) {
          char buf[128];
          snprintf(buf, sizeof(buf), "edge(%d,%d) L^2=%d is not a Loeschian number", u, v, L2i);
          fail("loeschian_form", ctx, buf);
          return false;
        }
      }
  return true;
}

bool IDTAudit::verify_no_multi_edges(const FulleroidDelaunay& D, const char* ctx) {
  n_checks++;
  for (node_t u = 0; u < D.N; u++) {
    auto row = D[u];
    for (size_t i = 0; i < row.size(); i++)
      for (size_t j = i+1; j < row.size(); j++)
        if (row[i] == row[j]) {
          char buf[128];
          snprintf(buf, sizeof(buf), "vertex %d has duplicate neighbour %d (positions %zu,%zu)",
                   u, row[i], i, j);
          fail("no_multi_edges", ctx, buf);
          return false;
        }
  }
  return true;
}

// ============================================================================
// 3D Embedding
// ============================================================================

matrix<double> FulleroidDelaunay::all_pairs_distances() const
{
  const double INF = 1e30;
  matrix<double> D(N, N, INF);

  for (node_t u = 0; u < N; u++) {
    D(u, u) = 0;
    for (node_t v : (*this)[u])
      D(u, v) = get_length(u, v);
  }

  // Floyd-Warshall
  for (node_t k = 0; k < N; k++)
    for (node_t i = 0; i < N; i++)
      for (node_t j = 0; j < N; j++)
        if (D(i, k) + D(k, j) < D(i, j))
          D(i, j) = D(i, k) + D(k, j);

  return D;
}

// ============================================================================
// 3D Embedding primitives
// ============================================================================

// Vector space over R^{3N}: element-wise operations on arrays of coord3d.
// Inherits from vector<coord3d> so it passes transparently to existing APIs.
struct V3 : vector<coord3d> {
  using vector::vector;

  void zero() { assign(size(), coord3d()); }
  double dot(const V3& b) const { double s = 0; for (size_t i = 0; i < size(); i++) s += (*this)[i].dot(b[i]); return s; }
  double norm() const { return sqrt(dot(*this)); }

  V3& operator+=(const V3& b) { for (size_t i = 0; i < size(); i++) (*this)[i] += b[i]; return *this; }
  V3 operator+(const V3& b) const { V3 r(*this); r += b; return r; }
  V3 operator-(const V3& b) const { V3 r(size()); for (size_t i = 0; i < size(); i++) r[i] = (*this)[i] - b[i]; return r; }
  V3 operator-() const { V3 r(size()); for (size_t i = 0; i < size(); i++) r[i] = -(*this)[i]; return r; }
  V3 operator*(double s) const { V3 r(*this); for (auto& v : r) v *= s; return r; }
  friend V3 operator*(double s, const V3& v) { return v * s; }
};

// Per-edge distance-matching geometry for E = sum (|xi-xj| - Lij)^2.
//   Hessian per edge: H = 2[(L/r) uu^T + (e/r)(I - uu^T)]
//   where u = diff/dist is the unit edge vector, e = dist - target.
struct EdgeStressData {
  int i, j;
  coord3d diff;                // x[i] - x[j]
  double dist, target, err;   // |diff|, L_ij, dist - L_ij

  static EdgeStressData compute(const vector<coord3d>& x, int i, int j, double L) {
    coord3d d = x[i] - x[j];
    double r = d.norm();
    return {i, j, d, r, L, r - L};
  }

  bool valid() const { return dist > 1e-15; }
  double energy() const { return err * err; }

  void scatter_gradient(vector<coord3d>& g) const {
    coord3d c = diff * (2.0 * err / dist);
    g[i] = g[i] + c;  g[j] = g[j] - c;
  }

  void scatter_hv(const vector<coord3d>& v, vector<coord3d>& Hv) const {
    coord3d dv = v[i] - v[j];
    double udv = diff.dot(dv) / (dist * dist);
    coord3d h = diff * (2.0 * target / dist * udv) + dv * (2.0 * err / dist);
    Hv[i] = Hv[i] + h;  Hv[j] = Hv[j] - h;
  }
};

// Apply f to each EdgeStressData in the triangulation.
template<typename F>
static void for_each_edge(int N, const neighbours_t& nbrs, const matrix<double>& el,
                          const vector<coord3d>& x, F&& f) {
  for (node_t u = 0; u < N; u++)
    for (node_t v : nbrs[u])
      if (u < v) {
        auto ed = EdgeStressData::compute(x, u, v, el(u, v));
        if (ed.valid()) f(ed);
      }
}

// Double-center a symmetric matrix: B = -0.5 * J * A * J
// where J = I - (1/N)*11^T is the centering matrix.
static matrix<double> double_center(const matrix<double>& A) {
  int N = A.m;
  vector<double> row_mean(N, 0), col_mean(N, 0);
  double grand_mean = 0;
  for (int i = 0; i < N; i++)
    for (int j = 0; j < N; j++) {
      row_mean[i] += A(i, j);
      col_mean[j] += A(i, j);
      grand_mean += A(i, j);
    }
  for (int i = 0; i < N; i++) { row_mean[i] /= N; col_mean[i] /= N; }
  grand_mean /= (N * N);

  matrix<double> B(N, N);
  for (int i = 0; i < N; i++)
    for (int j = 0; j < N; j++)
      B(i, j) = -0.5 * (A(i, j) - row_mean[i] - col_mean[j] + grand_mean);
  return B;
}

// Jacobi eigendecomposition of symmetric matrix A.
// Returns {eigenvalues, eigenvectors} where columns of V are eigenvectors.
static pair<vector<double>, matrix<double>> jacobi_eigen(matrix<double> A) {
  int N = A.m;
  matrix<double> V(N, N, 0.0);
  for (int i = 0; i < N; i++) V(i, i) = 1.0;

  for (int iter = 0; iter < 10000; iter++) {
    // Find largest off-diagonal element
    double max_val = 0; int p = 0, q = 1;
    for (int i = 0; i < N; i++)
      for (int j = i + 1; j < N; j++)
        if (fabs(A(i, j)) > max_val) { max_val = fabs(A(i, j)); p = i; q = j; }
    if (max_val < 1e-15) break;

    double app = A(p,p), aqq = A(q,q), apq = A(p,q);
    double theta = (fabs(app - aqq) < 1e-30) ? M_PI/4 : 0.5 * atan2(2*apq, app - aqq);
    double cs = cos(theta), sn = sin(theta);

    // Givens rotation on rows, columns, and eigenvector accumulator
    for (int j = 0; j < N; j++) {
      double bp = A(p,j), bq = A(q,j);
      A(p,j) = cs*bp + sn*bq;  A(q,j) = -sn*bp + cs*bq;
    }
    for (int i = 0; i < N; i++) {
      double bp = A(i,p), bq = A(i,q);
      A(i,p) = cs*bp + sn*bq;  A(i,q) = -sn*bp + cs*bq;
    }
    for (int i = 0; i < N; i++) {
      double vp = V(i,p), vq = V(i,q);
      V(i,p) = cs*vp + sn*vq;  V(i,q) = -sn*vp + cs*vq;
    }
  }

  vector<double> evals(N);
  for (int i = 0; i < N; i++) evals[i] = A(i, i);
  return {evals, V};
}

// Solve ||p + tau*d|| = Delta for the positive root tau.
static double trust_boundary_tau(const V3& p, const V3& d, double Delta) {
  double pd = p.dot(d), pp = p.dot(p), dd = d.dot(d);
  return (-pd + sqrt(max(0.0, pd*pd - dd*(pp - Delta*Delta)))) / dd;
}

// 3x3 determinant and adjugate for solving inertia-tensor systems.
static double det3(const matrix3d& M) {
  return M(0,0)*(M(1,1)*M(2,2)-M(1,2)*M(2,1))
       - M(0,1)*(M(1,0)*M(2,2)-M(1,2)*M(2,0))
       + M(0,2)*(M(1,0)*M(2,1)-M(1,1)*M(2,0));
}
static matrix3d adjugate(const matrix3d& M) {
  return matrix3d(
    M(1,1)*M(2,2)-M(1,2)*M(2,1), M(0,2)*M(2,1)-M(0,1)*M(2,2), M(0,1)*M(1,2)-M(0,2)*M(1,1),
    M(1,2)*M(2,0)-M(1,0)*M(2,2), M(0,0)*M(2,2)-M(0,2)*M(2,0), M(0,2)*M(1,0)-M(0,0)*M(1,2),
    M(1,0)*M(2,1)-M(1,1)*M(2,0), M(0,1)*M(2,0)-M(0,0)*M(2,1), M(0,0)*M(1,1)-M(0,1)*M(1,0));
}

// Project out 6 rigid-body DOFs (3 translation + 3 rotation) from gradient vector.
static void project_rigid_body(vector<coord3d>& g, const vector<coord3d>& x) {
  int N = g.size();

  // Translation: subtract mean
  coord3d mean;
  for (auto& gi : g) mean += gi;
  mean /= N;
  for (auto& gi : g) gi -= mean;

  // Rotation: omega = inertia^{-1} * torque,  then  g -= omega x r
  coord3d c;
  for (auto& xi : x) c += xi;
  c /= N;

  matrix3d I;
  coord3d tau;
  for (int k = 0; k < N; k++) {
    coord3d r = x[k] - c;
    I += matrix3d::unit_matrix() * r.dot(r) - r.outer(r);
    tau += r.cross(g[k]);
  }

  double d = det3(I);
  if (fabs(d) < 1e-30) return;
  coord3d omega = adjugate(I) * tau * (1.0 / d);

  for (int k = 0; k < N; k++)
    g[k] -= omega.cross(x[k] - c);
}

// Orient polyhedron so CCW face normals point outward.
static void orient_outward(vector<coord3d>& x, const neighbours_t& nbrs) {
  int N = x.size();
  coord3d c;
  for (auto& xi : x) c += xi;
  c /= N;

  double vol = 0;
  for (node_t u = 0; u < N; u++) {
    int deg = nbrs[u].size();
    for (int j = 0; j < deg; j++) {
      node_t v = nbrs[u][j], w = nbrs[u][(j+1)%deg];
      if (u < v && u < w)
        vol += (x[u]-c).dot((x[v]-c).cross(x[w]-c));
    }
  }
  if (vol < 0)
    for (auto& xi : x) xi = c*2.0 - xi;
}

// Trust-region Steihaug CG optimizer.
// Minimizes f(x) given callbacks for energy+gradient and Hessian-vector product.
template<typename Eval, typename HvProd>
static V3 steihaug_cg(V3 x, Eval eval, HvProd hv_prod,
                       int max_outer = 200, double tol = 1e-13) {
  int N = x.size();
  V3 g(N);
  double E = eval(x, g);
  double Delta = max(g.norm(), 1e-14);

  for (int outer = 0; outer < max_outer; outer++) {
    double gnorm = g.norm();
    if (gnorm < tol) break;

    // Steihaug CG subproblem: min g.p + 0.5 p.Hp  s.t. |p| <= Delta
    V3 p(N), r(g), d = -r;
    double rr = r.dot(r);
    bool hit_boundary = false;

    for (int cg = 0; cg < 3*N; cg++) {
      V3 Hd(N);
      hv_prod(x, d, Hd);
      double dHd = d.dot(Hd);

      if (dHd <= 1e-15 * rr) {
        p += d * trust_boundary_tau(p, d, Delta);
        hit_boundary = true; break;
      }

      double alpha = rr / dHd;
      V3 p_new = p + d * alpha;

      if (p_new.norm() >= Delta) {
        p += d * trust_boundary_tau(p, d, Delta);
        hit_boundary = true; break;
      }

      p = p_new;
      r += Hd * alpha;
      double rr_new = r.dot(r);
      if (sqrt(rr_new) < tol * gnorm) break;

      d = d * (rr_new / rr) - r;
      rr = rr_new;
    }

    // Evaluate trial point and update trust region
    double predicted = g.dot(p);
    { V3 Hp(N); hv_prod(x, p, Hp); predicted += 0.5 * p.dot(Hp); }

    V3 x_trial = x + p, g_trial(N);
    double E_trial = eval(x_trial, g_trial);

    double rho = (E - E_trial) / max(1e-30, -predicted);
    double pnorm = p.norm();
    if (rho < 0.25)                        Delta = 0.25 * pnorm;
    else if (rho > 0.75 && hit_boundary)   Delta = min(2.0 * Delta, 1e10);
    Delta = max(Delta, 1e-15);

    if (rho > 0.1) { x = x_trial; g = g_trial; E = E_trial; }
  }
  return x;
}

// Classical MDS: APSP distances -> double-center -> Jacobi eigen -> top-3 embedding.
// See initgeometry-BFSvMDS.tex for full description.
static V3 mds_placement(const matrix<double>& D) {
  int N = D.m;

  // Squared-distance matrix
  matrix<double> D_sq(N, N);
  for (int i = 0; i < N; i++)
    for (int j = 0; j < N; j++)
      D_sq(i, j) = D(i, j) * D(i, j);

  auto [evals, V] = jacobi_eigen(double_center(D_sq));

  // Sort eigenvalues descending, embed using top 3
  vector<int> order(N);
  iota(order.begin(), order.end(), 0);
  sort(order.begin(), order.end(), [&](int a, int b) { return evals[a] > evals[b]; });

  V3 x(N);
  for (int i = 0; i < N; i++)
    for (int dim = 0; dim < 3; dim++) {
      int col = order[dim];
      x[i][dim] = V(i, col) * sqrt(max(0.0, evals[col]));
    }

  // Flatness detection: if 3rd eigenvalue is degenerate, perturb
  double lam3 = max(0.0, evals[order[2]]);
  double lam1 = max(1e-30, evals[order[0]]);
  if (lam3 / lam1 < 1e-6)
    for (int i = 0; i < N; i++)
      x[i][2] += 0.1 * lam1 * sin(i * 2.0 * M_PI / N);

  return x;
}

vector<coord3d> FulleroidDelaunay::embed_3d() const
{
  // Step 1: Classical MDS initial placement (APSP + eigendecomposition)
  V3 x = mds_placement(all_pairs_distances());

  // Step 2: Stress minimization — E = sum_{edges} (|xi-xj| - Lij)^2
  auto eval = [&](const V3& pos, V3& g) -> double {
    double E = 0;
    g.zero();
    for_each_edge(N, *this, edge_lengths, pos, [&](const EdgeStressData& ed) {
      E += ed.energy();
      ed.scatter_gradient(g);
    });
    project_rigid_body(g, pos);
    return E;
  };

  auto hv_prod = [&](const V3& pos, const V3& dir, V3& Hv) {
    Hv.zero();
    for_each_edge(N, *this, edge_lengths, pos, [&](const EdgeStressData& ed) {
      ed.scatter_hv(dir, Hv);
    });
    project_rigid_body(Hv, pos);
  };

  x = steihaug_cg(std::move(x), eval, hv_prod);

  // Step 3: Orient so CCW face normals point outward
  orient_outward(x, *this);
  return x;
}

// ============================================================================
// OriginTracker — exact Eisenstein face-origin tracking
// ============================================================================
//
// Maintains a copy of the original equilateral triangulation so that during
// flips and vertex removals, the repartitioning of original faces among iDT
// faces can be done exactly using the Eisenstein integer turn predicate.
//
// The key operation is unfold_patch(): BFS-unfold a connected set of original
// equilateral faces into the Z[omega] grid.  Since all original edges have
// unit length, every vertex lands at an Eisenstein integer, and the turn()
// predicate gives exact orientation tests with zero floating-point error.

struct DelaunayTriangulation::OriginTracker {
  int N;                                    // vertex count in original triangulation
  vector<array<int,3>> face_verts;          // face_verts[fid] = {u, v, w} CCW
  unordered_map<int64_t, int> arc_to_face;  // directed arc (u,v) → face ID

  // Build from the sorted triangulation used by from_triangulation().
  // num_faces must match the number of faces assigned during construction.
  OriginTracker(const Triangulation& T, int num_faces);

  // BFS-unfold a connected patch of original equilateral faces into Z[omega].
  // Returns a map from original vertex ID to Eisenstein grid position.
  unordered_map<int, Eisenstein> unfold_patch(const vector<int>& face_ids) const;

  // Classify original faces across the directed line vtx_A → vtx_B.
  // Returns (left_of_line, right_of_line).  Faces whose centroid is exactly
  // on the line are assigned to both sides.
  pair<vector<int>, vector<int>>
  classify_across_line(const vector<int>& face_ids,
                       int vtx_A, int vtx_B) const;

  // Classify original faces into a set of triangles (for ear-clipping).
  // tri_verts[j] = {a, b, c} gives the CCW vertex IDs of triangle j.
  // Returns assignment[j] = list of original face IDs inside triangle j.
  vector<vector<int>>
  classify_into_triangles(const vector<int>& face_ids,
                          const vector<array<int,3>>& tri_verts) const;
};

DelaunayTriangulation::OriginTracker::OriginTracker(
    const Triangulation& T, int num_faces) : N(T.N)
{
  face_verts.resize(num_faces);
  int fid = 0;
  for (node_t u = 0; u < T.N; u++) {
    auto row = T[u];
    int deg = row.size();
    for (int j = 0; j < deg; j++) {
      node_t v = row[j], w = row[(j+1) % deg];
      if (u < v && u < w) {
        face_verts[fid] = {u, v, w};
        arc_to_face[int64_t(u) * N + v] = fid;
        arc_to_face[int64_t(v) * N + w] = fid;
        arc_to_face[int64_t(w) * N + u] = fid;
        fid++;
      }
    }
  }
  assert(fid == num_faces);
}

unordered_map<int, Eisenstein>
DelaunayTriangulation::OriginTracker::unfold_patch(
    const vector<int>& face_ids) const
{
  if (face_ids.empty()) return {};

  unordered_set<int> patch(face_ids.begin(), face_ids.end());
  unordered_map<int, Eisenstein> pos;
  unordered_set<int> placed;
  queue<int> Q;

  // Place first face at the standard triangle.
  int f0 = face_ids[0];
  auto& fv = face_verts[f0];
  pos[fv[0]] = Eisenstein(0, 0);
  pos[fv[1]] = Eisenstein(1, 0);
  pos[fv[2]] = Eisenstein(0, 1);
  placed.insert(f0);
  Q.push(f0);

  while (!Q.empty()) {
    int f = Q.front(); Q.pop();
    auto& v = face_verts[f];

    for (int i = 0; i < 3; i++) {
      int eu = v[i], ev = v[(i+1) % 3];

      // Adjacent face across edge eu→ev has the twin arc ev→eu.
      auto it = arc_to_face.find(int64_t(ev) * N + eu);
      if (it == arc_to_face.end()) continue;
      int adj = it->second;
      if (!patch.count(adj) || placed.count(adj)) continue;

      // Third vertex of adjacent face (ev, eu, third) in CCW order.
      auto& av = face_verts[adj];
      int third = av[0];
      if (third == eu || third == ev) third = av[1];
      if (third == eu || third == ev) third = av[2];

      // Place: third = ev + (eu - ev).nextCCW()  [standard equilateral unfolding]
      if (!pos.count(third))
        pos[third] = pos[ev] + (pos[eu] - pos[ev]).nextCCW();

      placed.insert(adj);
      Q.push(adj);
    }
  }
  return pos;
}

pair<vector<int>, vector<int>>
DelaunayTriangulation::OriginTracker::classify_across_line(
    const vector<int>& face_ids, int vtx_A, int vtx_B) const
{
  auto pos = unfold_patch(face_ids);

  // Scale line endpoints by 3 so we can test the centroid (p+q+r)
  // without dividing by 3.  sign(turn(3A, 3B, p+q+r)) = sign(turn(A, B, centroid)).
  Eisenstein A3 = pos.at(vtx_A) * 3;
  Eisenstein B3 = pos.at(vtx_B) * 3;

  vector<int> left, right;
  for (int fid : face_ids) {
    auto& fv = face_verts[fid];
    Eisenstein sum = pos.at(fv[0]) + pos.at(fv[1]) + pos.at(fv[2]);
    int t = Eisenstein::turn(A3, B3, sum);
    if (t > 0) left.push_back(fid);
    else if (t < 0) right.push_back(fid);
    else { left.push_back(fid); right.push_back(fid); }  // on the line: both
  }
  return {left, right};
}

vector<vector<int>>
DelaunayTriangulation::OriginTracker::classify_into_triangles(
    const vector<int>& face_ids,
    const vector<array<int,3>>& tri_verts) const
{
  auto pos = unfold_patch(face_ids);
  int nt = tri_verts.size();
  vector<vector<int>> assignment(nt);

  // Precompute scaled triangle corners.
  vector<array<Eisenstein,3>> corners(nt);
  for (int j = 0; j < nt; j++)
    for (int k = 0; k < 3; k++)
      corners[j][k] = pos.at(tri_verts[j][k]) * 3;

  for (int fid : face_ids) {
    auto& fv = face_verts[fid];
    Eisenstein sum = pos.at(fv[0]) + pos.at(fv[1]) + pos.at(fv[2]);

    for (int j = 0; j < nt; j++) {
      auto& c = corners[j];
      // Point inside CCW triangle iff all three turn tests are >= 0.
      if (Eisenstein::turn(c[0], c[1], sum) >= 0 &&
          Eisenstein::turn(c[1], c[2], sum) >= 0 &&
          Eisenstein::turn(c[2], c[0], sum) >= 0) {
        assignment[j].push_back(fid);
        break;  // each face goes to exactly one triangle
      }
    }
  }
  return assignment;
}

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
    f_origin.push_back({});
  }
  f_he[fid] = -1;
  f_origin[fid].clear();
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
  f_origin[f].clear();
  free_faces.push_back(f);
}

int DelaunayTriangulation::vertex_degree(int v) const
{
  if (v_out[v] < 0) return 0;
  int h0 = v_out[v], h = h0, deg = 0;
  do { deg++; h = cw(h); } while (h != h0);
  return deg;
}

// --- Construction ---

DelaunayTriangulation DelaunayTriangulation::from_triangulation(const Triangulation& T)
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
  D.he_length.assign(D.nh, 1.0);
  D.he_angle.assign(D.nh, M_PI / 3.0);

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
        D.f_origin.push_back({fid});
      }
    }
    if (deg > 0) D.v_out[u] = arc_to_he.at({u, row[0]});
  }

  // Phase 4: Vertex data.
  D.v_cone_angle.resize(D.nv);
  D.v_orig_degree.resize(D.nv);
  for (node_t v = 0; v < T.N; v++) {
    D.v_orig_degree[v] = T.degree(v);
    D.v_cone_angle[v] = T.degree(v) * M_PI / 3.0;
  }

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
  int h0 = f_he[f];
  int h1 = he_next[h0];
  int h2 = he_next[h1];
  double a = he_length[h0], b = he_length[h1], c = he_length[h2];
  // Angle at origin of h0 = angle opposite side h1 (length b)... no.
  // h0: u->v (length a), h1: v->w (length b), h2: w->u (length c).
  // Angle at u (origin of h0) is opposite side b (the side v->w).
  // Wait: the angle at u is between edges u->v (length a) and u->w (length c).
  // The opposite side is v->w (length b).
  // cos(angle_u) = (a^2 + c^2 - b^2) / (2*a*c)
  he_angle[h0] = acos(std::clamp((a*a + c*c - b*b) / (2*a*c), -1.0, 1.0));
  he_angle[h1] = acos(std::clamp((a*a + b*b - c*c) / (2*a*b), -1.0, 1.0));
  he_angle[h2] = acos(std::clamp((b*b + c*c - a*a) / (2*b*c), -1.0, 1.0));
}

void DelaunayTriangulation::recompute_all_angles()
{
  for (int f = 0; f < nf; f++)
    recompute_face_angles(f);
}

// --- Delaunay operations ---

bool DelaunayTriangulation::is_delaunay_edge(int h) const
{
  return diamond(h).is_delaunay();
}

bool DelaunayTriangulation::flip_edge(int h)
{
  int t = h ^ 1;
  int h1 = he_next[h], h2 = he_next[h1];   // face of h
  int h4 = he_next[t], h5 = he_next[h4];    // face of twin

  int u = he_origin[h], v = he_origin[t];
  int B = he_origin[h2], D = he_origin[h5];

  // Topological guard: B == D means the diamond is degenerate.
  if (B == D) return false;

  // Geometric guards.
  Diamond dm = diamond(h);
  if (!dm.is_convex()) return false;
  double f_len = dm.flipped_length();
  if (!std::isfinite(f_len) || f_len <= 0) return false;

  // Execute flip.
  // Before: face(h) = h -> h1 -> h2,  face(t) = t -> h4 -> h5
  // After:  face(h) = h -> h5 -> h1,  face(t) = t -> h2 -> h4

  int fh = he_face[h], ft = he_face[t];

  // Origins: h becomes B->D, t becomes D->B.
  he_origin[h] = B;
  he_origin[t] = D;

  // Next pointers (all 6 change).
  he_next[h]  = h5;  he_next[h5] = h1;  he_next[h1] = h;
  he_next[t]  = h2;  he_next[h2] = h4;  he_next[h4] = t;

  // Face assignments: h5 moves to face(h), h2 moves to face(t).
  he_face[h5] = fh;
  he_face[h2] = ft;
  // h, h1 stay in fh; t, h4 stay in ft.

  // Update face representatives.
  f_he[fh] = h;
  f_he[ft] = t;

  // Edge length.
  he_length[h] = f_len;
  he_length[t] = f_len;

  // Update v_out: u lost h (no longer leaves u), v lost t.
  if (v_out[u] == h) v_out[u] = h4;
  if (v_out[v] == t) v_out[v] = h1;

  // Recompute angles for both affected faces.
  recompute_face_angles(fh);
  recompute_face_angles(ft);

  // Repartition face origins across the new diagonal B→D.
  {
    vector<int> all;
    all.reserve(f_origin[fh].size() + f_origin[ft].size());
    all.insert(all.end(), f_origin[fh].begin(), f_origin[fh].end());
    all.insert(all.end(), f_origin[ft].begin(), f_origin[ft].end());
    sort(all.begin(), all.end());
    all.erase(unique(all.begin(), all.end()), all.end());

    if (origin_tracker) {
      // Exact: classify each original face by which side of B→D its centroid
      // falls on, using Eisenstein turn() in the Z[omega] grid.
      // After flip, face(h) = (B, D, v) is left of B→D;
      //             face(t) = (D, B, u) is right of B→D.
      auto [left, right] = origin_tracker->classify_across_line(all, B, D);
      f_origin[fh] = std::move(left);
      f_origin[ft] = std::move(right);
    } else {
      // Conservative: assign the full union to both faces.
      f_origin[fh] = all;
      f_origin[ft] = all;
    }
  }

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

  int budget = 200 * nv;
  while (!S.empty() && budget > 0) {
    int h = S.top(); S.pop();
    in_stack[h >> 1] = false;

    if (!alive(h) || is_delaunay_edge(h)) continue;

    // Record rim edges before flipping (they'll be checked next).
    int h1 = he_next[h], h2 = he_next[h1];
    int h4 = he_next[h ^ 1], h5 = he_next[h4];

    if (!flip_edge(h)) continue;
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
  int flips = lawson_sweep();
  if (!is_delaunay())
    flips += delaunay_resolve();
  return flips;
}

int DelaunayTriangulation::delaunay_resolve()
{
  int total_flips = 0;

  for (int round = 0; round < 20; round++) {
    int count0 = count_non_delaunay();
    if (count0 == 0) break;

    DelaunayTriangulation saved = *this;
    bool improved = false;

    // Collect all live edges.
    vector<int> edges;
    for (int h = 0; h < nh; h += 2)
      if (alive(h)) edges.push_back(h);

    // Depth 1: try flipping each edge, then Lawson sweep.
    for (size_t i = 0; i < edges.size() && !improved; i++) {
      *this = saved;
      if (!flip_edge(edges[i])) continue;
      int f = 1 + lawson_sweep();
      if (count_non_delaunay() < count0) {
        total_flips += f;
        improved = true;
      }
    }

    // Depth 2: try all pairs.
    for (size_t i = 0; i < edges.size() && !improved; i++) {
      for (size_t j = i + 1; j < edges.size() && !improved; j++) {
        *this = saved;
        if (!flip_edge(edges[i])) continue;
        if (!alive(edges[j]) || !flip_edge(edges[j])) continue;
        int f = 2 + lawson_sweep();
        if (count_non_delaunay() < count0) {
          total_flips += f;
          improved = true;
        }
      }
    }

    if (!improved) { *this = saved; break; }
  }
  return total_flips;
}

// --- Vertex removal ---

void DelaunayTriangulation::remove_flat_vertex(int v)
{
  int deg = vertex_degree(v);
  if (deg == 0) return;

  // Phase 1: Reduce degree by flipping incident edges.
  while (deg > 4) {
    bool progress = false;
    int h0 = v_out[v], h = h0;
    do {
      if (flip_edge(h)) { deg--; progress = true; break; }
      h = cw(h);
    } while (h != h0);
    if (!progress) break;
  }

  deg = vertex_degree(v);

  // Phase 2: Ear clipping for degree >= 4.
  if (deg >= 4) {
    int k = deg;

    // Collect spoke half-edges in CCW order around v.
    // CCW order: face of spoke_he[i] is (v, nb[i], nb[i+1]).
    //   spoke_he[i]: v -> nb[i]
    //   next(spoke_he[i]): nb[i] -> nb[i+1]  (inner rim)
    //   prev(spoke_he[i]): nb[i+1] -> v
    //   ccw(spoke_he[i]) = spoke_he[i+1]: v -> nb[i+1]
    vector<int> spoke_he(k);
    spoke_he[0] = v_out[v];
    for (int i = 1; i < k; i++)
      spoke_he[i] = ccw(spoke_he[i-1]);

    vector<int> nb(k);
    vector<double> spokes(k), rims(k);
    for (int i = 0; i < k; i++) {
      nb[i] = dest(spoke_he[i]);
      spokes[i] = he_length[spoke_he[i]];
      // Inner rim: next(spoke_he[i]) goes nb[i] -> nb[(i+1)%k].
      rims[i] = he_length[he_next[spoke_he[i]]];
    }

    // Cumulative fan angles.
    vector<double> cum(k + 1, 0);
    for (int i = 0; i < k; i++) {
      double si = spokes[i], sn = spokes[(i+1)%k], r = rims[i];
      double cos_a = std::clamp((si*si + sn*sn - r*r) / (2*si*sn), -1.0, 1.0);
      cum[i+1] = cum[i] + acos(cos_a);
    }

    auto diag_len = [&](int from, int to) -> double {
      double angle = (to > from) ? cum[to] - cum[from]
                                  : (cum[k] - cum[from]) + cum[to];
      double sf = spokes[from], st = spokes[to];
      double len2 = sf*sf + st*st - 2*sf*st*cos(angle);
      return (len2 > 0) ? sqrt(len2) : 0;
    };

    auto ear_signed_area = [&](int pp, int pi, int pn) -> double {
      double rp = spokes[pp], ri = spokes[pi], rn = spokes[pn];
      double tp = cum[pp], ti = cum[pi], tn = cum[pn];
      return rp*ri*sin(ti - tp) + ri*rn*sin(tn - ti) + rn*rp*sin(tp - tn);
    };

    // Ear clipping (same algorithm as FulleroidDelaunay).
    struct EarDiag { int prev_idx, ear_idx, next_idx; double len; };
    vector<EarDiag> diagonals;
    vector<int> poly(k);
    for (int i = 0; i < k; i++) poly[i] = i;

    while ((int)poly.size() > 3) {
      int n = poly.size();
      bool found = false;
      for (int j = 0; j < n; j++) {
        int jm = (j - 1 + n) % n, jp = (j + 1) % n;
        int pp = poly[jm], pi = poly[j], pn = poly[jp];

        if (nb[pp] == nb[pn]) continue;

        // Multi-edge check via vertex circulation.
        bool edge_exists = false;
        { int h0 = v_out[nb[pp]], hc = h0;
          if (h0 >= 0) do {
            if (dest(hc) == nb[pn]) { edge_exists = true; break; }
            hc = cw(hc);
          } while (hc != h0); }
        for (auto& d : diagonals)
          if ((nb[d.prev_idx] == nb[pp] && nb[d.next_idx] == nb[pn]) ||
              (nb[d.prev_idx] == nb[pn] && nb[d.next_idx] == nb[pp]))
            { edge_exists = true; break; }
        if (edge_exists) continue;

        if (ear_signed_area(pp, pi, pn) <= 1e-10) continue;

        double sub_angle = (pn > pp) ? cum[pn] - cum[pp]
                                      : (cum[k] - cum[pp]) + cum[pn];
        if (sub_angle > M_PI + 1e-10) continue;

        double len = diag_len(pp, pn);
        if (len <= 1e-15) continue;

        diagonals.push_back({pp, pi, pn, len});
        poly.erase(poly.begin() + j);
        found = true;
        break;
      }

      if (!found) {
        // Blocker-flip: try flipping blocking edges.
        bool flipped_blocker = false;
        int n2 = poly.size();
        for (int j = 0; j < n2 && !flipped_blocker; j++) {
          int jm = (j - 1 + n2) % n2, jp = (j + 1) % n2;
          int pp = poly[jm], pn = poly[jp];
          if (nb[pp] == nb[pn]) continue;
          int h0 = v_out[nb[pp]], hc = h0;
          if (h0 >= 0) do {
            if (dest(hc) == nb[pn]) {
              if (flip_edge(hc)) { flipped_blocker = true; break; }
              if (flip_edge(hc ^ 1)) { flipped_blocker = true; break; }
            }
            hc = cw(hc);
          } while (hc != h0);
        }
        if (!flipped_blocker) return;  // stuck
      }
    }

    // --- Splice the ear-clipped polygon into the DCEL ---

    // Collect face origins from the fan.
    vector<int> all_origins;
    { int h0 = v_out[v], h = h0;
      do {
        int f = he_face[h];
        all_origins.insert(all_origins.end(), f_origin[f].begin(), f_origin[f].end());
        h = ccw(h);
      } while (h != h0); }
    sort(all_origins.begin(), all_origins.end());
    all_origins.erase(unique(all_origins.begin(), all_origins.end()), all_origins.end());

    // Inner rim half-edges: next(spoke_he[i]) goes nb[i] -> nb[(i+1)%k].
    // These are in the fan faces and will be reassigned to new ear faces.
    // Outer rim half-edges (their twins) are in external faces — untouched.
    vector<int> inner_rim(k);
    for (int i = 0; i < k; i++)
      inner_rim[i] = he_next[spoke_he[i]];

    // Fix v_out for neighbors before deallocating spokes.
    for (int i = 0; i < k; i++) {
      int spoke_twin = spoke_he[i] ^ 1;  // nb[i] -> v
      if (v_out[nb[i]] == spoke_twin) {
        // Use outer rim half-edge: twin of inner_rim[(i-1+k)%k] goes
        // nb[i] -> nb[(i-1+k)%k], which is in an external face.
        // But we need an outgoing he from nb[i] that's NOT in the fan.
        // The outer rim twin of inner_rim[i] goes nb[(i+1)%k] -> nb[i] (incoming).
        // We need outgoing. Outer rim twin of inner_rim[(i-1+k)%k]:
        //   inner_rim[(i-1+k)%k] goes nb[(i-1+k)%k] -> nb[i], in fan.
        //   Its twin goes nb[i] -> nb[(i-1+k)%k], in external face. Outgoing from nb[i]!
        v_out[nb[i]] = inner_rim[(i-1+k)%k] ^ 1;
      }
    }

    // Deallocate fan faces and spoke edges.
    for (int i = 0; i < k; i++) {
      dealloc_face(he_face[spoke_he[i]]);
      dealloc_edge(spoke_he[i]);
    }
    v_out[v] = -1;

    // Build triangles from ear clipping.
    struct TriData { int v0, v1, v2; };
    vector<TriData> tris;
    { vector<int> rpoly(k);
      for (int i = 0; i < k; i++) rpoly[i] = i;
      for (auto& ear : diagonals) {
        tris.push_back({ear.prev_idx, ear.ear_idx, ear.next_idx});
        rpoly.erase(std::find(rpoly.begin(), rpoly.end(), ear.ear_idx));
      }
      assert(rpoly.size() == 3);
      tris.push_back({rpoly[0], rpoly[1], rpoly[2]});
    }

    // Build arc-to-halfedge map. Only inner rim half-edges are available for
    // new faces (outer rim stays with external faces).
    map<pair<int,int>, int> local_arc;
    map<pair<int,int>, double> edge_len_map;

    // Inner rim: goes nb[i] -> nb[(i+1)%k], available for new face wiring.
    for (int i = 0; i < k; i++) {
      local_arc[{i, (i+1)%k}] = inner_rim[i];
      edge_len_map[{i, (i+1)%k}] = edge_len_map[{(i+1)%k, i}] = rims[i];
    }

    // Diagonal edges: allocate new half-edge pairs.
    for (auto& d : diagonals) {
      int h_d = alloc_edge();
      he_origin[h_d] = nb[d.prev_idx];
      he_origin[h_d ^ 1] = nb[d.next_idx];
      he_length[h_d] = d.len;
      he_length[h_d ^ 1] = d.len;
      local_arc[{d.prev_idx, d.next_idx}] = h_d;
      local_arc[{d.next_idx, d.prev_idx}] = h_d ^ 1;
      edge_len_map[{d.prev_idx, d.next_idx}] = edge_len_map[{d.next_idx, d.prev_idx}] = d.len;
    }

    // Wire up each triangle.
    int origin_idx = 0;
    int origins_per_face = all_origins.size() / tris.size();
    for (size_t ti = 0; ti < tris.size(); ti++) {
      auto& tri = tris[ti];
      int fid = alloc_face();

      int h_01 = local_arc.at({tri.v0, tri.v1});
      int h_12 = local_arc.at({tri.v1, tri.v2});
      int h_20 = local_arc.at({tri.v2, tri.v0});

      he_next[h_01] = h_12;
      he_next[h_12] = h_20;
      he_next[h_20] = h_01;
      he_face[h_01] = fid;
      he_face[h_12] = fid;
      he_face[h_20] = fid;
      f_he[fid] = h_01;

      int n_assign = (ti + 1 < tris.size()) ? origins_per_face
                                              : (int)all_origins.size() - origin_idx;
      n_assign = std::max(n_assign, 0);
      f_origin[fid].assign(all_origins.begin() + origin_idx,
                            all_origins.begin() + origin_idx + n_assign);
      origin_idx += n_assign;
    }

    // Recompute angles for new faces.
    for (auto& tri : tris) {
      int h_01 = local_arc.at({tri.v0, tri.v1});
      recompute_face_angles(he_face[h_01]);
    }

    // Fix v_out for neighbors: ensure they point to live outgoing half-edges.
    for (int i = 0; i < k; i++) {
      int u = nb[i];
      if (v_out[u] < 0 || !alive(v_out[u]) || he_origin[v_out[u]] != u) {
        for (auto& [arc, hid] : local_arc)
          if (nb[arc.first] == u && alive(hid)) { v_out[u] = hid; break; }
      }
    }

    return;
  }

  // Phase 3: Degree 3 — the three fan faces merge into one triangle.
  assert(deg == 3);
  int h0 = v_out[v], h1 = ccw(h0), h2 = ccw(h1);

  // Collect face origins.
  vector<int> all_origins;
  for (int h : {h0, h1, h2}) {
    int f = he_face[h];
    all_origins.insert(all_origins.end(), f_origin[f].begin(), f_origin[f].end());
  }
  sort(all_origins.begin(), all_origins.end());
  all_origins.erase(unique(all_origins.begin(), all_origins.end()), all_origins.end());

  // Inner rim: next(spoke_he[i]) goes nb[i] -> nb[i+1].
  int inner0 = he_next[h0], inner1 = he_next[h1], inner2 = he_next[h2];

  // Fix v_out for neighbors before deallocating.
  for (int h : {h0, h1, h2}) {
    int nbr = dest(h);
    int spoke_twin = h ^ 1;
    if (v_out[nbr] == spoke_twin) {
      // Use outer rim: twin of the inner rim arriving at nbr.
      // Inner rim arriving at nbr = next(prev_spoke), where prev_spoke
      // is the spoke before h in CCW order.
      // Actually simpler: find any outgoing he from nbr not pointing at v.
      int h_s = cw(spoke_twin);
      while (dest(h_s) == v && h_s != spoke_twin) h_s = cw(h_s);
      v_out[nbr] = h_s;
    }
  }

  // Deallocate fan faces and spokes.
  dealloc_face(he_face[h0]);
  dealloc_face(he_face[h1]);
  dealloc_face(he_face[h2]);
  dealloc_edge(h0);
  dealloc_edge(h1);
  dealloc_edge(h2);
  v_out[v] = -1;

  // Wire the three inner rim half-edges into a single triangle face.
  // inner0: nb0 -> nb1, inner1: nb1 -> nb2, inner2: nb2 -> nb0.
  he_next[inner0] = inner1;
  he_next[inner1] = inner2;
  he_next[inner2] = inner0;

  int fid = alloc_face();
  he_face[inner0] = fid;
  he_face[inner1] = fid;
  he_face[inner2] = fid;
  f_he[fid] = inner0;
  f_origin[fid] = all_origins;

  recompute_face_angles(fid);
}

void DelaunayTriangulation::remove_flat_vertices()
{
  // Sort: flat (degree-6) vertices last.
  // We process them in reverse order.
  // The from_triangulation() constructor preserves the Triangulation's vertex
  // order (which was sorted by sort_flat_last in FulleroidDelaunay).

  while (nv > 0 && v_orig_degree[nv - 1] == 6 && v_out[nv - 1] >= 0) {
    int old_nv = nv;
    remove_flat_vertex(nv - 1);

    // Check if removal succeeded (vertex is now dead).
    if (v_out[nv - 1] >= 0) {
      // Stuck — try restructuring with Delaunay flips.
      bool removed = false;
      for (int retry = 0; retry < 5; retry++) {
        flip_to_delaunay();
        remove_flat_vertex(nv - 1);
        if (v_out[nv - 1] < 0) { removed = true; break; }
      }
      if (!removed) break;
    }

    // "Remove" vertex: decrement nv (dead vertices stay in arrays but are skipped).
    // Actually, we just mark it dead via v_out[v] = -1. Decrement nv.
    while (nv > 0 && v_out[nv - 1] < 0) nv--;

    lawson_sweep();
  }

  flip_to_delaunay();
}

// --- Full algorithm ---

DelaunayTriangulation DelaunayTriangulation::compute(const Triangulation& T)
{
  // Sort flat vertices last, then build DCEL and run the algorithm.
  Triangulation sorted = T.sort_flat_last();
  DelaunayTriangulation D = from_triangulation(sorted);
  D.remove_flat_vertices();
  return D;
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

  return true;
}

// ============================================================================
// DelaunayTriangulation::embed_3d()
//
// Embeds the reduced DCEL triangulation in 3D by minimizing:
//   E = E_edge + lambda * E_cone
//
// E_edge = sum_{shortest edge per vertex pair} (|xi-xj| - Lij)^2
//   Matches extrinsic distances to intrinsic geodesic lengths.
//   For multi-edges, only the shortest is used (longest are redundant).
//
// E_cone = sum_v (angle_sum(v) - cone_angle(v))^2
//   Matches the angle sum around each vertex to the intrinsic cone angle.
//   Provides 11 independent constraints (12 minus Gauss-Bonnet redundancy)
//   that compensate for dropped multi-edge constraints.
// ============================================================================

// Per-face-angle geometry for cone angle energy.
// Computes the 3D angle at a triangle corner and its gradients.
struct FaceAngleData {
  int a, b, d;        // arm1 end, apex, arm2 end (angle at b)
  double theta;       // angle at b
  coord3d p, q;       // dtheta/d(x[a]), dtheta/d(x[d])
  double ra, rc, S;   // |a-b|, |d-b|, sin(theta)
  coord3d ua, uc;     // unit arm vectors
  bool valid;

  static FaceAngleData compute(const vector<coord3d>& x, int a, int b, int d) {
    FaceAngleData fa;
    fa.a = a; fa.b = b; fa.d = d;

    coord3d va = x[a] - x[b], vc = x[d] - x[b];
    fa.ra = va.norm(); fa.rc = vc.norm();
    fa.valid = (fa.ra > 1e-15 && fa.rc > 1e-15);
    if (!fa.valid) return fa;

    fa.ua = va / fa.ra;
    fa.uc = vc / fa.rc;
    double C = max(-1.0, min(1.0, fa.ua.dot(fa.uc)));
    fa.theta = acos(C);
    fa.S = sin(fa.theta);
    if (fa.S < 1e-10) { fa.valid = false; return fa; }

    coord3d::dangle(va, vc, fa.p, fa.q);
    return fa;
  }

  void scatter_gradient(double w, vector<coord3d>& g) const {
    g[a] = g[a] + p * w;
    g[d] = g[d] + q * w;
    g[b] = g[b] - (p + q) * w;
  }

  // Hessian-vector product with weight k.
  // H = k * [alpha * (p⊗p, q⊗q, p⊗q) + correction terms]
  void scatter_hv(double k, const vector<coord3d>& v, vector<coord3d>& Hv) const {
    double C = ua.dot(uc);
    double dev = theta - M_PI / 3.0;  // deviation used for correction terms
    double alpha = 1.0 - dev * C / S;
    matrix3d sym_ac = ua.outer(uc) + uc.outer(ua);
    matrix3d I = matrix3d::unit_matrix();

    matrix3d Haa = p.outer(p) * (k * alpha)
      + (sym_ac + I * C - ua.outer(ua) * (3*C)) * (k * dev / (ra*ra * S));
    matrix3d Hcc = q.outer(q) * (k * alpha)
      + (sym_ac + I * C - uc.outer(uc) * (3*C)) * (k * dev / (rc*rc * S));
    matrix3d Hac = p.outer(q) * (k * alpha)
      + (ua.outer(ua) + uc.outer(uc) - ua.outer(uc) * C - I) * (k * dev / (ra*rc * S));

    matrix3d Hab = Haa + Hac;
    matrix3d Hcb = Hcc + Hac.transpose();

    Hv[a] = Hv[a] + Haa * v[a] + Hac * v[d] + (-Hab) * v[b];
    Hv[d] = Hv[d] + Hac.transpose() * v[a] + Hcc * v[d] + (-Hcb) * v[b];
    Hv[b] = Hv[b] + (-Hab.transpose()) * v[a] + (-Hcb.transpose()) * v[d] + (Hab + Hcb) * v[b];
  }
};

vector<coord3d> DelaunayTriangulation::embed_3d() const
{
  // --- Step 1: Extract shortest edge per vertex pair ---
  // Multi-edges between the same pair have different lengths (different geodesics).
  // For distance matching, use the shortest one (closest to extrinsic distance).
  struct EdgeInfo { int u, v; double L; };
  vector<EdgeInfo> edges;
  map<pair<int,int>, double> shortest;

  for (int h = 0; h < nh; h += 2) {
    if (!alive(h)) continue;
    int u = he_origin[h], v = dest(h);
    if (u > v) swap(u, v);
    double L = he_length[h];
    auto key = make_pair(u, v);
    auto it = shortest.find(key);
    if (it == shortest.end() || L < it->second)
      shortest[key] = L;
  }
  for (auto& [key, L] : shortest)
    edges.push_back({key.first, key.second, L});

  // --- Step 2: APSP via Floyd-Warshall on shortest edges ---
  matrix<double> D(nv, nv, 1e18);
  for (int i = 0; i < nv; i++) D(i, i) = 0;
  for (auto& e : edges) {
    D(e.u, e.v) = e.L;
    D(e.v, e.u) = e.L;
  }
  for (int k = 0; k < nv; k++)
    for (int i = 0; i < nv; i++)
      for (int j = 0; j < nv; j++)
        D(i, j) = min(D(i, j), D(i, k) + D(k, j));

  // --- Step 3: MDS initial placement ---
  V3 x = mds_placement(D);

  // Separate collapsed vertices.  When the DCEL has a symmetry that swaps
  // two vertices with identical APSP profiles, MDS places them at exactly
  // the same point.  Perturb such pairs along a direction perpendicular to
  // the centroid-to-vertex line, by half their target distance.
  for (auto& e : edges) {
    coord3d diff = x[e.u] - x[e.v];
    double d = diff.norm();
    if (d < 0.1 * e.L) {
      // Find a direction to perturb: use cross product with a non-parallel axis
      coord3d c;
      for (auto& xi : x) c += xi;
      c /= nv;
      coord3d radial = x[e.u] - c;
      if (radial.norm() < 1e-10) radial = coord3d(1, 0, 0);
      coord3d axis(0, 0, 1);
      if (fabs(radial.dot(axis)) / radial.norm() > 0.9) axis = coord3d(0, 1, 0);
      coord3d perp = radial.cross(axis);
      perp /= max(perp.norm(), 1e-15);
      x[e.u] = x[e.u] + perp * (0.5 * e.L);
      x[e.v] = x[e.v] - perp * (0.5 * e.L);
    }
  }

  // --- Step 4: Collect per-vertex fan structure for cone angle energy ---
  // For each vertex, store the CW-ordered list of half-edges leaving it.
  // Each consecutive pair (h, cw(h)) shares a face, giving a face angle.
  struct VertexFan {
    vector<int> out_halfedges;  // CW-ordered outgoing half-edges
  };
  vector<VertexFan> fans(nv);
  for (int v = 0; v < nv; v++) {
    if (v_out[v] < 0) continue;
    int h0 = v_out[v], h = h0;
    do {
      fans[v].out_halfedges.push_back(h);
      h = cw(h);
    } while (h != h0);
  }

  // Cone angle weight: scale so E_cone is comparable to E_edge.
  // E_edge ~ n_edges * L^2, E_cone ~ n_vertices * angle^2.
  // With typical L ~ O(1) and angle deficits ~ O(1), lambda ~ 1 is reasonable.
  double lambda = 1.0;

  // --- Step 5: Stress + cone angle optimization ---
  auto eval = [&](const V3& pos, V3& g) -> double {
    double E = 0;
    g.zero();

    // E_edge: sum (|xi-xj| - Lij)^2
    for (auto& e : edges) {
      auto ed = EdgeStressData::compute(pos, e.u, e.v, e.L);
      if (!ed.valid()) continue;
      E += ed.energy();
      ed.scatter_gradient(g);
    }

    // E_cone: sum_v (angle_sum(v) - cone_angle(v))^2
    for (int v = 0; v < nv; v++) {
      auto& fan = fans[v];
      int deg = fan.out_halfedges.size();
      if (deg < 2) continue;

      double angle_sum = 0;
      vector<FaceAngleData> fa(deg);
      for (int i = 0; i < deg; i++) {
        int h = fan.out_halfedges[i];
        int d_v = dest(h);
        int h_next = fan.out_halfedges[(i + 1) % deg];
        int d_next = dest(h_next);
        // Angle at v between edges v->d_v and v->d_next
        fa[i] = FaceAngleData::compute(pos, d_v, v, d_next);
        if (fa[i].valid) angle_sum += fa[i].theta;
      }

      double dev = angle_sum - v_cone_angle[v];
      E += lambda * dev * dev;

      double w = 2.0 * lambda * dev;
      for (int i = 0; i < deg; i++)
        if (fa[i].valid) fa[i].scatter_gradient(w, g);
    }

    project_rigid_body(g, pos);
    return E;
  };

  auto hv_prod = [&](const V3& pos, const V3& dir, V3& Hv) {
    Hv.zero();

    // E_edge Hv
    for (auto& e : edges) {
      auto ed = EdgeStressData::compute(pos, e.u, e.v, e.L);
      if (!ed.valid()) continue;
      ed.scatter_hv(dir, Hv);
    }

    // E_cone Hv: H = 2*lambda * [dA⊗dA + dev * d²A/dx²]
    // where A = angle_sum, dev = A - cone_angle
    for (int v = 0; v < nv; v++) {
      auto& fan = fans[v];
      int deg = fan.out_halfedges.size();
      if (deg < 2) continue;

      double angle_sum = 0;
      vector<FaceAngleData> fa(deg);
      for (int i = 0; i < deg; i++) {
        int h = fan.out_halfedges[i];
        int d_v = dest(h);
        int h_next = fan.out_halfedges[(i + 1) % deg];
        int d_next = dest(h_next);
        fa[i] = FaceAngleData::compute(pos, d_v, v, d_next);
        if (fa[i].valid) angle_sum += fa[i].theta;
      }
      double dev = angle_sum - v_cone_angle[v];

      // Rank-1 term: 2*lambda * (dA . dir) * dA
      double dA_dot_v = 0;
      for (int i = 0; i < deg; i++)
        if (fa[i].valid) {
          dA_dot_v += fa[i].p.dot(dir[fa[i].a] - dir[v])
                    + fa[i].q.dot(dir[fa[i].d] - dir[v]);
        }
      for (int i = 0; i < deg; i++)
        if (fa[i].valid) fa[i].scatter_gradient(2.0 * lambda * dA_dot_v, Hv);

      // Correction term: 2*lambda * dev * d²A/dx² . dir
      if (fabs(dev) > 1e-15) {
        double w2 = 2.0 * lambda * dev;
        for (int i = 0; i < deg; i++)
          if (fa[i].valid) fa[i].scatter_hv(w2, dir, Hv);
      }
    }

    project_rigid_body(Hv, pos);
  };

  x = steihaug_cg(std::move(x), eval, hv_prod);

  // --- Step 6: Orient outward using DCEL face iteration ---
  coord3d c;
  for (auto& xi : x) c += xi;
  c /= nv;

  double vol = 0;
  for (int f = 0; f < nf; f++) {
    if (f_he[f] < 0) continue;
    int h0 = f_he[f];
    int h1 = he_next[h0];
    int h2 = he_next[h1];
    int u = he_origin[h0], v = he_origin[h1], w = he_origin[h2];
    vol += (x[u] - c).dot((x[v] - c).cross(x[w] - c));
  }
  if (vol < 0)
    for (auto& xi : x) xi = c * 2.0 - xi;

  return x;
}
