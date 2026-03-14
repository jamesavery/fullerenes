#include "fullerenes/delaunay.hh"

#include <stack>
#include <cmath>
#include <cassert>
#include <algorithm>
#include <numeric>
#include <map>

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

  return true;
}

int FulleroidDelaunay::flip_to_delaunay()
{
  // Two-phase Lawson algorithm for intrinsic Delaunay triangulations.
  //
  // Phase 1: Standard Lawson — flip all unblocked non-Delaunay edges.
  // Phase 2: When only blocked edges remain (multi-edge prevents flip),
  //   resolve ONE blocked edge via support flip (flip a rim edge to change
  //   the diamond geometry), then return to Phase 1.

  int flips = 0;

  for (int round = 0; round < 50; round++) {
    // Phase 1: Exhaust all unblocked non-Delaunay flips.
    {
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
    }

    if (is_delaunay()) return flips;

    // Phase 2: Find one blocked non-Delaunay edge and try support flips.
    bool resolved = false;
    for (node_t u = 0; u < N && !resolved; u++)
      for (node_t v : (*this)[u]) {
        if (u >= v || is_delaunay_edge(u, v)) continue;
        node_t B = next(v, u), D = next(u, v);
        bool blocked = (B != D) && edge_exists(edge_t(std::min(B,D), std::max(B,D)));

        if (!blocked) continue;  // Phase 1 handles unblocked edges

        // Edge (u,v) is non-Delaunay, blocked by multi-edge (B,D).
        // Try flipping each rim edge to change the diamond geometry.
        edge_t e(u, v);
        for (auto [ra, rb] : {pair{(int)v,(int)D}, {(int)D,u}, {u,(int)B}, {(int)B,(int)v}}) {
          edge_t rim(std::min(ra,rb), std::max(ra,rb));
          if (!edge_exists(rim)) continue;
          node_t rB = next(rb, ra), rD = next(ra, rb);
          if (!flip_edge(ra, rb)) continue;
          flips++;

          // Check: did the support flip resolve the target AND not create
          // a new non-Delaunay edge?  If the created edge is non-Delaunay,
          // Phase 1 will flip it right back, causing an infinite cycle.
          edge_t created(std::min(rB,rD), std::max(rB,rD));
          bool target_ok = edge_exists(e) && is_delaunay_edge(u, v);
          bool created_ok = !edge_exists(created) || is_delaunay_edge(rB, rD);

          if (target_ok && created_ok) {
            resolved = true;
            break;
          }

          // Support flip didn't help or created new problem — undo.
          edge_t new_rim(std::min(rB,rD), std::max(rB,rD));
          if (edge_exists(new_rim) && flip_edge(rB, rD)) {
            flips--;
          } else {
            resolved = true;  // can't undo, state changed — re-run Phase 1
            break;
          }
        }
        if (resolved) break;

        // Atomic blocker resolution: flip blocker, then original.
        if (flip_edge(B, D)) {
          flips++;
          resolved = true;
          if (flip_edge(u, v)) flips++;
        }
        break;
      }

    if (!resolved) break;
  }

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
}

// Phase 1 only: standard Lawson flips, no blocker resolution.
// Safe to call between vertex removals — terminates on polyhedral surfaces.
int FulleroidDelaunay::flip_to_delaunay_phase1()
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
      bool removed = false;
      for (int retry = 0; retry < 5; retry++) {
        flip_to_delaunay();
        remove_flat_vertex(N - 1);
        if (N < old_N) { removed = true; break; }
      }
      if (!removed) break;
    }

    // Between removals: only do Phase 1 (standard Lawson) flips.
    // Phase 2 (blocker resolution) on graphs still containing flat vertices
    // can cycle because high-degree flat vertices create multi-edge blocker
    // patterns where support flips and Phase 1 undo each other.
    flip_to_delaunay_phase1();
  }

  // Full Delaunay flipping (with blocker resolution) on the final graph.
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
