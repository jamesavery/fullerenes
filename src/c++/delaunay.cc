#include "fullerenes/delaunay.hh"

#include <stack>
#include <cmath>
#include <cassert>
#include <algorithm>

// ============================================================================
// Constructor
// ============================================================================

FulleroidDelaunay::FulleroidDelaunay(const Triangulation& T)
  : Triangulation(T.sort_nodes()), edge_lengths(N, N, 0)
{
  // Edge flips during vertex removal can temporarily push vertex degrees
  // well above 6 (the max for fullerene duals). Restride to give headroom.
  if (dmax < 20) {
    restride_inplace(20);
  }
  for (node_t u = 0; u < N; u++)
    for (node_t v : (*this)[u])
      edge_lengths(u, v) = 1.0;
}

// ============================================================================
// Intrinsic geometry from edge lengths
// ============================================================================

double FulleroidDelaunay::cot_opposite(node_t vi, node_t vj, node_t vk) const
{
  // Cotangent of the angle at vk in triangle (vi, vj, vk), opposite edge vi-vj.
  double a = get_length(vi, vj);  // side opposite vk
  double b = get_length(vi, vk);
  double c = get_length(vj, vk);

  assert(a > 0 && b > 0 && c > 0);

  // Law of cosines: cos(alpha) = (b^2 + c^2 - a^2) / (2bc)
  double cos_a = (b*b + c*c - a*a) / (2.0 * b * c);

  // sin(alpha) from half-angle tangent (more stable than sqrt(1-cos^2)):
  //   tan(alpha/2) = sqrt((a-b+c)(a+b-c) / ((a+b+c)(-a+b+c)))
  //   sin(alpha) = 2*tan(alpha/2) / (1 + tan^2(alpha/2))
  double num = (a - b + c) * (a + b - c);
  double den = (a + b + c) * (-a + b + c);

  if (den <= 0 || num < 0) {
    // Degenerate triangle (violates triangle inequality)
    return (cos_a >= 0) ? 1e15 : -1e15;
  }

  double t = sqrt(num / den);  // tan(alpha/2)
  double t2 = t * t;
  double sin_a = 2.0 * t / (1.0 + t2);

  if (sin_a < 1e-15) {
    return (cos_a >= 0) ? 1e15 : -1e15;
  }

  return cos_a / sin_a;
}

// ============================================================================
// Delaunay operations
// ============================================================================

bool FulleroidDelaunay::is_delaunay_edge(node_t u, node_t v) const
{
  //    B
  //   / \     .
  //  u---v    (edge being checked)
  //   \ /
  //    D
  //
  // Delaunay criterion: cot(angle_B) + cot(angle_D) >= 0

  node_t B = next(v, u);  // successor of u in neighbours[v]
  node_t D = next(u, v);  // successor of v in neighbours[u]

  double cot_B = cot_opposite(u, v, B);
  double cot_D = cot_opposite(u, v, D);

  const double eps = 1e-10;
  return cot_B + cot_D >= -eps;
}

double FulleroidDelaunay::flipped_edge_length(node_t u, node_t v) const
{
  //    B
  //   / \     .         After flip, new edge is B-D.
  //  u---v              Lay out the diamond in 2D and compute |B-D|.
  //   \ /
  //    D

  node_t B = next(v, u);
  node_t D = next(u, v);

  double e_uv = get_length(u, v);
  double a_uB = get_length(u, B);
  double a_vB = get_length(v, B);
  double a_uD = get_length(u, D);
  double a_vD = get_length(v, D);

  assert(e_uv > 0 && a_uB > 0 && a_vB > 0 && a_uD > 0 && a_vD > 0);

  // Place u at origin, v at (e_uv, 0).
  // Place B above the line u-v using law of cosines for angle at u.
  double cos_alpha = (e_uv*e_uv + a_uB*a_uB - a_vB*a_vB) / (2.0 * e_uv * a_uB);
  cos_alpha = std::max(-1.0, std::min(1.0, cos_alpha));
  double sin_alpha = sqrt(std::max(0.0, 1.0 - cos_alpha * cos_alpha));
  double Bx = a_uB * cos_alpha;
  double By = a_uB * sin_alpha;

  // Place D below the line u-v using law of cosines for angle at u.
  double cos_delta = (e_uv*e_uv + a_uD*a_uD - a_vD*a_vD) / (2.0 * e_uv * a_uD);
  cos_delta = std::max(-1.0, std::min(1.0, cos_delta));
  double sin_delta = sqrt(std::max(0.0, 1.0 - cos_delta * cos_delta));
  double Dx = a_uD * cos_delta;
  double Dy = -a_uD * sin_delta;  // below the line

  double dx = Bx - Dx;
  double dy = By - Dy;
  return sqrt(dx*dx + dy*dy);
}

bool FulleroidDelaunay::flip_edge(node_t u, node_t v)
{
  //    B                    B
  //   / \                  /|\    .
  //  u---v    ---->       u | v
  //   \ /                  \|/
  //    D                    D

  node_t B = next(v, u);
  node_t D = next(u, v);

  // Guard: no self-loop
  if (B == D) return false;
  // Guard: no multi-edge
  if (edge_exists(edge_t(B, D))) return false;

  // Guard: diamond must be convex (total angle < π at both u and v).
  // A non-convex diamond creates overlapping triangles with wrong geometry.
  {
    double e_uv = get_length(u, v);
    double a_uB = get_length(u, B), a_vB = get_length(v, B);
    double a_uD = get_length(u, D), a_vD = get_length(v, D);

    // Angle at u in triangle (u, B, v)
    double cos1 = (e_uv*e_uv + a_uB*a_uB - a_vB*a_vB) / (2.0 * e_uv * a_uB);
    // Angle at u in triangle (u, v, D)
    double cos2 = (e_uv*e_uv + a_uD*a_uD - a_vD*a_vD) / (2.0 * e_uv * a_uD);
    // Diamond angle at u: reject if sum of angles >= π (non-convex).
    // sin(a+b) = sin(a)cos(b) + cos(a)sin(b). Non-convex iff sin(a+b) <= 0.
    cos1 = std::max(-1.0, std::min(1.0, cos1));
    cos2 = std::max(-1.0, std::min(1.0, cos2));
    double sin1 = sqrt(std::max(0.0, 1.0 - cos1*cos1));
    double sin2 = sqrt(std::max(0.0, 1.0 - cos2*cos2));
    double sin_total_u = sin1*cos2 + cos1*sin2;
    if (sin_total_u <= 1e-12) return false;  // angle at u >= π

    // Angle at v in triangle (v, B, u)
    double cos3 = (e_uv*e_uv + a_vB*a_vB - a_uB*a_uB) / (2.0 * e_uv * a_vB);
    // Angle at v in triangle (v, D, u)
    double cos4 = (e_uv*e_uv + a_vD*a_vD - a_uD*a_uD) / (2.0 * e_uv * a_vD);
    cos3 = std::max(-1.0, std::min(1.0, cos3));
    cos4 = std::max(-1.0, std::min(1.0, cos4));
    double sin3 = sqrt(std::max(0.0, 1.0 - cos3*cos3));
    double sin4 = sqrt(std::max(0.0, 1.0 - cos4*cos4));
    double sin_total_v = sin3*cos4 + cos3*sin4;
    if (sin_total_v <= 1e-12) return false;  // angle at v >= π
  }

  // Compute new edge length
  double f = flipped_edge_length(u, v);
  if (!std::isfinite(f) || f <= 0)
    return false;

  // Topological flip.
  //
  // Before:
  //   Triangle (u,B,v) CCW => neighbours[B]: ..., v, ..., u, ... with v before u
  //   Triangle (u,v,D) CCW => neighbours[D]: ..., u, ..., v, ... with u before v
  //   neighbours[u]: ..., B, v, D, ...
  //   neighbours[v]: ..., D, u, B, ...
  //
  // After removing u-v:
  //   neighbours[u]: ..., B, D, ...  (v removed)
  //   neighbours[v]: ..., D, B, ...  (u removed)
  //   neighbours[B]: unchanged  ..., v, u, ...
  //   neighbours[D]: unchanged  ..., u, v, ...
  //
  // After adding B-D:
  //   Triangle (u,B,D) CCW => neighbours[B] needs D before u
  //   Triangle (B,v,D) CCW => neighbours[D] needs B before v
  //
  //   So: insert D before u in neighbours[B]
  //       insert B before v in neighbours[D]
  //
  //   insert_edge(arc_t(B,D), u, v):
  //     D goes before u in neighbours[B]  ✓
  //     B goes before v in neighbours[D]  ✓

  this->remove_edge(edge_t(u, v));
  set_length(u, v, 0);

  this->insert_edge(arc_t(B, D), u, v);
  set_length(B, D, f);

  return true;
}

int FulleroidDelaunay::flip_to_delaunay()
{
  // Lawson's algorithm: flip non-Delaunay edges using a stack.
  map<edge_t, bool> in_stack;
  stack<edge_t> S;

  for (node_t u = 0; u < N; u++)
    for (node_t v : (*this)[u])
      if (u < v) {
        edge_t e(u, v);
        if (!in_stack[e]) {
          S.push(e);
          in_stack[e] = true;
        }
      }

  int flips = 0;
  while (!S.empty()) {
    edge_t e = S.top(); S.pop();
    in_stack[e] = false;

    node_t u = e.first, v = e.second;
    if (!edge_exists(e)) continue;

    if (!is_delaunay_edge(u, v)) {
      node_t B = next(v, u);
      node_t D = next(u, v);

      if (flip_edge(u, v)) {
        flips++;

        edge_t boundary[4] = {
          edge_t(u, B), edge_t(B, v), edge_t(v, D), edge_t(D, u)
        };
        for (auto& ec : boundary) {
          if (!in_stack[ec]) {
            S.push(ec);
            in_stack[ec] = true;
          }
        }
      }
    }
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
  // Reduce v's degree to 3 via edge flips, then remove.
  // If stuck at degree 4 (all diamonds degenerate), use direct star retriangulation.

  const int max_outer = 20;
  for (int outer = 0; degree(v) > 3; outer++) {
    if (outer >= max_outer) break;

    // Phase 1: try direct incident flips and blocker+incident combos.
    {
      bool made_progress = true;
      while (made_progress && degree(v) > 3) {
        made_progress = false;

        vector<node_t> vnbrs((*this)[v].begin(), (*this)[v].end());
        for (node_t u : vnbrs) {
          if (flip_edge(v, u)) { made_progress = true; break; }
        }
        if (made_progress) continue;

        vnbrs = vector<node_t>((*this)[v].begin(), (*this)[v].end());
        for (node_t u : vnbrs) {
          node_t B = next(u, v);
          node_t D = next(v, u);
          if (B == D) continue;
          if (edge_exists(edge_t(B, D)) && flip_edge(B, D)) {
            if (flip_edge(v, u)) { made_progress = true; break; }
          }
        }
      }
      if (degree(v) <= 3) break;
    }

    // If stuck at degree 4, use direct star retriangulation.
    if (degree(v) == 4) break;

    // Phase 2: flip a cross-edge to change local geometry, then retry.
    {
      vector<node_t> vnbrs((*this)[v].begin(), (*this)[v].end());
      int k = vnbrs.size();
      bool cross_ok = false;
      for (int i = 0; i < k; i++) {
        node_t a = vnbrs[i];
        node_t b = vnbrs[(i + 1) % k];
        if (flip_edge(a, b)) {
          cross_ok = true;
          bool any = true;
          while (any) {
            any = false;
            vector<node_t> cur((*this)[v].begin(), (*this)[v].end());
            for (node_t u : cur) {
              if (flip_edge(v, u)) { any = true; break; }
            }
          }
          break;
        }
      }
      if (!cross_ok) {
        flip_to_delaunay();
      }
    }
  }

  int deg = degree(v);

  if (deg == 4) {
    // Direct star retriangulation for the degree-4 stuck case.
    //
    // Unfold the 4 triangles around v into a flat polygon (cone angle = 2π),
    // compute both quad diagonal lengths from the 2D layout, pick a valid one,
    // and replace v with that single diagonal edge.
    //
    //   Before:          After (diagonal A-C):
    //     A                  A
    //    /|\                / \    .
    //   B-v-D     →       B   D
    //    \|/                \ /
    //     C                  C

    auto row = (*this)[v];
    node_t nb[4] = {row[0], row[1], row[2], row[3]};

    // Edge lengths from v to each neighbor and between consecutive neighbors.
    double lv[4], lx[4];
    for (int i = 0; i < 4; i++) {
      lv[i] = get_length(v, nb[i]);
      lx[i] = get_length(nb[i], nb[(i + 1) % 4]);
    }

    // Unfold: place v at origin, nb[0] along +x axis.
    double px[4], py[4];
    px[0] = lv[0]; py[0] = 0;
    double cum_angle = 0;
    for (int i = 1; i < 4; i++) {
      double a = lx[i - 1];  // side opposite v
      double b = lv[i - 1];  // from v to prev neighbor
      double c = lv[i];      // from v to current neighbor
      double cos_v = (b*b + c*c - a*a) / (2.0 * b * c);
      cos_v = std::max(-1.0, std::min(1.0, cos_v));
      cum_angle += acos(cos_v);
      px[i] = c * cos(cum_angle);
      py[i] = c * sin(cum_angle);
    }

    // Compute both diagonal lengths.
    double dx02 = px[0] - px[2], dy02 = py[0] - py[2];
    double len02 = sqrt(dx02*dx02 + dy02*dy02);  // nb[0]-nb[2]
    double dx13 = px[1] - px[3], dy13 = py[1] - py[3];
    double len13 = sqrt(dx13*dx13 + dy13*dy13);  // nb[1]-nb[3]

    // Check which diagonals are topologically valid (no multi-edge).
    bool can02 = !edge_exists(edge_t(nb[0], nb[2])) && len02 > 1e-15;
    bool can13 = !edge_exists(edge_t(nb[1], nb[3])) && len13 > 1e-15;

    // Check geometric validity: diagonal must be interior to the quad.
    // Points on opposite sides of the diagonal line ⟹ valid.
    auto cross2d = [](double ax, double ay, double bx, double by) {
      return ax * by - ay * bx;
    };
    if (can02) {
      double s1 = cross2d(px[2]-px[0], py[2]-py[0], px[1]-px[0], py[1]-py[0]);
      double s3 = cross2d(px[2]-px[0], py[2]-py[0], px[3]-px[0], py[3]-py[0]);
      if (s1 * s3 >= 0) can02 = false;
    }
    if (can13) {
      double s0 = cross2d(px[3]-px[1], py[3]-py[1], px[0]-px[1], py[0]-py[1]);
      double s2 = cross2d(px[3]-px[1], py[3]-py[1], px[2]-px[1], py[2]-py[1]);
      if (s0 * s2 >= 0) can13 = false;
    }

    assert((can02 || can13) && "remove_deg4: neither diagonal is valid");

    // Choose the Delaunay diagonal: the one where cot-sum >= 0 in the
    // resulting two triangles.  Fall back to whichever is valid.
    bool use02 = can02;
    if (can02 && can13) {
      // Pick the one with smaller circumradius ratio (more Delaunay).
      // Simple heuristic: pick the shorter diagonal.
      use02 = (len02 <= len13);
    }

    // Save neighbor list, then remove v's 4 edges.
    for (int i = 0; i < 4; i++)
      set_length(v, nb[i], 0);
    // Remove in reverse to keep indices stable (v is last vertex).
    for (int i = 3; i >= 0; i--)
      this->remove_edge(edge_t(v, nb[i]));

    // After removing v's edges, neighbor lists become:
    //   nb[0]: ..., nb[1], nb[3], ...   (v was between nb[1] and nb[3])
    //   nb[1]: ..., nb[2], nb[0], ...
    //   nb[2]: ..., nb[3], nb[1], ...
    //   nb[3]: ..., nb[0], nb[2], ...
    //
    // For diagonal nb[0]-nb[2]: insert_edge(arc_t(nb[0],nb[2]), nb[3], nb[1])
    // For diagonal nb[1]-nb[3]: insert_edge(arc_t(nb[1],nb[3]), nb[0], nb[2])

    if (use02) {
      this->insert_edge(arc_t(nb[0], nb[2]), nb[3], nb[1]);
      set_length(nb[0], nb[2], len02);
    } else {
      this->insert_edge(arc_t(nb[1], nb[3]), nb[0], nb[2]);
      set_length(nb[1], nb[3], len13);
    }

    // Remove v from the graph.
    assert(v == N - 1);
    assert(degree(v) == 0);
    pop_back();
    return;
  }

  // v has degree 3. The three surrounding triangles collapse into one.
  assert(deg == 3);
  auto vrow = (*this)[v];
  node_t a = vrow[0], b = vrow[1], c = vrow[2];

  set_length(v, a, 0);
  set_length(v, b, 0);
  set_length(v, c, 0);
  this->remove_edge(edge_t(v, a));
  this->remove_edge(edge_t(v, b));
  this->remove_edge(edge_t(v, c));

  assert(v == N - 1);
  pop_back();
}

void FulleroidDelaunay::remove_flat_vertices()
{
  // Record original degrees before any modifications.
  vector<int> original_degrees(N);
  for (node_t v = 0; v < N; v++)
    original_degrees[v] = degree(v);

  // Remove flat vertices from the back, with local Delaunay flipping
  // after each removal to keep the triangulation close to Delaunay.
  while (N > 0 && original_degrees[N - 1] == 6) {
    remove_flat_vertex(N - 1);
    flip_to_delaunay();
  }
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
