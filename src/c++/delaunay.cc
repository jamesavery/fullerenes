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
  for (node_t u = 0; u < N; u++)
    for (node_t v : neighbours[u])
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

  // Compute new edge length
  double f = flipped_edge_length(u, v);
  if (!std::isfinite(f) || f <= 0) return false;

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

  Graph::remove_edge(edge_t(u, v));
  set_length(u, v, 0);

  Graph::insert_edge(arc_t(B, D), u, v);
  set_length(B, D, f);

  return true;
}

int FulleroidDelaunay::flip_to_delaunay()
{
  // Lawson's algorithm: flip non-Delaunay edges using a stack.
  map<edge_t, bool> in_stack;
  stack<edge_t> S;

  for (node_t u = 0; u < N; u++)
    for (node_t v : neighbours[u])
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
    for (node_t v : neighbours[u])
      if (u < v && !is_delaunay_edge(u, v))
        return false;
  return true;
}

// ============================================================================
// Vertex removal
// ============================================================================

vector<coord2d> FulleroidDelaunay::layout_fan(node_t v) const
{
  // Lay out the fan of triangles around v in 2D.
  // Place v at origin, neighbours[v][0] along positive x-axis.
  // For each subsequent neighbor, accumulate the angle at v.

  const auto& nbrs = neighbours[v];
  int k = nbrs.size();
  vector<coord2d> pos(k);

  double d0 = get_length(v, nbrs[0]);
  pos[0] = coord2d(d0, 0);

  double angle_acc = 0;

  for (int i = 1; i < k; i++) {
    // Angle at v in triangle (v, nbrs[i-1], nbrs[i]).
    // Sides: a = length(nbrs[i-1], nbrs[i])  (opposite v)
    //        b = length(v, nbrs[i-1])
    //        c = length(v, nbrs[i])
    double a = get_length(nbrs[i-1], nbrs[i]);
    double b = get_length(v, nbrs[i-1]);
    double c = get_length(v, nbrs[i]);

    assert(a > 0 && b > 0 && c > 0);

    // cos(theta) at v by law of cosines
    double cos_theta = (b*b + c*c - a*a) / (2.0 * b * c);
    cos_theta = std::max(-1.0, std::min(1.0, cos_theta));
    double theta = acos(cos_theta);

    angle_acc += theta;

    pos[i] = coord2d(c * cos(angle_acc), c * sin(angle_acc));
  }

  return pos;
}

void FulleroidDelaunay::remove_flat_vertex(node_t v)
{
  vector<node_t> hole(neighbours[v]);
  int k = hole.size();
  assert(k >= 3);

  // 1. Lay out the fan in 2D to get distances between hole vertices.
  vector<coord2d> fan_pos = layout_fan(v);

  // 2. Find a starting vertex for the fan triangulation that is NOT already
  //    connected to non-adjacent hole vertices (to avoid multi-edges).
  int start = -1;
  for (int attempt = 0; attempt < k; attempt++) {
    bool ok = true;
    for (int i = 2; i < k - 1; i++) {
      int idx = (attempt + i) % k;
      if (edge_exists(edge_t(hole[attempt], hole[idx]))) {
        ok = false;
        break;
      }
    }
    if (ok) { start = attempt; break; }
  }

  if (start == -1) {
    // Fallback: shouldn't happen for reasonable triangulations
    start = 0;
  }

  // Rotate hole and fan positions to start at the chosen vertex.
  if (start > 0) {
    std::rotate(hole.begin(), hole.begin() + start, hole.end());
    std::rotate(fan_pos.begin(), fan_pos.begin() + start, fan_pos.end());
  }

  // 3. Compute fan edge distances: from hole[0] to hole[i].
  vector<double> fan_distances(k);
  for (int i = 0; i < k; i++) {
    double dx = fan_pos[0].first - fan_pos[i].first;
    double dy = fan_pos[0].second - fan_pos[i].second;
    fan_distances[i] = sqrt(dx*dx + dy*dy);
  }

  // 4. Remove all edges from v.
  for (int i = 0; i < k; i++) {
    Graph::remove_edge(edge_t(v, hole[i]));
    set_length(v, hole[i], 0);
  }

  // 5. Remove v from the graph (it must be the last vertex).
  assert(v == N - 1);
  neighbours.pop_back();
  N--;

  // 6. Fan-triangulate the hole from hole[0].
  //    Add edges: hole[0]-hole[2], hole[0]-hole[3], ..., hole[0]-hole[k-2].
  //
  //    For each new edge (hole[0], hole[i]):
  //      In neighbours[hole[0]]: hole[i] goes before hole[k-1]
  //      In neighbours[hole[i]]: hole[0] goes before hole[i-1]
  //
  //    This creates CCW triangles:
  //      (hole[0], hole[1], hole[2]), (hole[0], hole[2], hole[3]), ...
  //
  vector<edge_t> new_edges;
  for (int i = 2; i < k - 1; i++) {
    Graph::insert_edge(arc_t(hole[0], hole[i]), hole[k-1], hole[i-1]);
    set_length(hole[0], hole[i], fan_distances[i]);
    new_edges.push_back(edge_t(hole[0], hole[i]));
  }

  // 7. Locally Delaunayify: check all hole edges (new fan + boundary).
  map<edge_t, bool> in_stack;
  stack<edge_t> S;

  for (auto& e : new_edges) {
    S.push(e);
    in_stack[e] = true;
  }
  for (int i = 0; i < k; i++) {
    edge_t e(hole[i], hole[(i+1) % k]);
    if (!in_stack[e]) {
      S.push(e);
      in_stack[e] = true;
    }
  }

  while (!S.empty()) {
    edge_t e = S.top(); S.pop();
    in_stack[e] = false;

    node_t a = e.first, b = e.second;
    if (!edge_exists(e)) continue;

    if (!is_delaunay_edge(a, b)) {
      node_t B = next(b, a);
      node_t D = next(a, b);

      if (flip_edge(a, b)) {
        edge_t boundary[4] = {
          edge_t(a, B), edge_t(B, b), edge_t(b, D), edge_t(D, a)
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
}

void FulleroidDelaunay::remove_flat_vertices()
{
  // Record original degrees before any modifications.
  vector<int> original_degrees(N);
  for (node_t v = 0; v < N; v++)
    original_degrees[v] = neighbours[v].size();

  // Remove flat vertices from the back.
  while (N > 0 && original_degrees[N - 1] == 6) {
    remove_flat_vertex(N - 1);
  }

  // Final global Delaunay sweep for correctness.
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
