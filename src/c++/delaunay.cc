#include "fullerenes/delaunay.hh"
#include "fullerenes/symmetry.hh"
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

// Old FulleroidDelaunay + IDTAudit implementation moved to delaunay_old.cc.

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

  // BFS-unfold original equilateral faces into Z[omega].
  // Two-phase BFS: first through target faces only (to avoid wrapping on
  // closed surfaces), then through all faces if required vertices are still
  // unreached.  If anchor_vertex >= 0, starts from a target face containing
  // that vertex (ensures the unfolding is centered on the fan).
  // Returns a map from original vertex ID to Eisenstein grid position.
  unordered_map<int, Eisenstein> unfold_patch(
      const vector<int>& face_ids,
      const vector<int>& required_vertices = {},
      int anchor_vertex = -1) const;

  // Classify original faces across the directed line vtx_A → vtx_B.
  // Returns (left_of_line, right_of_line).  Faces whose centroid is exactly
  // on the line are assigned to both sides.
  pair<vector<int>, vector<int>>
  classify_across_line(const vector<int>& face_ids,
                       int vtx_A, int vtx_B) const;

  // Classify original faces into a set of triangles (for ear-clipping).
  // tri_verts[j] = {a, b, c} gives the CCW vertex IDs of triangle j.
  // anchor_vertex: the removed flat vertex (ensures unfolding is fan-centered).
  // Returns assignment[j] = list of original face IDs inside triangle j.
  vector<vector<int>>
  classify_into_triangles(const vector<int>& face_ids,
                          const vector<array<int,3>>& tri_verts,
                          int anchor_vertex = -1) const;
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
    const vector<int>& face_ids,
    const vector<int>& required_vertices,
    int anchor_vertex) const
{
  if (face_ids.empty()) return {};

  // Two-phase BFS unfolding.
  //
  // On a closed surface, naive BFS can reach a vertex via a path that wraps
  // around the surface, giving it coordinates inconsistent with the local
  // geometry of the target patch.  To avoid this:
  //
  // Phase 1: expand only through TARGET faces.  This keeps the unfolding
  //   within the patch (e.g., a fan polygon around a removed vertex).
  //   All target face vertices get consistent local coordinates.
  //
  // Phase 2: if required vertices are still unreached, expand through
  //   ALL adjacent faces to find them.  Positions set in phase 1 are
  //   never overwritten.

  unordered_set<int> target_faces(face_ids.begin(), face_ids.end());
  unordered_set<int> target_verts(required_vertices.begin(),
                                  required_vertices.end());
  unordered_map<int, Eisenstein> pos;
  unordered_set<int> placed;
  queue<int> Q;

  int faces_remaining = target_faces.size();
  int verts_remaining = target_verts.size();

  auto done = [&]() { return faces_remaining <= 0 && verts_remaining <= 0; };

  auto place_vertex = [&](int vtx, Eisenstein p) {
    if (!pos.count(vtx)) {
      pos[vtx] = p;
      if (target_verts.count(vtx)) verts_remaining--;
    }
  };

  // Helper: expand from face f, placing the third vertex of each adjacent
  // face.  If target_only, skip non-target adjacent faces.
  auto expand = [&](int f, bool target_only) {
    auto& v = face_verts[f];
    for (int i = 0; i < 3; i++) {
      int eu = v[i], ev = v[(i+1) % 3];
      auto it = arc_to_face.find(int64_t(ev) * N + eu);
      if (it == arc_to_face.end()) continue;
      int adj = it->second;
      if (placed.count(adj)) continue;
      if (target_only && !target_faces.count(adj)) continue;

      auto& av = face_verts[adj];
      int third = av[0];
      if (third == eu || third == ev) third = av[1];
      if (third == eu || third == ev) third = av[2];

      place_vertex(third, pos[ev] + (pos[eu] - pos[ev]).nextCCW());
      placed.insert(adj);
      if (target_faces.count(adj)) faces_remaining--;
      Q.push(adj);
    }
  };

  // Choose starting face: prefer one containing anchor_vertex (ensures the
  // unfolding starts at the fan center for vertex removal).
  int f0 = face_ids[0];
  if (anchor_vertex >= 0) {
    for (int fid : face_ids) {
      auto& fv2 = face_verts[fid];
      if (fv2[0] == anchor_vertex || fv2[1] == anchor_vertex || fv2[2] == anchor_vertex) {
        f0 = fid;
        break;
      }
    }
  }
  auto& fv = face_verts[f0];
  place_vertex(fv[0], Eisenstein(0, 0));
  place_vertex(fv[1], Eisenstein(1, 0));
  place_vertex(fv[2], Eisenstein(0, 1));
  placed.insert(f0);
  faces_remaining--;
  Q.push(f0);

  // Phase 1: BFS through target faces only.
  while (!Q.empty() && !done()) {
    int f = Q.front(); Q.pop();
    expand(f, /*target_only=*/true);
  }

  // If target faces are disconnected, seed each unvisited component.
  if (faces_remaining > 0) {
    for (int fid : face_ids) {
      if (placed.count(fid)) continue;
      // Try to find an adjacent placed face to seed from.
      auto& fv2 = face_verts[fid];
      bool seeded = false;
      for (int i = 0; i < 3 && !seeded; i++) {
        int eu = fv2[i], ev = fv2[(i+1) % 3];
        // Check if the adjacent face across this edge is already placed.
        auto it = arc_to_face.find(int64_t(ev) * N + eu);
        if (it != arc_to_face.end() && placed.count(it->second)) {
          // Place this face's third vertex.
          int third = fv2[0];
          if (third == eu || third == ev) third = fv2[1];
          if (third == eu || third == ev) third = fv2[2];
          place_vertex(third, pos[ev] + (pos[eu] - pos[ev]).nextCCW());
          placed.insert(fid);
          faces_remaining--;
          Q.push(fid);
          seeded = true;
        }
      }
      if (!seeded) {
        // Bridge: expand through non-target faces to reach this component.
        // Re-enqueue all placed faces and expand without target restriction
        // until we reach this unvisited target face.
        queue<int> bridge_Q;
        for (int p : placed) bridge_Q.push(p);
        while (!bridge_Q.empty() && !placed.count(fid)) {
          int f = bridge_Q.front(); bridge_Q.pop();
          auto& v = face_verts[f];
          for (int i = 0; i < 3; i++) {
            int eu = v[i], ev = v[(i+1) % 3];
            auto it = arc_to_face.find(int64_t(ev) * N + eu);
            if (it == arc_to_face.end()) continue;
            int adj = it->second;
            if (placed.count(adj)) continue;
            auto& av = face_verts[adj];
            int third = av[0];
            if (third == eu || third == ev) third = av[1];
            if (third == eu || third == ev) third = av[2];
            place_vertex(third, pos[ev] + (pos[eu] - pos[ev]).nextCCW());
            placed.insert(adj);
            if (target_faces.count(adj)) faces_remaining--;
            bridge_Q.push(adj);
          }
        }
      }
      // Continue phase 1 BFS from any newly added faces.
      while (!Q.empty() && !done()) {
        int f = Q.front(); Q.pop();
        expand(f, /*target_only=*/true);
      }
    }
  }

  // Phase 2: expand through all faces if required vertices still unreached.
  if (verts_remaining > 0) {
    // Re-seed from all placed faces.
    queue<int>().swap(Q);  // clear
    for (int f : placed) Q.push(f);
    while (!Q.empty() && !done()) {
      int f = Q.front(); Q.pop();
      expand(f, /*target_only=*/false);
    }
  }

  return pos;
}

pair<vector<int>, vector<int>>
DelaunayTriangulation::OriginTracker::classify_across_line(
    const vector<int>& face_ids, int vtx_A, int vtx_B) const
{
  if (face_ids.empty()) return {{}, {}};
  auto pos = unfold_patch(face_ids, {vtx_A, vtx_B});

  // Scale line endpoints by 3 so we can test the centroid (p+q+r)
  // without dividing by 3.  sign(turn(3A, 3B, p+q+r)) = sign(turn(A, B, centroid)).
  Eisenstein A3 = pos.at(vtx_A) * 3;
  Eisenstein B3 = pos.at(vtx_B) * 3;

  vector<int> left, right;
  for (int fid : face_ids) {
    auto& fv = face_verts[fid];
    Eisenstein sum = pos.at(fv[0]) + pos.at(fv[1]) + pos.at(fv[2]);
    int t = Eisenstein::turn(A3, B3, sum);
    if (t >= 0) left.push_back(fid);     // on the line → left (by convention)
    else right.push_back(fid);
  }
  return {left, right};
}

vector<vector<int>>
DelaunayTriangulation::OriginTracker::classify_into_triangles(
    const vector<int>& face_ids,
    const vector<array<int,3>>& tri_verts,
    int anchor_vertex) const
{
  // Collect all triangle corner vertices as required.
  vector<int> req;
  for (auto& tv : tri_verts)
    for (int v : tv)
      req.push_back(v);

  auto pos = unfold_patch(face_ids, req, anchor_vertex);
  int nt = tri_verts.size();
  vector<vector<int>> assignment(nt);

  // Precompute scaled triangle corners (3x to avoid centroid division).
  vector<array<Eisenstein,3>> corners(nt);
  for (int j = 0; j < nt; j++)
    for (int k = 0; k < 3; k++)
      corners[j][k] = pos.at(tri_verts[j][k]) * 3;

  for (int fid : face_ids) {
    auto& fv = face_verts[fid];
    Eisenstein sum = pos.at(fv[0]) + pos.at(fv[1]) + pos.at(fv[2]);

    for (int j = 0; j < nt; j++) {
      auto& c = corners[j];
      if (Eisenstein::turn(c[0], c[1], sum) >= 0 &&
          Eisenstein::turn(c[1], c[2], sum) >= 0 &&
          Eisenstein::turn(c[2], c[0], sum) >= 0) {
        assignment[j].push_back(fid);
        break;
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

  // Topological guard: B == D means the flip would create a self-loop.
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
  if (origin_tracker) {
    vector<int> all;
    all.reserve(f_origin[fh].size() + f_origin[ft].size());
    all.insert(all.end(), f_origin[fh].begin(), f_origin[fh].end());
    all.insert(all.end(), f_origin[ft].begin(), f_origin[ft].end());
    sort(all.begin(), all.end());
    all.erase(unique(all.begin(), all.end()), all.end());

    // Classify each original face by which side of B→D its centroid
    // falls on, using Eisenstein turn() in the Z[omega] grid.
    auto [left, right] = origin_tracker->classify_across_line(all, B, D);
    f_origin[fh] = std::move(left);
    f_origin[ft] = std::move(right);
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
  bool complete = false;        // true if all ears were successfully clipped
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
  for (int i = 0; i < fan.k; i++) {
    double si = fan.spokes[i], sn = fan.spokes[(i+1)%fan.k], r = fan.rims[i];
    double cos_a = std::clamp((si*si + sn*sn - r*r) / (2*si*sn), -1.0, 1.0);
    fan.cum[i+1] = fan.cum[i] + acos(cos_a);
  }

  return fan;
}

// Ear-clip the fan polygon into triangles.
// May mutate D via blocker-flips when no valid ear exists.
static FanTriangulation ear_clip_fan(DelaunayTriangulation& D,
                                      const FanPolygon& fan) {
  int k = fan.k;
  FanTriangulation tri;

  vector<int> poly(k);
  for (int i = 0; i < k; i++) poly[i] = i;

  while ((int)poly.size() > 3) {
    int n = poly.size();
    bool found = false;

    for (int j = 0; j < n; j++) {
      int jm = (j - 1 + n) % n, jp = (j + 1) % n;
      int pp = poly[jm], pi = poly[j], pn = poly[jp];

      // Self-loop: diagonal endpoints must be distinct global vertices.
      if (fan.nb[pp] == fan.nb[pn]) continue;

      // Multi-edge: diagonal must not duplicate an existing edge.
      bool blocked = false;
      { int h0 = D.v_out[fan.nb[pp]], hc = h0;
        if (h0 >= 0) do {
          if (D.dest(hc) == fan.nb[pn]) { blocked = true; break; }
          hc = D.cw(hc);
        } while (hc != h0); }
      for (auto& d : tri.diagonals)
        if ((fan.nb[d.from] == fan.nb[pp] && fan.nb[d.to] == fan.nb[pn]) ||
            (fan.nb[d.from] == fan.nb[pn] && fan.nb[d.to] == fan.nb[pp]))
          { blocked = true; break; }
      if (blocked) continue;

      // Ear triangle must have positive area.
      if (fan.ear_area(pp, pi, pn) <= 1e-10) continue;

      // Subtended angle < pi (diagonal is interior to the polygon).
      double sub_angle = (pn > pp) ? fan.cum[pn] - fan.cum[pp]
                                    : (fan.cum[k] - fan.cum[pp]) + fan.cum[pn];
      if (sub_angle > M_PI + 1e-10) continue;

      // Diagonal must have positive length.
      double len = fan.diag_length(pp, pn);
      if (len <= 1e-15) continue;

      tri.diagonals.push_back({pp, pi, pn, len});
      poly.erase(poly.begin() + j);
      found = true;
      break;
    }

    if (!found) {
      // --- Blocker-flip: try flipping blocking edges ---
      // (Ugly imperative code — to be cleaned up separately.)
      bool flipped = false;
      int n2 = poly.size();
      for (int j = 0; j < n2 && !flipped; j++) {
        int jm = (j - 1 + n2) % n2, jp = (j + 1) % n2;
        int pp = poly[jm], pn = poly[jp];
        if (fan.nb[pp] == fan.nb[pn]) continue;
        int h0 = D.v_out[fan.nb[pp]], hc = h0;
        if (h0 >= 0) do {
          if (D.dest(hc) == fan.nb[pn]) {
            if (D.flip_edge(hc)) { flipped = true; break; }
            if (D.flip_edge(hc ^ 1)) { flipped = true; break; }
          }
          hc = D.cw(hc);
        } while (hc != h0);
      }
      if (!flipped) return tri;  // stuck — return incomplete
    }
  }

  // Build triangle list from diagonals.
  { vector<int> rpoly(k);
    for (int i = 0; i < k; i++) rpoly[i] = i;
    for (auto& ear : tri.diagonals) {
      tri.triangles.push_back({ear.from, ear.ear, ear.to});
      rpoly.erase(std::find(rpoly.begin(), rpoly.end(), ear.ear));
    }
    assert(rpoly.size() == 3);
    tri.triangles.push_back({rpoly[0], rpoly[1], rpoly[2]});
  }

  tri.complete = true;
  return tri;
}

// Compute per-triangle origin assignments via Eisenstein unfolding.
// Each original face in the fan is mapped to a 2D point in the fan's
// coordinate system, then classified into one of the ear triangles.
static vector<vector<int>> compute_origin_assignment(
    const DelaunayTriangulation& D, int v,
    const FanPolygon& fan, const FanTriangulation& tri,
    const vector<vector<int>>& sector_origins)
{
  int k = fan.k;
  int nt = tri.triangles.size();

  // Fan 2D positions of boundary vertices.
  vector<double> fan_x(k), fan_y(k);
  for (int i = 0; i < k; i++) {
    fan_x[i] = fan.x(i);
    fan_y[i] = fan.y(i);
  }

  // Ear triangle vertices in fan 2D.
  vector<array<double,2>> ear_v0(nt), ear_v1(nt), ear_v2(nt);
  for (int ti = 0; ti < nt; ti++) {
    ear_v0[ti] = {fan_x[tri.triangles[ti].v0], fan_y[tri.triangles[ti].v0]};
    ear_v1[ti] = {fan_x[tri.triangles[ti].v1], fan_y[tri.triangles[ti].v1]};
    ear_v2[ti] = {fan_x[tri.triangles[ti].v2], fan_y[tri.triangles[ti].v2]};
  }

  // 2D cross product sign for point-in-triangle test.
  auto cross2d_sign = [](double ax, double ay, double bx, double by,
                          double px, double py) -> int {
    double c = (bx - ax) * (py - ay) - (by - ay) * (px - ax);
    return (c > 1e-12) ? 1 : (c < -1e-12) ? -1 : 0;
  };

  // Classify a 2D point into an ear triangle.
  auto classify_point = [&](double px, double py) -> int {
    for (int ti = 0; ti < nt; ti++) {
      int t0 = cross2d_sign(ear_v0[ti][0], ear_v0[ti][1],
                            ear_v1[ti][0], ear_v1[ti][1], px, py);
      int t1 = cross2d_sign(ear_v1[ti][0], ear_v1[ti][1],
                            ear_v2[ti][0], ear_v2[ti][1], px, py);
      int t2 = cross2d_sign(ear_v2[ti][0], ear_v2[ti][1],
                            ear_v0[ti][0], ear_v0[ti][1], px, py);
      if (t0 >= 0 && t1 >= 0 && t2 >= 0) return ti;
      if (t0 <= 0 && t1 <= 0 && t2 <= 0) return ti;  // CW winding
    }
    return -1;
  };

  // Find the ear triangle adjacent to sector sec (contains rim edge sec -> sec+1).
  auto sector_fallback = [&](int sec) -> int {
    int sn = (sec + 1) % k;
    for (int ti = 0; ti < nt; ti++) {
      auto& t = tri.triangles[ti];
      if ((t.v0 == sec || t.v1 == sec || t.v2 == sec) &&
          (t.v0 == sn  || t.v1 == sn  || t.v2 == sn))
        return ti;
    }
    return 0;
  };

  vector<vector<int>> assignment(nt);

  for (int sec = 0; sec < k; sec++) {
    auto& sec_origs = sector_origins[sec];
    if (sec_origs.empty()) continue;

    // Unfold this sector's original faces in Eisenstein coordinates.
    auto pos = D.origin_tracker->unfold_patch(
        sec_origs, {v, fan.nb[sec], fan.nb[(sec+1)%k]}, v);

    auto it_v = pos.find(v);
    auto it_a = pos.find(fan.nb[sec]);
    auto it_b = pos.find(fan.nb[(sec+1)%k]);

    if (it_v == pos.end() || it_a == pos.end() || it_b == pos.end()) {
      // Couldn't unfold sector — fallback to rim assignment.
      int fb = sector_fallback(sec);
      for (int fid : sec_origs) assignment[fb].push_back(fid);
      continue;
    }

    // Affine transform from Eisenstein to fan 2D via barycentric coords
    // in the sector triangle (v, nb[sec], nb[(sec+1)%k]).
    double eVx = it_v->second.first, eVy = it_v->second.second;
    double eAx = it_a->second.first, eAy = it_a->second.second;
    double eBx = it_b->second.first, eBy = it_b->second.second;
    double fAx = fan_x[sec], fAy = fan_y[sec];
    double fBx = fan_x[(sec+1)%k], fBy = fan_y[(sec+1)%k];
    double det = (eAx - eVx) * (eBy - eVy) - (eBx - eVx) * (eAy - eVy);

    for (int fid : sec_origs) {
      auto& fverts = D.origin_tracker->face_verts[fid];

      // Centroid in Eisenstein coordinates.
      double cx = 0, cy = 0;
      bool all_found = true;
      for (int vi = 0; vi < 3; vi++) {
        auto it = pos.find(fverts[vi]);
        if (it == pos.end()) { all_found = false; break; }
        cx += it->second.first;
        cy += it->second.second;
      }

      if (!all_found || std::abs(det) < 1e-15) {
        assignment[sector_fallback(sec)].push_back(fid);
        continue;
      }
      cx /= 3.0; cy /= 3.0;

      // Barycentric coords in Eisenstein triangle -> fan 2D position.
      double ba = ((cx - eVx)*(eBy - eVy) - (cy - eVy)*(eBx - eVx)) / det;
      double bb = ((cy - eVy)*(eAx - eVx) - (cx - eVx)*(eAy - eVy)) / det;
      double fan_px = ba * fAx + bb * fBx;
      double fan_py = ba * fAy + bb * fBy;

      int ti = classify_point(fan_px, fan_py);
      assignment[ti >= 0 ? ti : sector_fallback(sec)].push_back(fid);
    }
  }

  return assignment;
}

// Replace the star of v with the fan triangulation (DCEL surgery).
// Deallocates spokes and fan faces, allocates diagonal edges and ear faces,
// wires the inner rim and diagonals into triangular face boundaries.
static void splice_fan(DelaunayTriangulation& D, int v,
                        const FanPolygon& fan, const FanTriangulation& tri) {
  int k = fan.k;

  // --- Collect face origins before deallocation ---
  vector<vector<int>> sector_origins(k);
  if (D.origin_tracker) {
    int h0 = D.v_out[v], h = h0;
    int sec = 0;
    unordered_set<int> seen;
    do {
      int f = D.he_face[h];
      for (int fid : D.f_origin[f])
        if (seen.insert(fid).second)
          sector_origins[sec].push_back(fid);
      h = D.ccw(h);
      sec++;
    } while (h != h0);
  }

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

  // --- Build arc-to-halfedge map ---
  map<pair<int,int>, int> local_arc;

  for (int i = 0; i < k; i++)
    local_arc[{i, (i + 1) % k}] = fan.inner_rim[i];

  for (auto& d : tri.diagonals) {
    int h_d = D.alloc_edge();
    D.he_origin[h_d]     = fan.nb[d.from];
    D.he_origin[h_d ^ 1] = fan.nb[d.to];
    D.he_length[h_d]     = d.length;
    D.he_length[h_d ^ 1] = d.length;
    local_arc[{d.from, d.to}] = h_d;
    local_arc[{d.to, d.from}] = h_d ^ 1;
  }

  // --- Compute origin assignments ---
  vector<vector<int>> origin_assignment;
  if (D.origin_tracker)
    origin_assignment = compute_origin_assignment(D, v, fan, tri, sector_origins);

  // --- Wire up each triangle ---
  for (size_t ti = 0; ti < tri.triangles.size(); ti++) {
    auto& t = tri.triangles[ti];
    int fid = D.alloc_face();

    int h_01 = local_arc.at({t.v0, t.v1});
    int h_12 = local_arc.at({t.v1, t.v2});
    int h_20 = local_arc.at({t.v2, t.v0});

    D.he_next[h_01] = h_12;
    D.he_next[h_12] = h_20;
    D.he_next[h_20] = h_01;
    D.he_face[h_01] = fid;
    D.he_face[h_12] = fid;
    D.he_face[h_20] = fid;
    D.f_he[fid] = h_01;

    if (D.origin_tracker)
      D.f_origin[fid] = std::move(origin_assignment[ti]);
  }

  // --- Recompute angles for new faces ---
  for (auto& t : tri.triangles) {
    int h_01 = local_arc.at({t.v0, t.v1});
    D.recompute_face_angles(D.he_face[h_01]);
  }

  // --- Safety: fix v_out if still broken ---
  for (int i = 0; i < k; i++) {
    int u = fan.nb[i];
    if (D.v_out[u] < 0 || !D.alive(D.v_out[u]) || D.he_origin[D.v_out[u]] != u) {
      for (auto& [arc, hid] : local_arc)
        if (fan.nb[arc.first] == u && D.alive(hid)) { D.v_out[u] = hid; break; }
    }
  }
}

// Reduce the star degree of vertex v by flipping incident edges.
static void reduce_star_degree(DelaunayTriangulation& D, int v, int target) {
  int deg = D.vertex_degree(v);
  while (deg > target) {
    bool progress = false;
    int h0 = D.v_out[v], h = h0;
    do {
      if (D.flip_edge(h)) { deg--; progress = true; break; }
      h = D.cw(h);
    } while (h != h0);
    if (!progress) break;
  }
}

// Remove a degree-3 vertex: three fan faces merge into one triangle.
static void remove_degree3(DelaunayTriangulation& D, int v) {
  int h0 = D.v_out[v], h1 = D.ccw(h0), h2 = D.ccw(h1);

  // Collect face origins.
  vector<int> all_origins;
  if (D.origin_tracker) {
    for (int h : {h0, h1, h2}) {
      int f = D.he_face[h];
      all_origins.insert(all_origins.end(),
                         D.f_origin[f].begin(), D.f_origin[f].end());
    }
    sort(all_origins.begin(), all_origins.end());
    all_origins.erase(unique(all_origins.begin(), all_origins.end()),
                      all_origins.end());
  }

  // Inner rim half-edges connect consecutive neighbors.
  int inner0 = D.he_next[h0], inner1 = D.he_next[h1], inner2 = D.he_next[h2];

  // Fix v_out for neighbors before deallocating.
  for (int h : {h0, h1, h2}) {
    int nbr = D.dest(h);
    int spoke_twin = h ^ 1;
    if (D.v_out[nbr] == spoke_twin) {
      int h_s = D.cw(spoke_twin);
      while (D.dest(h_s) == v && h_s != spoke_twin) h_s = D.cw(h_s);
      D.v_out[nbr] = h_s;
    }
  }

  // Deallocate fan faces and spokes.
  D.dealloc_face(D.he_face[h0]);
  D.dealloc_face(D.he_face[h1]);
  D.dealloc_face(D.he_face[h2]);
  D.dealloc_edge(h0);
  D.dealloc_edge(h1);
  D.dealloc_edge(h2);
  D.v_out[v] = -1;

  // Wire the three inner rim half-edges into a single triangle face.
  D.he_next[inner0] = inner1;
  D.he_next[inner1] = inner2;
  D.he_next[inner2] = inner0;

  int fid = D.alloc_face();
  D.he_face[inner0] = fid;
  D.he_face[inner1] = fid;
  D.he_face[inner2] = fid;
  D.f_he[fid] = inner0;
  if (D.origin_tracker)
    D.f_origin[fid] = all_origins;

  D.recompute_face_angles(fid);
}

// ============================================================================
// Vertex removal: main entry point
// ============================================================================

void DelaunayTriangulation::remove_flat_vertex(int v)
{
  int deg = vertex_degree(v);
  if (deg == 0) return;

  // Phase 1: Reduce degree by flipping incident edges.
  reduce_star_degree(*this, v, 4);
  deg = vertex_degree(v);

  if (deg >= 4) {
    // Phase 2: Ear clipping + DCEL surgery.
    //   extract_fan:  read star geometry -> FanPolygon (isometric 2D embedding)
    //   ear_clip_fan: triangulate the fan polygon (with blocker-flips if stuck)
    //   splice_fan:   replace the star with the triangulation (DCEL surgery)
    FanPolygon fan = extract_fan(*this, v);
    FanTriangulation tri = ear_clip_fan(*this, fan);
    if (!tri.complete) return;  // stuck
    splice_fan(*this, v, fan, tri);
  } else if (deg == 3) {
    // Phase 3: Direct removal (three faces merge into one triangle).
    remove_degree3(*this, v);
  }
}

void DelaunayTriangulation::remove_flat_vertices()
{
  // Sort: flat (degree-6) vertices last.
  // We process them in reverse order.
  // The from_triangulation() constructor preserves the Triangulation's vertex
  // order (which was sorted by sort_flat_last in FulleroidDelaunay).

  while (nv > 0 && v_orig_degree[nv - 1] == 6 && v_out[nv - 1] >= 0) {
    int old_nv = nv;
    int target = nv - 1;

    // Flip away self-loops at the target vertex before removal.
    // Self-loops arise from self-loop ear diagonals in previous removals;
    // they can't be handled by remove_flat_vertex because the inner rim
    // half-edge at the self-loop position has origin = target (becomes dead).
    if (v_out[target] >= 0) {
      int h0 = v_out[target], h = h0;
      bool flipped_any = true;
      while (flipped_any) {
        flipped_any = false;
        h0 = v_out[target];
        if (h0 < 0) break;
        h = h0;
        do {
          if (dest(h) == target) {
            // Self-loop: try to flip it away.
            if (flip_edge(h)) {
              flipped_any = true;
              break;  // Restart scan (circulation changed).
            }
            // Can't flip — try the twin.
            if (flip_edge(h ^ 1)) {
              flipped_any = true;
              break;
            }
          }
          h = cw(h);
        } while (h != h0);
      }
    }

    remove_flat_vertex(target);

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

DelaunayTriangulation DelaunayTriangulation::compute(
    const Triangulation& T, bool track_origins)
{
  // Sort flat vertices last, then build DCEL and run the algorithm.
  Triangulation sorted = T.sort_flat_last();
  DelaunayTriangulation D = from_triangulation(sorted);
  if (track_origins)
    D.origin_tracker = std::make_shared<OriginTracker>(sorted, D.nf);
  D.remove_flat_vertices();
  return D;
}

// --- Symmetric iDT ---

int DelaunayTriangulation::check_edge_symmetry(const vector<vector<int>>& cone_perms) const
{
  set<pair<int,int>> edges;
  for (int h = 0; h < nh; h += 2) {
    if (!alive(h)) continue;
    int u = he_origin[h], v = dest(h);
    edges.insert({min(u,v), max(u,v)});
  }
  int violations = 0;
  for (auto& perm : cone_perms)
    for (auto& e : edges) {
      int pu = perm[e.first], pv = perm[e.second];
      if (!edges.count({min(pu,pv), max(pu,pv)})) violations++;
    }
  return violations;
}

// Find half-edge from u to v, or -1.
static int find_half_edge(const DelaunayTriangulation& D, int u, int v) {
  if (D.v_out[u] < 0) return -1;
  int h0 = D.v_out[u], h = h0;
  do { if (D.dest(h) == v) return h; h = D.cw(h); } while (h != h0);
  return -1;
}

// Insert a Steiner vertex at the center of a co-circular quad, splitting its
// 2 triangles into 4.  The new vertex has zero curvature (flat).
// h_diag: half-edge of the current diagonal of the quad.
// Returns the index of the new vertex.
static int insert_steiner_at_quad(DelaunayTriangulation& D, int h_diag) {
  // The diamond around h_diag:
  //      B
  //     / \        upper face: h_diag(u->v), h1(v->B), h2(B->u)
  //    /   \       lower face: h_twin(v->u), h4(u->D), h5(D->v)
  //   u-----v
  //    \   /
  //     \ /
  //      D
  int h = h_diag, t = h ^ 1;
  int u = D.he_origin[h], v = D.dest(h);
  int h1 = D.he_next[h], h2 = D.he_next[h1];
  int h4 = D.he_next[t], h5 = D.he_next[h4];
  int B = D.he_origin[h2];  // = dest(h1)
  int Dv = D.he_origin[h5]; // = dest(h4)

  // Edge lengths of the quad.
  double uB = D.he_length[h2], Bv = D.he_length[h1];
  double vD = D.he_length[h5], Du = D.he_length[h4];
  double uv = D.he_length[h];  // current diagonal

  // Compute the alternative diagonal length (B-D).
  // For a co-circular quad: use the law of cosines on the two triangles.
  // In triangle (u,v,B): angle at B
  double cos_B = (uB*uB + Bv*Bv - uv*uv) / (2*uB*Bv);
  // In triangle (u,v,D): angle at D
  double cos_D = (Du*Du + vD*vD - uv*uv) / (2*Du*vD);
  // For co-circular: angle_B + angle_D = pi, so cos_D = -cos_B (verified by cot-sum=0).
  // Diagonal BD via triangle (B,u,D): angle at u = angle_Bux + angle_xuD
  // Simpler: use the 2D unfolding.
  // Place u at origin, v along x-axis.
  double ux = 0, uy = 0, vx = uv, vy = 0;
  double ang_uB = acos(clamp((uv*uv + uB*uB - Bv*Bv)/(2*uv*uB), -1.0, 1.0));
  double Bx = uB * cos(ang_uB), By = uB * sin(ang_uB);
  double ang_uD = acos(clamp((uv*uv + Du*Du - vD*vD)/(2*uv*Du), -1.0, 1.0));
  double Dx = Du * cos(-ang_uD), Dy = Du * sin(-ang_uD);  // below the u-v line

  // Intersection of diagonals u-v and B-D.
  // Line u-v: parametric (t*vx, 0) for t in [0,1].
  // Line B-D: parametric B + s*(D-B) for s in [0,1].
  // Solve: t*vx = Bx + s*(Dx-Bx), 0 = By + s*(Dy-By).
  double s = -By / (Dy - By);
  double cx = Bx + s * (Dx - Bx);
  double cy = By + s * (Dy - By);

  // Distances from center to quad vertices.
  double cu = sqrt(cx*cx + cy*cy);
  double cv = sqrt((cx-vx)*(cx-vx) + cy*cy);
  double cB = sqrt((cx-Bx)*(cx-Bx) + (cy-By)*(cy-By));
  double cD = sqrt((cx-Dx)*(cx-Dx) + (cy-Dy)*(cy-Dy));

  // Allocate the new Steiner vertex.
  int sv = D.nv++;
  if (sv >= (int)D.v_out.size()) {
    D.v_out.push_back(-1);
    D.v_cone_angle.push_back(2.0 * M_PI);  // flat (zero curvature)
    D.v_orig_degree.push_back(4);  // degree 4 in the new triangulation
  } else {
    D.v_out[sv] = -1;
    D.v_cone_angle[sv] = 2.0 * M_PI;
    D.v_orig_degree[sv] = 4;
  }

  // Delete the diagonal edge and both faces.
  int fu = D.he_face[h], fl = D.he_face[t];
  D.dealloc_face(fu);
  D.dealloc_face(fl);
  D.dealloc_edge(h);

  // Create 4 new spoke edges: sv-u, sv-v, sv-B, sv-D.
  int su = D.alloc_edge();  // su: sv->u, su^1: u->sv
  int svv = D.alloc_edge(); // svv: sv->v, svv^1: v->sv
  int sB = D.alloc_edge();  // sB: sv->B, sB^1: B->sv
  int sD = D.alloc_edge();  // sD: sv->D, sD^1: D->sv

  D.he_origin[su] = sv;  D.he_origin[su^1] = u;
  D.he_origin[svv] = sv; D.he_origin[svv^1] = v;
  D.he_origin[sB] = sv;  D.he_origin[sB^1] = B;
  D.he_origin[sD] = sv;  D.he_origin[sD^1] = Dv;

  D.he_length[su] = D.he_length[su^1] = cu;
  D.he_length[svv] = D.he_length[svv^1] = cv;
  D.he_length[sB] = D.he_length[sB^1] = cB;
  D.he_length[sD] = D.he_length[sD^1] = cD;

  D.v_out[sv] = su;

  // Create 4 faces using the INNER rim half-edges (orphaned from deleted faces).
  // Inner rim directions: h2(B->u), h1(v->B), h5(D->v), h4(u->D).
  // Faces (CCW, matching inner rim direction):
  //   Face A: (sv, B, u) — arcs sB(sv->B), h2(B->u), su^1(u->sv)
  //   Face B: (sv, v, B) — arcs svv(sv->v), h1(v->B), sB^1(B->sv)
  //   Face C: (sv, D, v) — arcs sD(sv->D), h5(D->v), svv^1(v->sv)
  //   Face D: (sv, u, D) — arcs su(sv->u), h4(u->D), sD^1(D->sv)
  auto make_face = [&](int a, int b, int c, double La, double Lb, double Lc) {
    int fid = D.alloc_face();
    D.he_next[a] = b; D.he_next[b] = c; D.he_next[c] = a;
    D.he_face[a] = D.he_face[b] = D.he_face[c] = fid;
    D.f_he[fid] = a;
    auto ang = [](double s1, double s2, double opp) {
      return acos(clamp((s1*s1+s2*s2-opp*opp)/(2*s1*s2), -1.0, 1.0));
    };
    D.he_angle[a] = ang(Lc, La, Lb);  // angle at origin(a)
    D.he_angle[b] = ang(La, Lb, Lc);  // angle at origin(b)
    D.he_angle[c] = ang(Lb, Lc, La);  // angle at origin(c)
  };

  make_face(sB, h2, su^1,  cB, uB, cu);   // Face A: sv->B->u->sv
  make_face(svv, h1, sB^1,  cv, Bv, cB);   // Face B: sv->v->B->sv
  make_face(sD, h5, svv^1,  cD, vD, cv);   // Face C: sv->D->v->sv
  make_face(su, h4, sD^1,  cu, Du, cD);   // Face D: sv->u->D->sv

  // Fix v_out for quad vertices if they pointed to the deleted diagonal.
  if (D.v_out[u] == h) D.v_out[u] = su ^ 1;   // u->sv
  if (D.v_out[v] == t) D.v_out[v] = svv ^ 1;  // v->sv
  // B and Dv's v_out can't point to the diagonal (they're not endpoints of it).

  return sv;
}

// Compute cot-sum for an edge.
static double edge_cot_sum(const DelaunayTriangulation& D, int h) {
  Diamond dm = D.diamond(h);
  auto cot_opp = [](double e, double a, double b) {
    double cos_C = (a*a + b*b - e*e) / (2*a*b);
    double sin_C = sqrt(max(0.0, 1.0 - cos_C*cos_C));
    return sin_C > 1e-15 ? cos_C / sin_C : 1e15;
  };
  return cot_opp(dm.e, dm.a, dm.b) + cot_opp(dm.e, dm.c, dm.d);
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

  // Compute |w-B| via Stewart's theorem or law of cosines.
  // In triangle (u, v, B) with sides L (u-v), a (B-u), b (v-B):
  // Angle at u: cos(A_u) = (L² + a² - b²) / (2*L*a)
  // |w-B|² = (L/2)² + a² - 2*(L/2)*a*cos(A_u)
  //        = L²/4 + a² - a*(L² + a² - b²)/L   (simplified from Stewart's)
  //        = (a² + b²)/2 - L²/4                  (Stewart's median formula... not quite)
  // Actually, Stewart's theorem for a cevian of length m from B to midpoint w of uv:
  //   b²*(L/2) + a²*(L/2) - m²*L = (L/2)*(L/2)*L  ... no.
  // Stewart's: if cevian from B to point w on uv at distance d_uw from u, d_wv from v:
  //   b²*d_uw + a²*d_wv - m²*L = d_uw * d_wv * L
  // For midpoint: d_uw = d_wv = L/2:
  //   b²*(L/2) + a²*(L/2) - m²*L = (L/2)*(L/2)*L
  //   (a² + b²)/2 - m² = L²/4
  //   m² = (a² + b²)/2 - L²/4
  double wB2 = (a*a + b*b) / 2.0 - L*L / 4.0;
  double wB = wB2 > 0 ? sqrt(wB2) : 0;

  // Similarly for |w-D|.
  double wD2 = (c*c + d*d) / 2.0 - L*L / 4.0;
  double wD = wD2 > 0 ? sqrt(wD2) : 0;

  // Allocate new vertex w.
  int w = D.nv++;
  if (w >= (int)D.v_out.size()) {
    D.v_out.push_back(-1);
    D.v_cone_angle.push_back(2.0 * M_PI);  // flat
    D.v_orig_degree.push_back(4);
  } else {
    D.v_out[w] = -1;
    D.v_cone_angle[w] = 2.0 * M_PI;
    D.v_orig_degree[w] = 4;
  }

  // Delete original edge and its two faces.
  int fu = D.he_face[h], ft = D.he_face[t];
  D.dealloc_face(fu);
  D.dealloc_face(ft);
  D.dealloc_edge(h);

  // Allocate new edges: u-w, w-v, w-B, w-D.
  int uw = D.alloc_edge();  // uw: u->w, uw^1: w->u
  int wv = D.alloc_edge();  // wv: w->v, wv^1: v->w
  int wB_he = D.alloc_edge();  // wB: w->B, wB^1: B->w
  int wD_he = D.alloc_edge();  // wD: w->D, wD^1: D->w

  D.he_origin[uw] = u;    D.he_origin[uw^1] = w;
  D.he_origin[wv] = w;    D.he_origin[wv^1] = v;
  D.he_origin[wB_he] = w; D.he_origin[wB_he^1] = B;
  D.he_origin[wD_he] = w; D.he_origin[wD_he^1] = Dv;

  D.he_length[uw] = D.he_length[uw^1] = half;
  D.he_length[wv] = D.he_length[wv^1] = half;
  D.he_length[wB_he] = D.he_length[wB_he^1] = wB;
  D.he_length[wD_he] = D.he_length[wD_he^1] = wD;

  D.v_out[w] = uw ^ 1;  // w -> u

  // Create 4 faces using the surviving rim half-edges (h_vB, h_Bu, h_uD, h_Dv)
  // and the new edges.
  //
  // Face 1: (u, w, B) — arcs: uw(u->w), wB(w->B), h_Bu(B->u)
  // Face 2: (w, v, B) — arcs: wv(w->v), h_vB(v->B), wB^1(B->w)
  // Face 3: (v, w, D) — arcs: wv^1(v->w), wD(w->D), h_Dv(D->v)  ... wait.
  //
  // The original face of t was: v -> u -> D -> v (half-edges t, h_uD, h_Dv).
  // After removing t (u->v side), the surviving rim is h_uD (u->D) and h_Dv (D->v).
  // But we split u-v into u-w and w-v. So:
  // Face 3: (w, u, D) — arcs: uw^1(w->u), h_uD(u->D), wD^1(D->w)
  // Face 4: (v, w, D) — arcs: wv^1(v->w), wD(w->D)... no, D->v is h_Dv.
  //
  // Let me be careful with orientation. The original faces were CCW:
  // Face of h: u -> v -> B -> u  (CCW)
  // Face of t: v -> u -> D -> v  (CCW)
  //
  // After bisection, the 4 new faces (all CCW):
  // Face 1: u -> w -> B -> u    uses uw, wB_he, h_Bu
  // Face 2: w -> v -> B -> w    uses wv, h_vB, wB_he^1
  // Face 3: v -> w -> D -> v    uses wv^1, wD_he, h_Dv
  // Face 4: w -> u -> D -> w    uses uw^1, h_uD, wD_he^1

  auto safe_angle = [](double s1, double s2, double opp) {
    return acos(std::clamp((s1*s1 + s2*s2 - opp*opp) / (2*s1*s2), -1.0, 1.0));
  };

  auto make_face = [&](int e1, int e2, int e3, double L1, double L2, double L3) {
    int fid = D.alloc_face();
    D.he_next[e1] = e2; D.he_next[e2] = e3; D.he_next[e3] = e1;
    D.he_face[e1] = D.he_face[e2] = D.he_face[e3] = fid;
    D.f_he[fid] = e1;
    D.he_angle[e1] = safe_angle(L3, L1, L2);
    D.he_angle[e2] = safe_angle(L1, L2, L3);
    D.he_angle[e3] = safe_angle(L2, L3, L1);
  };

  // Face 1: u -> w -> B -> u  (sides: uw=half, wB, Bu=a)
  make_face(uw, wB_he, h_Bu,  half, wB, a);
  // Face 2: w -> v -> B -> w  (sides: wv=half, vB=b, Bw=wB)
  make_face(wv, h_vB, wB_he^1,  half, b, wB);
  // Face 3: v -> w -> D -> v  (sides: vw=half, wD, Dv=d)
  make_face(wv^1, wD_he, h_Dv,  half, wD, d);
  // Face 4: w -> u -> D -> w  (sides: wu=half, uD=c, Dw=wD)
  make_face(uw^1, h_uD, wD_he^1,  half, c, wD);

  // Fix v_out for u and v if they pointed to the deleted edge.
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

DelaunayTriangulation DelaunayTriangulation::compute_symmetric(
    const Triangulation& T, const Symmetry& S)
{
  // Symmetric iDT via Steiner vertex insertion at self-dual co-circular quads.
  //
  // 1. Compute non-symmetric iDT, correct edge lengths to geodesic distances.
  // 2. Identify self-dual co-circular edge orbits (quads with strictly Delaunay
  //    rims where neither diagonal choice is G-invariant).
  // 3. Insert a Steiner vertex at each quad center, splitting 2 triangles -> 4.
  //    The Steiner vertices are flat (zero curvature) and related by the group
  //    action, so the result is G-invariant by construction.

  // Step 1: Non-symmetric iDT with geodesic distances.
  DelaunayTriangulation D = compute(T);
  if (D.nv != 12) return D;

  vector<vector<int>> G_full(S.G.begin(), S.G.end());
  auto cone_perms = restrict_to_cone_points(G_full, T);
  if (cone_perms.empty()) return D;

  // Correct edge lengths.
  vector<int> cone_pts;
  for (int v = 0; v < T.N; v++)
    if (T.degree(v) != 6) cone_pts.push_back(v);
  auto dist2 = T.surface_distances(cone_pts);
  for (int h = 0; h < D.nh; h += 2) {
    if (!D.alive(h)) continue;
    D.he_length[h] = D.he_length[h^1] = sqrt(dist2(D.he_origin[h], D.dest(h)));
  }
  D.recompute_all_angles();
  D.flip_to_delaunay();

  if (D.check_edge_symmetry(cone_perms) == 0) return D;

  // Step 2: Find self-dual co-circular quads.
  // A co-circular edge is in a self-dual orbit if its orbit image is not present.
  set<pair<int,int>> edge_set;
  for (int h = 0; h < D.nh; h += 2)
    if (D.alive(h))
      edge_set.insert({min(D.he_origin[h],D.dest(h)), max(D.he_origin[h],D.dest(h))});

  // Collect co-circular diagonal half-edges that have missing orbit images.
  // Deduplicate by edge ID (each quad has one diagonal).
  set<int> quads_to_split;  // edge IDs of co-circular diagonals to split
  for (auto& perm : cone_perms)
    for (auto& e : edge_set) {
      auto pe = make_pair(min(perm[e.first],perm[e.second]),
                          max(perm[e.first],perm[e.second]));
      if (edge_set.count(pe)) continue;
      // e is present but pe is missing -> e is in a self-dual orbit.
      // Find the half-edge for e and check it's co-circular.
      for (int h = 0; h < D.nh; h += 2) {
        if (!D.alive(h)) continue;
        int u = D.he_origin[h], v = D.dest(h);
        if (make_pair(min(u,v),max(u,v)) == e) {
          if (fabs(edge_cot_sum(D, h)) < 1e-8)
            quads_to_split.insert(h >> 1);
          break;
        }
      }
    }

  if (quads_to_split.empty()) return D;

  // Step 3: Insert Steiner vertices at quad centers.
  int n_inserted = 0;
  for (int eid : quads_to_split) {
    int h = eid * 2;
    if (!D.alive(h)) continue;
    insert_steiner_at_quad(D, h);
    n_inserted++;
  }

  if (!D.check_consistency()) {
    fprintf(stderr, "  WARNING: DCEL inconsistent after %d Steiner insertions\n", n_inserted);
    return compute(T);  // fallback
  }

  // Build extended permutations covering Steiner vertices.
  vector<set<int>> steiner_nbs;
  for (int sv = 12; sv < D.nv; sv++) {
    set<int> nbs;
    int h0 = D.v_out[sv], h = h0;
    if (h0 >= 0) do { nbs.insert(D.dest(h)); h = D.cw(h); } while (h != h0);
    steiner_nbs.push_back(nbs);
  }

  auto build_extended_perms = [&]() {
    // Recompute Steiner neighbor sets (they change after flips).
    steiner_nbs.clear();
    for (int sv = 12; sv < D.nv; sv++) {
      set<int> nbs;
      int h0 = D.v_out[sv], h = h0;
      if (h0 >= 0) do { nbs.insert(D.dest(h)); h = D.cw(h); } while (h != h0);
      steiner_nbs.push_back(nbs);
    }
    vector<vector<int>> ext_perms;
    for (auto& cp : cone_perms) {
      vector<int> fp(D.nv);
      for (int v = 0; v < 12; v++) fp[v] = cp[v];
      for (int si = 0; si < (int)steiner_nbs.size(); si++) {
        set<int> mapped;
        for (int v : steiner_nbs[si]) if (v < 12) mapped.insert(cp[v]);
        int mapped_sv = 12 + si;
        for (int sj = 0; sj < (int)steiner_nbs.size(); sj++) {
          set<int> sj_cone;
          for (int v : steiner_nbs[sj]) if (v < 12) sj_cone.insert(v);
          if (sj_cone == mapped) { mapped_sv = 12 + sj; break; }
        }
        fp[12 + si] = mapped_sv;
      }
      ext_perms.push_back(fp);
    }
    return ext_perms;
  };

  // Orbit-aware Lawson: flip non-Delaunay edges, but flip ALL orbit members
  // simultaneously to preserve G-invariance.
  for (int pass = 0; pass < 200; pass++) {
    auto ext_perms = build_extended_perms();

    // Find a non-Delaunay edge.
    int flip_h = -1;
    for (int h = 0; h < D.nh; h += 2) {
      if (!D.alive(h)) continue;
      if (!D.diamond(h).is_delaunay() && D.diamond(h).is_convex()) {
        flip_h = h; break;
      }
    }
    if (flip_h < 0) break;  // all Delaunay

    // Find ALL orbit members of this edge and flip them all.
    int fu = D.he_origin[flip_h], fv = D.dest(flip_h);
    set<int> to_flip;
    to_flip.insert(flip_h);
    for (auto& perm : ext_perms) {
      int pu = perm[fu], pv = perm[fv];
      // Find the half-edge from pu to pv.
      if (D.v_out[pu] >= 0) {
        int h0 = D.v_out[pu], h = h0;
        do {
          if (D.dest(h) == pv) { to_flip.insert(h & ~1); break; }
          h = D.cw(h);
        } while (h != h0);
      }
    }

    for (int h : to_flip) {
      if (D.alive(h) && !D.diamond(h).is_delaunay() && D.diamond(h).is_convex())
        D.flip_edge(h);
    }
  }

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
// Symmetry utilities
// ============================================================================

vector<vector<int>> restrict_to_cone_points(
    const vector<vector<int>>& G, const Triangulation& T)
{
  if (G.empty()) return {};

  // Cone points: vertices with degree != 6, in original index order.
  // sort_flat_last() places these first in the same relative order,
  // so iDT vertex i corresponds to the i-th cone point by original index.
  int N = T.N;
  vector<int> cone;
  vector<int> orig_to_idt(N, -1);
  for (int v = 0; v < N; v++)
    if (T.degree(v) != 6) {
      orig_to_idt[v] = cone.size();
      cone.push_back(v);
    }
  int nc = cone.size();

  vector<vector<int>> result;
  for (auto& pi : G) {
    vector<int> r(nc);
    bool valid = true;
    for (int i = 0; i < nc; i++) {
      int img = orig_to_idt[pi[cone[i]]];
      if (img < 0) { valid = false; break; }
      r[i] = img;
    }
    if (valid) result.push_back(r);
  }

  // Remove duplicate permutations
  sort(result.begin(), result.end());
  result.erase(unique(result.begin(), result.end()), result.end());
  return result;
}
