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

// Old FulleroidDelaunay + IDTAudit implementation moved to delaunay_old.cc.

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

    // Collect face origins from the fan (both combined and per-sector).
    // Faces on the boundary between iDT faces (from classify_across_line's
    // "both sides" case) may appear in multiple sectors.  Deduplicate so
    // each face appears in exactly one sector (the first one that claims it).
    vector<int> all_origins;
    vector<vector<int>> sector_origins(k);
    { int h0 = v_out[v], h = h0;
      int sec = 0;
      unordered_set<int> seen;
      do {
        int f = he_face[h];
        for (int fid : f_origin[f]) {
          all_origins.push_back(fid);
          if (seen.insert(fid).second)
            sector_origins[sec].push_back(fid);
        }
        h = ccw(h);
        sec++;
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

    // Compute per-triangle origin assignment.
    vector<vector<int>> origin_assignment;
    if (origin_tracker) {
      // Per-sector local unfolding + barycentric mapping to fan coordinates.
      //
      // Global Eisenstein unfolding fails for vertex removal because the fan
      // can contain cone points, and BFS wraps around them giving inconsistent
      // coordinates.  Instead, unfold each sector (iDT fan face) separately
      // and map each original face's centroid to fan 2D coordinates via
      // barycentric interpolation within the iDT triangle.
      //
      // Fan 2D positions of boundary vertices:
      vector<double> fan_x(k), fan_y(k);
      for (int i = 0; i < k; i++) {
        fan_x[i] = spokes[i] * cos(cum[i]);
        fan_y[i] = spokes[i] * sin(cum[i]);
      }

      // Ear triangle vertices in fan 2D.
      int nt = tris.size();
      vector<array<double,2>> ear_v0(nt), ear_v1(nt), ear_v2(nt);
      for (int ti = 0; ti < nt; ti++) {
        ear_v0[ti] = {fan_x[tris[ti].v0], fan_y[tris[ti].v0]};
        ear_v1[ti] = {fan_x[tris[ti].v1], fan_y[tris[ti].v1]};
        ear_v2[ti] = {fan_x[tris[ti].v2], fan_y[tris[ti].v2]};
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
        return -1;  // not found (shouldn't happen)
      };

      origin_assignment.resize(nt);

      // Process each sector separately.
      for (int sec = 0; sec < k; sec++) {
        auto& sec_origins = sector_origins[sec];

        if (!sec_origins.empty()) {
          // Local Eisenstein unfolding of this sector's original faces.
          auto pos = origin_tracker->unfold_patch(sec_origins, {v, nb[sec], nb[(sec+1)%k]}, v);

          // Get Eisenstein positions of sector triangle vertices.
          auto it_v = pos.find(v);
          auto it_a = pos.find(nb[sec]);
          auto it_b = pos.find(nb[(sec+1)%k]);

          if (it_v != pos.end() && it_a != pos.end() && it_b != pos.end()) {
            // Compute the affine transform from Eisenstein to fan 2D.
            // Eisenstein triangle: V, A, B (as doubles).
            double eVx = it_v->second.first, eVy = it_v->second.second;
            double eAx = it_a->second.first, eAy = it_a->second.second;
            double eBx = it_b->second.first, eBy = it_b->second.second;

            // Fan 2D triangle: (0,0), (fan_x[sec], fan_y[sec]), (fan_x[(sec+1)%k], fan_y[(sec+1)%k]).
            double fAx = fan_x[sec], fAy = fan_y[sec];
            double fBx = fan_x[(sec+1)%k], fBy = fan_y[(sec+1)%k];

            // Determinant for barycentric coords in Eisenstein triangle.
            double det = (eAx - eVx) * (eBy - eVy) - (eBx - eVx) * (eAy - eVy);

            for (int fid : sec_origins) {
              auto& fverts = origin_tracker->face_verts[fid];
              // Centroid in Eisenstein coordinates.
              double cx = 0, cy = 0;
              bool all_found = true;
              for (int vi = 0; vi < 3; vi++) {
                auto it = pos.find(fverts[vi]);
                if (it == pos.end()) { all_found = false; break; }
                cx += it->second.first;
                cy += it->second.second;
              }
              if (!all_found) {
                // Fallback: assign to ear triangle containing the sector's rim edge.
                for (int ti = 0; ti < nt; ti++) {
                  auto& t = tris[ti];
                  if ((t.v0 == sec || t.v1 == sec || t.v2 == sec) &&
                      (t.v0 == (sec+1)%k || t.v1 == (sec+1)%k || t.v2 == (sec+1)%k)) {
                    origin_assignment[ti].push_back(fid);
                    break;
                  }
                }
                continue;
              }
              cx /= 3.0; cy /= 3.0;

              // Barycentric coords in Eisenstein triangle.
              if (std::abs(det) < 1e-15) {
                // Degenerate sector — assign to first ear triangle touching this sector.
                for (int ti = 0; ti < nt; ti++) {
                  auto& t = tris[ti];
                  if ((t.v0 == sec || t.v1 == sec || t.v2 == sec) &&
                      (t.v0 == (sec+1)%k || t.v1 == (sec+1)%k || t.v2 == (sec+1)%k)) {
                    origin_assignment[ti].push_back(fid);
                    break;
                  }
                }
                continue;
              }
              double bary_a = ((cx - eVx) * (eBy - eVy) - (cy - eVy) * (eBx - eVx)) / det;
              double bary_b = ((cy - eVy) * (eAx - eVx) - (cx - eVx) * (eAy - eVy)) / det;

              // Map to fan 2D.
              double fan_px = bary_a * fAx + bary_b * fBx;
              double fan_py = bary_a * fAy + bary_b * fBy;

              int ti = classify_point(fan_px, fan_py);
              if (ti >= 0) {
                origin_assignment[ti].push_back(fid);
              } else {
                // Fallback: assign to nearest ear triangle.
                for (int j = 0; j < nt; j++) {
                  auto& t = tris[j];
                  if ((t.v0 == sec || t.v1 == sec || t.v2 == sec) &&
                      (t.v0 == (sec+1)%k || t.v1 == (sec+1)%k || t.v2 == (sec+1)%k)) {
                    origin_assignment[j].push_back(fid);
                    break;
                  }
                }
              }
            }
          } else {
            // Couldn't unfold sector — fallback to rim assignment.
            for (int fid : sec_origins) {
              for (int ti = 0; ti < nt; ti++) {
                auto& t = tris[ti];
                if ((t.v0 == sec || t.v1 == sec || t.v2 == sec) &&
                    (t.v0 == (sec+1)%k || t.v1 == (sec+1)%k || t.v2 == (sec+1)%k)) {
                  origin_assignment[ti].push_back(fid);
                  break;
                }
              }
            }
          }
        }

      }
    } else {
      // Approximate: distribute origins evenly across triangles.
      origin_assignment.resize(tris.size());
      int per = all_origins.size() / tris.size();
      int idx = 0;
      for (size_t i = 0; i < tris.size(); i++) {
        int n = (i + 1 < tris.size()) ? per
                                       : (int)all_origins.size() - idx;
        n = std::max(n, 0);
        origin_assignment[i].assign(all_origins.begin() + idx,
                                    all_origins.begin() + idx + n);
        idx += n;
      }
    }

    // Wire up each triangle.
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

      f_origin[fid] = std::move(origin_assignment[ti]);
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

DelaunayTriangulation DelaunayTriangulation::compute(
    const Triangulation& T, bool exact_origins)
{
  // Sort flat vertices last, then build DCEL and run the algorithm.
  Triangulation sorted = T.sort_flat_last();
  DelaunayTriangulation D = from_triangulation(sorted);

  if (exact_origins)
    D.origin_tracker = std::make_shared<OriginTracker>(sorted, D.nf);

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
// Symmetry utilities for embed_3d
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

// Compute vertex orbits from a group of permutations (union-find).
vector<vector<int>> compute_orbits(int n, const vector<vector<int>>& G) {
  vector<int> parent(n);
  iota(parent.begin(), parent.end(), 0);
  std::function<int(int)> find = [&](int x) {
    return parent[x] == x ? x : parent[x] = find(parent[x]);
  };
  for (auto& pi : G)
    for (int v = 0; v < n; v++) {
      int a = find(v), b = find(pi[v]);
      if (a != b) parent[a] = b;
    }

  map<int, vector<int>> m;
  for (int v = 0; v < n; v++) m[find(v)].push_back(v);
  vector<vector<int>> orbits;
  for (auto& [_, orbit] : m) orbits.push_back(orbit);
  return orbits;
}

// 3x3 Jacobi eigendecomposition of symmetric matrix A.
// On return: A is diagonal (eigenvalues), V columns are eigenvectors.
static void jacobi_eigen_3x3(matrix3d& A, matrix3d& V) {
  V = matrix3d::unit_matrix();
  for (int iter = 0; iter < 100; iter++) {
    double mx = 0; int p = 0, q = 1;
    for (int i = 0; i < 3; i++)
      for (int j = i+1; j < 3; j++)
        if (fabs(A(i,j)) > mx) { mx = fabs(A(i,j)); p = i; q = j; }
    if (mx < 1e-15) break;

    double app = A(p,p), aqq = A(q,q), apq = A(p,q);
    double theta = (fabs(app-aqq) < 1e-30) ? M_PI/4 : 0.5*atan2(2*apq, app-aqq);
    double cs = cos(theta), sn = sin(theta);

    for (int j = 0; j < 3; j++) {
      double bp = A(p,j), bq = A(q,j);
      A(p,j) = cs*bp + sn*bq; A(q,j) = -sn*bp + cs*bq;
    }
    for (int i = 0; i < 3; i++) {
      double bp = A(i,p), bq = A(i,q);
      A(i,p) = cs*bp + sn*bq; A(i,q) = -sn*bp + cs*bq;
    }
    for (int i = 0; i < 3; i++) {
      double vp = V(i,p), vq = V(i,q);
      V(i,p) = cs*vp + sn*vq; V(i,q) = -sn*vp + cs*vq;
    }
  }
}

// Matrix-vector product for matrix3d * coord3d.
static coord3d matvec(const matrix3d& M, const coord3d& v) {
  return coord3d(M(0,0)*v[0]+M(0,1)*v[1]+M(0,2)*v[2],
                 M(1,0)*v[0]+M(1,1)*v[1]+M(1,2)*v[2],
                 M(2,0)*v[0]+M(2,1)*v[1]+M(2,2)*v[2]);
}

// Procrustes rotation: find the orthogonal R minimizing sum |R*src[i] - dst[i]|^2.
// Allows improper rotations (reflections) for orientation-reversing automorphisms.
static matrix3d procrustes_rotation(const vector<coord3d>& src,
                                    const vector<coord3d>& dst) {
  // Cross-covariance M = sum dst[i] src[i]^T
  matrix3d M;
  for (size_t i = 0; i < src.size(); i++)
    M += dst[i].outer(src[i]);

  // SVD via eigendecomposition of M^T M
  matrix3d MtM = M.transpose() * M;
  matrix3d V;
  jacobi_eigen_3x3(MtM, V);

  // U = M V Sigma^{-1}, processing columns in decreasing singular value order
  double svals[3] = { sqrt(max(0.0, MtM(0,0))),
                      sqrt(max(0.0, MtM(1,1))),
                      sqrt(max(0.0, MtM(2,2))) };

  matrix3d U;
  int computed = 0;
  for (int col = 0; col < 3; col++) {
    if (svals[col] > 1e-12) {
      coord3d v_col(V(0,col), V(1,col), V(2,col));
      coord3d u_col = matvec(M, v_col) * (1.0/svals[col]);
      U(0,col) = u_col[0]; U(1,col) = u_col[1]; U(2,col) = u_col[2];
      computed++;
    }
  }

  // Complete U for zero singular values via cross product
  if (computed == 2) {
    int z = -1;
    for (int i = 0; i < 3; i++) if (svals[i] <= 1e-12) { z = i; break; }
    int a = (z+1)%3, b = (z+2)%3;
    coord3d ua(U(0,a),U(1,a),U(2,a)), ub(U(0,b),U(1,b),U(2,b));
    coord3d uc = ua.cross(ub);
    double n = uc.norm();
    if (n > 1e-15) uc /= n;
    U(0,z) = uc[0]; U(1,z) = uc[1]; U(2,z) = uc[2];
  }

  return U * V.transpose();
}

// Symmetrize coordinates via the Reynolds operator:
//   x_sym[v] = (1/|G|) sum_{pi in G} R_pi^T x[pi(v)]
// where R_pi is the Procrustes rotation for group element pi.
static void symmetrize_embedding(V3& x, const vector<vector<int>>& G) {
  if (G.empty()) return;
  int n = x.size();

  // Center
  coord3d c;
  for (auto& xi : x) c += xi;
  c /= n;
  for (auto& xi : x) xi -= c;

  // Compute Procrustes rotation for each group element
  vector<matrix3d> rotations(G.size());
  for (size_t g = 0; g < G.size(); g++) {
    auto& pi = G[g];
    vector<coord3d> src(n), dst(n);
    for (int v = 0; v < n; v++) {
      src[v] = x[v];
      dst[v] = x[pi[v]];
    }
    rotations[g] = procrustes_rotation(src, dst);
  }

  // Reynolds average: x_sym[v] = (1/|G|) sum_pi R_pi^T x[pi(v)]
  V3 x_sym(n);
  for (size_t g = 0; g < G.size(); g++) {
    auto& pi = G[g];
    matrix3d Rt = rotations[g].transpose();
    for (int v = 0; v < n; v++)
      x_sym[v] += matvec(Rt, x[pi[v]]);
  }
  double inv_g = 1.0 / G.size();
  for (int v = 0; v < n; v++)
    x[v] = x_sym[v] * inv_g + c;
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

vector<coord3d> DelaunayTriangulation::embed_3d(const vector<vector<int>>& sym) const
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

  // Resolve MDS degeneracies from symmetry.  When vertices in the same
  // orbit have identical APSP profiles to all other vertices, the full
  // MDS places them at the same point.  Fix this by running a secondary
  // MDS on the intra-orbit distances: orbit members have different
  // distances to EACH OTHER, so this local MDS is non-degenerate.
  auto orbits = compute_orbits(nv, sym);  // trivial singletons if sym={}
  for (auto& orbit : orbits) {
    if (orbit.size() < 2) continue;

    // Orbit centroid and spread
    coord3d center;
    for (int v : orbit) center += x[v];
    center /= orbit.size();

    double max_spread = 0;
    for (int v : orbit)
      max_spread = max(max_spread, (x[v] - center).norm());

    // Minimum intra-orbit APSP distance
    double min_intra = 1e18;
    for (size_t i = 0; i < orbit.size(); i++)
      for (size_t j = i+1; j < orbit.size(); j++)
        min_intra = min(min_intra, D(orbit[i], orbit[j]));

    if (max_spread > 0.1 * min_intra) continue;  // already separated

    // Collapsed: local MDS on intra-orbit distances
    int m = orbit.size();
    matrix<double> D_local(m, m);
    for (int i = 0; i < m; i++)
      for (int j = 0; j < m; j++)
        D_local(i, j) = D(orbit[i], orbit[j]);

    V3 local = mds_placement(D_local);
    for (int i = 0; i < m; i++)
      x[orbit[i]] = center + local[i];
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

  if (sym.empty()) {
    x = steihaug_cg(std::move(x), eval, hv_prod);
  } else {
    // Symmetry-constrained optimization in two phases:
    //
    // Phase 1: Optimize without symmetry constraint.  The heuristic
    //   perturbation (above) breaks MDS degeneracy for collapsed orbits.
    //   The optimizer converges to a non-degenerate embedding that's
    //   approximately but not exactly symmetric.
    //
    // Phase 2: Symmetrize via Procrustes + Reynolds, then polish.
    //   Now that the vertices are at well-separated positions, the
    //   Procrustes rotation for each group element is well-determined.
    //   Symmetrization projects onto the symmetric subspace, and a
    //   short optimization pass restores edge-length accuracy.
    x = steihaug_cg(std::move(x), eval, hv_prod);
    for (int pass = 0; pass < 5; pass++) {
      V3 x_before(x);
      symmetrize_embedding(x, sym);
      double shift = (x - x_before).norm() / max(x.norm(), 1e-15);
      if (shift < 1e-12) break;
      x = steihaug_cg(std::move(x), eval, hv_prod, 100, 1e-13);
    }
  }

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
