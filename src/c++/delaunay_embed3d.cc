#include "fullerenes/delaunay_embed3d.hh"
#include "fullerenes/union_find.hh"
#include "fullerenes/dense_linalg.hh"

#include <cmath>
#include <cassert>
#include <algorithm>
#include <numeric>
#include <map>
#include <stdexcept>
#include <vector>

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
//
// Collapse-aware eigenvector selection: for non-Euclidean distance matrices
// (geodesic distances), the 3rd and 4th Gram eigenvalues can be near-degenerate.
// When the top-3 eigenvectors all have the same symmetry character (e.g. all
// symmetric under a Cs mirror), symmetry-related vertex pairs collapse to the
// same point.  When edges are provided and collapse is detected, the 3rd
// eigenvector is replaced with the best alternative from lower positive
// eigenvalues.
struct MDSEdge { int u, v; double L; };

static V3 mds_placement(const matrix<double>& D,
                         const vector<MDSEdge>& edges = {}) {
  int N = D.m;

  // Squared-distance matrix
  matrix<double> D_sq(N, N);
  for (int i = 0; i < N; i++)
    for (int j = 0; j < N; j++)
      D_sq(i, j) = D(i, j) * D(i, j);

  // Gram-matrix eigendecomposition (BLAS-free cyclic Jacobi).  double_center is
  // symmetric by construction, so the guard trip (empty result) is unreachable.
  LinAlg::SymEigen::Decomp ed = LinAlg::SymEigen::decompose(double_center(D_sq));
  if (ed.lambda.empty())
    throw std::runtime_error("mds_placement: SymEigen::decompose returned no "
                             "spectrum for a symmetric double-centered matrix");
  const vector<double>&  evals = ed.lambda;   // ascending
  const matrix<double>&  V     = ed.Q;        // eigenvectors as columns

  // Sort eigenvalues descending
  vector<int> order(N);
  iota(order.begin(), order.end(), 0);
  sort(order.begin(), order.end(), [&](int a, int b) { return evals[a] > evals[b]; });

  // Collect indices of positive eigenvalues (candidates for MDS axes).
  vector<int> pos_idx;
  for (int k = 0; k < N; k++)
    if (evals[order[k]] > 1e-10) pos_idx.push_back(k);

  // Helpers for collapse-aware eigenvector selection.
  auto build_embedding = [&](int i0, int i1, int i2) -> V3 {
    int cols[3] = {order[i0], order[i1], order[i2]};
    V3 x(N);
    for (int i = 0; i < N; i++)
      for (int dim = 0; dim < 3; dim++)
        x[i][dim] = V(i, cols[dim]) * sqrt(max(0.0, evals[cols[dim]]));
    return x;
  };

  auto edge_stress = [&](const V3& x) -> double {
    double stress = 0;
    for (auto& e : edges) {
      double d = (x[e.u] - x[e.v]).norm();
      double err = d - e.L;
      stress += err * err;
    }
    return stress;
  };

  // Default: top 3 eigenvalues.
  V3 x = build_embedding(0, 1, min(2, N-1));

  if (!edges.empty() && pos_idx.size() > 3) {
    // Check for near-coincident vertex pairs: if any edge (u,v) has
    // MDS distance < 10% of its target, the top-3 eigenvectors missed
    // the component that separates them (typically due to symmetry-related
    // near-degenerate eigenvalues).
    bool has_collapse = false;
    for (auto& e : edges) {
      double d = (x[e.u] - x[e.v]).norm();
      if (d < 0.1 * e.L) { has_collapse = true; break; }
    }

    if (has_collapse) {
      // Keep top-2 eigenvalues (large-scale shape) and search over
      // which eigenvalue to use for the 3rd axis.
      double best_stress = edge_stress(x);
      int K = min((int)pos_idx.size(), 8);
      for (int c = 2; c < K; c++) {
        V3 candidate = build_embedding(pos_idx[0], pos_idx[1], pos_idx[c]);
        double s = edge_stress(candidate);
        if (s < best_stress) { best_stress = s; x = std::move(candidate); }
      }
    }
  }

  // Flatness detection: if 3rd eigenvalue is degenerate, perturb
  double lam3 = max(0.0, evals[order[min(2, N-1)]]);
  double lam1 = max(1e-30, evals[order[0]]);
  if (lam3 / lam1 < 1e-6)
    for (int i = 0; i < N; i++)
      x[i][2] += 0.1 * lam1 * sin(i * 2.0 * M_PI / N);

  return x;
}

// ============================================================================
// Symmetry utilities for embed_3d
// ============================================================================

SymmetryConstraint restrict_symmetry_to_cone_points(
    const vector<vector<int>>& G, const vector<matrix3d>& R,
    const Triangulation& T)
{
  if (G.empty() || G.size() != R.size()) return {};

  int N = T.N;
  vector<int> cone;
  vector<int> orig_to_idt(N, -1);
  for (int v = 0; v < N; v++)
    if (T.degree(v) != 6) {
      orig_to_idt[v] = cone.size();
      cone.push_back(v);
    }
  int nc = cone.size();

  // Build paired (perm, matrix), deduplicating by perm.
  SymmetryConstraint result;
  map<vector<int>, int> seen;  // perm -> index in result

  for (size_t g = 0; g < G.size(); g++) {
    vector<int> r(nc);
    bool valid = true;
    for (int i = 0; i < nc; i++) {
      int img = orig_to_idt[G[g][cone[i]]];
      if (img < 0) { valid = false; break; }
      r[i] = img;
    }
    if (!valid) continue;

    if (seen.find(r) == seen.end()) {
      seen[r] = result.perms.size();
      result.perms.push_back(r);
      result.matrices.push_back(R[g]);
    }
  }
  return result;
}

// Vertex orbits from a group of permutations -- the connected components of
// u ~ pi[u]. Members are ascending, so orbit[0] is each orbit's smallest member
// (the representative callers use); orbit-vector order is an internal relabel.
vector<vector<int>> compute_orbits(int n, const vector<vector<int>>& G) {
  return orbits(n, G);
}

// ============================================================================
// Orbit structure for symmetry-constrained embedding
// ============================================================================

// Per-orbit data for reduced-parameter optimization.
struct OrbitInfo {
  int rep;                    // orbit representative (smallest index)
  vector<int> members;        // all orbit members (including rep)
  vector<int> gen_element;    // gen_element[k]: group index g such that perms[g][rep] == members[k]
  int subspace_dim;           // dimension of stabilizer fixed-point subspace (1, 2, or 3)
  matrix3d basis;             // columns 0..subspace_dim-1 are ON basis for the subspace
};

// Compute the fixed-point subspace of a set of O(3) matrices.
// Returns (dimension, basis) where basis columns 0..dim-1 span the subspace.
static pair<int, matrix3d> fixed_point_subspace(const vector<matrix3d>& stabilizer) {
  // Form A = sum (R_h - I)^T (R_h - I).  Null space of A = fixed-point subspace.
  matrix3d I3 = matrix3d::unit_matrix();
  matrix3d A;
  for (auto& R : stabilizer) {
    matrix3d D = R - I3;
    A += D.transpose() * D;
  }

  // Eigendecompose A (symmetric PSD).  eigensystem() sorts eigenvalues by
  // absolute value ascending -- for a PSD matrix that is value-ascending, so
  // the null-space eigenvectors (eigenvalue ~ 0) come first -- and returns each
  // unit eigenvector as row k of C.  The fixed-point subspace is the null space.
  auto [evals, C] = A.eigensystem();

  int dim = 0;
  matrix3d basis;
  for (int k = 0; k < 3; k++) {
    if (evals[k] < 1e-10) {
      for (int r = 0; r < 3; r++) basis(r, dim) = C(k, r);   // C row k -> basis column
      dim++;
    }
  }

  // Dimension 0 means the origin is the only fixed point -- shouldn't happen
  // for a vertex of a convex polyhedron.  Fall back to full R^3.
  if (dim == 0) {
    dim = 3;
    basis = matrix3d::unit_matrix();
  }

  return {dim, basis};
}

// Build orbit structure from a SymmetryConstraint.
static vector<OrbitInfo> compute_orbit_structure(
    int nv, const SymmetryConstraint& sym)
{
  auto orbits = compute_orbits(nv, sym.perms);
  int n_g = sym.perms.size();

  vector<OrbitInfo> result;
  for (auto& orbit : orbits) {
    OrbitInfo oi;
    oi.rep = orbit[0];  // compute_orbits uses union-find; just use first
    oi.members = orbit;

    // For each member, find a group element mapping rep to it.
    oi.gen_element.resize(orbit.size());
    for (size_t k = 0; k < orbit.size(); k++) {
      oi.gen_element[k] = -1;
      for (int g = 0; g < n_g; g++) {
        if (sym.perms[g][oi.rep] == orbit[k]) {
          oi.gen_element[k] = g;
          break;
        }
      }
      assert(oi.gen_element[k] >= 0);
    }

    // Compute stabilizer (elements fixing the rep)
    vector<matrix3d> stabilizer;
    for (int g = 0; g < n_g; g++)
      if (sym.perms[g][oi.rep] == oi.rep)
        stabilizer.push_back(sym.matrices[g]);

    auto [dim, basis] = fixed_point_subspace(stabilizer);
    oi.subspace_dim = dim;
    oi.basis = basis;

    result.push_back(oi);
  }
  return result;
}

// Expand orbit-rep coordinates to full vertex set via group action.
// reduced[k] is the coord3d for orbit k's representative.
// Returns V3 of size nv.
static V3 expand_to_full(const V3& reduced, const vector<OrbitInfo>& orbits,
                          int nv, const SymmetryConstraint& sym) {
  V3 full(nv);
  for (size_t k = 0; k < orbits.size(); k++) {
    auto& oi = orbits[k];
    // Project reduced[k] into the orbit's fixed-point subspace
    coord3d x_rep;
    for (int d = 0; d < oi.subspace_dim; d++) {
      coord3d b(oi.basis(0,d), oi.basis(1,d), oi.basis(2,d));
      x_rep += b * reduced[k].dot(b);
    }
    for (size_t m = 0; m < oi.members.size(); m++) {
      int g = oi.gen_element[m];
      full[oi.members[m]] = sym.matrices[g] * x_rep;
    }
  }
  return full;
}

// Restrict full gradient to orbit-rep space via Reynolds averaging + subspace projection.
// For orbit rep v:  g_reduced[v] = P_sub * (1/|orbit|) sum_{m in orbit} R_g^T * g_full[m]
// where g maps rep to m, and P_sub projects onto the fixed-point subspace.
static V3 restrict_gradient(const V3& g_full, const vector<OrbitInfo>& orbits,
                             const SymmetryConstraint& sym) {
  V3 g_red(orbits.size());
  for (size_t k = 0; k < orbits.size(); k++) {
    auto& oi = orbits[k];
    // Average back-rotated gradients from all orbit members
    coord3d avg;
    for (size_t m = 0; m < oi.members.size(); m++) {
      int g = oi.gen_element[m];
      avg += sym.matrices[g].transpose() * g_full[oi.members[m]];
    }
    avg /= oi.members.size();

    // Project onto fixed-point subspace
    coord3d proj;
    for (int d = 0; d < oi.subspace_dim; d++) {
      coord3d b(oi.basis(0,d), oi.basis(1,d), oi.basis(2,d));
      proj += b * avg.dot(b);
    }
    g_red[k] = proj;
  }
  return g_red;
}

// Symmetry-aware initialization from MDS result.
// Uses the Reynolds operator to project MDS positions into the standard
// crystallographic frame defined by the symmetry matrices, then extracts
// reduced parameters in each orbit rep's fixed-point subspace.
static V3 symmetry_aware_init(const V3& mds_full, const matrix<double>& D,
                               const vector<OrbitInfo>& orbits,
                               const SymmetryConstraint& sym) {
  V3 reduced(orbits.size());
  for (size_t k = 0; k < orbits.size(); k++) {
    auto& oi = orbits[k];

    // Reynolds average: pull all orbit members' MDS positions back to the
    // standard frame via R_g^T, then average.  This extracts the component
    // of the MDS that is consistent with the standard-frame group action.
    coord3d reynolds;
    for (size_t m = 0; m < oi.members.size(); m++) {
      int g = oi.gen_element[m];
      reynolds += sym.matrices[g].transpose() * mds_full[oi.members[m]];
    }
    reynolds /= oi.members.size();

    // Project onto fixed-point subspace
    coord3d proj;
    for (int d = 0; d < oi.subspace_dim; d++) {
      coord3d b(oi.basis(0,d), oi.basis(1,d), oi.basis(2,d));
      proj += b * reynolds.dot(b);
    }

    // If projection is too small (collapsed MDS or unlucky cancellation),
    // perturb along subspace basis using APSP distances as scale.
    double scale = 0;
    for (size_t i = 0; i < oi.members.size(); i++)
      for (size_t j = i+1; j < oi.members.size(); j++)
        scale = max(scale, D(oi.members[i], oi.members[j]));
    if (oi.members.size() == 1) {
      for (int v = 0; v < (int)mds_full.size(); v++)
        if (v != oi.rep) scale = max(scale, D(oi.rep, v));
    }

    if (proj.norm() < 0.01 * scale) {
      for (int d = 0; d < oi.subspace_dim; d++) {
        coord3d b(oi.basis(0,d), oi.basis(1,d), oi.basis(2,d));
        proj += b * (0.3 * scale * (d == 0 ? 1.0 : 0.3));
      }
    }
    reduced[k] = proj;
  }
  return reduced;
}

// Polar decomposition: extract the orthogonal factor Q from M = Q*S.
// Uses SVD (M = U Sigma W^T) to return Q = U W^T, with W the right singular
// vectors (eigenvectors of MtM = M^T M).  eigensystem() returns each unit
// eigenvector as a ROW of C, so W = C^T and Q = U W^T = U C.
static matrix3d polar_orthogonal(const matrix3d& M) {
  matrix3d MtM = M.transpose() * M;
  auto [evals, C] = MtM.eigensystem();

  double svals[3] = { sqrt(max(0.0, evals[0])),
                      sqrt(max(0.0, evals[1])),
                      sqrt(max(0.0, evals[2])) };

  matrix3d U;
  int computed = 0;
  for (int col = 0; col < 3; col++) {
    if (svals[col] > 1e-12) {
      coord3d w_col(C(col,0), C(col,1), C(col,2));   // right singular vector = C row
      coord3d u_col = M * w_col * (1.0 / svals[col]);
      U(0,col) = u_col[0]; U(1,col) = u_col[1]; U(2,col) = u_col[2];
      computed++;
    }
  }

  if (computed == 2) {
    int z = -1;
    for (int i = 0; i < 3; i++) if (svals[i] <= 1e-12) { z = i; break; }
    if (z >= 0) {
      int a = (z+1)%3, b = (z+2)%3;
      coord3d ua(U(0,a),U(1,a),U(2,a)), ub(U(0,b),U(1,b),U(2,b));
      coord3d uc = ua.cross(ub);
      double n = uc.norm();
      if (n > 1e-15) uc /= n;
      U(0,z) = uc[0]; U(1,z) = uc[1]; U(2,z) = uc[2];
    }
  } else if (computed < 2) {
    return matrix3d::unit_matrix();
  }

  return U * C;
}

// Align symmetry matrices to the embedding frame via generalized intertwiner.
//
// The standard-frame matrices R_g from representation_3d() are related to the
// embedding-frame action by conjugation: M_g = Q R_g Q^T.  We find Q by:
//   1. Compute per-element Procrustes P_g ~ M_g from positions.
//   2. Form generalized intertwiner T_A = sum_g P_g * A * R_g^T.
//      For irreducible representations, T_A = (|G|/3) tr(Q^T A) Q.
//      The standard choice A=I gives coefficient tr(Q), which vanishes
//      when Q is a 120 degree rotation.  We try A = e_i e_j^T (coefficient Q_{ji})
//      and pick the probe with largest ||T_A||, guaranteeing a robust signal.
//   3. Extract Q via polar decomposition of T_A.
//   4. Return conjugated matrices Q R_g Q^T (exact group structure).
static SymmetryConstraint align_symmetry_to_embedding(
    const V3& x, const SymmetryConstraint& sym)
{
  int nv = x.size();
  int ng = sym.perms.size();

  // Step 1: Compute all Procrustes rotations P_g
  vector<matrix3d> P(ng);
  for (int g = 0; g < ng; g++) {
    matrix3d M_g;
    for (int v = 0; v < nv; v++)
      M_g += x[sym.perms[g][v]].outer(x[v]);
    P[g] = polar_orthogonal(M_g);
  }

  // Step 2: Try generalized intertwiners T_A for different seed matrices A.
  // A=I: T = sum P_g R_g^T, coefficient = (|G|/3) tr(Q)
  // A=e_i*e_j^T: T = sum (P_g e_i)(R_g e_j)^T, coefficient = (|G|/3) Q_{ji}
  // Pick the one with largest Frobenius norm.
  matrix3d best_T;
  double best_norm_sq = -1;

  auto try_probe = [&](const matrix3d& T_candidate) {
    double nsq = 0;
    for (int a = 0; a < 3; a++)
      for (int b = 0; b < 3; b++)
        nsq += T_candidate(a,b) * T_candidate(a,b);
    if (nsq > best_norm_sq) {
      best_norm_sq = nsq;
      best_T = T_candidate;
    }
  };

  // Probe 0: A = I (standard intertwiner)
  {
    matrix3d T;
    for (int g = 0; g < ng; g++)
      T += P[g] * sym.matrices[g].transpose();
    try_probe(T);
  }

  // Probes 1-9: A = e_i e_j^T
  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      matrix3d T;
      for (int g = 0; g < ng; g++) {
        // P_g e_i = column i of P_g
        // R_g e_j = column j of R_g
        // outer product (P_g e_i)(R_g e_j)^T
        coord3d pe(P[g](0,i), P[g](1,i), P[g](2,i));
        coord3d re(sym.matrices[g](0,j), sym.matrices[g](1,j), sym.matrices[g](2,j));
        T += pe.outer(re);
      }
      try_probe(T);
    }
  }

  // Step 3: Extract frame rotation Q
  matrix3d Q = polar_orthogonal(best_T);

  // Step 4: Conjugate standard-frame matrices: R'_g = Q R_g Q^T
  SymmetryConstraint aligned;
  aligned.perms = sym.perms;
  aligned.matrices.resize(ng);
  matrix3d Qt = Q.transpose();
  for (int g = 0; g < ng; g++)
    aligned.matrices[g] = Q * sym.matrices[g] * Qt;

  return aligned;
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
  // H = k * [alpha * (p@p, q@q, p@q) + correction terms]
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

vector<coord3d> embed_delaunay_3d(const DelaunayTriangulation& DT,
                                   const SymmetryConstraint& sym)
{
  // Aliases for DCEL accessors used throughout.
  const int nv = DT.nv, nh = DT.nh, nf = DT.nf;
  const auto& he_next   = DT.he_next;
  const auto& he_origin = DT.he_origin;
  const auto& he_length = DT.he_length;
  const auto& v_out     = DT.v_out;
  const auto& v_cone_angle = DT.v_cone_angle;
  const auto& f_he      = DT.f_he;
  auto alive = [&](int h) { return DT.alive(h); };
  auto dest  = [&](int h) { return DT.dest(h); };
  auto cw    = [&](int h) { return DT.cw(h); };
  // --- Step 1: Extract shortest edge per vertex pair ---
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
  vector<MDSEdge> mds_edges;
  for (auto& e : edges) mds_edges.push_back({e.u, e.v, e.L});
  V3 x_mds = mds_placement(D, mds_edges);

  // Center at origin (needed for Procrustes alignment in symmetry path)
  {
    coord3d center;
    for (int v = 0; v < nv; v++) center += x_mds[v];
    center /= nv;
    for (int v = 0; v < nv; v++) x_mds[v] = x_mds[v] - center;
  }

  // Orient MDS outward BEFORE optimization.  MDS eigenvector signs are
  // arbitrary, so ~50% of the time the embedding is inside-out.  If the
  // optimizer starts from an inside-out configuration it may converge to
  // a non-convex local minimum with inverted faces.
  // Signed volume is translation-invariant for a closed surface, so we
  // can use the origin directly (the MDS is already centered there).
  {
    double vol = 0;
    for (int f = 0; f < nf; f++) {
      if (f_he[f] < 0) continue;
      int h0 = f_he[f], h1 = he_next[h0], h2 = he_next[h1];
      int u = he_origin[h0], v = he_origin[h1], w = he_origin[h2];
      vol += x_mds[u].dot(x_mds[v].cross(x_mds[w]));
    }
    if (vol < 0)
      for (int v = 0; v < nv; v++) x_mds[v] = x_mds[v] * (-1.0);
  }

  // --- Step 4: Collect per-vertex fan structure for cone angle energy ---
  struct VertexFan {
    vector<int> out_halfedges;
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

  double lambda = 1.0;
  // No butterfly-fold penalty — empirically it increases self-intersections
  // by distorting the energy landscape (the penalty fires on concave edges,
  // which are normal geometry, and pushes the optimizer into worse basins).

  // --- Collect face vertex triples from DCEL ---
  struct FaceVerts { int u, v, w; };
  vector<FaceVerts> face_verts;
  for (int f = 0; f < nf; f++) {
    if (f_he[f] < 0) continue;
    int h0 = f_he[f], h1 = he_next[h0], h2 = he_next[h1];
    face_verts.push_back({he_origin[h0], he_origin[h1], he_origin[h2]});
  }


  // --- Step 5: Full-space energy and Hessian-vector product ---
  // These operate on V3(nv) -- the full vertex set.
  auto eval_full = [&](const V3& pos, V3& g) -> double {
    double E = 0;
    g.zero();

    for (auto& e : edges) {
      auto ed = EdgeStressData::compute(pos, e.u, e.v, e.L);
      if (!ed.valid()) continue;
      E += ed.energy();
      ed.scatter_gradient(g);
    }

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
      E += lambda * dev * dev;

      double w = 2.0 * lambda * dev;
      for (int i = 0; i < deg; i++)
        if (fa[i].valid) fa[i].scatter_gradient(w, g);
    }


    return E;
  };

  auto hv_full = [&](const V3& pos, const V3& dir, V3& Hv) {
    Hv.zero();

    for (auto& e : edges) {
      auto ed = EdgeStressData::compute(pos, e.u, e.v, e.L);
      if (!ed.valid()) continue;
      ed.scatter_hv(dir, Hv);
    }

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

      double dA_dot_v = 0;
      for (int i = 0; i < deg; i++)
        if (fa[i].valid) {
          dA_dot_v += fa[i].p.dot(dir[fa[i].a] - dir[v])
                    + fa[i].q.dot(dir[fa[i].d] - dir[v]);
        }
      for (int i = 0; i < deg; i++)
        if (fa[i].valid) fa[i].scatter_gradient(2.0 * lambda * dA_dot_v, Hv);

      if (fabs(dev) > 1e-15) {
        double w2 = 2.0 * lambda * dev;
        for (int i = 0; i < deg; i++)
          if (fa[i].valid) fa[i].scatter_hv(w2, dir, Hv);
      }
    }

  };

  V3 x;

  if (sym.empty()) {
    // --- No symmetry: optimize in full space with rigid-body projection ---
    x = x_mds;

    auto eval = [&](const V3& pos, V3& g) -> double {
      double E = eval_full(pos, g);
      project_rigid_body(g, pos);
      return E;
    };
    auto hv = [&](const V3& pos, const V3& dir, V3& Hv) {
      hv_full(pos, dir, Hv);
      project_rigid_body(Hv, pos);
    };
    x = steihaug_cg(std::move(x), eval, hv);

  } else {
    // --- Symmetric initialization + full-space optimization ---
    //
    // The iDT edge set is NOT G-invariant (co-circular obstructions make a
    // symmetric iDT impossible for some isomers). Optimizing in the reduced
    // symmetric subspace loses gradient information from the asymmetric edges,
    // creating phantom stationary points.
    //
    // Instead: use the symmetry to compute a good INITIALIZATION (Reynolds-
    // projected MDS positions), then optimize in the FULL space (same as the
    // nosym path). This gives the symmetric starting basin advantage without
    // losing gradient information.

    SymmetryConstraint aligned = align_symmetry_to_embedding(x_mds, sym);
    auto orbit_info = compute_orbit_structure(nv, aligned);

    // Symmetric initialization: Reynolds-project MDS into the symmetric subspace,
    // expand to full vertex set. This is a better starting point than raw MDS
    // because it respects the approximate symmetry of the polyhedron.
    V3 x_red = symmetry_aware_init(x_mds, D, orbit_info, aligned);
    V3 x_sym = expand_to_full(x_red, orbit_info, nv, aligned);

    // Full-space optimizer (same as nosym path).
    auto eval = [&](const V3& pos, V3& g) -> double {
      double E = eval_full(pos, g);
      project_rigid_body(g, pos);
      return E;
    };
    auto hv = [&](const V3& pos, const V3& dir, V3& Hv) {
      hv_full(pos, dir, Hv);
      project_rigid_body(Hv, pos);
    };

    // Try symmetric init first; if it's worse than plain MDS, use MDS.
    V3 g_sym(nv), g_mds(nv);
    double E_sym = eval(x_sym, g_sym);
    double E_mds = eval(x_mds, g_mds);

    x = (E_sym < E_mds) ? x_sym : x_mds;
    x = steihaug_cg(std::move(x), eval, hv);
  }

  // --- Orient outward using signed volume (translation-invariant) ---
  // Point inversion (negate all coords) reverses all face orientations
  // without changing the shape.
  auto orient_outward = [&](V3& pos) {
    double vol = 0;
    for (auto& fv : face_verts)
      vol += pos[fv.u].dot(pos[fv.v].cross(pos[fv.w]));
    if (vol < 0)
      for (auto& xi : pos) xi = xi * (-1.0);
  };

  orient_outward(x);

  // --- Reflect-optimize loop: fix concave vertices ---
  // For each vertex, compute its height h above the neighbor centroid plane.
  // If h < 0, the vertex is concave (dipped inward).  Reflect it to h = +|h|
  // (basin-switching), then re-optimize.  Repeat until no concave vertices.
  {
    auto reflect_concave = [&](V3& pos) -> int {
      int count = 0;
      for (int v = 0; v < nv; v++) {
        if (v_out[v] < 0) continue;
        auto& fan = fans[v];
        int deg = fan.out_halfedges.size();
        if (deg < 3) continue;

        // Neighbor centroid
        coord3d centroid;
        for (int i = 0; i < deg; i++)
          centroid += pos[dest(fan.out_halfedges[i])];
        centroid /= deg;

        // Fan normal (average of consecutive edge cross products)
        coord3d normal;
        for (int i = 0; i < deg; i++) {
          coord3d a = pos[dest(fan.out_halfedges[i])] - centroid;
          coord3d b = pos[dest(fan.out_halfedges[(i+1)%deg])] - centroid;
          normal += a.cross(b);
        }
        double nn = normal.norm();
        if (nn < 1e-15) continue;
        normal /= nn;

        double h = (pos[v] - centroid).dot(normal);
        // normal points inward (CW fan order), so h > 0 = concave.
        if (h > 0) {
          pos[v] = centroid + normal * (-h);
          count++;
        }
      }
      return count;
    };

    auto eval = [&](const V3& pos, V3& g) -> double {
      double E = eval_full(pos, g);
      project_rigid_body(g, pos);
      return E;
    };
    auto hv = [&](const V3& pos, const V3& dir, V3& Hv) {
      hv_full(pos, dir, Hv);
      project_rigid_body(Hv, pos);
    };

    for (int pass = 0; pass < 10; pass++) {
      int n_refl = reflect_concave(x);
      if (n_refl == 0) break;
      x = steihaug_cg(std::move(x), eval, hv);
      orient_outward(x);
    }
  }

  return x;
}
