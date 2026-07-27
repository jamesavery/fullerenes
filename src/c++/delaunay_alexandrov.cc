// Alexandrov embeddings of fullerene metrics (Bobenko-Izmestiev).
//
// Given an intrinsic Delaunay triangulation T of a polyhedral metric on S²
// with geodesic edge lengths, find the unique convex polyhedron P ⊂ R³
// whose boundary is isometric to T.  Two metrics are provided:
//
//   AlexandrovSolver   — any cone iDT (the fullerene dual's 12-cone metric
//                        in production use; the solver is n-generic).
//   AlexandrovIDTCubic — the CUBIC polyhedral metric (flat regular unit
//                        pentagons/hexagons; 20..60 cones), a thin wrapper
//                        that builds the kis metric and feeds the solver.
//
// The solver parameterizes P by the radii r = (|a−v_i|), where a is an
// interior apex point.  The curvature κ_i = 2π − ω_i (angle deficit at
// the radial edge a−v_i) satisfies κ = 0 iff the GCP is a genuine polytope.
// We trace the homotopy κ(r) = t·κ₁ from (t=1, large R) to (t→0, r*)
// by natural t-continuation (BI eq. 38), then extrapolate and polish.
// (The legacy pseudo-arc-length path and its experiment apparatus live in
// attic/delaunay_alexandrov_palc.cc.attic.)
//
// Layers:
//   GCP          — observables: κ(T,r), J(T,r), θ(T,r,h)
//   TrustRegion  — LM subproblem solver, accept/reject rule
//   Topology     — weighted Delaunay flip maintenance
//   Continuation — natural-t predictor-corrector homotopy tracking
//   Newton       — trust-region Newton polish for κ(r)=0
//   Reconstruct  — BFS vertex placement from Gram matrix entries
//   AlexandrovSolver / AlexandrovIDTCubic — the public entry points

#include "fullerenes/delaunay_alexandrov.hh"
#include "fullerenes/dense_linalg.hh"
#include "fullerenes/matrix.hh"
#include <cstdio>
#include <cmath>
#include <map>
#include <numeric>
#include <set>
#include <stdexcept>
#include <string>
#include <vector>
#include <algorithm>

using namespace std;

// All dense linear algebra goes through fullerenes/dense_linalg.hh —
// deliberately NOT BLAS/LAPACK (see the constraint note there: at least
// one deployed OpenBLAS silently corrupts dgesv from n ≳ 60).

namespace {

// Forward declaration for Reconstruct::from_radii — used by make_traj
// before Reconstruct is defined.  Implementation in Layer 6b below.
namespace Reconstruct {
  std::vector<coord3d> from_radii(const DelaunayTriangulation& T,
                                   const std::vector<double>& r);
}

// ============================================================================
// Layer 1: Generalized Convex Polytope observables
//
// Pure functions of (DCEL T, radii r).  The pyramid over face f with
// apex a has base triangle from T and lateral edges of length r[v].
// ============================================================================

namespace GCP {

// Angle at vertex i in the Euclidean triangle (apex, i, j).
// ρ_e in the Bobenko-Izmestiev notation, for oriented edge e: i→j.
double rho(const vector<double>& r, double L, int i, int j) {
  if (r[i] < 1e-15 || L < 1e-15) return NAN;
  double cs = (r[i]*r[i] + L*L - r[j]*r[j]) / (2*r[i]*L);
  return acos(clamp(cs, -1.0, 1.0));
}

// Angle at apex in the Euclidean triangle (apex, i, j).
// φ_e in the paper, for oriented edge e: i→j.
double phi(const vector<double>& r, double L, int i, int j) {
  if (r[i] < 1e-15 || r[j] < 1e-15) return NAN;
  double cs = (r[i]*r[i] + r[j]*r[j] - L*L) / (2*r[i]*r[j]);
  return acos(clamp(cs, -1.0, 1.0));
}

// Development of the pyramid over face(h) into the base plane.
//
// The pyramid has apex a and base triangle (u,v,w) with |a-u|=r_u, etc.
// Develop the base into the plane: u at origin, v on the x-axis, so
// w = (wx, wy); project a onto this plane: p = (px, py) satisfies
// |p-u|²+h²=r_u², |p-v|²+h²=r_v², |p-w|²+h²=r_w², with h_sq the squared
// pyramid height.  base_ok = false iff the base triangle inequality
// fails (wy_sq < −1e-10·lwu²); h_sq < 0 (beyond the callers' guards)
// means the pyramid fails to close (Cayley-Menger violation).
// Single source for alpha() and the min-h² diagnostic (pyramid_h_sq_at).
struct PyramidDev { double wy, py, h_sq; bool base_ok; };

PyramidDev develop_pyramid(const DelaunayTriangulation& T,
                            const vector<double>& r, int h) {
  int u = T.he_origin[h], v = T.he_origin[T.he_next[h]], w = T.he_origin[T.prev(h)];
  double luv = T.he_length[h], lvw = T.he_length[T.he_next[h]], lwu = T.he_length[T.prev(h)];

  double wx = (luv*luv + lwu*lwu - lvw*lvw) / (2*luv);  // cosine rule
  double wy_sq = lwu*lwu - wx*wx;
  if (wy_sq < -1e-10 * lwu*lwu) return {0, 0, 0, false};
  double wy = sqrt(max(0.0, wy_sq));

  double px = (r[u]*r[u] - r[v]*r[v] + luv*luv) / (2*luv);
  double py = (wy > 1e-15) ? (r[u]*r[u] - r[w]*r[w] + wx*wx + wy*wy - 2*px*wx) / (2*wy) : 0;
  double h_sq = r[u]*r[u] - px*px - py*py;
  return {wy, py, h_sq, true};
}

// Dihedral angle α at base edge h in the pyramid over face(h):
// α = atan2(h, py) in the development above, where h is the pyramid
// height and py the signed distance from the apex projection to edge uv.
double alpha(const DelaunayTriangulation& T, const vector<double>& r, int h) {
  auto d = develop_pyramid(T, r, h);
  if (!d.base_ok) return NAN;                    // triangle inequality violation
  int u = T.he_origin[h];
  if (d.h_sq < -1e-10 * r[u]*r[u]) return NAN;   // pyramid doesn't close
  double height = sqrt(max(0.0, d.h_sq));

  return atan2(height, d.py);
}

// Total dihedral angle θ at edge h: sum of α from both adjacent pyramids.
double theta(const DelaunayTriangulation& T, const vector<double>& r, int h) {
  return alpha(T, r, h) + alpha(T, r, T.twin(h));
}

// Curvature κ_v = 2π − ω_v, where ω_v is the total solid angle around
// the radial edge a−v.  For each incident face (v,j,k), the spherical
// section at v is a spherical triangle with sides ρ_vj, ρ_vk and
// included angle = intrinsic face angle at v.  By the spherical cosine
// rule, the opposite angle (= dihedral at a−v in that pyramid) is:
//
//   cos(ω_face) = (cos(face_angle) − cos ρ_vj cos ρ_vk) / (sin ρ_vj sin ρ_vk)
//
double curvature_at(const DelaunayTriangulation& T, const vector<double>& r, int v) {
  if (T.v_out[v] < 0) return 0;
  double omega = 0;
  for (int h : T.incident(v)) {
    double rho_vj = rho(r, T.he_length[h], v, T.dest(h));
    double rho_vk = rho(r, T.he_length[T.prev(h)], v, T.dest(T.he_next[h]));
    double sj = sin(rho_vj), sk = sin(rho_vk);
    if (sj < 1e-15 || sk < 1e-15) { omega = NAN; break; }
    double cos_omega = (cos(T.he_angle[h]) - cos(rho_vj)*cos(rho_vk)) / (sj*sk);
    omega += acos(clamp(cos_omega, -1.0, 1.0));
  }
  return 2*M_PI - omega;
}

// Curvature vector κ(T, r).
vector<double> kappa(const DelaunayTriangulation& T, const vector<double>& r) {
  vector<double> k(T.nv);
  for (int v = 0; v < T.nv; v++) k[v] = curvature_at(T, r, v);
  return k;
}

// Bobenko–Izmestiev total scalar curvature (BI 2008, Definition 3.1):
//   H(T, r) = Σ_v r_v · κ_v(r) + Σ_e ℓ_e · (π − θ_e(r))
// where θ_e = α_e + α_{e^1} is the full dihedral at base edge e.
// Per BI Proposition 5 eq. (13), ∂H/∂r_v = κ_v exactly.
// Note: H is NOT strictly concave — the Hessian has Lorentzian signature
// (1, n−1) by BI Theorem 4 + Lemma 3.4.  H is exposed for diagnostics
// and verification, not as a merit function for direct optimisation.
double H(const DelaunayTriangulation& T, const vector<double>& r) {
  double s = 0;
  auto kv = kappa(T, r);
  for (int v = 0; v < T.nv; v++) s += r[v] * kv[v];
  for (int h = 0; h < T.nh; h += 2) {
    if (!T.alive(h)) continue;
    s += T.he_length[h] * (M_PI - theta(T, r, h));
  }
  return s;
}

// Feasibility: r ∈ F(T) iff every incident pyramid closes.  Any guard in
// alpha() (wy_sq<0 or h_sq<0) yields NAN; a single non-finite return
// short-circuits the scan.
bool feasible(const DelaunayTriangulation& T, const vector<double>& r) {
  for (int h = 0; h < T.nh; h++) {
    if (!T.alive(h)) continue;
    if (!isfinite(alpha(T, r, h))) return false;
  }
  return true;
}

// Largest s ∈ [0, 1] such that r_from + s·(r_to − r_from) ∈ F(T).
// Assumes r_from ∈ F.  Returns 1 if the entire segment is feasible;
// otherwise bisects in s to ~2⁻ⁿ precision.  Used by the endgame guard
// (Tier 1) and prospectively by any other line search constrained to F.
double feasibility_max_step(const DelaunayTriangulation& T,
                             const vector<double>& r_from,
                             const vector<double>& r_to,
                             int n_iter = 40) {
  if (feasible(T, r_to)) return 1.0;
  vector<double> rs(r_from.size());
  double lo = 0, hi = 1;
  for (int it = 0; it < n_iter; it++) {
    double mid = 0.5 * (lo + hi);
    for (size_t i = 0; i < r_from.size(); i++)
      rs[i] = r_from[i] + mid * (r_to[i] - r_from[i]);
    if (feasible(T, rs)) lo = mid; else hi = mid;
  }
  return lo;
}

// Safety factor for landing strictly inside F(T) when a step would otherwise
// reach the feasibility boundary.  Matches B-I numerical practice.
constexpr double FEAS_SAFETY = 0.95;

// Clip step δ so r_from + δ ∈ F(T).  Returns the (possibly scaled) step δ'.
//   - if r_from + δ ∈ F(T):  δ' = δ                 (no clip)
//   - else:                  δ' = (FEAS_SAFETY · s_max) · δ
// where s_max ∈ [0,1] is the largest feasible step (feasibility_max_step).
//
// Pre: r_from ∈ F(T).
// Post: r_from + result ∈ F(T) strictly.
//
// Single source of truth for the F(T)-feasible step rule.  Used by:
//   - solve() endgame extrapolation
//   - Newton::polish trust-region step
// `clipped` (if non-null) is set true iff a scale was applied.
vector<double> feasible_step(const DelaunayTriangulation& T,
                              const vector<double>& r_from,
                              const vector<double>& delta,
                              bool* clipped = nullptr) {
  vector<double> r_to(r_from.size());
  for (size_t i = 0; i < r_from.size(); i++) r_to[i] = r_from[i] + delta[i];
  double s_max = feasibility_max_step(T, r_from, r_to);
  bool clip = (s_max < 1.0);
  double s = clip ? FEAS_SAFETY * s_max : 1.0;
  if (clipped) *clipped = clip;
  vector<double> out(delta.size());
  for (size_t i = 0; i < delta.size(); i++) out[i] = s * delta[i];
  return out;
}

// Per-oriented-edge Jacobian contribution.
// For half-edge h (oriented edge e: i→j):
//   J_e = (cot α_e + cot α_{−e}) / (ℓ_e sin ρ_e sin ρ_{−e})
// Returns NAN if degenerate (propagates failure to linear solve).
double J_edge(const DelaunayTriangulation& T, const vector<double>& r, int h) {
  int i = T.he_origin[h], j = T.dest(h);
  double L = T.he_length[h];

  double rho_e  = rho(r, L, i, j);
  double rho_me = rho(r, L, j, i);
  double alpha_e  = alpha(T, r, h);
  double alpha_me = alpha(T, r, T.twin(h));

  double sr = sin(rho_e), srm = sin(rho_me);
  double sa = sin(alpha_e), sam = sin(alpha_me);
  if (sr < 1e-15 || srm < 1e-15 || sa < 1e-15 || sam < 1e-15) return NAN;

  return (cos(alpha_e)/sa + cos(alpha_me)/sam) / (L * sr * srm);
}

// Jacobian J(T, r) = ∂κ/∂r.
// Off-diagonal: J(i,j) = Σ_{edges i→j} J_e(h)
// Diagonal:     J(i,i) = −Σ_{e: a(e)=i} cos(φ_e) · J_e(h)
// (per-oriented-edge formula, correct for multi-edges)
// @post symmetric BITWISE: each undirected edge's J_e is computed once
//       and the same double is added to J(i,j) and J(j,i), so Jᵀ = J
//       exactly in floating point (not merely in exact arithmetic).
//       TrustRegion::solve's Gauss-Newton form J² = JᵀJ relies on this;
//       if the construction ever changes to independent per-entry
//       evaluations, restore an explicit transpose there.
matrix<double> jacobian(const DelaunayTriangulation& T, const vector<double>& r) {
  int n = T.nv;
  matrix<double> J(n, n, 0.0);

  // Off-diagonal
  for (int h = 0; h < T.nh; h += 2) {
    if (!T.alive(h)) continue;
    double Je = J_edge(T, r, h);
    J(T.he_origin[h], T.dest(h)) += Je;
    J(T.dest(h), T.he_origin[h]) += Je;
  }

  // Diagonal (per-oriented-edge, not per-vertex-pair)
  for (int i = 0; i < n; i++) {
    if (T.v_out[i] < 0) continue;
    double diag = 0;
    for (int h : T.incident(i)) {
      double Je = J_edge(T, r, h);
      double phi_e = phi(r, T.he_length[h], i, T.dest(h));
      diag -= cos(phi_e) * Je;
    }
    J(i, i) = diag;
  }
  return J;
}

} // namespace GCP

// ============================================================================
// Layer 2: Linear algebra
// ============================================================================

// The linear-algebra primitives — vector reductions, the cyclic-Jacobi
// eigensolver, and the LU with determinant sign — live in
// fullerenes/dense_linalg.hh (namespace LinAlg); the vector arithmetic
// operators come from auxiliary.hh's global templates.
using LinAlg::V;

// Pack one continuation- or Newton-step diagnostic into a TraceEntry.  Caller
// decides whether to record (avoids the cost of κ and J spectrum when
// trace_jacobian is off); this just bundles the fields.
static AlexandrovSolver::TraceEntry make_trace(
    char phase, int step, double t, double ds, int nit,
    const V& kappa, const matrix<double>& J) {
  return {phase, step, t, ds, nit,
          LinAlg::max_abs(kappa), LinAlg::norm(kappa), LinAlg::sym_eigvals(J)};
}

// Pyramid height squared at half-edge h, via GCP::develop_pyramid.
// Negative result means the pyramid fails to close (Cayley-Menger
// violation); −1 when the base triangle itself is degenerate.
static double pyramid_h_sq_at(const DelaunayTriangulation& T,
                                const V& r, int h) {
  auto d = GCP::develop_pyramid(T, r, h);
  return d.base_ok ? d.h_sq : -1.0;
}

// Mean edge length ⟨ℓ⟩ over alive edges — the length scale used by the
// volume normalisation and the convexity tolerance.  1 on an edgeless T.
static double mean_edge_length(const DelaunayTriangulation& T) {
  double sum_l = 0; int n_e = 0;
  for (int h = 0; h < T.nh; h += 2) {
    if (!T.alive(h)) continue;
    sum_l += T.he_length[h]; n_e++;
  }
  return (n_e > 0) ? sum_l / n_e : 1.0;
}

// Six times the signed volume enclosed by pos under T's face structure:
// Σ_f a·(b×c) over live faces (signed tetrahedra from the origin) —
// positive iff the CCW half-edge convention has outward normals.  The
// same divergence-theorem quantity as PolyhedronView::volume_tetra, on
// (DCEL, pos) data.  Used by the outward-orientation flip, the volume-
// degeneracy gate, and the convexity precondition.
static double signed_volume6(const DelaunayTriangulation& T,
                              const vector<coord3d>& pos) {
  double vol6 = 0;
  for (int f = 0; f < T.nf; f++) {
    if (T.f_he[f] < 0) continue;
    const auto v = T.face_vertices(f);
    vol6 += pos[v[0]].dot(pos[v[1]].cross(pos[v[2]]));
  }
  return vol6;
}

// Pack one continuation- or Newton-step trajectory-diagnostic record.  Cheap to
// compute (O(nh) for theta + h_sq scans, plus one LU for det sign);
// gated by AlexandrovSolver::record_diag at call site.
static AlexandrovSolver::DiagEntry make_diag(
    char phase, int step, double t, double ds, int nit,
    const DelaunayTriangulation& T, const V& r,
    const V& kappa, const matrix<double>& J, int n_flips_cum) {
  AlexandrovSolver::DiagEntry e;
  e.phase = phase; e.step = step; e.t = t; e.ds = ds; e.nit = nit;
  e.kappa_max = LinAlg::max_abs(kappa);
  e.n_flips_cum = n_flips_cum;

  // θ stats over non-bigon alive edges.
  double min_dist = M_PI;
  int n01 = 0, n001 = 0, n0001 = 0, n_alive = 0;
  for (int h = 0; h < T.nh; h += 2) {
    if (!T.alive(h)) continue;
    if (T.is_bigon(h)) continue;                     // θ undefined on bigons
    double th = GCP::theta(T, r, h);
    if (!std::isfinite(th)) continue;
    n_alive++;
    double d = M_PI - th;
    if (d < min_dist) min_dist = d;
    if (d < 0.1)   n01++;
    if (d < 0.01)  n001++;
    if (d < 1e-3)  n0001++;
  }
  e.theta_min_dist_to_pi = min_dist;
  e.n_near_pi_01 = n01;
  e.n_near_pi_001 = n001;
  e.n_near_pi_0001 = n0001;
  e.n_non_bigon_alive = n_alive;

  // F(T) margin: smallest pyramid h_sq across all alive half-edges.
  double mhs = std::numeric_limits<double>::infinity();
  for (int h = 0; h < T.nh; h++) {
    if (!T.alive(h)) continue;
    double hs = pyramid_h_sq_at(T, r, h);
    if (hs < mhs) mhs = hs;
  }
  e.min_h_sq = mhs;

  // r coefficient of variation.
  double mean = 0;
  for (int i = 0; i < (int)r.size(); i++) mean += r[i];
  mean /= r.size();
  double var = 0;
  for (int i = 0; i < (int)r.size(); i++) var += (r[i] - mean) * (r[i] - mean);
  var /= r.size();
  e.r_cv = (mean > 0) ? std::sqrt(var) / mean : 0;

  // sign(det J): one extra LU.
  V dummy(J.m, 0.0);
  auto sol = LinAlg::solve_with_sign(J, dummy);
  e.det_J_sign = sol ? sol->det_sign : 0;

  return e;
}

// Pack one per-step homotopy-trajectory entry for visualization (gated by
// AlexandrovSolver::record_trajectory at the call site).  Reconstructs the GCP
// cone positions (apex at origin) and lists the triangular faces of T, so the
// deforming pyramids and the flipping base triangulation can be drawn directly.
static AlexandrovSolver::TrajEntry make_traj(
    char phase, int step, double t,
    const DelaunayTriangulation& T, const V& r, const V& kappa) {
  AlexandrovSolver::TrajEntry e;
  e.phase = phase; e.step = step; e.t = t;
  e.kappa_max = LinAlg::max_abs(kappa);
  e.kappa = kappa;
  e.positions = Reconstruct::from_radii(T, r);   // empty/NaN if Gram-BFS fails
  for (int f = 0; f < T.nf; f++) {
    if (T.f_he[f] < 0 || !T.alive(T.f_he[f])) continue;   // dead face slot
    const auto h = T.face_halfedges(f);
    e.faces.push_back({T.he_origin[h[0]], T.he_origin[h[1]], T.he_origin[h[2]]});
    e.face_len.push_back({T.he_length[h[0]], T.he_length[h[1]], T.he_length[h[2]]});
  }
  return e;
}

// ============================================================================
// Layer 3: Trust-region subproblem
// ============================================================================

namespace TrustRegion {

// Solve the trust-region subproblem for E = ½||κ||²: try the pure
// Newton root step δ = −J⁻¹κ first (λ=0); outside the radius (or with
// singular J), fall back to Gauss-Newton Levenberg-Marquardt on the
// NORMAL equations, (JᵀJ + λI)δ = −Jᵀκ, bisecting on λ ≥ 0 for
// ||δ|| ≈ Δ.
//
// The normal-equations form is essential, not cosmetic: J is the
// Lorentzian Hessian of the B-I functional (one positive eigenvalue,
// the rest negative), so the former shifted-J system (J+λI)δ = −κ was
// singular at every λ in the negative spectrum, and its large-λ limit
// −κ/λ has directional derivative −κᵀJκ > 0 for that signature — a
// systematic ASCENT direction for E that made every damped recovery
// step reject (the C134 hard stall: trust radius collapsing 13 orders
// with κ frozen).  JᵀJ + λI is SPD for every λ > 0 — no poles,
// ||δ(λ)|| monotone (well-posed bisection) — and its large-λ limit
// −Jᵀκ/λ = −∇E/λ is steepest descent on E, so damped steps always
// make progress and reject cascades terminate.
V solve(const matrix<double>& J, const V& kappa, double Delta) {
  // Try pure Newton (λ=0) — identical to the pre-GN behaviour.
  auto delta = LinAlg::solve(J, -kappa);
  if (LinAlg::is_usable_step(delta) && LinAlg::norm(delta) <= Delta)
    return delta;

  // J is symmetric bitwise (see jacobian() @post), so Jᵀ = J exactly
  // and the Gauss-Newton objects need no transpose:
  // JᵀJ = J² (SPD), via THE view-level product (dense_linalg_view.hh's
  // matmul -- the one bit-pinned i-j-k body, shared with the batch port;
  // matrix.hh's generic operator* is no longer on any solver path).
  const matrix<double> JtJ = J * J;   // JtJ = J^2 (SPD); operator* delegates
                                      // to LinAlg::matmul (@ref matmul-ijk-order)
  const V minus_Jtk = -LinAlg::matvec(J, kappa);         // −Jᵀκ = −∇E

  // Bisect on λ to find (JᵀJ+λI)⁻¹(−Jᵀκ) with ||δ|| ≈ Δ
  double lo = 0, hi = LinAlg::max_abs(minus_Jtk) / Delta + 1.0;
  for (int probe = 0; probe < 10; probe++) {
    delta = LinAlg::solve_shifted(JtJ, minus_Jtk, hi);
    if (LinAlg::is_usable_step(delta) && LinAlg::norm(delta) <= Delta) break;
    hi *= 4;
  }
  for (int bis = 0; bis < 20; bis++) {
    double mid = 0.5 * (lo + hi);
    delta = LinAlg::solve_shifted(JtJ, minus_Jtk, mid);
    if (!LinAlg::is_usable_step(delta) || LinAlg::norm(delta) > Delta) lo = mid;
    else hi = mid;
  }
  return LinAlg::solve_shifted(JtJ, minus_Jtk, hi);
}

// Predicted reduction: E(κ) − E(κ + J·δ) where E = ½||·||².
double predicted_reduction(const matrix<double>& J, const V& kappa, const V& delta) {
  return LinAlg::energy(kappa) - LinAlg::energy(kappa + LinAlg::matvec(J, delta));
}

// Trust-region accept/reject.  Returns (accepted, new_Delta).
pair<bool, double> update(double actual, double predicted, double dnorm,
                           double Delta, double Delta_max) {
  double rho = (predicted > 1e-30) ? actual / predicted : -1;
  if (rho > 0.1) {
    double new_Delta = (rho > 0.75 && dnorm > 0.5 * Delta)
                       ? min(2.0 * Delta, Delta_max) : Delta;
    return {true, new_Delta};
  } else {
    return {false, max(Delta * 0.25, 1e-14)};
  }
}

} // namespace TrustRegion

// ============================================================================
// Layer 4: Topology operations
// ============================================================================

namespace Topology {

// B-I 2008 §3.4 (lines 614–640) define an edge h as "bad" iff the function
// q̃_T fails Q-concavity across h.  They give two cases:
//
//   (ConcQuadr) — h has two distinct adjacent triangles forming a strictly-
//     convex quadrilateral, with angle at the opposite vertex ≥ π.  The bad
//     condition is θ_h > π (the GCP dihedral exceeds π).
//
//   (CloGeod)   — h is the i–j edge of an "iji" bigon-face (a face with both
//     half-edges of h, two i-corners, and one j-corner; B-I lines 562–563).
//     The bad condition is q_j < q_i − ℓ²_ij, i.e. r_j² < r_i² − ℓ²_ij.
//
// Our previous needs_flip implemented only (ConcQuadr).  The (CloGeod) clause
// is needed because we start the homotopy from a Δ-complex iDT — bigon
// faces may be present and need flipping during the path to maintain
// Q-concavity.
// Closure is checked with a −1e-10 numerical-noise buffer: flip when
// θ > π + 1e-10, resp. q_j < q_i − ℓ²_ij − 1e-10, so that legitimate
// flat-face diagonals (θ → π at the polytope limit) do not trigger flips.
int needs_flip(const DelaunayTriangulation& T, const vector<double>& r) {
  constexpr double margin = -1e-10;   // closure buffer
  // (ConcQuadr): edge with two distinct adjacent triangles, θ > π − margin.
  // For a bigon edge, GCP::theta is not meaningful; skip and let the
  // (CloGeod) pass below handle it.
  // NaN θ — returned by alpha() when the abstract pyramid is degenerate
  // (h_sq < 0, i.e. Cayley-Menger fails) — also counts as bad: the
  // configuration is outside P(M) and a flip is required.  IEEE
  // comparisons with NaN are false, so we test isnan() explicitly.
  double theta_threshold = M_PI - margin;
  for (int h = 0; h < T.nh; h += 2) {
    if (!T.alive(h)) continue;
    if (T.is_bigon(h)) continue;                           // handled below
    double theta = GCP::theta(T, r, h);
    if (std::isnan(theta) || theta > theta_threshold) return h;
  }
  // (CloGeod): bigon i–j edge of an iji-face, q_j < q_i − ℓ²_ij + margin.
  for (int h = 0; h < T.nh; h += 2) {
    if (!T.alive(h)) continue;
    if (!T.is_bigon(h)) continue;
    int u = T.he_origin[h], v = T.dest(h);
    if (u == v) continue;                                  // pure self-loop, not i–j edge
    // Identify i (the doubled vertex) and j (the single vertex) from the
    // bigon face's third half-edge: it must be a self-loop with origin = i.
    int f = T.he_face[h];
    int hh = T.f_he[f], h_self = -1;
    for (int s = 0; s < 3; s++, hh = T.he_next[hh]) {
      if (T.he_origin[hh] == T.dest(hh)) { h_self = hh; break; }
    }
    if (h_self < 0) continue;                              // no self-loop in this face — not iji-shape
    int i = T.he_origin[h_self];
    int j = (i == u) ? v : u;
    if (i != u && i != v) continue;                        // i must be one of h's endpoints
    double q_i = r[i] * r[i], q_j = r[j] * r[j];
    double ell_ij_sq = T.he_length[h] * T.he_length[h];
    if (q_j < q_i - ell_ij_sq + margin) return h;
  }
  return -1;
}

// Flip until needs_flip(T, r) returns -1 (or 20 iter cap).
// Returns number of flips performed.
int flip_to_weighted_delaunay(DelaunayTriangulation& T, const vector<double>& r) {
  int total = 0;
  for (int iter = 0; iter < 20; iter++) {
    int h = needs_flip(T, r);
    if (h < 0) break;
    if (T.flip_edge(h)) total++;
    else break;
  }
  return total;
}

} // namespace Topology

// ============================================================================
// Layer 5: Natural t-continuation of the BI homotopy
//
// Traces the homotopy curve F(t,r) = κ(r) − t·κ₁ = 0 from t=1 toward t=0
// by predictor-corrector steps in t.  (The legacy pseudo-arc-length
// tracker and its experiment apparatus are retired to
// attic/delaunay_alexandrov_palc.cc.attic.)
// ============================================================================

namespace Continuation {

// Homotopy velocity dr/dt = J⁻¹κ₁ along κ(r)=t·κ₁ (the un-normalized tangent
// direction). Returns an invalid vector (LinAlg::is_usable_step == false) if J is
// singular.
V dr_dt(const matrix<double>& J, const V& kappa1) {
  return LinAlg::solve(J, kappa1);
}

// Adapt step size based on corrector iterations (AUTO strategy).
double adapt_dt(double dt, int nit, int max_nit, double dt_max) {
  if (nit <= 1) dt *= 2.0;
  else if (nit == 2) dt *= 1.5;
  else if (nit <= max_nit / 2) dt *= 1.1;
  else if (nit >= max_nit) dt *= 0.5;
  return min(fabs(dt), dt_max);
}

constexpr int CORRECTOR_MAX_ITER = 8;
constexpr double CORRECTOR_TOL = 1e-12;

// One predictor-corrector step.  Returns the corrector iteration count
// alongside the (t, r) iterate — always present; "accepted" iff
// nit ∈ [0, CORRECTOR_MAX_ITER).  t, r are valid only if accepted.
struct StepResult {
  int nit;
  double t; V r;
  bool accepted() const { return nit >= 0 && nit < CORRECTOR_MAX_ITER; }
};

struct TrackStats {
  int steps = 0, flips = 0, newton_total = 0;
};

struct TrackResult {
  double t_final;
  V r_final;
  vector<pair<double, V>> history;
  TrackStats stats;
};

// Initial radii for the BI homotopy: r = 2·R_max·1.
// Pre:  T is a valid iDT with positive edge lengths.
// Post: r ∈ F(T) and 0 < κ_i(T, r) < δ_i (BI admissibility).
V initial_radii(const DelaunayTriangulation& T) {
  double R = 0;
  for (int h = 0; h < T.nh; h += 2)
    if (T.alive(h)) R = max(R, T.he_length[h]);
  return V(T.nv, 2 * R);
}

// Plain Newton driving κ(r) → target: r −= J⁻¹(κ(r) − target).
// Returns the iteration count: nit ∈ [0, max_iter) converged to tol, max_iter
// if not, −1 on a failed linear solve.  (Newton::polish is the trust-region
// κ→0 analogue.)
int newton_correct(const DelaunayTriangulation& T, V& r, const V& target,
                   double tol, int max_iter) {
  for (int nit = 0; nit < max_iter; nit++) {
    V F = GCP::kappa(T, r) - target;
    if (LinAlg::max_abs(F) < tol) return nit;
    auto dr = LinAlg::solve(GCP::jacobian(T, r), F);
    if (!LinAlg::is_usable_step(dr)) return -1;
    r = r - dr;
  }
  return max_iter;
}

// ---- Natural (t-parameterized) continuation ----
//
// Bobenko–Izmestiev 2008 (eq. 38) lift the curvature path κ(r)=t·κ₁ as a
// function of t:  dr/dt = J⁻¹κ₁.  By Thm 5 / Lemma 3.4, J = ∂κ/∂r is
// non-degenerate with constant Lorentzian signature (1,n−1) while 0<κᵢ<δᵢ —
// which the homotopy preserves for all t∈(0,1] (Lemma 4.2) — so the path is
// monotone in t with NO turning points.  Plain t-continuation therefore needs
// neither arclength nor the bordered pseudo-arclength system.  It is
// scale-invariant by construction: t is dimensionless, dr/dt=J⁻¹κ₁ is degree
// +1, and the corrector residual κ(r)−t·κ₁ is degree 0 — nothing mixes the
// dimensionless homotopy parameter with the length-dimensioned radii the way
// the PALC arclength ds²=dt²+‖dr‖² does (the root cause of PALC's scale bug).

// Step bounds on dt, the dimensionless homotopy parameter t∈[0,1]:
// DT_MAX caps it at 10% of the t-range per step; DT_MIN is the give-up floor.
constexpr double DT_MIN = 1e-9;
constexpr double DT_MAX = 0.1;

// One natural-continuation step to t1 = t0 − dt: Euler predictor r0 − dt·(dr/dt)
// then a fixed-t1 Newton correction of κ(r) onto t1·κ₁.
StepResult natural_step(const DelaunayTriangulation& T,
                        double t0, const V& r0, const V& kappa1,
                        const matrix<double>& J, double dt) {
  double t1 = t0 - dt;
  auto v = dr_dt(J, kappa1);
  if (!LinAlg::is_usable_step(v)) return {-1, t1, r0};
  V r = r0 - v * dt;
  int nit = newton_correct(T, r, kappa1 * t1, CORRECTOR_TOL, CORRECTOR_MAX_ITER);
  return nit < 0 ? StepResult{-1, t1, r0} : StepResult{nit, t1, std::move(r)};
}

// Record (t, r) in the extrapolation history, collapsing same-t entries
// (a re-accepted step at unchanged t replaces its predecessor).
void record_history(vector<pair<double, V>>& history, double t, const V& r) {
  if (!history.empty() && fabs(t - history.back().first) < 1e-14)
    history.back() = {t, r};
  else
    history.push_back({t, r});
}

// Natural t-continuation: trace κ(r)=t·κ₁ from t=1 toward t_target by
// natural_step predictor-corrector steps, re-Delaunay-flipping after each
// accepted step and adapting dt (clamped against overshooting t_target).
//   Pre:  r ∈ F(T), 0 < κᵢ(T, r) < δᵢ, t_target > 0.
//   Post: t_final ≤ t_target if the continuation reached it, else it stalled.
TrackResult natural_track(DelaunayTriangulation& T, V r, const V& kappa1,
                          double t_target, double dt_init,
                          vector<AlexandrovSolver::TraceEntry>* trace,
                          vector<AlexandrovSolver::DiagEntry>* diag,
                          vector<AlexandrovSolver::TrajEntry>* traj = nullptr) {
  double t = 1.0, dt = dt_init;
  vector<pair<double, V>> history;
  TrackStats stats;

  for (int step_i = 0; step_i < 500 && t > t_target; step_i++) {
    auto J = GCP::jacobian(T, r);
    auto result = natural_step(T, t, r, kappa1, J, min(dt, t - t_target));

    if (result.accepted()) {
      t = result.t; r = std::move(result.r);
      stats.flips += Topology::flip_to_weighted_delaunay(T, r);
      record_history(history, t, r);
      dt = adapt_dt(dt, result.nit, CORRECTOR_MAX_ITER, DT_MAX);
    } else {
      dt *= 0.5;
      if (dt < DT_MIN) break;
    }
    if (trace)
      trace->push_back(make_trace('T', step_i, t, dt, result.nit,
                                   GCP::kappa(T, r), J));
    if (diag) {
      auto Jd = GCP::jacobian(T, r);  // J on the (possibly post-flip) state
      diag->push_back(make_diag('T', step_i, t, dt, result.nit,
                                  T, r, GCP::kappa(T, r), Jd, stats.flips));
    }
    if (traj)
      traj->push_back(make_traj('T', step_i, t, T, r, GCP::kappa(T, r)));
    stats.steps++;
    stats.newton_total += max(result.nit, 0);
  }
  return {t, std::move(r), std::move(history), stats};
}

// Polynomial extrapolation: given (t_i, r_i) pairs approaching t=0,
// fit a polynomial r(t) and evaluate at t=0.
vector<double> extrapolate(const vector<pair<double, vector<double>>>& history) {
  int k = history.size();
  if (k == 0) return {};
  if (k == 1) return history[0].second;

  int n = history[0].second.size();

  // Use the last min(k, 4) points for polynomial fit
  int m = min(k, 4);
  int start = k - m;

  // Lagrange interpolation at t = 0
  V result(n, 0.0);
  for (int j = 0; j < m; j++) {
    double tj = history[start + j].first;
    double basis = 1.0;
    for (int l = 0; l < m; l++) {
      if (l != j) basis *= -history[start + l].first / (tj - history[start + l].first);
    }
    result = result + history[start + j].second * basis;
  }
  return result;
}

} // namespace Continuation

// ============================================================================
// Layer 6a: Trust-region Newton for κ(r) = 0
// ============================================================================

namespace Newton {

// Minimize E = ½||κ||² using LM trust-region Newton.
// Returns (converged, final_kappa).
// If `out_trace` is non-null, records one TraceEntry per iteration.
//
// κ is only piecewise-smooth: crossing a flip boundary (θ_e = π)
// changes the weighted-Delaunay cell, and κ/J evaluated on the stale
// triangulation disagree with the true energy.  Two invariants keep
// every evaluation inside its correct cell:
//   - ENTRY: the endgame extrapolation moves r without flipping, so the
//     iterate can arrive past a boundary (θ > π); restore
//     weighted-Delaunay before the first κ/J evaluation.  (Without this
//     the C134 stall: 21 straight rejects with Δ shrinking 13 orders,
//     κ frozen, because flips only ran after ACCEPTED steps.)
//   - TRIAL: evaluate each trial on its own flipped copy (T_trial,
//     r_trial), so acceptance compares true energies across cells; on
//     acceptance adopt the copy atomically, on rejection discard it.
//     (Without this the C124 soft stall: all-reject collapse at a
//     feasible point near a boundary.)
pair<bool, double> polish(DelaunayTriangulation& T, V& r,
                           double tol = 1e-10, int max_iter = 50,
                           std::vector<AlexandrovSolver::TraceEntry>* out_trace = nullptr,
                           std::vector<AlexandrovSolver::DiagEntry>* out_diag = nullptr,
                           std::vector<AlexandrovSolver::TrajEntry>* out_traj = nullptr,
                           int* flips_cum = nullptr) {
  using LinAlg::energy; using LinAlg::norm; using LinAlg::max_abs;

  double r_avg = LinAlg::dot(r, V(r.size(), 1.0)) / r.size();
  double Delta = 0.5 * r_avg, Delta_max = 2.0 * r_avg;
  int rejects = 0;
  int flips_local = flips_cum ? *flips_cum : 0;

  flips_local += Topology::flip_to_weighted_delaunay(T, r);   // ENTRY invariant

  for (int iter = 0; iter < max_iter; iter++) {
    auto kappa = GCP::kappa(T, r);
    if (max_abs(kappa) < tol) {
      if (flips_cum) *flips_cum = flips_local;
      return {true, max_abs(kappa)};
    }
    if (rejects > 20) break;

    double E         = energy(kappa);
    auto   J         = GCP::jacobian(T, r);
    // Levenberg-Marquardt-bisected trust-region step.
    V delta_raw = TrustRegion::solve(J, kappa, Delta);
    // Clip step to F(T): keep r_trial inside the feasibility region so
    // pyramids stay non-degenerate and κ(r_trial) is finite (no NaN θ on
    // multi-edges, no false |κ|=0 convergence outside F).  Same helper
    // used by the endgame extrapolation — single source of truth.
    bool   clipped;
    auto   delta     = GCP::feasible_step(T, r, delta_raw, &clipped);
    double pred      = TrustRegion::predicted_reduction(J, kappa, delta);
    V      r_trial   = r + delta;
    // TRIAL invariant: κ(r_trial) on r_trial's own weighted-Delaunay
    // cell.  The copy snapshots any active point tracker; adopting it on
    // acceptance keeps transport commit-or-nothing.
    DelaunayTriangulation T_trial = T;
    int    trial_flips = Topology::flip_to_weighted_delaunay(T_trial, r_trial);
    double E_trial   = energy(GCP::kappa(T_trial, r_trial));

    auto [ok, D2] = TrustRegion::update(E - E_trial, pred, norm(delta), Delta, Delta_max);
    // If we clipped, cap Δ at the step we actually took so the next
    // subproblem doesn't keep proposing the same infeasible direction.
    Delta = clipped ? min(D2, norm(delta)) : D2;
    if (out_trace)
      out_trace->push_back(make_trace('N', iter, 0.0, Delta, ok ? 1 : 0, kappa, J));
    if (ok) {
      r = r_trial;
      T = std::move(T_trial);
      flips_local += trial_flips;
      rejects = 0;
    } else rejects++;

    if (out_diag) {
      auto Jd = GCP::jacobian(T, r);
      out_diag->push_back(make_diag('N', iter, 0.0, Delta, ok ? 1 : 0,
                                      T, r, GCP::kappa(T, r), Jd, flips_local));
    }
    if (out_traj)
      out_traj->push_back(make_traj('N', iter, 0.0, T, r, GCP::kappa(T, r)));
  }
  if (flips_cum) *flips_cum = flips_local;
  return {false, max_abs(GCP::kappa(T, r))};
}

} // namespace Newton

// ============================================================================
// Layer 6b: 3D reconstruction from converged (T, r)
//
// With κ ≈ 0, the apex sits at the origin and |pos[v]| = r[v].
// The Gram matrix entry pos[u]·pos[v] = (r[u]² + r[v]² − L²_{uv})/2
// is known for each edge.  We place vertices face-by-face via BFS,
// solving 3 inner-product constraints per new vertex.
// ============================================================================

namespace Reconstruct {

// Gram entry: pos[u]·pos[v] = (r_u² + r_v² − L²_uv) / 2.
double gram(double r_u, double r_v, double L_uv) {
  return (r_u*r_u + r_v*r_v - L_uv*L_uv) / 2;
}

// Place vertex w given: |w| = r_w, w·pos[u] = g_wu, w·pos[v] = g_wv.
coord3d place_vertex(coord3d pu, coord3d pv, double g_wu, double g_wv,
                      double r_w, coord3d p_old) {
  double uu = pu.dot(pu), uv = pu.dot(pv), vv = pv.dot(pv);
  double det = uu*vv - uv*uv;
  double a = 0, b = 0;
  if (fabs(det) > 1e-20) {
    a = (g_wu*vv - g_wv*uv) / det;
    b = (g_wv*uu - g_wu*uv) / det;
  }

  coord3d proj  = pu*a + pv*b;
  coord3d n     = pu.cross(pv);
  double  nl    = n.norm();
  double  gamma_sq = r_w*r_w - proj.dot(proj);
  if (gamma_sq < -1e-10 * r_w*r_w) return coord3d(NAN, NAN, NAN);
  double  gamma = sqrt(max(0.0, gamma_sq));
  double  scale = (nl > 1e-15) ? gamma / nl : 0;

  double side = (p_old - proj).dot(n);
  return (side > 0) ? proj - n*scale : proj + n*scale;
}

// 3D reconstruction from converged (T, r) with κ ≈ 0.
vector<coord3d> from_radii(const DelaunayTriangulation& T, const V& r) {
  int n = T.nv;
  vector<coord3d> pos(n, coord3d(0,0,0));
  vector<bool> placed(n, false), face_done(T.nf, false);

  // Seed face: place vertex i on x-axis, j in xy-plane, k with z ≥ 0.
  int f0 = -1;
  for (int f = 0; f < T.nf; f++) if (T.f_he[f] >= 0) { f0 = f; break; }
  if (f0 < 0) return pos;

  const auto [h0, h1, h2] = T.face_halfedges(f0);
  const auto [i, j, k]    = T.face_vertices(f0);

  pos[i] = coord3d(r[i], 0, 0);

  double g_ij = gram(r[i], r[j], T.he_length[h0]);
  double jx   = g_ij / r[i];
  double jy2  = r[j]*r[j] - jx*jx;
  if (jy2 < -1e-10 * r[j]*r[j]) return {};
  pos[j] = coord3d(jx, sqrt(max(0.0, jy2)), 0);

  double g_ik = gram(r[i], r[k], T.he_length[h2]);
  double g_jk = gram(r[j], r[k], T.he_length[h1]);
  double kx   = g_ik / r[i];
  double ky   = (pos[j][1] > 1e-15) ? (g_jk - kx*pos[j][0]) / pos[j][1] : 0;
  double kz2  = r[k]*r[k] - kx*kx - ky*ky;
  if (kz2 < -1e-10 * r[k]*r[k]) return {};
  pos[k] = coord3d(kx, ky, sqrt(max(0.0, kz2)));
  placed[i] = placed[j] = placed[k] = true;
  face_done[f0] = true;

  vector<int> queue = {f0};
  int head = 0;
  while (head < (int)queue.size()) {
    int f = queue[head++];
    int hf = T.f_he[f];
    for (int s = 0; s < 3; s++, hf = T.he_next[hf]) {
      int ht = hf ^ 1, fa = T.he_face[ht];
      if (fa < 0 || face_done[fa]) continue;

      int u = T.he_origin[hf], v = T.dest(hf), w = T.he_origin[T.prev(ht)];
      if (!placed[w]) {
        double g_wu = gram(r[w], r[u], T.he_length[T.he_next[ht]]);
        double g_wv = gram(r[w], r[v], T.he_length[T.prev(ht)]);

        int old_w = -1;
        int hf2 = T.f_he[f];
        for (int s2 = 0; s2 < 3; s2++, hf2 = T.he_next[hf2])
          if (T.he_origin[hf2] != u && T.he_origin[hf2] != v) { old_w = T.he_origin[hf2]; break; }

        pos[w] = place_vertex(pos[u], pos[v], g_wu, g_wv, r[w],
                               old_w >= 0 ? pos[old_w] : coord3d(0,0,0));
        placed[w] = true;
      }
      face_done[fa] = true;
      queue.push_back(fa);
    }
  }

  // Orient outward
  if (signed_volume6(T, pos) < 0)
    for (auto& p : pos) p = p * (-1.0);

  return pos;
}

} // namespace Reconstruct

} // anonymous namespace

// ============================================================================
// AlexandrovSolver: top-level 5-step algorithm
// ============================================================================

const char* AlexandrovSolver::status_str(ValidationStatus s) {
  switch (s) {
    case ValidationStatus::OK:                       return "OK";
    case ValidationStatus::FAIL_KAPPA_NOT_CONVERGED: return "FAIL_KAPPA_NOT_CONVERGED";
    case ValidationStatus::FAIL_NOT_SIMPLE:          return "FAIL_NOT_SIMPLE";
    case ValidationStatus::FAIL_RECONSTRUCT:         return "FAIL_RECONSTRUCT";
    case ValidationStatus::FAIL_VOLUME_DEGENERATE:   return "FAIL_VOLUME_DEGENERATE";
    case ValidationStatus::FAIL_SELF_INTERSECTING:   return "FAIL_SELF_INTERSECTING";
    case ValidationStatus::FAIL_NOT_CONVEX:          return "FAIL_NOT_CONVEX";
  }
  return "UNKNOWN";
}

// Validation gate: a returned polytope must satisfy three named
// properties, ALL required for valid output.  Fills S's stats_*
// diagnostics and returns the verdict; called by solve() after
// κ-convergence and reconstruction have already succeeded.
//
//   SIMPLICITY    — T̄(0) is a simple polygonal tesselation: every
//                   polygon has ≥ 3 distinct vertex labels and no
//                   repeated label.  At κ = 0 the inessential
//                   collapse reduces T(0) (which may carry redundant
//                   multi-edges per refined I-1) to T̄(0); simplicity
//                   is enforced on T̄(0).  Plus F ≥ 3 (no drum-cap,
//                   which would force all V = 12 vertices coplanar by
//                   Euler).
//
//   WELL-FORMEDNESS — reconstruct() returns a closed manifold, the
//                     polytope has non-degenerate volume
//                     (vol_norm > 0.01, well below the 0.12 healthy
//                     floor and well above the ~1e-6 degenerate
//                     ceiling observed across 1.03M scan), and no two
//                     non-adjacent triangles intersect in 3D.
//
//   CONVEXITY     — every non-face vertex sits on the inside of every
//                   face plane (Alexandrov's theorem requires this).
//                   Outward normal taken from the half-edge CCW
//                   convention; defensive precondition rejects
//                   inverted-volume positions.
static AlexandrovSolver::ValidationStatus validate_polytope(
    AlexandrovSolver& S, const vector<coord3d>& pos) {
  using VS = AlexandrovSolver::ValidationStatus;
  const DelaunayTriangulation& D = S.D;

  // ── SIMPLICITY ──────────────────────────────────────────────────────
  //   T̄(0) must be a simple polygonal tesselation with F ≥ 3.
  //   inessential_eps stays at the strict default (1e-7).  Loosening it
  //   risks falsely collapsing borderline simple edges with |θ − π|
  //   slightly above 1e-7 but well below "almost flat", which can fold
  //   a non-degenerate polytope into a drum-cap.  Multi-edges that
  //   legitimately collapsed to one polytope edge typically reach
  //   |θ − π| ≲ 1e-7 once Newton converges; if the continuation stalls
  //   before that (κ residual > tol), the F ≥ 3 check catches the
  //   resulting degeneracy.
  S.stats_t0_simplicial = AlexandrovSolver::is_simplicial(D);  // diagnostic only
  vector<int> labels(D.nv);
  std::iota(labels.begin(), labels.end(), 0);
  auto tbar = AlexandrovSolver::polytope_tesselation(D, S.r, labels);
  S.stats_tbar_n_cells = tbar.n_cells();
  S.stats_tbar_simple_polygonal = AlexandrovSolver::is_simple_polygonal(tbar);

  bool simple = S.stats_tbar_simple_polygonal && S.stats_tbar_n_cells >= 3;
  if (!simple) {
    if (S.verbose)
      printf("  VALIDATION (simplicity) failed: T̄(0) simple_polygonal=%d, "
             "n_cells=%d (need F≥3), T(0) simplicial=%d (diagnostic).\n",
             S.stats_tbar_simple_polygonal, S.stats_tbar_n_cells,
             S.stats_t0_simplicial);
    return VS::FAIL_NOT_SIMPLE;
  }

  // Compute volume_norm = |V| / ⟨ℓ⟩³ for the well-formedness gate and
  // diagnostic stat.
  double volume = std::abs(signed_volume6(D, pos)) / 6.0;
  double mean_l = mean_edge_length(D);
  S.stats_volume_norm = volume / (mean_l * mean_l * mean_l);

  // ── WELL-FORMEDNESS ─────────────────────────────────────────────────
  //   Non-degenerate volume AND embedded surface (no 3D self-intersection).
  //   Both are central definitional requirements for a valid polytope.
  //   - Volume: empirical 1.03M scan shows a clean ~5-orders-of-magnitude
  //     gap between healthy polytopes (vol_norm ≥ 0.12, median 1.01) and
  //     degenerate ones (vol_norm ≤ 1e-6, includes drum-caps and
  //     Newton-stalled cases).  vol_norm > 0.01 sits in the empty gap.
  //   - Self-intersection: two non-adjacent triangles crossing in 3D
  //     means the surface is not embedded — it's not a valid polytope
  //     at all.  Convexity does imply non-self-intersection, so on
  //     convex outputs the test is informationally subsumed; but this
  //     does NOT make the test optional.  It is part of the definition
  //     of "well-formed polytope" and must be enforced independently
  //     so that any failure (numerical or otherwise) is flagged with
  //     its true cause rather than swept into a different bucket.
  constexpr double VOLUME_NORM_DEGENERATE = 0.01;
  bool well_formed_volume = std::isfinite(S.stats_volume_norm) &&
                              S.stats_volume_norm >= VOLUME_NORM_DEGENERATE;
  S.stats_polytope_no_self_intersect =
      !AlexandrovSolver::has_self_intersection(D, pos);
  // Volume gate first, then self-intersection gate.  Order is the order
  // of stats_status when multiple fail.
  if (!well_formed_volume) {
    if (S.verbose)
      printf("  VALIDATION (well-formedness/volume) failed: vol_norm=%.3e "
             "(threshold %.2e).\n",
             S.stats_volume_norm, VOLUME_NORM_DEGENERATE);
    return VS::FAIL_VOLUME_DEGENERATE;
  }
  if (!S.stats_polytope_no_self_intersect) {
    if (S.verbose)
      printf("  VALIDATION (well-formedness/self-intersection) failed: "
             "two non-adjacent triangles cross in 3D.\n");
    return VS::FAIL_SELF_INTERSECTING;
  }

  // ── CONVEXITY ───────────────────────────────────────────────────────
  //   Every non-face vertex on the inside of every face plane.
  //   Alexandrov's theorem requires the polytope to be convex; failure
  //   here indicates the reconstruction landed in a geometrically
  //   non-convex configuration.  Outward normal taken from the half-
  //   edge CCW convention; defensive precondition inside is_convex
  //   verifies signed volume is strictly positive.
  S.stats_polytope_convex = AlexandrovSolver::is_convex(D, pos);
  if (!S.stats_polytope_convex) {
    if (S.verbose)
      printf("  VALIDATION (convexity) failed: some vertex sticks out "
             "beyond a face plane.\n");
    return VS::FAIL_NOT_CONVEX;
  }

  return VS::OK;
}

// The 5-step B-I algorithm: initial radii → continuation (κ(r)=t·κ₁, t:1→0)
// → endgame extrapolation → Newton polish → reconstruct → validate.
vector<coord3d> AlexandrovSolver::solve() {
  stats_steps = stats_flips = stats_newton_total = 0;
  trace.clear();
  diag_trace.clear();
  trajectory.clear();
  stats_status = ValidationStatus::FAIL_KAPPA_NOT_CONVERGED;

  // 1. Initialize: uniform radii (or caller-provided override).
  if ((int)r_init_override.size() == D.nv) {
    r = r_init_override;
  } else {
    r = Continuation::initial_radii(D);
  }
  auto kappa1 = GCP::kappa(D, r);

  // 2. Continuation: natural t-continuation (BI eq. 38) of κ(r) = t·κ₁
  //    from t=1 toward t_target — scale-invariant, no arclength.
  auto track = Continuation::natural_track(
      D, r, kappa1, /*t_target=*/0.1, /*dt_init=*/0.05,
      trace_jacobian ? &trace : nullptr,
      record_diag ? &diag_trace : nullptr,
      record_trajectory ? &trajectory : nullptr);
  r = std::move(track.r_final);
  stats_steps        = track.stats.steps;
  stats_flips        = track.stats.flips;
  stats_newton_total = track.stats.newton_total;

  // 3. Endgame: guarded extrapolation (Tier 1)
  stats_extrap_kappa = 0;
  r_before_extrap = r;
  if (!track.history.empty()) {
    auto r_ext = Continuation::extrapolate(track.history);
    if (!r_ext.empty()) r = r + GCP::feasible_step(D, r, r_ext - r);
    stats_extrap_kappa = LinAlg::max_abs(GCP::kappa(D, r));
  }

  // 4. Polish: trust-region Newton on κ(r) = 0.  An installed
  // polish_override (incubation seam, see header) replaces the internal
  // polish on the identical post-extrapolation state; the internal
  // trace/diag recorders are not populated on that path.
  constexpr double KAPPA_POLISH_TOL = 1e-10;
  bool ok;
  double mk;
  if (polish_override) {
    ok = polish_override(D, r);
    mk = LinAlg::max_abs(GCP::kappa(D, r));
  } else {
    int polish_flips = stats_flips;
    std::tie(ok, mk) = Newton::polish(D, r, KAPPA_POLISH_TOL, 50,
                                      trace_jacobian ? &trace : nullptr,
                                      record_diag ? &diag_trace : nullptr,
                                      record_trajectory ? &trajectory : nullptr,
                                      &polish_flips);
    stats_flips = polish_flips;
  }
  stats_final_kappa = mk;

  if (verbose)
    printf("  %d continuation steps, %d flips, max|κ|=%.2e (%s)\n",
           stats_steps, stats_flips, mk, ok ? "converged" : "FAILED");

  // Reconstruct positions whether or not validation passes — we always
  // return SOMETHING inspectable, with stats_status communicating
  // validity.
  auto pos = Reconstruct::from_radii(D, r);
  if (pos.empty()) {
    stats_status = ValidationStatus::FAIL_RECONSTRUCT;
    if (verbose)
      printf("  VALIDATION (reconstruct) failed: Gram-BFS yielded "
             "negative perpendicular squared distance.\n");
    return pos;   // empty
  }

  // Acceptance = the polish's converged verdict AND the κ residual
  // itself below the polish target.  For the internal polish the two
  // are equivalent (ok ⟺ max|κ| < tol); the residual check is the
  // numerical backstop for the polish_override seam, whose ok verdict
  // would otherwise be trusted unchecked.  The former 0.01 band
  // accepted stalled solves with κ up to 1e-2 as OK, leaking κ-scale
  // edge errors into the reconstruction (delaunay-fillin failing
  // suite, 2026-07-25); any residual stall now fails loudly here.
  // NB max_abs poisons NaN to +inf, so a NaN κ cannot pass.
  if (!ok || !(mk < KAPPA_POLISH_TOL)) {
    stats_status = ValidationStatus::FAIL_KAPPA_NOT_CONVERGED;
    return pos;   // failed-but-inspectable positions
  }

  // 5. Validation: the three-property gate (see validate_polytope above).
  stats_status = validate_polytope(*this, pos);
  return pos;
}

vector<coord3d> AlexandrovSolver::reconstruct(const DelaunayTriangulation& T,
                                              const vector<double>& r) {
  return Reconstruct::from_radii(T, r);
}

AlexandrovSolver::AlexandrovPolytope
AlexandrovSolver::solve_polytope(const vector<int>& vertex_labels) {
  AlexandrovPolytope out;
  out.positions = solve();
  out.status = stats_status;
  if (out.positions.empty()) return out;   // FAIL_RECONSTRUCT
  // T(0) is now in D; r in this->r.  Build labels (identity if not given).
  vector<int> labels = vertex_labels;
  if (labels.empty()) {
    labels.resize(D.nv);
    std::iota(labels.begin(), labels.end(), 0);
  }
  out.tesselation = polytope_tesselation(D, r, labels);
  return out;
}

vector<double> AlexandrovSolver::kappa(const DelaunayTriangulation& T,
                                        const vector<double>& r) {
  return GCP::kappa(T, r);
}

double AlexandrovSolver::H(const DelaunayTriangulation& T,
                            const vector<double>& r) {
  return GCP::H(T, r);
}

matrix<double> AlexandrovSolver::jacobian(const DelaunayTriangulation& T,
                                            const vector<double>& r) {
  return GCP::jacobian(T, r);
}

vector<double> AlexandrovSolver::jacobian_eigvals(const DelaunayTriangulation& T,
                                                    const vector<double>& r) {
  return LinAlg::sym_eigvals(GCP::jacobian(T, r));
}

int AlexandrovSolver::jacobian_det_sign(const DelaunayTriangulation& T,
                                          const vector<double>& r) {
  auto J = GCP::jacobian(T, r);
  V dummy(J.m, 0.0);
  auto sol = LinAlg::solve_with_sign(J, dummy);
  return sol ? sol->det_sign : 0;
}

bool AlexandrovSolver::feasible(const DelaunayTriangulation& T,
                                 const vector<double>& r) {
  return GCP::feasible(T, r);
}

vector<double> AlexandrovSolver::feasible_step(const DelaunayTriangulation& T,
                                                const vector<double>& r,
                                                const vector<double>& delta,
                                                bool* clipped) {
  return GCP::feasible_step(T, r, delta, clipped);
}

int AlexandrovSolver::flip_to_weighted_delaunay(DelaunayTriangulation& T,
                                                 const vector<double>& r) {
  return Topology::flip_to_weighted_delaunay(T, r);
}

double AlexandrovSolver::theta(const DelaunayTriangulation& T,
                                const vector<double>& r, int h) {
  return GCP::theta(T, r, h);
}

vector<bool> AlexandrovSolver::inessential_edges(const DelaunayTriangulation& T,
                                                  const vector<double>& r,
                                                  double eps) {
  // B-I §3.4 (line 798): an edge h is inessential iff q̃_T is Q on a
  // neighborhood of any interior point of h.  At κ=0 this equates to
  // θ_h = π exactly — the two adjacent pyramids over h are coplanar in
  // 3D, so h is interior to a flat 2-face of P.  We implement the
  // numerical version with a tolerance.
  vector<bool> tight(T.nh, false);
  for (int h = 0; h < T.nh; h += 2) {
    if (!T.alive(h)) continue;
    // For bigon edges, GCP::theta is not meaningful and the edge is by
    // definition not a flat-face diagonal of any 2-face of P (it's part
    // of a degenerate iji-bigon face).  Mark non-inessential.
    if (T.is_bigon(h)) continue;
    double theta = GCP::theta(T, r, h);
    if (std::isfinite(theta) && std::fabs(theta - M_PI) < eps) {
      tight[h]         = true;
      tight[T.twin(h)] = true;
    }
  }
  return tight;
}

CanonicalTesselation AlexandrovSolver::polytope_tesselation(
    const DelaunayTriangulation& T,
    const vector<double>& r,
    const vector<int>& vertex_labels,
    double inessential_eps) {
  return T.canonical_tesselation(vertex_labels,
                                  inessential_edges(T, r, inessential_eps));
}

bool AlexandrovSolver::is_simplicial(const DelaunayTriangulation& T) {
  // Per invariant I-1 (CLAUDE.md): any non-simple feature in T contradicts
  // an isometric R³ embedding for a non-degenerate polytope.  Delegates to
  // the DCEL's own predicate (arc-map injectivity).  That predicate also
  // subsumes the same-face-both-sides ("bigon") test on well-formed
  // triangulated complexes: a triangle containing both h and h^1 forces
  // its third half-edge to be a self-loop, which injectivity rejects.
  return T.is_simplicial();
}

bool AlexandrovSolver::is_simple_polygonal(const CanonicalTesselation& tess) {
  // Each polygon: ≥ 3 entries, all distinct labels.
  for (const auto& poly : tess.cells) {
    if (poly.size() < 3) return false;
    std::set<int> seen;
    for (const auto& [label, L] : poly) {
      if (!seen.insert(label).second) return false;            // repeated label
    }
  }
  return true;
}

bool AlexandrovSolver::is_convex(const DelaunayTriangulation& T,
                                   const vector<coord3d>& pos,
                                   double tol) {
  if ((int)pos.size() != T.nv) return false;

  double mean_l  = mean_edge_length(T);
  double abs_tol = tol * mean_l;   // relative tolerance

  // Defensive precondition: signed volume from the CCW half-edge convention
  // must be strictly positive.  The vertex test below uses (b−a) × (c−a)
  // for three consecutive vertices in he_next order as the outward normal;
  // that direction is outward only if the global signed volume comes out
  // positive in this same convention.  Reconstruct::from_radii enforces
  // this by flipping positions when needed, but other callers (diagnostics
  // bypassing the solver pipeline) might not — so we double-check here.
  //
  // Threshold strictly > 0 (with a numerical-noise margin scaled by
  // mean_edge_length³): rejects (a) flat polytopes (vol = 0, e.g. drum-cap)
  // since the vertex test would vacuously pass on coplanar configurations
  // and (b) globally inverted positions (vol < 0) where the CCW normal
  // points inward.
  double vol6 = signed_volume6(T, pos);
  double vol_threshold6 = 1e-6 * mean_l * mean_l * mean_l;  // 1e-7 vol_norm
  if (!std::isfinite(vol6) || vol6 < vol_threshold6) return false;

  // Outward normal is defined locally by the CCW half-edge order: walking
  // he_next around any face traverses its boundary CCW as seen from outside,
  // so (b−a) × (c−a) for three consecutive vertices points outward.  This
  // assumes Reconstruct::from_radii has already flipped positions to enforce
  // positive signed volume — verified by the precondition above.  No
  // spherical-approximation assumption: works for nanotubes, oblate
  // polytopes, etc.
  for (int f = 0; f < T.nf; f++) {
    if (T.f_he[f] < 0) continue;
    const auto [va, vb, vc] = T.face_vertices(f);
    coord3d a = pos[va], b = pos[vb], c = pos[vc];
    coord3d nf_raw = (b - a).cross(c - a);
    double nlen = sqrt(nf_raw.dot(nf_raw));
    if (nlen < 1e-15) continue;                 // degenerate triangle, skip
    coord3d nf = nf_raw * (1.0 / nlen);

    // Every other vertex must be on inside (signed dist ≤ tol·mean_edge).
    for (int v = 0; v < T.nv; v++) {
      if (v == va || v == vb || v == vc) continue;
      double d = (pos[v] - a).dot(nf);
      if (d > abs_tol) return false;
    }
  }
  return true;
}

namespace {
// Möller-Trumbore-style triangle-triangle intersection in 3D.
// Returns true if triangles (a0,b0,c0) and (a1,b1,c1) intersect with
// non-empty common interior, with `tol` slack on signed-distance tests.
// Sharing a vertex/edge is allowed (returns false in that case — the
// caller filters adjacent triangles up front).
//
// Algorithm: for each triangle, classify the other triangle's vertices
// by signed distance to its plane.  If all three are on one side
// (strictly), no intersection.  Otherwise, compute the line of
// intersection of the two planes, project both triangles' edge crossings
// onto it, check parameter intervals overlap.
bool tri_tri_intersect(const coord3d& a0, const coord3d& b0, const coord3d& c0,
                        const coord3d& a1, const coord3d& b1, const coord3d& c1,
                        double tol) {
  auto plane_dist = [&](const coord3d& a, const coord3d& b, const coord3d& c,
                          const coord3d& p) {
    coord3d n = (b - a).cross(c - a);
    return n.dot(p - a);   // unscaled signed distance · 2·area
  };
  double da0 = plane_dist(a1, b1, c1, a0);
  double db0 = plane_dist(a1, b1, c1, b0);
  double dc0 = plane_dist(a1, b1, c1, c0);
  if ((da0 > tol && db0 > tol && dc0 > tol) ||
      (da0 < -tol && db0 < -tol && dc0 < -tol)) return false;
  double da1 = plane_dist(a0, b0, c0, a1);
  double db1 = plane_dist(a0, b0, c0, b1);
  double dc1 = plane_dist(a0, b0, c0, c1);
  if ((da1 > tol && db1 > tol && dc1 > tol) ||
      (da1 < -tol && db1 < -tol && dc1 < -tol)) return false;

  // Both triangles' planes interleave.  Compute their line of
  // intersection direction L and project each triangle to it; check
  // overlap of parameter intervals.
  coord3d n0 = (b0 - a0).cross(c0 - a0);
  coord3d n1 = (b1 - a1).cross(c1 - a1);
  coord3d L  = n0.cross(n1);
  double L2 = L.dot(L);
  // |L|² = |n0|²·|n1|² · sin²(angle).  Treat as coplanar when
  // sin² < 1e-12 (≈ 1 micro-radian misalignment).  Healthy polytopes
  // cannot have two distinct flat 2-faces in the same plane (that
  // would be drum-cap, caught by the volume gate), so coplanar pairs
  // are sub-triangulations of the same T̄(0) face and never represent
  // real surface self-intersection.
  double n0_2 = n0.dot(n0), n1_2 = n1.dot(n1);
  if (L2 < 1e-12 * n0_2 * n1_2) {
    // Coplanar non-adjacent triangles.  On a valid convex polytope this
    // happens for sub-triangulations of the same flat 2-face: the iDT
    // triangulates a flat polygonal face into multiple coplanar
    // triangles, which then "share coplanar interior" without that
    // being a true 3D self-intersection of the polytope's surface.
    // The other case — two distinct coplanar 2-faces in the same plane
    // — only happens for degenerate drum-cap polytopes, which the
    // upstream volume gate already rejects.  So coplanar pairs are
    // either benign (same face) or already caught (drum-cap); reporting
    // intersection here would only false-positive on legitimate
    // multi-triangle flat faces (e.g. C24-D6d, C36-D6h hexagonal caps).
    return false;
  }
  // For a triangle, project vertices to L (as scalar t = (p−origin)·L̂).
  // Compute t-interval where triangle crosses the other plane.
  auto interval = [&](const coord3d& a, const coord3d& b, const coord3d& c,
                        double da, double db, double dc) {
    auto edge_t = [&](const coord3d& p, const coord3d& q,
                        double dp, double dq) {
      double s = dp / (dp - dq);
      coord3d x = p + (q - p) * s;
      return x.dot(L);
    };
    pair<double,double> r{1e300, -1e300};
    auto upd = [&](double t) {
      if (t < r.first)  r.first  = t;
      if (t > r.second) r.second = t;
    };
    if ((da > 0) != (db > 0)) upd(edge_t(a, b, da, db));
    if ((db > 0) != (dc > 0)) upd(edge_t(b, c, db, dc));
    if ((dc > 0) != (da > 0)) upd(edge_t(c, a, dc, da));
    // Vertices exactly on plane (within tol) included.
    if (fabs(da) <= tol) upd(a.dot(L));
    if (fabs(db) <= tol) upd(b.dot(L));
    if (fabs(dc) <= tol) upd(c.dot(L));
    return r;
  };
  auto i0 = interval(a0, b0, c0, da0, db0, dc0);
  auto i1 = interval(a1, b1, c1, da1, db1, dc1);
  if (i0.first  > i0.second) return false;        // triangle 0 doesn't cross plane 1
  if (i1.first  > i1.second) return false;        // triangle 1 doesn't cross plane 0
  return !(i0.second < i1.first - tol || i1.second < i0.first - tol);
}
} // anonymous namespace

bool AlexandrovSolver::has_self_intersection(const DelaunayTriangulation& T,
                                                const vector<coord3d>& pos,
                                                double tol) {
  if ((int)pos.size() != T.nv) return false;

  // Collect each face's three vertex indices once.
  struct Tri { int a, b, c; };
  vector<Tri> tris;
  tris.reserve(T.nf);
  for (int f = 0; f < T.nf; f++) {
    if (T.f_he[f] < 0) continue;
    const auto v = T.face_vertices(f);
    tris.push_back({v[0], v[1], v[2]});
  }

  for (size_t i = 0; i < tris.size(); i++) {
    for (size_t j = i + 1; j < tris.size(); j++) {
      const auto& A = tris[i];
      const auto& B = tris[j];
      // Skip adjacent triangles (any shared vertex).
      bool share = (A.a == B.a || A.a == B.b || A.a == B.c ||
                    A.b == B.a || A.b == B.b || A.b == B.c ||
                    A.c == B.a || A.c == B.b || A.c == B.c);
      if (share) continue;
      if (tri_tri_intersect(pos[A.a], pos[A.b], pos[A.c],
                              pos[B.a], pos[B.b], pos[B.c], tol)) {
        return true;
      }
    }
  }
  return false;
}

// ============================================================================
// AlexandrovIDTCubic: the cubic polyhedral metric
// ============================================================================

AlexandrovIDTCubic::KisMetric AlexandrovIDTCubic::kis_metric(const TriangulationView& T)
{
  KisMetric M;
  M.Nv = T.N;
  M.fdeg.resize(M.Nv);
  for (node_t u = 0; u < T.N; u++) {
    M.fdeg[u] = T.degree(u);
    if (M.fdeg[u] != 5 && M.fdeg[u] != 6)
      throw logic_error("AlexandrovIDTCubic: dual vertex " + to_string(u) +
                        " has degree " + to_string(M.fdeg[u]) +
                        " (input is not a fullerene dual)");
  }

  // The oriented faces come from the lib: triangles() emits each face once,
  // CCW, with each directed arc a->b in exactly one face (the face left of
  // it, per the compute_faces_oriented convention).  Index them by arc, so
  // face_left(a, b) is a lookup and across(a, b) is face_left(b, a).
  M.triangle = T.triangles();
  const int ntri = M.triangle.size();  // = 2*Nv - 4

  map<arc_t, int> face_of_arc;
  for (int t = 0; t < ntri; t++) {
    const tri_t& f = M.triangle[t];
    for (int j = 0; j < 3; j++)
      face_of_arc[{f[j], f[(j + 1) % 3]}] = t;
  }
  auto face_left = [&face_of_arc](node_t a, node_t b) -> int {
    auto it = face_of_arc.find({a, b});
    if (it == face_of_arc.end())
      throw logic_error("AlexandrovIDTCubic::kis_metric: arc " + to_string(a) +
                        "->" + to_string(b) + " has no registered face");
    return it->second;
  };

  // Adjacency of the kis complex, CCW rings built directly (orientation
  // invariant: never orient after the fact).
  vector<vector<node_t>> kis(M.Nv + ntri);
  // Around face center u: the incident triangle barycenters, in u's ring order.
  for (node_t u = 0; u < T.N; u++)
    for (node_t v : T[u])
      kis[u].push_back(M.Nv + face_left(u, v));
  // Around the barycenter of (u,v,w) (CCW): interleave corners with the
  // across-edge barycenters, [u, across(u,v), v, across(v,w), w, across(w,u)],
  // which is CCW around the barycenter when (u,v,w) is CCW.
  for (int t = 0; t < ntri; t++) {
    const tri_t& f = M.triangle[t];
    for (int j = 0; j < 3; j++) {
      node_t a = f[j], b = f[(j + 1) % 3];
      kis[M.Nv + t].push_back(a);
      kis[M.Nv + t].push_back(M.Nv + face_left(b, a));
    }
  }
  M.K = Triangulation(Graph(Spanify::OwnedDenseGraph<node_t>(kis)));
  return M;
}

DelaunayTriangulation::EdgeLengthFn AlexandrovIDTCubic::KisMetric::edge_length_fn() const
{
  return [Nv = Nv, fdeg = fdeg](node_t a, node_t b) -> double {
    const bool fa = a < Nv, fb = b < Nv;
    if (fa && fb)
      throw logic_error("AlexandrovIDTCubic: two face centers cannot be adjacent");
    if (!fa && !fb) return 1.0;                // cubic edge (barycenter-barycenter)
    return fdeg[fa ? a : b] == 5 ? R5 : 1.0;   // pentagon / hexagon spoke
  };
}

void AlexandrovIDTCubic::build(const TriangulationView& T)
{
  KisMetric M = kis_metric(T);

  vector<int> new_to_old;
  DelaunayTriangulation D = DelaunayTriangulation::compute(
      M.K, M.edge_length_fn(), FLAT_TOL, &new_to_old, track_removed);

  cone_triangle.clear();
  cone_npent.clear();
  cone_kis_vertex = new_to_old;
  double total = 0;
  for (int i = 0; i < D.nv; i++) {
    const int old = new_to_old[i];
    if (old < M.Nv)
      throw logic_error("AlexandrovIDTCubic: face center " + to_string(old) +
                        " survived flat removal");
    const tri_t& f = M.triangle[old - M.Nv];
    int k = 0;
    for (int j = 0; j < 3; j++) k += (M.fdeg[f[j]] == 5);
    const double kappa = 2 * M_PI - D.v_cone_angle[i];
    if (k < 1 || fabs(kappa - k * M_PI / 15) > KAPPA_TOL)
      throw logic_error("AlexandrovIDTCubic: cone " + to_string(i) +
                        " has kappa = " + to_string(kappa) + ", expected " +
                        to_string(k) + "*pi/15");
    total += kappa;
    cone_triangle.push_back(f);
    cone_npent.push_back(k);
  }
  if (fabs(total - 4 * M_PI) > TOTAL_KAPPA_TOL)
    throw logic_error("AlexandrovIDTCubic: total curvature " +
                      to_string(total) + " != 4*pi");

  solver.D = std::move(D);
}

std::vector<coord3d> AlexandrovIDTCubic::solve(const TriangulationView& T)
{
  build(T);
  return solver.solve();
}

AlexandrovSolver::AlexandrovPolytope
AlexandrovIDTCubic::solve_polytope(const TriangulationView& T)
{
  build(T);
  return solver.solve_polytope();   // default labels = cone index (identity)
}

AlexandrovIDTCubic::FlatFaceCensus
AlexandrovIDTCubic::flat_face_census(const TriangulationView& T,
                                     const CanonicalTesselation& tess) const
{
  // Cells as sorted cone-label sets.
  set<vector<int>> cells;
  for (const auto& cell : tess.cells) {
    vector<int> vs;
    for (const auto& [label, len2] : cell) vs.push_back(label);
    sort(vs.begin(), vs.end());
    cells.insert(vs);
  }

  // Cone ring of each face u: invert cone_triangle once.  Rings come out
  // ascending because cone indices are pushed in increasing order.
  vector<vector<int>> ring(T.N);
  for (int i = 0; i < (int)cone_triangle.size(); i++)
    for (int j = 0; j < 3; j++) ring[cone_triangle[i][j]].push_back(i);

  FlatFaceCensus census;
  census.n_cells = tess.n_cells();
  for (node_t u = 0; u < T.N; u++) {
    const bool is_pent = (T.degree(u) == 5);
    if (!is_pent) census.n_hex++;
    const bool flat = ((int)ring[u].size() == T.degree(u)) &&
                      cells.count(ring[u]);
    if (is_pent) census.pent_flat += flat;
    else         census.hex_flat  += flat;
  }
  return census;
}
