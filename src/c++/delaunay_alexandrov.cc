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
//   Flat         — the doubling witness: the exact certificate of a flat
//                  realization (a convex polygon doubled along its boundary)
//   AlexandrovSolver / AlexandrovIDTCubic — the public entry points

#include "fullerenes/delaunay_alexandrov.hh"
#include "fullerenes/dense_linalg.hh"
#include "fullerenes/matrix.hh"
#include <cstdio>
#include <cmath>
#include <limits>
#include <map>
#include <numeric>
#include <optional>
#include <set>
#include <span>
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

PyramidDev develop_pyramid(const DelaunayView& T, span<const double> r, int h) {
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
double alpha(const DelaunayView& T, span<const double> r, int h) {
  auto d = develop_pyramid(T, r, h);
  if (!d.base_ok) return NAN;                    // triangle inequality violation
  int u = T.he_origin[h];
  if (d.h_sq < -1e-10 * r[u]*r[u]) return NAN;   // pyramid doesn't close
  double height = sqrt(max(0.0, d.h_sq));

  return atan2(height, d.py);
}

// Total dihedral angle θ at edge h: sum of α from both adjacent pyramids.
double theta(const DelaunayView& T, span<const double> r, int h) {
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
// The scale s of that rule on its own: 1 if r_from + δ ∈ F(T), else
// FEAS_SAFETY · s_max.  feasible_step is s · δ; the optimizer
// framework's step-clip hook consumes s directly (the paradigm applies
// the scaling, in the same order: out[i] = s * delta[i]).
double feasible_fraction(const DelaunayTriangulation& T,
                         const vector<double>& r_from,
                         const vector<double>& delta,
                         bool* clipped = nullptr) {
  vector<double> r_to(r_from.size());
  for (size_t i = 0; i < r_from.size(); i++) r_to[i] = r_from[i] + delta[i];
  double s_max = feasibility_max_step(T, r_from, r_to);
  bool clip = (s_max < 1.0);
  if (clipped) *clipped = clip;
  return clip ? FEAS_SAFETY * s_max : 1.0;
}

vector<double> feasible_step(const DelaunayTriangulation& T,
                              const vector<double>& r_from,
                              const vector<double>& delta,
                              bool* clipped = nullptr) {
  const double s = feasible_fraction(T, r_from, delta, clipped);
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

// The mean edge length and the enclosed signed volume are the acceptance
// gate's own bodies (delaunay_polytope.hh), shared with the device
// validator; named here for the call sites below.
using polytope::mean_edge_length;
using polytope::signed_volume6;

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
  e.r = r;
  e.positions = Reconstruct::from_radii(T, r);   // empty/NaN if Gram-BFS fails
  // Live faces in compact order; compact[f] maps a DCEL face id to its index in
  // e.faces so the gluing below refers to compact indices.
  std::vector<int> compact(T.nf, -1);
  for (int f = 0; f < T.nf; f++) {
    if (T.f_he[f] < 0 || !T.alive(T.f_he[f])) continue;   // dead face slot
    compact[f] = (int)e.faces.size();
    const auto h = T.face_halfedges(f);
    e.faces.push_back({T.he_origin[h[0]], T.he_origin[h[1]], T.he_origin[h[2]]});
    e.face_len.push_back({T.he_length[h[0]], T.he_length[h[1]], T.he_length[h[2]]});
  }
  // The gluing across each base edge, per slot (see TrajEntry): the twin half-edge's
  // face and its cycle slot there, plus the GCP dihedral at the edge.
  for (int f = 0; f < T.nf; f++) {
    if (compact[f] < 0) continue;
    const auto h = T.face_halfedges(f);
    std::array<int,3>    tw;
    std::array<double,3> th;
    for (int s = 0; s < 3; s++) {
      const int t = T.twin(h[s]);
      const int g = T.he_face[t];
      tw[s] = (g >= 0 && compact[g] >= 0) ? 3 * compact[g] + T.cycle_slot(t) : -1;
      th[s] = GCP::theta(T, r, h[s]);
    }
    e.face_twin.push_back(tw);
    e.face_theta.push_back(th);
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

using Repair = AlexandrovSolver::Repair;   // the outcome of a repair, read by every layer below

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

// The unconditional flip loop (AlexandrovSolver::flip_to_weighted_delaunay's
// contract): every bad edge is flipped, a non-finite dihedral included, up
// to the cap.  Not on any solve path.
int flip_to_weighted_delaunay(DelaunayTriangulation& T, const vector<double>& r) {
  int total = 0;
  for (int iter = 0, cap = AlexandrovSolver::flip_cap(T.nh); iter < cap; iter++) {
    int h = needs_flip(T, r);
    if (h < 0) break;
    if (T.flip_edge(h)) total++;
    else break;
  }
  return total;
}

// Repair T to the weighted-Delaunay complex of r by legal flips
// (AlexandrovSolver::repair's contract).  Feasibility is tested before the
// first flip and after every flip, so needs_flip's NaN clause cannot fire
// here: a pyramid that does not close is Infeasible, never a bad edge.
// @variant the B-I piecewise-quadratic extension, which every legal flip
//          strictly increases over a finite set of triangulations: the loop
//          terminates without the cap, which guards a broken predicate only
// @post Delaunay: T weighted-Delaunay for r and r ∈ F(T); Budget: T's
//       status latch reads BudgetExceeded
Repair repair(DelaunayTriangulation& T, const vector<double>& r, int& flips) {
  const int cap = AlexandrovSolver::flip_cap(T.nh);
  for (int iter = 0; iter <= cap; iter++) {
    if (!GCP::feasible(T, r)) return Repair::Infeasible;
    int h = needs_flip(T, r);
    if (h < 0) return Repair::Delaunay;
    if (iter == cap) break;
    if (!T.flip_edge(h)) return Repair::Unflippable;
    flips++;
  }
  T.trip(DelaunayView::Status::BudgetExceeded,
         "Alexandrov repair: flip cap reached with a bad edge left", cap);
  return Repair::Budget;
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
// alongside the (t, r, T) iterate — always present; "accepted" iff
// nit ∈ [0, CORRECTOR_MAX_ITER).  t, r, T are valid only if accepted: T is
// then the weighted-Delaunay complex of r (the trial copy the corrector
// repaired), and `flips` counts the flips that repair applied.
// predictor_refused: the tangent system J·ṙ = κ₁ had no usable solution.
// That system does not contain dt, so halving dt re-runs an identical
// computation; the track stops instead.
struct StepResult {
  int nit;
  double t; V r;
  DelaunayTriangulation T;
  int flips;
  bool predictor_refused;
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

// Newton driving κ(r) → target: r −= J⁻¹(κ(r) − target), where κ is the
// curvature of the generalized convex polyhedron of r.  Returns the
// iteration count: nit ∈ [0, max_iter) converged to tol, max_iter if not,
// −1 when an iterate left the admissible set (the repair did not reach
// Delaunay) or a linear solve failed.  `flips` counts the flips the repairs
// applied.  (Newton::polish is the trust-region κ→0 analogue; its trials
// keep the same invariant.)
// @inv  every κ and J is evaluated on the weighted-Delaunay complex of the
//       iterate it is evaluated at: T is repaired in place before each
//       evaluation (Topology::repair)
// @post nit ∈ [0, max_iter): κ(T, r) = target to tol, T weighted-Delaunay
//       for r and r ∈ F(T)
//
// Evaluating on a stale complex is not a harmless approximation: past a flip
// boundary the stale κ is not the curvature of any convex polyhedron, and a
// corrector that converges there can land where the flipped complex has no
// closing pyramid for the radii at all, with no step size that recovers it.
int newton_correct(DelaunayTriangulation& T, V& r, const V& target,
                   double tol, int max_iter, int& flips) {
  for (int nit = 0; nit < max_iter; nit++) {
    if (Topology::repair(T, r, flips) != Repair::Delaunay) return -1;
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

// One natural-continuation step to t1 = t0 − dt on a TRIAL copy of T: Euler
// predictor r0 − dt·(dr/dt), then the fixed-t1 Newton correction of κ(r) onto
// t1·κ₁ with every iterate on its own weighted-Delaunay complex.  The copy is
// the step's commit-or-nothing unit (the polish's TRIAL invariant): an
// accepted step's (r, T) are adopted together, a rejected step leaves the
// caller's T untouched.  The copy carries any active point tracker with it.
StepResult natural_step(const DelaunayTriangulation& T,
                        double t0, const V& r0, const V& kappa1,
                        const matrix<double>& J, double dt) {
  double t1 = t0 - dt;
  auto v = dr_dt(J, kappa1);
  if (!LinAlg::is_usable_step(v)) return {-1, t1, r0, T, 0, true};
  StepResult s{0, t1, r0 - v * dt, T, 0, false};
  s.nit = newton_correct(s.T, s.r, kappa1 * t1, CORRECTOR_TOL, CORRECTOR_MAX_ITER, s.flips);
  return s;
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
// natural_step predictor-corrector steps, adopting each accepted step's
// repaired complex and adapting dt (clamped against overshooting t_target).
//   Pre:  r ∈ F(T), T weighted-Delaunay for r, 0 < κᵢ(T, r) < δᵢ, t_target > 0.
//   Post: t_final ≤ t_target if the continuation reached it, else it stalled;
//         T is weighted-Delaunay for r and r ∈ F(T) at every accepted step.
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

    // A refused predictor ends the track after this step's records: dt is
    // not in the tangent system, so no halving can change the outcome.
    const bool stop = !result.accepted() && result.predictor_refused;
    if (result.accepted()) {
      t = result.t; r = std::move(result.r);
      T = std::move(result.T);                  // weighted-Delaunay for r
      stats.flips += result.flips;
      record_history(history, t, r);
      dt = adapt_dt(dt, result.nit, CORRECTOR_MAX_ITER, DT_MAX);
    } else if (!stop) {
      dt *= 0.5;
      if (dt < DT_MIN) break;
    }
    if (trace)
      trace->push_back(make_trace('T', step_i, t, dt, result.nit,
                                   GCP::kappa(T, r), J));
    if (diag) {
      auto Jd = GCP::jacobian(T, r);  // J on the adopted (repaired) state
      diag->push_back(make_diag('T', step_i, t, dt, result.nit,
                                  T, r, GCP::kappa(T, r), Jd, stats.flips));
    }
    if (traj)
      traj->push_back(make_traj('T', step_i, t, T, r, GCP::kappa(T, r)));
    stats.steps++;
    stats.newton_total += max(result.nit, 0);
    if (stop) break;
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

  // ENTRY invariant, as a repair: the state must be feasible on its own
  // weighted-Delaunay complex, or there is nothing to polish.
  if (Topology::repair(T, r, flips_local) != Repair::Delaunay) {
    if (flips_cum) *flips_cum = flips_local;
    return {false, max_abs(GCP::kappa(T, r))};
  }

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
    int    trial_flips = 0;
    // A trial outside the admissible set (its repair does not reach
    // Delaunay) has infinite energy: the update below rejects it, and
    // nothing is evaluated there.
    double E_trial   = Topology::repair(T_trial, r_trial, trial_flips) == Repair::Delaunay
                           ? energy(GCP::kappa(T_trial, r_trial))
                           : std::numeric_limits<double>::infinity();

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

// 3D reconstruction from converged (T, r) with κ ≈ 0.  Returns the empty
// vector when a placement refuses: the seed face's two closures (the jy2 and
// kz2 tests below), a later pyramid's closure (place_vertex's NaN), or a
// complex with no live face.
vector<coord3d> from_radii(const DelaunayTriangulation& T, const V& r) {
  int n = T.nv;
  vector<coord3d> pos(n, coord3d(0,0,0));
  vector<bool> placed(n, false), face_done(T.nf, false);

  // Seed face: place vertex i on x-axis, j in xy-plane, k with z ≥ 0.
  int f0 = -1;
  for (int f = 0; f < T.nf; f++) if (T.f_he[f] >= 0) { f0 = f; break; }
  if (f0 < 0) return {};

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
        if (std::isnan(pos[w][0])) return {};   // the placement's own refusal
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

// ============================================================================
// Flat: the doubling witness (the header's alexandrov-doubled-polygon
// paragraph) -- the exact verification that twelve edges of a twelve-cone
// complex of the DUAL metric bound a convex polygon doubled along its
// boundary.  Lattice arithmetic only (eisenstein.hh); the radii enter only
// through the proposal, the twelve edges of smallest dihedral angle.
// ============================================================================

namespace Flat {

using Verdict = AlexandrovSolver::DoublingVerdict;
using Witness = AlexandrovSolver::DoubledPolygon;

// The integer squared length of a live half-edge: the dual metric's edge
// lengths are square roots of Eisenstein norms, so the double rounds to its
// norm (the canonical tesselation rounds the same way).
long long lattice_norm(const DelaunayView& D, int h) {
  return llround(D.he_length[h] * D.he_length[h]);
}

// The apex c on the LEFT of the placed base a -> b of a positively oriented
// lattice triangle with integer squared sides |ab|^2, |ac|^2, |bc|^2, or
// nothing when no lattice triangle has these sides on this base.
optional<Eisenstein> apex(Eisenstein a, Eisenstein b, long long ab2, long long ac2,
                          long long bc2) {
  return place_third_eis_total(a, b, (int)ab2, (int)ac2, (int)bc2, +1);
}

// The development of one sheet from one root: the positions of its faces'
// corner occurrences (indexed by half-edge, the occurrence of the origin)
// in the chart of the sheet's first face, whose first edge is laid along
// `root`, a lattice vector of that edge's norm.  Faces are placed one after
// another across the sheet's interior edges, each apex from the face's own
// integer squared sides; a face reached twice must agree, and every
// occurrence of a cone must land at one position (a cone inside the sheet
// would show as a holonomy).
struct Sheet {
  Verdict verdict = Verdict::NoLatticeDevelopment;
  vector<Eisenstein> pos;     // per half-edge
  vector<char> placed;        // per half-edge
};

Sheet develop_sheet(const DelaunayView& D, span<const int> faces,
                    const vector<char>& on_cycle, Eisenstein root) {
  Sheet S;
  S.pos.assign(D.nh, Eisenstein(0, 0));
  S.placed.assign(D.nh, 0);
  vector<char> in_sheet(D.nf, 0), face_placed(D.nf, 0);
  for (int f : faces) in_sheet[f] = 1;
  {
    // The root face: corner 0 at the origin, corner 1 along the root, the
    // third corner on their left.
    const array<int, 3> hs = D.face_halfedges(faces[0]);
    const Eisenstein p0(0, 0), p1 = root;
    const auto p2 = apex(p0, p1, lattice_norm(D, hs[0]), lattice_norm(D, hs[2]),
                         lattice_norm(D, hs[1]));
    if (!p2) return S;
    S.pos[hs[0]] = p0; S.pos[hs[1]] = p1; S.pos[hs[2]] = *p2;
    for (int h : hs) S.placed[h] = 1;
    face_placed[faces[0]] = 1;
  }
  vector<int> queue{faces[0]};
  for (size_t q = 0; q < queue.size(); q++) {
    const int f = queue[q];
    for (int h : D.face_halfedges(f)) {
      if (on_cycle[D.edge(h)]) continue;
      const int t = D.twin(h), g = D.he_face[t];
      if (!in_sheet[g]) { S.verdict = Verdict::NotTwoSheets; return S; }
      // The twin face's corners: origin(t) = dest(h) at p_j, dest(t) =
      // origin(h) at p_i, and its third corner l on the left of j -> i.
      const Eisenstein p_i = S.pos[h], p_j = S.pos[D.he_next[h]];
      const int nt = D.he_next[t], pt = D.prev(t);
      const auto p_l = apex(p_j, p_i, lattice_norm(D, t), lattice_norm(D, pt),
                            lattice_norm(D, nt));
      if (!p_l) { S.verdict = Verdict::NoLatticeDevelopment; return S; }
      if (face_placed[g]) {
        if (S.pos[t] != p_j || S.pos[nt] != p_i || S.pos[pt] != *p_l) {
          S.verdict = Verdict::SheetNotFlat;   // two routes to a face disagree
          return S;
        }
        continue;
      }
      S.pos[t] = p_j; S.pos[nt] = p_i; S.pos[pt] = *p_l;
      S.placed[t] = S.placed[nt] = S.placed[pt] = 1;
      face_placed[g] = 1;
      queue.push_back(g);
    }
  }
  for (int f : faces)
    if (!face_placed[f]) { S.verdict = Verdict::NotTwoSheets; return S; }   // not connected
  // Every occurrence of a cone in the sheet at one position.
  vector<char> seen(D.nv, 0);
  vector<Eisenstein> at(D.nv, Eisenstein(0, 0));
  for (int h = 0; h < D.nh; h++) {
    if (!S.placed[h]) continue;
    const int v = D.he_origin[h];
    if (!seen[v]) { seen[v] = 1; at[v] = S.pos[h]; }
    else if (at[v] != S.pos[h]) { S.verdict = Verdict::SheetNotFlat; return S; }   // a cone inside
  }
  S.verdict = Verdict::Witness;
  return S;
}

// The developed boundary turns left by exactly 30 degrees at every corner:
// with u, v the consecutive edge vectors, wedge(u, v) > 0, dot2(u, v) > 0
// and dot2(u, v)^2 = 3 |u|^2 |v|^2 (dot2 being twice the inner product, the
// last is cos^2 = 3/4).
bool turns_by_thirty_degrees(const array<Eisenstein, 12>& p) {
  for (int k = 0; k < 12; k++) {
    const Eisenstein u = p[(k + 1) % 12] - p[k], v = p[(k + 2) % 12] - p[(k + 1) % 12];
    const long long w = wedge(u, v), d = dot2(u, v);
    const long long uu = u.norm2(), vv = v.norm2();
    if (!(w > 0) || !(d > 0) || d * d != 3 * uu * vv) return false;
  }
  return true;
}

// The verification of a proposal: twelve edge indices (edge e is the pair
// of half-edges 2e, 2e+1).
Witness verify(const DelaunayView& D, span<const int, 12> edges) {
  Witness out;
  if (D.nv != 12) { out.verdict = Verdict::NotTwelveCones; return out; }
  const int ne = D.nh / 2;
  vector<char> on_cycle(ne, 0);
  vector<int> degree(D.nv, 0);
  for (int e : edges) {
    if (e < 0 || e >= ne || on_cycle[e]) { out.verdict = Verdict::NotTwelveEdges; return out; }
    const int h = 2 * e;
    if (!D.alive(h) || D.he_origin[h] == D.dest(h)) {
      out.verdict = Verdict::NotTwelveEdges;   // dead, or a loop cannot bound a polygon
      return out;
    }
    on_cycle[e] = 1;
    degree[D.he_origin[h]]++;
    degree[D.dest(h)]++;
  }
  for (int v = 0; v < D.nv; v++)
    if (degree[v] != 2) { out.verdict = Verdict::NotOneCycle; return out; }
  // Every live edge length a lattice norm: the dual metric.
  for (int h = 0; h < D.nh; h += 2) {
    if (!D.alive(h)) continue;
    const long long n = lattice_norm(D, h);
    if (n <= 0 || n > numeric_limits<int>::max() || !first_rep_of_norm((int)n)) {
      out.verdict = Verdict::LengthNotANorm;
      return out;
    }
  }

  // The two sheets: the components of the live faces across the non-cycle
  // edges.
  vector<int> component(D.nf, -1);
  int components = 0;
  for (int f0 = 0; f0 < D.nf; f0++) {
    if (D.f_he[f0] < 0 || component[f0] >= 0) continue;
    if (components == 2) { out.verdict = Verdict::NotTwoSheets; return out; }
    vector<int> queue{f0};
    component[f0] = components;
    for (size_t q = 0; q < queue.size(); q++)
      for (int h : D.face_halfedges(queue[q])) {
        if (on_cycle[D.edge(h)]) continue;
        const int g = D.he_face[D.twin(h)];
        if (component[g] < 0) { component[g] = components; queue.push_back(g); }
      }
    components++;
  }
  if (components != 2) { out.verdict = Verdict::NotTwoSheets; return out; }
  vector<int> a, b;
  for (int f = 0; f < D.nf; f++)
    if (D.f_he[f] >= 0) (component[f] == 0 ? a : b).push_back(f);
  if (a.size() != 10 || b.size() != 10) { out.verdict = Verdict::NotTwoSheets; return out; }
  copy(a.begin(), a.end(), out.sheet_a.begin());
  copy(b.begin(), b.end(), out.sheet_b.begin());

  // The boundary of the first sheet, oriented: the cycle half-edges whose
  // face lies in it (a face is on the left of its half-edges), chained.
  // Twelve cones of degree two on twelve edges form one cycle iff the walk
  // below closes after twelve steps through twelve distinct cones.
  vector<int> leaving(D.nv, -1);
  for (int e : edges)
    for (int h : {2 * e, 2 * e + 1})
      if (component[D.he_face[h]] == 0) leaving[D.he_origin[h]] = h;
  int h = leaving[D.he_origin[2 * edges[0]]];
  if (h < 0) h = leaving[D.dest(2 * edges[0])];
  for (int k = 0; k < 12; k++) {
    if (h < 0) { out.verdict = Verdict::NotOneCycle; return out; }
    out.boundary[k] = h;
    out.corner[k] = D.he_origin[h];
    h = leaving[D.dest(h)];
  }
  if (h != out.boundary[0]) { out.verdict = Verdict::NotOneCycle; return out; }
  for (int k = 0; k < 12; k++)
    for (int l = 0; l < k; l++)
      if (out.corner[k] == out.corner[l]) { out.verdict = Verdict::NotOneCycle; return out; }

  // Both sheets develop as flat disks, the root edge's unit-orbit
  // representatives tried in turn: a representative of the wrong orbit
  // fails to place some apex on the lattice and is passed over, while a
  // holonomy is a property of the sheet and ends the search.
  auto develop = [&](span<const int> faces) {
    const int h0 = D.face_halfedges(faces[0])[0];
    Sheet S;
    for (Eisenstein root : Sector0Reps((int)lattice_norm(D, h0))) {
      S = develop_sheet(D, faces, on_cycle, root);
      if (S.verdict != Verdict::NoLatticeDevelopment) break;
    }
    return S;
  };
  const Sheet SA = develop(out.sheet_a);
  if (SA.verdict != Verdict::Witness) { out.verdict = SA.verdict; return out; }
  const Sheet SB = develop(out.sheet_b);
  if (SB.verdict != Verdict::Witness) { out.verdict = SB.verdict; return out; }

  // Both boundaries turn by 30 degrees at every corner.  The second sheet's
  // boundary runs the cycle the other way: read it with that sheet on its
  // left, from the twins in reverse order.
  for (int k = 0; k < 12; k++) out.polygon[k] = SA.pos[out.boundary[k]] - SA.pos[out.boundary[0]];
  if (!turns_by_thirty_degrees(out.polygon)) { out.verdict = Verdict::NotThirtyDegrees; return out; }
  array<Eisenstein, 12> pb{};
  for (int k = 0; k < 12; k++) pb[k] = SB.pos[D.twin(out.boundary[(12 - k) % 12])];
  const Eisenstein origin_b = pb[0];
  for (int k = 0; k < 12; k++) pb[k] = pb[k] - origin_b;
  if (!turns_by_thirty_degrees(pb)) { out.verdict = Verdict::NotThirtyDegrees; return out; }
  // The area, in unit triangles: the sum of wedges around the boundary.
  for (int k = 0; k < 12; k++) out.area += wedge(out.polygon[k], out.polygon[(k + 1) % 12]);
  out.verdict = Verdict::Witness;
  return out;
}

// The dihedral at a base edge, for the PROPOSAL only: the flat limit's
// pyramids have no height, and a state reached at the floor of its
// arithmetic misses zero by that arithmetic's resolution, so a pyramid
// that fails to close by rounding counts as flat, its dihedral 0 or pi by
// the side of the apex's projection (GCP::alpha would refuse it).  A base
// that fails the triangle inequality has no dihedral (non-finite, sorted
// last).  The verdict never reads this.
double proposal_alpha(const DelaunayView& D, span<const double> r, int h) {
  const GCP::PyramidDev d = GCP::develop_pyramid(D, r, h);
  if (!d.base_ok) return numeric_limits<double>::quiet_NaN();
  return atan2(sqrt(max(0.0, d.h_sq)), d.py);
}

// ── The SEARCH: every fold the complex contains, verified.  Flatness is a
//    property of the metric (the header's alexandrov-doubled-polygon-search
//    paragraph), so no radii enter.  Backtracking over the cones' incident
//    edges from cone 0; a delta-complex may join two cones by several edges,
//    so the walk carries EDGE indices, never vertex pairs.  Each cone is
//    entered once, so the walk is bounded by the twelve cones' incidences.
// ──
struct FoldSearch {
  const DelaunayView& D;
  struct Inc { int edge, to; };
  vector<vector<Inc>> inc;
  vector<char> seen;
  array<int, 12> edges{};
  Witness out;
  bool done = false;

  explicit FoldSearch(const DelaunayView& D_) : D(D_), inc(D_.nv), seen(D_.nv, 0) {
    for (int h = 0; h < D.nh; h += 2) {
      if (!D.alive(h)) continue;
      const int u = D.he_origin[h], v = D.dest(h);
      if (u == v) continue;                       // a loop cannot bound a polygon
      inc[u].push_back({h / 2, v});
      inc[v].push_back({h / 2, u});
    }
  }

  // Depth is the number of cones already on the path (cone 0 included).
  void walk(int v, int depth) {
    if (done) return;
    for (const Inc& c : inc[v]) {
      if (done) return;
      if (depth == 12) {                          // the twelfth edge closes at cone 0
        if (c.to != 0 || c.edge == edges[0]) continue;
        edges[11] = c.edge;
        const Witness w = verify(D, span<const int, 12>(edges.data(), 12));
        if (w.ok()) { out = w; done = true; }
        continue;
      }
      if (seen[c.to]) continue;
      edges[depth - 1] = c.edge;
      seen[c.to] = 1;
      walk(c.to, depth + 1);
      seen[c.to] = 0;
    }
  }

  Witness run() {
    out.verdict = Verdict::NoFoldFound;
    if (D.nv != 12) { out.verdict = Verdict::NotTwelveCones; return out; }
    seen[0] = 1;
    walk(0, 1);
    return out;
  }
};

// The proposal from a state, verified: the twelve live non-loop edges of
// smallest dihedral angle, a non-finite angle sorting last.
Witness proposal(const DelaunayView& D, span<const double> r) {
  vector<pair<double, int>> order;
  for (int h = 0; h < D.nh; h += 2) {
    if (!D.alive(h) || D.he_origin[h] == D.dest(h)) continue;
    const double th = proposal_alpha(D, r, h) + proposal_alpha(D, r, D.twin(h));
    order.push_back({isfinite(th) ? th : numeric_limits<double>::infinity(), D.edge(h)});
  }
  stable_sort(order.begin(), order.end(),
              [](const pair<double, int>& x, const pair<double, int>& y) { return x.first < y.first; });
  Witness out;
  if (order.size() < 12) { out.verdict = Verdict::NotTwelveEdges; return out; }
  array<int, 12> edges{};
  for (int k = 0; k < 12; k++) edges[k] = order[k].second;
  return verify(D, edges);
}

} // namespace Flat

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
AlexandrovSolver::ValidationStatus AlexandrovSolver::validate_polytope(
    const DelaunayTriangulation& D, const vector<double>& r,
    const vector<coord3d>& pos, PolytopeValidation* out, bool verbose) {
  using VS = AlexandrovSolver::ValidationStatus;
  // Public entry: r/pos are caller-supplied (the instance path's were
  // solver-sized by construction), and the ladder indexes both by D.nv —
  // enforce the @pre up front.
  if ((int)r.size() < D.nv || (int)pos.size() != D.nv)
    throw std::invalid_argument(
        "AlexandrovSolver::validate_polytope: r.size() >= T.nv and "
        "pos.size() == T.nv required (r " + std::to_string(r.size()) +
        ", pos " + std::to_string(pos.size()) + ", nv " +
        std::to_string(D.nv) + ")");

  // THE GATE IS polytope::validate (delaunay_polytope.hh) — the same bodies
  // a device batch runs over one isomer's workspace.  All this wrapper adds
  // is the interior-edge predicate in the host's arithmetic (an edge of the
  // 2-skeleton's interior is one whose dihedral is π, bigons excepted: the
  // inessential mask of bi-inessential-edges, here as a predicate so no
  // per-edge vector is built), the record's shape, and the diagnostics.
  const auto tight = [&](int h) {
    if (!D.alive(h) || D.is_bigon(h)) return false;
    const double th = GCP::theta(D, r, h);
    return std::isfinite(th) && std::fabs(th - M_PI) < 1e-7;
  };
  polytope::Record rec;
  const polytope::Verdict verdict = polytope::validate(D, pos, tight, &rec);

  PolytopeValidation v;
  v.t0_simplicial         = AlexandrovSolver::is_simplicial(D);  // diagnostic only
  v.tbar_simple_polygonal = rec.simple;
  v.tbar_n_cells          = rec.n_cells;
  v.volume_norm           = rec.volume_norm;
  v.no_self_intersect     = rec.no_self_cross;
  v.convex                = rec.convex;
  if (out) *out = v;

  switch (verdict) {
    case polytope::Verdict::WalkUnclosed:
      // A cell boundary that does not close is a corrupt complex, not a
      // verdict on a polytope: loud, as it has always been (a silently
      // empty tesselation would compare equal to a legitimately empty one).
      throw std::runtime_error(
          "AlexandrovSolver::validate_polytope: the cell-boundary walk from "
          "half-edge " + std::to_string(rec.witness) + " did not close — a "
          "well-formed tesselation always does");
    case polytope::Verdict::NotSimple:
      if (verbose)
        printf("  VALIDATION (simplicity) failed: T̄(0) simple_polygonal=%d, "
               "n_cells=%d (need F≥3), T(0) simplicial=%d (diagnostic).\n",
               v.tbar_simple_polygonal, v.tbar_n_cells, v.t0_simplicial);
      return VS::FAIL_NOT_SIMPLE;
    case polytope::Verdict::VolumeDegenerate:
      if (verbose)
        printf("  VALIDATION (well-formedness/volume) failed: vol_norm=%.3e "
               "(threshold %.2e).\n", v.volume_norm, polytope::kVolumeFloor);
      return VS::FAIL_VOLUME_DEGENERATE;
    case polytope::Verdict::SelfIntersecting:
      if (verbose)
        printf("  VALIDATION (well-formedness/self-intersection) failed: "
               "two non-adjacent triangles cross in 3D.\n");
      return VS::FAIL_SELF_INTERSECTING;
    case polytope::Verdict::NotConvex:
      if (verbose)
        printf("  VALIDATION (convexity) failed: some vertex sticks out "
               "beyond a face plane.\n");
      return VS::FAIL_NOT_CONVEX;
    case polytope::Verdict::Ok:
      return VS::OK;
  }
  return VS::FAIL_NOT_SIMPLE;   // unreachable: the switch is closed
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

  // 5. Validation: the three-property gate (the public static — one
  //    ladder for this instance path and for external (T, r, pos) callers),
  //    with the per-check record copied into the stats_* diagnostics.
  {
    PolytopeValidation pv;
    stats_status = validate_polytope(D, r, pos, &pv, verbose);
    stats_t0_simplicial            = pv.t0_simplicial;
    stats_tbar_n_cells             = pv.tbar_n_cells;
    stats_tbar_simple_polygonal    = pv.tbar_simple_polygonal;
    stats_volume_norm              = pv.volume_norm;
    stats_polytope_no_self_intersect = pv.no_self_intersect;
    stats_polytope_convex          = pv.convex;
  }
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

double AlexandrovSolver::feasible_fraction(const DelaunayTriangulation& T,
                                           const vector<double>& r,
                                           const vector<double>& delta,
                                           bool* clipped) {
  return GCP::feasible_fraction(T, r, delta, clipped);
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

AlexandrovSolver::Repair AlexandrovSolver::repair(DelaunayTriangulation& T,
                                                  const vector<double>& r, int& flips) {
  return Topology::repair(T, r, flips);
}

double AlexandrovSolver::theta(const DelaunayTriangulation& T,
                                const vector<double>& r, int h) {
  return GCP::theta(T, r, h);
}

const char* AlexandrovSolver::doubling_verdict_str(DoublingVerdict v) {
  switch (v) {
    case DoublingVerdict::Witness:              return "doubled polygon";
    case DoublingVerdict::NotTwelveCones:       return "not a twelve-cone complex";
    case DoublingVerdict::NotTwelveEdges:       return "not twelve distinct live non-loop edges";
    case DoublingVerdict::NotOneCycle:          return "not one cycle through the twelve cones";
    case DoublingVerdict::NotTwoSheets:         return "the cycle does not cut two sheets of ten faces";
    case DoublingVerdict::LengthNotANorm:       return "an edge length is not a lattice norm";
    case DoublingVerdict::NoLatticeDevelopment: return "a sheet has no lattice development";
    case DoublingVerdict::SheetNotFlat:         return "a sheet is not a flat disk";
    case DoublingVerdict::NotThirtyDegrees:     return "a corner does not turn by 30 degrees";
    case DoublingVerdict::NoFoldFound:          return "no fold of this complex verifies";
  }
  return "?";
}

AlexandrovSolver::DoubledPolygon AlexandrovSolver::doubled_polygon(const DelaunayView& D) {
  return Flat::FoldSearch(D).run();
}

AlexandrovSolver::DoubledPolygon AlexandrovSolver::doubled_polygon(const DelaunayView& D,
                                                                   span<const double> r) {
  return Flat::proposal(D, r);
}

AlexandrovSolver::DoubledPolygon AlexandrovSolver::verify_doubled_polygon(
    const DelaunayView& D, span<const int, 12> edges) {
  return Flat::verify(D, edges);
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
  // The gate's own body (delaunay_polytope.hh), shared with the device
  // validator; the size check is this entry point's (the body indexes pos
  // by T's vertex count).
  if ((int)pos.size() != T.nv) return false;
  return polytope::is_convex(T, pos, tol);
}

bool AlexandrovSolver::has_self_intersection(const DelaunayTriangulation& T,
                                                const vector<coord3d>& pos,
                                                double tol) {
  if ((int)pos.size() != T.nv) return false;
  return polytope::has_self_intersection(T, pos, tol);
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

AlexandrovIDTCubic::ConeLabels
AlexandrovIDTCubic::cone_labels(const KisMetric& M, const DelaunayTriangulation& D,
                                const vector<int>& new_to_old)
{
  ConeLabels L;
  L.kis_vertex = new_to_old;
  for (int i = 0; i < D.nv; i++) {
    const int old = new_to_old[i];
    if (old < M.Nv)
      throw logic_error("AlexandrovIDTCubic: face center " + to_string(old) +
                        " survived flat removal");
    const tri_t& f = M.triangle[old - M.Nv];
    int k = 0;
    for (int j = 0; j < 3; j++) k += (M.fdeg[f[j]] == 5);
    L.triangle.push_back(f);
    L.npent.push_back(k);
  }
  return L;
}

void AlexandrovIDTCubic::check_curvature_quantum(const DelaunayTriangulation& D,
                                                 const ConeLabels& L)
{
  double total = 0;
  for (int i = 0; i < D.nv; i++) {
    const int k = L.npent[i];
    const double kappa = 2 * M_PI - D.v_cone_angle[i];
    if (k < 1 || fabs(kappa - k * M_PI / 15) > KAPPA_TOL)
      throw logic_error("AlexandrovIDTCubic: cone " + to_string(i) +
                        " has kappa = " + to_string(kappa) + ", expected " +
                        to_string(k) + "*pi/15");
    total += kappa;
  }
  if (fabs(total - 4 * M_PI) > TOTAL_KAPPA_TOL)
    throw logic_error("AlexandrovIDTCubic: total curvature " +
                      to_string(total) + " != 4*pi");
}

void AlexandrovIDTCubic::set_cones(ConeLabels L)
{
  cone_triangle   = std::move(L.triangle);
  cone_npent      = std::move(L.npent);
  cone_kis_vertex = std::move(L.kis_vertex);
}

void AlexandrovIDTCubic::build_banded(const TriangulationView& T)
{
  KisMetric M = kis_metric(T);
  vector<int> new_to_old;
  DelaunayTriangulation D = DelaunayTriangulation::compute(
      M.K, M.edge_length_fn(), FLAT_TOL, &new_to_old, track_removed);
  ConeLabels L = cone_labels(M, D, new_to_old);
  check_curvature_quantum(D, L);
  set_cones(std::move(L));
  solver.D = std::move(D);
}

DelaunayView::CompletionStats AlexandrovIDTCubic::build(const TriangulationView& T)
{
  KisMetric M = kis_metric(T);
  DelaunayTriangulation D =
      DelaunayTriangulation::from_intrinsic_metric(M.K, M.edge_length_fn());
  // Tracking, when requested, precedes the reduction (the compute(...,
  // track_removed) contract); remove_and_complete_cyclotomic_kis passes the
  // point tracker through every flip and every star retriangulation, the
  // completion included.
  if (track_removed) D.enable_point_tracking();
  const DelaunayView::CompletionStats completion =
      D.remove_and_complete_cyclotomic_kis(M.Nv, std::span<const int>(M.fdeg));
  // Compaction after the completion: the exact lengths live only inside
  // that call, and compaction is a monotone relabelling (compact_vertices,
  // @post monotone), under which the completion's corner order is unchanged.
  const vector<int> new_to_old = D.compact_vertices();
  ConeLabels L = cone_labels(M, D, new_to_old);
  check_curvature_quantum(D, L);
  set_cones(std::move(L));
  solver.D = std::move(D);
  return completion;
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
