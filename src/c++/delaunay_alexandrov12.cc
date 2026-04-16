// Bobenko-Izmestiev algorithm for the Alexandrov embedding of a
// polyhedral metric on S² with 12 cone points (fullerene iDT).
//
// Given an intrinsic Delaunay triangulation T with geodesic edge lengths,
// find the unique convex polyhedron P ⊂ R³ whose boundary is isometric to T.
//
// The algorithm parameterizes P by the radii r = (|a−v_i|), where a is an
// interior apex point.  The curvature κ_i = 2π − ω_i (angle deficit at
// the radial edge a−v_i) satisfies κ = 0 iff the GCP is a genuine polytope.
//
// We trace the homotopy κ(r) = t·κ₁ from (t=1, large R) to (t→0, r*)
// using pseudo-arc-length continuation (PALC), then extrapolate and polish.
//
// Layers:
//   GCP          — observables: κ(T,r), J(T,r), θ(T,r,h)
//   LinAlg       — vector arithmetic, linear solvers
//   TrustRegion  — LM subproblem solver, accept/reject rule
//   Topology     — weighted Delaunay flip maintenance
//   PALC         — predictor-corrector arc-length continuation
//   Newton       — trust-region Newton polish for κ(r)=0
//   Reconstruct  — BFS vertex placement from Gram matrix entries
//   AlexandrovSolver — the 5-step algorithm

#include "fullerenes/delaunay_alexandrov12.hh"
#include "fullerenes/matrix.hh"
#include <cstdio>
#include <cmath>
#include <vector>
#include <algorithm>
#include <expected>

using namespace std;

extern "C" {
  void dgesv_(const int* n, const int* nrhs, double* A, const int* lda,
              int* ipiv, double* b, const int* ldb, int* info);
  void dsyev_(const char* jobz, const char* uplo, const int* n,
              double* A, const int* lda, double* w,
              double* work, const int* lwork, int* info);
}

namespace {

// Full symmetric eigenvalue decomposition of a 12×12-ish matrix.
// Returns eigenvalues sorted ascending.  On LAPACK failure, returns empty.
// Used by both PALC step tracing and Newton polish tracing.
static std::vector<double> sym_eigvals(const matrix<double>& A) {
  int n = A.m;
  std::vector<double> Af(n*n);
  for (int i = 0; i < n; i++)
    for (int j = 0; j < n; j++)
      Af[j*n + i] = 0.5 * (A(i,j) + A(j,i));
  std::vector<double> w(n), work(1);
  int lwork = -1, info;
  dsyev_("N", "U", &n, Af.data(), &n, w.data(), work.data(), &lwork, &info);
  if (info != 0) return {};
  lwork = (int)work[0]; work.assign(lwork, 0.0);
  dsyev_("N", "U", &n, Af.data(), &n, w.data(), work.data(), &lwork, &info);
  if (info != 0) return {};
  return w;
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

// Dihedral angle α at base edge h in the pyramid over face(h).
//
// The pyramid has apex a and base triangle (u,v,w) with |a-u|=r_u, etc.
// Develop the base into the plane: u at origin, v on the x-axis.
// Project a onto this plane: p = (px, py) satisfies |p-u|²+h²=r_u², etc.
// Then α = atan2(h, py), where h is the pyramid height and py is the
// signed distance from the projection to edge uv (the x-axis).
double alpha(const DelaunayTriangulation& T, const vector<double>& r, int h) {
  int u = T.he_origin[h], v = T.he_origin[T.he_next[h]], w = T.he_origin[T.prev(h)];
  double luv = T.he_length[h], lvw = T.he_length[T.he_next[h]], lwu = T.he_length[T.prev(h)];

  // w in the base plane (u at origin, v at (luv,0)):  w = (wx, wy)
  double wx = (luv*luv + lwu*lwu - lvw*lvw) / (2*luv);  // cosine rule
  double wy_sq = lwu*lwu - wx*wx;
  if (wy_sq < -1e-10 * lwu*lwu) return NAN;  // triangle inequality violation
  double wy = sqrt(max(0.0, wy_sq));

  // Apex projection p = (px, py):  |p-u|²+h²=r_u², |p-v|²+h²=r_v², |p-w|²+h²=r_w²
  double px = (r[u]*r[u] - r[v]*r[v] + luv*luv) / (2*luv);
  double py = (wy > 1e-15) ? (r[u]*r[u] - r[w]*r[w] + wx*wx + wy*wy - 2*px*wx) / (2*wy) : 0;
  double h_sq = r[u]*r[u] - px*px - py*py;
  if (h_sq < -1e-10 * r[u]*r[u]) return NAN;  // pyramid doesn't close
  double height = sqrt(max(0.0, h_sq));

  return atan2(height, py);
}

// Total dihedral angle θ at edge h: sum of α from both adjacent pyramids.
double theta(const DelaunayTriangulation& T, const vector<double>& r, int h) {
  return alpha(T, r, h) + alpha(T, r, h ^ 1);
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
  int h0 = T.v_out[v], h = h0;
  do {
    double rho_vj = rho(r, T.he_length[h], v, T.dest(h));
    double rho_vk = rho(r, T.he_length[T.prev(h)], v, T.dest(T.he_next[h]));
    double sj = sin(rho_vj), sk = sin(rho_vk);
    if (sj < 1e-15 || sk < 1e-15) { omega = NAN; break; }
    double cos_omega = (cos(T.he_angle[h]) - cos(rho_vj)*cos(rho_vk)) / (sj*sk);
    omega += acos(clamp(cos_omega, -1.0, 1.0));
    h = T.cw(h);
  } while (h != h0);
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
  double alpha_me = alpha(T, r, h ^ 1);

  double sr = sin(rho_e), srm = sin(rho_me);
  double sa = sin(alpha_e), sam = sin(alpha_me);
  if (sr < 1e-15 || srm < 1e-15 || sa < 1e-15 || sam < 1e-15) return NAN;

  return (cos(alpha_e)/sa + cos(alpha_me)/sam) / (L * sr * srm);
}

// Jacobian J(T, r) = ∂κ/∂r.
// Off-diagonal: J(i,j) = Σ_{edges i→j} J_e(h)
// Diagonal:     J(i,i) = −Σ_{e: a(e)=i} cos(φ_e) · J_e(h)
// (per-oriented-edge formula, correct for multi-edges)
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
    int h0 = T.v_out[i], h = h0;
    do {
      double Je = J_edge(T, r, h);
      double phi_e = phi(r, T.he_length[h], i, T.dest(h));
      diag -= cos(phi_e) * Je;
      h = T.cw(h);
    } while (h != h0);
    J(i, i) = diag;
  }
  return J;
}

} // namespace GCP

// ============================================================================
// Layer 2: Linear algebra
// ============================================================================

namespace LinAlg {

// --- Vector arithmetic ---
using V = vector<double>;

V operator+(const V& a, const V& b)   { V c(a.size()); for (size_t i=0;i<a.size();i++) c[i]=a[i]+b[i]; return c; }
V operator-(const V& a, const V& b)   { V c(a.size()); for (size_t i=0;i<a.size();i++) c[i]=a[i]-b[i]; return c; }
V operator*(const V& a, double s)     { V c(a.size()); for (size_t i=0;i<a.size();i++) c[i]=a[i]*s;     return c; }
V operator*(double s, const V& a)     { return a * s; }
V operator-(const V& a)               { return a * (-1.0); }

double dot(const V& a, const V& b) { double s=0; for (size_t i=0;i<a.size();i++) s+=a[i]*b[i]; return s; }
double norm(const V& v)            { return sqrt(dot(v, v)); }
double max_abs(const V& v)         { double m=0; for (double x:v) { if (!(x==x)) return HUGE_VAL; m=max(m,fabs(x)); } return m; }
double sum_sq(const V& v)          { return dot(v, v); }
bool   is_valid(const V& v)        { double n=sum_sq(v); return isfinite(n) && n > 1e-30; }

// Residual energy E = ½||v||²
double energy(const V& v) { return 0.5 * sum_sq(v); }

// --- Linear solvers ---

// LU-with-sign solve: returns (x, det_sign) on success.
// det_sign ∈ {−1, +1} is sign(det A), computed as the product of:
//   (−1)^{# row swaps from ipiv}  ×  product of signs of diag(U).
// Failures are tagged: Singular for a zero pivot, LapackError for bad info.
struct LuSolved { V x; int det_sign; };
enum class LuFail { Singular, LapackError };

std::expected<LuSolved, LuFail>
solve_with_sign(const matrix<double>& A, const V& b) {
  int n = A.m;
  vector<double> Af(n*n);
  V x(b);
  for (int i = 0; i < n; i++)
    for (int j = 0; j < n; j++)
      Af[j*n + i] = A(i, j);
  vector<int> ipiv(n);
  int nrhs = 1, info;
  dgesv_(&n, &nrhs, Af.data(), &n, ipiv.data(), x.data(), &n, &info);
  if (info < 0) return std::unexpected(LuFail::LapackError);
  if (info > 0) return std::unexpected(LuFail::Singular);
  int sign = 1;
  for (int i = 0; i < n; i++) {
    if (ipiv[i] != i + 1) sign = -sign;       // row swap parity
    double d = Af[i*n + i];                   // diagonal of U
    if (d < 0)       sign = -sign;
    else if (d == 0) return std::unexpected(LuFail::Singular);
  }
  return LuSolved{std::move(x), sign};
}

// Solve A·x = b via LAPACK.  Returns zero vector on failure.
// (Backward-compat wrapper around solve_with_sign that discards the sign.)
V solve(const matrix<double>& A, const V& b) {
  auto r = solve_with_sign(A, b);
  return r ? std::move(r->x) : V(A.m, 0.0);
}

// Solve (A + λI)·x = b.
V solve_shifted(const matrix<double>& A, const V& b, double lambda) {
  matrix<double> Al = A;
  for (int i = 0; i < A.m; i++) Al(i,i) += lambda;
  return solve(Al, b);
}

// Matrix-vector product A·v (for predicted reduction).
V matvec(const matrix<double>& A, const V& v) {
  int n = A.m;
  V r(n, 0.0);
  for (int i = 0; i < n; i++)
    for (int j = 0; j < n; j++)
      r[i] += A(i,j) * v[j];
  return r;
}

} // namespace LinAlg

// Bring vector ops into scope for readability
using LinAlg::V;

// Pack one PALC- or Newton-step diagnostic into a TraceEntry.  Caller
// decides whether to record (avoids the cost of κ and J spectrum when
// trace_jacobian is off); this just bundles the fields.
static AlexandrovSolver::TraceEntry make_trace(
    char phase, int step, double t, double ds, int nit,
    const V& kappa, const matrix<double>& J) {
  return {phase, step, t, ds, nit,
          LinAlg::max_abs(kappa), LinAlg::norm(kappa), sym_eigvals(J)};
}

// ============================================================================
// Layer 3: Trust-region subproblem
// ============================================================================

namespace TrustRegion {

// Solve the LM subproblem: find δ such that (J+λI)δ = −κ and ||δ|| ≤ Δ.
// Bisects on λ ≥ 0.
V solve(const matrix<double>& J, const V& kappa, double Delta) {
  // Try pure Newton (λ=0)
  auto delta = LinAlg::solve(J, -kappa);
  if (LinAlg::is_valid(delta) && LinAlg::norm(delta) <= Delta)
    return delta;

  // Bisect on λ to find (J+λI)⁻¹(-κ) with ||δ|| ≈ Δ
  double lo = 0, hi = LinAlg::max_abs(kappa) / Delta + 1.0;
  for (int probe = 0; probe < 10; probe++) {
    delta = LinAlg::solve_shifted(J, -kappa, hi);
    if (LinAlg::is_valid(delta) && LinAlg::norm(delta) <= Delta) break;
    hi *= 4;
  }
  for (int bis = 0; bis < 20; bis++) {
    double mid = 0.5 * (lo + hi);
    delta = LinAlg::solve_shifted(J, -kappa, mid);
    if (!LinAlg::is_valid(delta) || LinAlg::norm(delta) > Delta) lo = mid;
    else hi = mid;
  }
  return LinAlg::solve_shifted(J, -kappa, hi);
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

// Find first edge with θ > π that can be flipped (skip B==D multi-edges).
int needs_flip(const DelaunayTriangulation& T, const vector<double>& r) {
  for (int h = 0; h < T.nh; h += 2) {
    if (!T.alive(h)) continue;
    if (T.he_origin[T.prev(h)] == T.he_origin[T.prev(h ^ 1)]) continue;
    if (GCP::theta(T, r, h) > M_PI + 1e-10) return h;
  }
  return -1;
}

// Flip all edges where θ > π.  Returns number of flips performed.
int flip_to_delaunay(DelaunayTriangulation& T, const vector<double>& r) {
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
// Layer 5: Pseudo-arc-length continuation (Keller 1977)
//
// Traces the homotopy curve F(t,r) = κ(r) − t·κ₁ = 0 using predictor-
// corrector steps parameterized by arc length.  The bordering algorithm
// solves the augmented (n+1)×(n+1) system using only n×n solves.
// ============================================================================

namespace PALC {

struct Tangent { double t_dot; vector<double> r_dot; };

// Tangent to the homotopy curve at (t, r).
// From J·ṙ = κ₁·ṫ and ||(ṫ, ṙ)|| = 1, ṫ < 0.
Tangent compute_tangent(const matrix<double>& J, const vector<double>& kappa1) {
  auto v = LinAlg::solve(J, kappa1);  // v = J⁻¹ κ₁
  if (!LinAlg::is_valid(v))
    return {-1.0, vector<double>(kappa1.size(), 0.0)};
  double vn = LinAlg::norm(v);
  double td = -1.0 / sqrt(1.0 + vn * vn);
  vector<double> rd(v.size());
  for (size_t i = 0; i < v.size(); i++) rd[i] = v[i] * td;
  return {td, rd};
}

// Predictor: Euler step along tangent by arc length ds.
pair<double, V> predict(double t, const V& r, const Tangent& tau, double ds) {
  return { t + tau.t_dot * ds,  r + tau.r_dot * ds };
}

// Corrector: Newton on the bordered system using the bordering algorithm.
//
//   [ J    −κ₁ ] [δr]   [−F                        ]
//   [ ṙ₀ᵀ   ṫ₀] [δt] = [ds − ṙ₀ᵀ(r−r₀) − ṫ₀(t−t₀)]
//
// Solved via two n×n solves per iteration:
//   z₁ = J⁻¹(−F),  z₂ = J⁻¹(κ₁)
//   δt = (g − ṙ₀·z₁) / (ṫ₀ + ṙ₀·z₂)     where g = ds − ṙ₀·(r−r₀) − ṫ₀(t−t₀)
//   δr = z₁ + δt·z₂
//
// Returns number of iterations (−1 on failure).
int correct(const DelaunayTriangulation& T, double& t, V& r,
            const V& kappa1, const Tangent& tau0,
            double t0, const V& r0, double ds,
            int max_iter, double tol) {
  using LinAlg::dot;
  for (int nit = 0; nit < max_iter; nit++) {
    V F = GCP::kappa(T, r) - kappa1 * t;                    // residual
    if (LinAlg::max_abs(F) < tol) return nit;

    auto J = GCP::jacobian(T, r);
    auto z1 = LinAlg::solve(J, -F);                          // J·z₁ = −F
    auto z2 = LinAlg::solve(J, kappa1);                      // J·z₂ = κ₁
    if (!LinAlg::is_valid(z1) || !LinAlg::is_valid(z2)) return -1;

    double g   = ds - dot(tau0.r_dot, r - r0) - tau0.t_dot * (t - t0);
    double den = tau0.t_dot + dot(tau0.r_dot, z2);
    if (fabs(den) < 1e-30) return -1;

    double dt = (g - dot(tau0.r_dot, z1)) / den;
    V dr = z1 + z2 * dt;

    t += dt;
    r = r + dr;
  }
  return max_iter;
}

// Adapt step size based on corrector iterations (AUTO strategy).
double adapt_ds(double ds, int nit, int max_nit, double ds_max) {
  if (nit <= 1) ds *= 2.0;
  else if (nit == 2) ds *= 1.5;
  else if (nit <= max_nit / 2) ds *= 1.1;
  else if (nit >= max_nit) ds *= 0.5;
  return min(fabs(ds), ds_max);
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

} // namespace PALC

// ============================================================================
// Layer 6a: Trust-region Newton for κ(r) = 0
// ============================================================================

namespace Newton {

// Minimize E = ½||κ||² using LM trust-region Newton.
// Returns (converged, final_kappa).
// If `out_trace` is non-null, records one TraceEntry per iteration.
pair<bool, double> polish(DelaunayTriangulation& T, V& r,
                           double tol = 1e-10, int max_iter = 50,
                           std::vector<AlexandrovSolver::TraceEntry>* out_trace = nullptr) {
  using LinAlg::energy; using LinAlg::norm; using LinAlg::max_abs;

  double r_avg = LinAlg::dot(r, V(r.size(), 1.0)) / r.size();
  double Delta = 0.5 * r_avg, Delta_max = 2.0 * r_avg;
  int rejects = 0;

  for (int iter = 0; iter < max_iter; iter++) {
    auto kappa = GCP::kappa(T, r);
    if (max_abs(kappa) < tol) return {true, max_abs(kappa)};
    if (rejects > 20) break;

    double E    = energy(kappa);
    auto   J    = GCP::jacobian(T, r);
    auto   delta = TrustRegion::solve(J, kappa, Delta);
    double pred = TrustRegion::predicted_reduction(J, kappa, delta);
    V      r_trial = r + delta;
    double E_trial = energy(GCP::kappa(T, r_trial));

    auto [ok, D2] = TrustRegion::update(E - E_trial, pred, norm(delta), Delta, Delta_max);
    Delta = D2;
    if (out_trace)
      out_trace->push_back(make_trace('N', iter, 0.0, Delta, ok ? 1 : 0, kappa, J));
    if (ok) { r = r_trial; Topology::flip_to_delaunay(T, r); rejects = 0; }
    else rejects++;
  }
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

  int h0 = T.f_he[f0], h1 = T.he_next[h0], h2 = T.he_next[h1];
  int i = T.he_origin[h0], j = T.he_origin[h1], k = T.he_origin[h2];

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
  double vol = 0;
  for (int f = 0; f < T.nf; f++) {
    if (T.f_he[f] < 0) continue;
    int ha = T.f_he[f], hb = T.he_next[ha], hc = T.he_next[hb];
    vol += pos[T.he_origin[ha]].dot(pos[T.he_origin[hb]].cross(pos[T.he_origin[hc]]));
  }
  if (vol < 0) for (auto& p : pos) p = p * (-1.0);

  return pos;
}

} // namespace Reconstruct

} // anonymous namespace

// ============================================================================
// AlexandrovSolver: top-level 5-step algorithm
// ============================================================================

vector<coord3d> AlexandrovSolver::solve() {
  int n = D.nv;
  stats_steps = stats_flips = stats_newton_total = 0;
  trace.clear();

  // 1. Initialize: uniform radii
  double R = 0;
  for (int h = 0; h < D.nh; h += 2)
    if (D.alive(h)) R = max(R, D.he_length[h]);
  r.assign(n, 2 * R);

  // 2. PALC: trace κ(r) = t·κ₁ from t=1 toward the endgame zone
  auto kappa1 = GCP::kappa(D, r);
  double t = 1.0, ds = 0.1;
  vector<pair<double, vector<double>>> history;

  for (int step = 0; step < 500 && t > 0.1; step++) {
    auto J0 = GCP::jacobian(D, r);
    auto tau = PALC::compute_tangent(J0, kappa1);
    auto [tp, rp] = PALC::predict(t, r, tau, ds);
    double tc = tp; auto rc = rp;
    int nit = PALC::correct(D, tc, rc, kappa1, tau, t, r, ds, 8, 1e-12);

    if (nit >= 0 && nit < 8) {
      t = tc; r = rc;
      stats_flips += Topology::flip_to_delaunay(D, r);
      if (!history.empty() && fabs(t - history.back().first) < 1e-14)
        history.back() = {t, r};
      else
        history.push_back({t, r});
      ds = PALC::adapt_ds(ds, nit, 8, 0.5);
    } else {
      ds *= 0.5;
      if (ds < 1e-15) break;
    }
    if (trace_jacobian)
      trace.push_back(make_trace('P', step, t, ds, nit, GCP::kappa(D, r), J0));
    stats_steps++;
    stats_newton_total += max(nit, 0);
  }

  // 3. Endgame: extrapolate to t = 0, with Tier-1 feasibility guard.
  // Lagrange extrapolation can overstep the pyramid-closure boundary when
  // the PALC path bends near a fold-fold singularity (observed on
  // C120 #1061350 and C140 #5207982).  When r_ext leaves F, we land on
  // the line from r_last to r_ext at 0.95·s_max for a conditioning margin.
  stats_extrap_kappa = 0;
  r_before_extrap = r;
  if (!history.empty()) {
    auto r_ext = PALC::extrapolate(history);
    if (!r_ext.empty()) {
      double s_max = GCP::feasibility_max_step(D, r, r_ext);
      double s = (s_max < 1.0) ? 0.95 * s_max : 1.0;
      for (size_t i = 0; i < r.size(); i++) r[i] += s * (r_ext[i] - r[i]);
    }
    stats_extrap_kappa = LinAlg::max_abs(GCP::kappa(D, r));
  }

  // 4. Polish: trust-region Newton on κ(r) = 0
  auto [ok, mk] = Newton::polish(D, r, 1e-10, 50, trace_jacobian ? &trace : nullptr);
  stats_final_kappa = mk;

  if (verbose)
    printf("  %d PALC steps, %d flips, max|κ|=%.2e (%s)\n",
           stats_steps, stats_flips, mk, ok ? "converged" : "FAILED");

  if (mk > 0.01) return {};

  // 5. Reconstruct 3D positions
  return Reconstruct::from_radii(D, r);
}

vector<coord3d> AlexandrovSolver::reconstruct(const DelaunayTriangulation& T,
                                              const vector<double>& r) {
  return Reconstruct::from_radii(T, r);
}

vector<double> AlexandrovSolver::kappa(const DelaunayTriangulation& T,
                                        const vector<double>& r) {
  return GCP::kappa(T, r);
}

double AlexandrovSolver::H(const DelaunayTriangulation& T,
                            const vector<double>& r) {
  return GCP::H(T, r);
}

vector<double> AlexandrovSolver::jacobian_eigvals(const DelaunayTriangulation& T,
                                                    const vector<double>& r) {
  return sym_eigvals(GCP::jacobian(T, r));
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
