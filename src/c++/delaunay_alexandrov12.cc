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
#include <map>
#include <random>
#include <set>
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

// Forward declaration for Reconstruct::from_radii — used by Gauge::snap
// before Reconstruct is defined.  Implementation in Layer 6b below.
namespace Reconstruct {
  std::vector<coord3d> from_radii(const DelaunayTriangulation& T,
                                   const std::vector<double>& r);
}
namespace GCP {
  std::vector<double> feasible_step(const DelaunayTriangulation& T,
                                     const std::vector<double>& r_from,
                                     const std::vector<double>& delta,
                                     bool* clipped);
}

// ============================================================================
// Translation-gauge fixing.
//
// At any solution r* of κ(r) = 0, dκ/dr has a 3-D kernel from translating
// the apex of the B-I star-pyramid construction (apex shift α gives
// δr_v = −α · p̂_v).  Off-manifold this becomes a "soft kernel" that
// makes Newton ill-conditioned.  Centroid-snap factors out the gauge by
// canonicalising r to the apex-at-centroid representative of its
// translation orbit at each step:
//   1. reconstruct positions p_v from current r via Gram-BFS,
//   2. compute centroid c = (1/n) Σ p_v,
//   3. propose r' with r'_v = ‖p_v − c‖,
//   4. apply r ← r + GCP::feasible_step(T, r, r' − r) to keep r ∈ F(T).
// Step (4) is essential: at off-manifold iterates the reconstructed
// positions are a B-I generalized polytope (not a real polytope), so
// the proposed r' may itself violate pyramid closure.  Clipping via
// feasible_step keeps the iterate inside F(T).
// Idempotent at κ = 0; no-op if Gram-BFS fails.
// ============================================================================
namespace Gauge {

// Returns true if r was updated; false if Gram-BFS failed (r unchanged).
bool snap(const DelaunayTriangulation& T, std::vector<double>& r) {
  auto pos = Reconstruct::from_radii(T, r);
  if ((int)pos.size() != T.nv) return false;
  for (const auto& p : pos) {
    if (std::isnan(p[0]) || std::isnan(p[1]) || std::isnan(p[2])) return false;
  }
  coord3d c(0, 0, 0);
  for (const auto& p : pos) c = c + p;
  c = c * (1.0 / pos.size());
  std::vector<double> delta(T.nv);
  for (int v = 0; v < T.nv; v++) delta[v] = (pos[v] - c).norm() - r[v];
  // Clip the gauge-snap update to F(T).  At large-κ iterates the snap
  // can push r outside F(T) and the next pyramid closure check will
  // produce NaN.
  auto clipped_delta = GCP::feasible_step(T, r, delta, nullptr);
  for (int v = 0; v < T.nv; v++) r[v] += clipped_delta[v];
  return true;
}

// Build the 3-D gauge basis Q (12×3 orthonormal columns) at the current
// iterate r.  The k-th gauge mode in r-space is δr^{(k)}_v = −ê_k · p̂_v
// with p̂_v = p_v / r_v (apex at origin), the change in r induced by
// translating the apex by ê_k.  Returns Q as 3 row vectors of length n,
// or empty if Gram-BFS fails / basis is rank-deficient.
std::vector<std::vector<double>> basis(const DelaunayTriangulation& T,
                                         const std::vector<double>& r) {
  auto pos = Reconstruct::from_radii(T, r);
  if ((int)pos.size() != T.nv) return {};
  for (const auto& p : pos)
    if (std::isnan(p[0]) || std::isnan(p[1]) || std::isnan(p[2])) return {};

  int n = T.nv;
  std::vector<std::vector<double>> Q(3, std::vector<double>(n));
  for (int v = 0; v < n; v++) {
    double inv_r = (r[v] > 1e-15) ? 1.0 / r[v] : 0.0;
    for (int k = 0; k < 3; k++)
      Q[k][v] = -pos[v][k] * inv_r;
  }
  // Modified Gram-Schmidt
  for (int k = 0; k < 3; k++) {
    for (int j = 0; j < k; j++) {
      double dot = 0;
      for (int v = 0; v < n; v++) dot += Q[k][v] * Q[j][v];
      for (int v = 0; v < n; v++) Q[k][v] -= dot * Q[j][v];
    }
    double norm2 = 0;
    for (int v = 0; v < n; v++) norm2 += Q[k][v] * Q[k][v];
    double norm = std::sqrt(norm2);
    if (norm < 1e-12) return {};   // basis rank-deficient
    for (int v = 0; v < n; v++) Q[k][v] /= norm;
  }
  return Q;
}

// Project step Δr onto gauge-orthogonal subspace: Δr ← Δr − Q Qᵀ Δr.
// No-op if Q is empty.
void project_step(std::vector<double>& delta,
                   const std::vector<std::vector<double>>& Q) {
  for (const auto& q : Q) {
    double dot = 0;
    for (size_t i = 0; i < delta.size(); i++) dot += q[i] * delta[i];
    for (size_t i = 0; i < delta.size(); i++) delta[i] -= dot * q[i];
  }
}

} // namespace Gauge

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

// Pseudo-arc-length bordered system, treated as a first-class object:
//
//   B = [   J         −κ₁     ]    (n × (n+1))
//       [   r_dotᵀ     t_dot  ]    (1 × (n+1))
//
// .solve(neg_F, g) returns (δr, δt, sign(det B)) by the bordering
// algorithm: two n-solves with J share an LU, and
//   det(B) = det(J) · (t_dot + r_dotᵀ · J⁻¹ · κ₁)
// by Schur complement, so the sign is free.
struct Bordered {
  const matrix<double>& J;
  const V&              kappa1;
  const V&              r_dot;
  double                t_dot;

  struct Solution { V dr; double dt; int det_sign; };

  // Pre: J is n×n; kappa1, r_dot are n-vectors; neg_F is n-vector; g is scalar.
  // Post on success: B · [dr; dt] = [neg_F; g] within LAPACK precision,
  //                  and det_sign ∈ {−1, +1} is sign(det B).
  // Post on failure: LuFail::Singular if J near-singular or denominator near 0.
  std::expected<Solution, LuFail> solve(const V& neg_F, double g) const {
    auto z1 = solve_with_sign(J, neg_F);
    if (!z1) return std::unexpected(z1.error());
    auto z2 = solve_with_sign(J, kappa1);
    if (!z2) return std::unexpected(z2.error());
    double den = t_dot + dot(r_dot, z2->x);
    if (fabs(den) < 1e-30) return std::unexpected(LuFail::Singular);
    double dt = (g - dot(r_dot, z1->x)) / den;
    V      dr = z1->x + z2->x * dt;
    int sign_den = (den > 0) - (den < 0);   // {−1, +1}
    return Solution{std::move(dr), dt, z1->det_sign * sign_den};
  }
};

// Solve (A + λI)·x = b.
V solve_shifted(const matrix<double>& A, const V& b, double lambda) {
  matrix<double> Al = A;
  for (int i = 0; i < A.m; i++) Al(i,i) += lambda;
  return solve(Al, b);
}

// --- Block constructors / destructors for the bordered system ---
// These are the named primitives that turn the corrector's imperative
// matrix-and-vector assembly into a sequence of pure functions.

// Bordered matrix:    [ J     c ]    (n+1) × (n+1)
//                     [ rᵀ    t ]
matrix<double> bordered(const matrix<double>& J,
                         const V& c, const V& r, double t) {
  int n = J.m;
  matrix<double> B(n + 1, n + 1, 0.0);
  for (int i = 0; i < n; i++) {
    for (int j = 0; j < n; j++) B(i, j) = J(i, j);
    B(i, n) = c[i];
    B(n, i) = r[i];
  }
  B(n, n) = t;
  return B;
}

// Concatenate vector and scalar:  [v₀, …, v_{n-1}, s]
V concat(const V& v, double s) {
  V out = v;
  out.push_back(s);
  return out;
}

// Inverse of concat: split [x₀, …, x_{n-1}, x_n] ↦ ([x₀, …, x_{n-1}], x_n)
std::pair<V, double> split(const V& x) {
  V head(x.begin(), x.end() - 1);
  return {std::move(head), x.back()};
}

// --- Symmetric truncated-pseudoinverse solve ---
//
// J is the Hessian of the BI total scalar curvature H, hence symmetric.
// Its spectral decomposition  J = Q · diag(λ) · Qᵀ  yields the pseudo-
// inverse  J⁺ = Q · diag(λ⁺) · Qᵀ  with (λ⁺)ᵢ = 1/λᵢ when |λᵢ| > rcond·
// max|λ|, else 0.
//
// For the bordered system   B = [J -κ₁; τᵀ τ_t]   used in the PALC
// corrector, we apply the symmetric pseudoinverse twice via the standard
// bordering algorithm rather than directly factoring the (non-symmetric)
// 13×13 B.  This has three benefits:
//
//   1. Half the work: one 12×12 dsyev cached, applied twice.
//   2. The bordered solution exactly satisfies the arc-length constraint
//      and lies in image(J) = ker(J)^⊥ — the gauge-orthogonal subspace.
//      Residual is confined to ker(J), the translation gauge of κ, where
//      it is harmless to homotopy convergence by translation invariance.
//   3. The eigenvalues yield the Lorentzian signature directly, giving
//      a per-step diagnostic that distinguishes the BI failure classes
//      (drum-cap, soft-kernel, fold) empirically.
namespace SymEigen {

// Symmetric eigendecomposition: J = Q·diag(λ)·Qᵀ.
// Q has eigenvectors as columns (column-major); λ ascending.
// Empty on LAPACK failure.
struct Decomp {
  matrix<double> Q;
  V              lambda;
};

Decomp decompose(const matrix<double>& A) {
  int n = A.m;
  std::vector<double> Af(n * n);
  for (int i = 0; i < n; i++)                       // column-major + symmetrize
    for (int j = 0; j < n; j++)
      Af[j * n + i] = 0.5 * (A(i, j) + A(j, i));
  V lambda(n);
  // dsyev with JOBZ='V' needs LWORK ≥ max(1, 1 + 6n + 2n²) for the
  // optimum, but any larger value works.  For our typical n=12 the
  // requirement is ≤ 361; we allocate generously and skip the
  // workspace-query call entirely.  This saves one dsyev invocation
  // per decomposition — a measurable speedup at small n.
  int lwork = std::max(64, 1 + 6 * n + 2 * n * n);
  std::vector<double> work(lwork);
  int info = 0;
  dsyev_("V", "U", &n, Af.data(), &n, lambda.data(),
         work.data(), &lwork, &info);
  if (info != 0) return Decomp{matrix<double>(0, 0, 0.0), {}};
  matrix<double> Q(n, n, 0.0);
  for (int i = 0; i < n; i++)
    for (int j = 0; j < n; j++)
      Q(j, i) = Af[i * n + j];                      // eigvec i is column i
  return {std::move(Q), std::move(lambda)};
}

// Lorentzian signature of a symmetric operator at a given numerical
// "zero" tolerance:  (#λ > tol, #λ < -tol, #|λ| ≤ tol).
struct Signature { int pos, neg, zero; };

Signature signature(const V& lambda, double tol) {
  Signature s{0, 0, 0};
  for (double l : lambda) {
    if      (l >  tol) s.pos++;
    else if (l < -tol) s.neg++;
    else               s.zero++;
  }
  return s;
}

// Apply truncated pseudoinverse to a vector:
//   x = Σ_{|λᵢ|>cutoff} (qᵢᵀ b / λᵢ) qᵢ
// Returns rank = #non-truncated modes; lambda_min_kept = smallest |λ|
// among those modes (0 if all truncated).
struct ApplyResult { V x; int rank; double lambda_min_kept; };

ApplyResult apply(const Decomp& d, const V& b, double cutoff) {
  int n = (int)d.lambda.size();
  V   x(n, 0.0);
  int rank = 0;
  double lam_min = 0;
  for (int i = 0; i < n; i++) {
    double l = d.lambda[i];
    if (std::fabs(l) <= cutoff) continue;
    rank++;
    if (lam_min == 0 || std::fabs(l) < lam_min) lam_min = std::fabs(l);
    double bq = 0;
    for (int j = 0; j < n; j++) bq += d.Q(j, i) * b[j];
    double coeff = bq / l;
    for (int j = 0; j < n; j++) x[j] += coeff * d.Q(j, i);
  }
  return {std::move(x), rank, lam_min};
}

struct Solution {
  V         x;
  int       rank          = 0;
  double    lambda_abs_min = 0;
  double    lambda_abs_max = 0;
  Signature sig            = {0, 0, 0};
};

// Standalone symmetric pseudoinverse solve: x = J⁺·b for symmetric J.
Solution solve(const matrix<double>& A, const V& b, double rcond) {
  auto d = decompose(A);
  if (d.lambda.empty()) return {};
  double lam_max = 0;
  for (double l : d.lambda) lam_max = std::max(lam_max, std::fabs(l));
  double cutoff = rcond * lam_max;
  auto   r      = apply(d, b, cutoff);
  return {std::move(r.x), r.rank, r.lambda_min_kept, lam_max,
          signature(d.lambda, cutoff)};
}

// Bordered solve via bordering algorithm + ONE symmetric eigendecomp of J.
//
//   B·[Δr; Δt] = [-F; g]    where    B = [J -κ₁; τ_rᵀ τ_t]
//
// Decompose J once, apply the truncated pseudoinverse twice (cached):
//   z₁ = J⁺·(-F)
//   z₂ = J⁺·κ₁
//   Δt = (g − τ_r·z₁) / (τ_t + τ_r·z₂)
//   Δr = z₁ + z₂·Δt
// Bottom row (arc-length closure) is satisfied EXACTLY by construction.
// Top row is satisfied exactly within image(J); residual lives in ker(J),
// which is the translation gauge of κ — a direction along which κ is
// invariant, so the residual is harmless to homotopy convergence.
struct BorderedSolution {
  V         dr;
  double    dt             = 0;
  int       rank           = 0;
  double    lambda_abs_min = 0;
  double    lambda_abs_max = 0;
  Signature sig            = {0, 0, 0};
};

BorderedSolution bordered_solve(const matrix<double>& J, const V& kappa1,
                                 const V& tau_r, double tau_t,
                                 const V& neg_F, double g, double rcond) {
  auto d = decompose(J);
  if (d.lambda.empty()) return {};
  double lam_max = 0;
  for (double l : d.lambda) lam_max = std::max(lam_max, std::fabs(l));
  double cutoff = rcond * lam_max;
  auto   r1     = apply(d, neg_F, cutoff);
  auto   r2     = apply(d, kappa1, cutoff);
  double den    = tau_t + LinAlg::dot(tau_r, r2.x);
  if (std::fabs(den) < 1e-30) return {};
  double dt     = (g - LinAlg::dot(tau_r, r1.x)) / den;
  V      dr     = r1.x;
  for (size_t i = 0; i < dr.size(); i++) dr[i] += r2.x[i] * dt;
  // r1.rank == r2.rank since both share (Q, λ, cutoff); take r1's.
  double lam_min = (r1.lambda_min_kept > 0) ? r1.lambda_min_kept : 0;
  return {std::move(dr), dt, r1.rank, lam_min, lam_max,
          signature(d.lambda, cutoff)};
}

}  // namespace SymEigen

// Predicate: TSVD-solve outcome is acceptable.
// Up to 3 truncations is the gauge regime; more indicates non-gauge
// collapse and the caller should treat as failure.
bool tsvd_ok(int rank, int n) {
  return rank >= n - 3;
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

// Pyramid height squared at half-edge h.  Inlined replica of the
// computation inside GCP::alpha — exposes h_sq as a diagnostic without
// changing alpha's public contract.  Negative result means the pyramid
// fails to close (Cayley-Menger violation).
static double pyramid_h_sq_at(const DelaunayTriangulation& T,
                                const V& r, int h) {
  int u = T.he_origin[h];
  int v = T.he_origin[T.he_next[h]];
  int w = T.he_origin[T.prev(h)];
  double luv = T.he_length[h];
  double lvw = T.he_length[T.he_next[h]];
  double lwu = T.he_length[T.prev(h)];
  double wx = (luv*luv + lwu*lwu - lvw*lvw) / (2*luv);
  double wy_sq = lwu*lwu - wx*wx;
  if (wy_sq < -1e-10 * lwu*lwu) return -1.0;
  double wy = std::sqrt(std::max(0.0, wy_sq));
  double px = (r[u]*r[u] - r[v]*r[v] + luv*luv) / (2*luv);
  double py = (wy > 1e-15)
                ? (r[u]*r[u] - r[w]*r[w] + wx*wx + wy*wy - 2*px*wx) / (2*wy)
                : 0;
  return r[u]*r[u] - px*px - py*py;
}

// Pack one PALC- or Newton-step trajectory-diagnostic record.  Cheap to
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
    if (T.he_face[h] == T.he_face[h ^ 1]) continue;  // bigon — θ undefined
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
// `margin` parameterizes how strictly the iterate is held inside P(M):
//   margin = -1e-10  closure with numerical noise buffer (existing behavior;
//                    flip when θ > π + 1e-10 or q_j < q_i − ℓ²_ij − 1e-10)
//   margin = 0       exact closure (boundary triggers a flip)
//   margin = c > 0   strict-interior with safety distance c from ∂P(M):
//                    flip when θ > π − c or q_j < q_i − ℓ²_ij + c
// PALC uses a t-dependent schedule c·t so the margin shrinks to 0 as t → 0,
// allowing legitimate flat-face diagonals (θ → π) at the polytope limit
// while pushing the path strictly inside at intermediate t.  At t = 0
// (Newton polish) the default margin = -1e-10 reproduces the closure check.
int needs_flip(const DelaunayTriangulation& T, const vector<double>& r,
                double margin = -1e-10) {
  // (ConcQuadr): edge with two distinct adjacent triangles, θ > π − margin.
  // For a bigon edge (he_face[h] == he_face[h^1]), GCP::theta is not
  // meaningful; skip and let the (CloGeod) pass below handle it.
  // NaN θ — returned by alpha() when the abstract pyramid is degenerate
  // (h_sq < 0, i.e. Cayley-Menger fails) — also counts as bad: the
  // configuration is outside P(M) and a flip is required.  IEEE
  // comparisons with NaN are false, so we test isnan() explicitly.
  double theta_threshold = M_PI - margin;
  for (int h = 0; h < T.nh; h += 2) {
    if (!T.alive(h)) continue;
    if (T.he_face[h] == T.he_face[h ^ 1]) continue;       // bigon — handled below
    double theta = GCP::theta(T, r, h);
    if (std::isnan(theta) || theta > theta_threshold) return h;
  }
  // (CloGeod): bigon i–j edge of an iji-face, q_j < q_i − ℓ²_ij + margin.
  for (int h = 0; h < T.nh; h += 2) {
    if (!T.alive(h)) continue;
    if (T.he_face[h] != T.he_face[h ^ 1]) continue;       // not a bigon
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

// Flip until needs_flip(T, r, margin) returns -1 (or 20 iter cap).
// Returns number of flips performed.
int flip_to_delaunay(DelaunayTriangulation& T, const vector<double>& r,
                      double margin = -1e-10) {
  int total = 0;
  for (int iter = 0; iter < 20; iter++) {
    int h = needs_flip(T, r, margin);
    if (h < 0) break;
    if (T.flip_edge(h)) total++;
    else break;
  }
  return total;
}

// Synchronized batch multi-flip (Direction 2).  Identify all alive
// non-bigon edges with θ > π − threshold, then flip them all in one
// pass.  Single-edge flip cycles (which defeat sequential flip_to_-
// delaunay on the symmetric drum-cap metrics) cannot occur here
// because no θ re-evaluation happens between flips in the batch.
// Returns number of flips actually performed (some may have been
// skipped if they share vertices with an earlier flip in the batch
// and are no longer alive by the time we reach them).
int batch_multi_flip(DelaunayTriangulation& T, const vector<double>& r,
                      double threshold) {
  if (threshold <= 0) return 0;
  vector<int> targets;
  double cutoff = M_PI - threshold;
  for (int h = 0; h < T.nh; h += 2) {
    if (!T.alive(h)) continue;
    if (T.he_face[h] == T.he_face[h ^ 1]) continue;       // bigon
    double theta = GCP::theta(T, r, h);
    if (std::isnan(theta) || theta > cutoff) targets.push_back(h);
  }
  int total = 0;
  for (int h : targets) {
    if (!T.alive(h)) continue;                             // killed by earlier flip
    if (T.flip_edge(h)) total++;
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

// Deflation context for Direction 4: augments the homotopy residual by
// a repulsive term centered at `target` with magnitude `strength`.
//   F(r, t) → F(r, t) + (1 − t) · α · (r − r*) / ‖r − r*‖²
// active() ⇒ apply.  inactive Deflation is the identity transformation
// — F, J, and κ₁ are returned unchanged.
struct Deflation {
  vector<double> target;
  double         strength = 0.0;

  bool active(int n) const {
    return strength > 0 && (int)target.size() == n;
  }

  // D(r) = α · (r − r*) / ‖r − r*‖²
  vector<double> D(const vector<double>& r) const {
    int n = (int)r.size();
    vector<double> delta(n), out(n);
    double dn2 = 0;
    for (int i = 0; i < n; i++) { delta[i] = r[i] - target[i]; dn2 += delta[i]*delta[i]; }
    if (dn2 < 1e-30) return out;   // singular at target — caller should avoid
    double s = strength / dn2;
    for (int i = 0; i < n; i++) out[i] = s * delta[i];
    return out;
  }

  // ∂D_i/∂r_j = α/‖δ‖² · (δ_ij − 2·δ_i·δ_j/‖δ‖²)
  matrix<double> dDdr(const vector<double>& r) const {
    int n = (int)r.size();
    matrix<double> dD(n, n, 0.0);
    vector<double> delta(n);
    double dn2 = 0;
    for (int i = 0; i < n; i++) { delta[i] = r[i] - target[i]; dn2 += delta[i]*delta[i]; }
    if (dn2 < 1e-30) return dD;
    double s = strength / dn2;
    double inv_dn2 = 1.0 / dn2;
    for (int i = 0; i < n; i++)
      for (int j = 0; j < n; j++)
        dD(i, j) = s * ((i == j ? 1.0 : 0.0) - 2.0 * delta[i] * delta[j] * inv_dn2);
    return dD;
  }
};

struct Tangent { double t_dot; vector<double> r_dot; };

// Tangent to the homotopy curve at (t, r).
// From J·ṙ = κ₁·ṫ and ||(ṫ, ṙ)|| = 1, ṫ < 0.
// With deflation:  J_def·ṙ = κ₁_eff·ṫ where J_def = J + (1−t)·∂D/∂r
// and κ₁_eff = κ₁ + D(r).
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

// Arc-length residual:  g = ds − τ_r·(r − r₀) − τ_t·(t − t₀)
// (Closure of the predictor-corrector arc to length ds.)
double arc_residual(double t, const V& r, double t0, const V& r0,
                     const Tangent& tau, double ds) {
  return ds - LinAlg::dot(tau.r_dot, r - r0) - tau.t_dot * (t - t0);
}

// Cap step magnitude to a trust radius.  No-op if already inside.
V scale_to_radius(const V& dr, double radius) {
  double n = LinAlg::norm(dr);
  return (n > radius) ? dr * (radius / n) : dr;
}

// Apply the deflation modification (in place) to F, J, and κ₁_eff.
// Models F_def = F + (1−t)·D(r);  J_def = J + (1−t)·∂D/∂r;  κ₁_eff = κ₁ + D(r).
// No-op if defl is null or inactive.
void apply_deflation(const Deflation* defl, double t, const V& r,
                      V& F, matrix<double>& J, V& kappa1_eff) {
  if (!defl || !defl->active((int)r.size())) return;
  auto Dvec = defl->D(r);
  auto dD   = defl->dDdr(r);
  double w  = 1.0 - t;
  for (size_t i = 0; i < F.size(); i++)         F[i]          += w * Dvec[i];
  for (int    i = 0; i < J.m;       i++)
    for (int  j = 0; j < J.n;       j++)        J(i, j)       += w * dD(i, j);
  for (size_t i = 0; i < kappa1_eff.size(); i++) kappa1_eff[i] += Dvec[i];
}

// Corrector: Newton on the bordered system at each iterate.
//
// Two modes:
//   use_tsvd = false  →  classic Bordered LU solve.  Standard Newton.
//   use_tsvd = true   →  TSVD pseudoinverse on the bordered matrix.
//                         Truncates singular values below rcond·σ_max,
//                         producing the unique gauge-orthogonal min-norm
//                         step.  See tsvd-design.md.
//
// Returns number of iterations (−1 on failure).
int correct(const DelaunayTriangulation& T, double& t, V& r,
            const V& kappa1, const Tangent& tau0,
            double t0, const V& r0, double ds,
            int max_iter, double tol,
            const Deflation* defl = nullptr,
            bool /*gauge_project*/ = false,
            bool use_tsvd = false,
            double rcond = 5e-3) {
  for (int nit = 0; nit < max_iter; nit++) {
    auto F          = GCP::kappa(T, r) - kappa1 * t;
    auto J          = GCP::jacobian(T, r);
    auto kappa1_eff = kappa1;
    apply_deflation(defl, t, r, F, J, kappa1_eff);

    if (LinAlg::max_abs(F) < tol) return nit;

    double g = arc_residual(t, r, t0, r0, tau0, ds);

    if (use_tsvd) {
      // Symmetric pseudoinverse on J + bordering algorithm.
      auto sol = LinAlg::SymEigen::bordered_solve(
                    J, kappa1_eff, tau0.r_dot, tau0.t_dot, -F, g, rcond);
      if (sol.dr.empty() || !LinAlg::tsvd_ok(sol.rank, J.m)) return -1;
      // Step-norm termination: TSVD truncates 2–3 directions where κ
      // is invariant, so the corrector cannot drive |F| to zero there
      // — there is an irreducible kernel residual at O(rcond·σ_max).
      // When the proposed step is tiny, no further progress is possible;
      // accept.  Scaled by ‖r‖ to auto-adapt to problem size.
      double step  = LinAlg::norm(sol.dr) + std::fabs(sol.dt);
      double scale = LinAlg::norm(r) + 1.0;
      if (step < tol * scale) return nit;
      t += sol.dt;
      r  = r + sol.dr;
    } else {
      LinAlg::Bordered B{J, kappa1_eff, tau0.r_dot, tau0.t_dot};
      auto sol = B.solve(-F, g);
      if (!sol) return -1;
      t += sol->dt;
      r  = r + sol->dr;
    }
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

// ---- Structured PALC primitives (step 3 of Tier-2 plan) ----

constexpr int CORRECTOR_MAX_ITER = 8;
constexpr double CORRECTOR_TOL = 1e-12;
constexpr double DS_MIN = 1e-15;
constexpr double DS_MAX = 0.5;

// One predictor-corrector step.  Returns the corrector iteration count
// alongside the (t, r) iterate — always present; "accepted" iff
// nit ∈ [0, CORRECTOR_MAX_ITER).  t, r are valid only if accepted.
struct StepResult {
  int nit;
  double t; V r;
  bool accepted() const { return nit >= 0 && nit < CORRECTOR_MAX_ITER; }
};

// Pre:  (t₀, r₀) ∈ F(T) with 0 < κᵢ(t₀, r₀) for the BI homotopy.
// Post: if accepted(), (t, r) is on the homotopy path at arc-length
//       distance ds from (t₀, r₀), to within CORRECTOR_TOL.
StepResult palc_step(const DelaunayTriangulation& T,
                      double t0, const V& r0, const V& kappa1,
                      const matrix<double>& J, double ds,
                      const Deflation* defl = nullptr,
                      bool gauge_project = false,
                      bool use_tsvd = false,
                      double rcond = 5e-3) {
  // Tangent on the deflated curve uses J_def and κ₁_eff.
  matrix<double> J_eff = J;
  V kappa1_eff = kappa1;
  if (defl && defl->active((int)r0.size())) {
    V F_dummy(r0.size(), 0.0);
    apply_deflation(defl, t0, r0, F_dummy, J_eff, kappa1_eff);
  }
  auto tau = compute_tangent(J_eff, kappa1_eff);
  auto [tp, rp] = predict(t0, r0, tau, ds);
  double tc = tp; V rc = rp;
  int nit = correct(T, tc, rc, kappa1, tau, t0, r0, ds,
                     CORRECTOR_MAX_ITER, CORRECTOR_TOL, defl,
                     gauge_project, use_tsvd, rcond);
  return {nit, tc, std::move(rc)};
}

struct PALCStats {
  int steps = 0, flips = 0, newton_total = 0;
};

struct TrackResult {
  double t_final;
  V r_final;
  vector<pair<double, V>> history;
  PALCStats stats;
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

// Track the homotopy κ(r) = t·κ₁ from t=1 toward t_target.
//
// Pre:  r ∈ F(T), 0 < κᵢ(T, r) < δᵢ, t_target > 0.
// Post: t_final ≤ t_target (if PALC converges) or t_final > t_target
//       (if PALC stalls — caller should fall back to Tier-2).
TrackResult palc_track(DelaunayTriangulation& T, V r, const V& kappa1,
                        double t_target, double ds_init,
                        vector<AlexandrovSolver::TraceEntry>* trace,
                        vector<AlexandrovSolver::DiagEntry>* diag,
                        double interior_margin_coeff = 0.0,
                        double stochastic_eps = 0.0,
                        uint32_t stochastic_seed = 1,
                        double batch_multiflip_threshold = 0.0,
                        const Deflation* deflation = nullptr,
                        bool gauge_snap = false,
                        bool gauge_project = false,
                        bool use_tsvd = false,
                        double rcond = 5e-3) {
  double t = 1.0, ds = ds_init;
  vector<pair<double, V>> history;
  PALCStats stats;
  std::mt19937 rng(stochastic_seed);
  std::uniform_real_distribution<double> jitter(-1.0, 1.0);

  for (int step = 0; step < 500 && t > t_target; step++) {
    auto J = GCP::jacobian(T, r);
    auto result = palc_step(T, t, r, kappa1, J, ds, deflation,
                             gauge_project, use_tsvd, rcond);

    if (result.accepted()) {
      t = result.t; r = std::move(result.r);
      // Strict-interior margin schedule: c·t at intermediate t shrinks to
      // 0 at t = 0 (where flat-face diagonals legitimately have θ → π).
      // c = 0 (default) reproduces pre-#28 behaviour (closure check via
      // the negative -1e-10 numerical buffer in needs_flip).
      double margin = (interior_margin_coeff > 0)
                        ? interior_margin_coeff * t
                        : -1e-10;
      stats.flips += Topology::flip_to_delaunay(T, r, margin);

      // Batch multi-flip (Direction 2): after the regular post-accept
      // flip cleanup, do a single synchronized batch flip of any
      // remaining near-π edges to escape sequential-flip cycles.
      // Then a final flip_to_delaunay cleans up the post-batch state.
      if (batch_multiflip_threshold > 0) {
        int batch_flips = Topology::batch_multi_flip(T, r,
                                                       batch_multiflip_threshold);
        if (batch_flips > 0) {
          stats.flips += batch_flips;
          stats.flips += Topology::flip_to_delaunay(T, r, margin);
        }
      }

      // Stochastic perturbation (Direction 1): r ← r·(1 + ε·u_v) per
      // vertex, with u_v ~ uniform[−1,+1].  Applied after the post-
      // accept flip phase so the corrector starts the next step from
      // a perturbed but still in-F(T)-near r (small eps preserves
      // F(T)).  Breaks symmetry continuously rather than only at
      // t = 1.  ε = 0 disables.
      if (stochastic_eps > 0) {
        for (size_t i = 0; i < r.size(); i++) {
          r[i] *= 1.0 + stochastic_eps * jitter(rng);
        }
      }

      // Translation-gauge fixing (Direction 5): canonicalise r to the
      // apex-at-centroid representative of its translation orbit.
      // Gated on small-κ iterates: the snap reconstructs positions via
      // Gram-BFS, but at high-κ those are a B-I generalized polytope
      // (not a real polytope), so its centroid is geometrically
      // meaningless and snapping adds noise.  Threshold 0.5 catches
      // the moderate-κ PALC-stall regime where the soft kernel is
      // dominant, while excluding initial-uniform-r iterates where
      // |κ|_max ≈ 1 (B-I admissibility κ_i < δ_i ≈ 1.05).
      if (gauge_snap && LinAlg::max_abs(GCP::kappa(T, r)) < 0.5)
        Gauge::snap(T, r);
      if (!history.empty() && fabs(t - history.back().first) < 1e-14)
        history.back() = {t, r};
      else
        history.push_back({t, r});
      ds = adapt_ds(ds, result.nit, CORRECTOR_MAX_ITER, DS_MAX);
    } else {
      ds *= 0.5;
      if (ds < DS_MIN) break;
    }
    if (trace)
      trace->push_back(make_trace('P', step, t, ds, result.nit,
                                   GCP::kappa(T, r), J));
    if (diag) {
      auto Jd = GCP::jacobian(T, r);  // J on the (possibly post-flip) state
      diag->push_back(make_diag('P', step, t, ds, result.nit,
                                  T, r, GCP::kappa(T, r), Jd, stats.flips));
    }
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
                           std::vector<AlexandrovSolver::TraceEntry>* out_trace = nullptr,
                           std::vector<AlexandrovSolver::DiagEntry>* out_diag = nullptr,
                           int* flips_cum = nullptr,
                           bool gauge_snap = false,
                           bool gauge_project = false,
                           bool use_tsvd = false,
                           double rcond = 5e-3) {
  using LinAlg::energy; using LinAlg::norm; using LinAlg::max_abs;

  double r_avg = LinAlg::dot(r, V(r.size(), 1.0)) / r.size();
  double Delta = 0.5 * r_avg, Delta_max = 2.0 * r_avg;
  int rejects = 0;
  int flips_local = flips_cum ? *flips_cum : 0;

  for (int iter = 0; iter < max_iter; iter++) {
    auto kappa = GCP::kappa(T, r);
    if (max_abs(kappa) < tol) {
      if (flips_cum) *flips_cum = flips_local;
      return {true, max_abs(kappa)};
    }
    if (rejects > 20) break;

    double E         = energy(kappa);
    auto   J         = GCP::jacobian(T, r);
    // TSVD or trust-region Newton — both produce a candidate step.
    // TSVD: min-norm pseudoinverse step capped to trust radius.
    //       Naturally gauge-orthogonal at the current iterate.
    // LM:   classic Levenberg-Marquardt-bisected trust-region step.
    V delta_raw;
    if (use_tsvd) {
      // Symmetric pseudoinverse: J = Q diag(λ) Qᵀ, J⁺ truncates the
      // soft kernel, Δr = -J⁺·κ is the gauge-orthogonal min-norm step.
      auto sol = LinAlg::SymEigen::solve(J, -kappa, rcond);
      if (sol.x.empty() || !LinAlg::tsvd_ok(sol.rank, J.m)) {
        rejects++;
        if (out_trace)
          out_trace->push_back(make_trace('N', iter, 0.0, Delta, 0, kappa, J));
        continue;
      }
      delta_raw = PALC::scale_to_radius(sol.x, Delta);
    } else {
      delta_raw = TrustRegion::solve(J, kappa, Delta);
    }
    // Gauge step-projection: remove the 3-D translation-gauge component
    // from the trust-region step before clipping to F(T).  Gated on
    // |κ|_max < 0.5 (high-κ Gram-BFS positions are meaningless).
    if (gauge_project && max_abs(kappa) < 0.5) {
      auto Q = Gauge::basis(T, r);
      if (!Q.empty()) Gauge::project_step(delta_raw, Q);
    }
    // Clip step to F(T): keep r_trial inside the feasibility region so
    // pyramids stay non-degenerate and κ(r_trial) is finite (no NaN θ on
    // multi-edges, no false |κ|=0 convergence outside F).  Same helper
    // used by the endgame extrapolation — single source of truth.
    bool   clipped;
    auto   delta     = GCP::feasible_step(T, r, delta_raw, &clipped);
    double pred      = TrustRegion::predicted_reduction(J, kappa, delta);
    V      r_trial   = r + delta;
    double E_trial   = energy(GCP::kappa(T, r_trial));

    auto [ok, D2] = TrustRegion::update(E - E_trial, pred, norm(delta), Delta, Delta_max);
    // If we clipped, cap Δ at the step we actually took so the next
    // subproblem doesn't keep proposing the same infeasible direction.
    Delta = clipped ? min(D2, norm(delta)) : D2;
    if (out_trace)
      out_trace->push_back(make_trace('N', iter, 0.0, Delta, ok ? 1 : 0, kappa, J));
    if (ok) {
      r = r_trial;
      flips_local += Topology::flip_to_delaunay(T, r);
      if (gauge_snap) Gauge::snap(T, r);
      rejects = 0;
    } else rejects++;

    if (out_diag) {
      auto Jd = GCP::jacobian(T, r);
      out_diag->push_back(make_diag('N', iter, 0.0, Delta, ok ? 1 : 0,
                                      T, r, GCP::kappa(T, r), Jd, flips_local));
    }
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

vector<coord3d> AlexandrovSolver::solve() {
  stats_steps = stats_flips = stats_newton_total = 0;
  trace.clear();
  diag_trace.clear();
  stats_status = ValidationStatus::FAIL_KAPPA_NOT_CONVERGED;

  // 1. Initialize: uniform radii (or caller-provided override).
  if ((int)r_init_override.size() == D.nv) {
    r = r_init_override;
  } else {
    r = PALC::initial_radii(D);
  }
  auto kappa1 = GCP::kappa(D, r);

  // 2. PALC: trace κ(r) = t·κ₁ from t=1 toward t_target
  PALC::Deflation defl;
  defl.target = r_deflate_target;
  defl.strength = deflate_strength;
  auto track = PALC::palc_track(D, r, kappa1, /*t_target=*/0.1, /*ds_init=*/0.1,
                                 trace_jacobian ? &trace : nullptr,
                                 record_diag ? &diag_trace : nullptr,
                                 palc_interior_margin,
                                 stochastic_perturbation_eps,
                                 stochastic_seed,
                                 palc_batch_multiflip_threshold,
                                 defl.active(D.nv) ? &defl : nullptr,
                                 palc_gauge_snap,
                                 palc_gauge_project,
                                 palc_tsvd,
                                 palc_tsvd_rcond);
  r = std::move(track.r_final);
  stats_steps        = track.stats.steps;
  stats_flips        = track.stats.flips;
  stats_newton_total = track.stats.newton_total;

  // 3. Endgame: guarded extrapolation (Tier 1)
  stats_extrap_kappa = 0;
  r_before_extrap = r;
  if (!track.history.empty()) {
    auto r_ext = PALC::extrapolate(track.history);
    if (!r_ext.empty()) r = r + GCP::feasible_step(D, r, r_ext - r);
    stats_extrap_kappa = LinAlg::max_abs(GCP::kappa(D, r));
  }

  // 4. Polish: trust-region Newton on κ(r) = 0
  int polish_flips = stats_flips;
  auto [ok, mk] = Newton::polish(D, r, 1e-10, 50,
                                   trace_jacobian ? &trace : nullptr,
                                   record_diag ? &diag_trace : nullptr,
                                   record_diag ? &polish_flips : nullptr,
                                   palc_gauge_snap,
                                   palc_gauge_project,
                                   palc_tsvd,
                                   palc_tsvd_rcond);
  if (record_diag) stats_flips = polish_flips;
  stats_final_kappa = mk;

  if (verbose)
    printf("  %d PALC steps, %d flips, max|κ|=%.2e (%s)\n",
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

  if (mk > 0.01) {
    stats_status = ValidationStatus::FAIL_KAPPA_NOT_CONVERGED;
    return pos;   // failed-but-inspectable positions
  }

  // 5–9. Validation: a returned polytope must satisfy three named properties,
  // ALL of which are required for valid output:
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
  //
  // ── SIMPLICITY ──────────────────────────────────────────────────────
  //   T̄(0) must be a simple polygonal tesselation with F ≥ 3.
  //   inessential_eps stays at the strict default (1e-7).  Loosening it
  //   risks falsely collapsing borderline simple edges with |θ − π|
  //   slightly above 1e-7 but well below "almost flat", which can fold
  //   a non-degenerate polytope into a drum-cap.  Multi-edges that
  //   legitimately collapsed to one polytope edge typically reach
  //   |θ − π| ≲ 1e-7 once Newton converges; if PALC stalls before that
  //   (κ residual > tol), the F ≥ 3 check catches the resulting
  //   degeneracy.
  stats_t0_simplicial = is_simplicial(D);   // diagnostic only
  vector<int> labels(D.nv);
  for (int v = 0; v < D.nv; v++) labels[v] = v;
  auto tbar = polytope_tesselation(D, r, labels);
  stats_tbar_n_cells = tbar.n_cells();
  stats_tbar_simple_polygonal = is_simple_polygonal(tbar);

  bool simple = stats_tbar_simple_polygonal && stats_tbar_n_cells >= 3;
  if (!simple) {
    stats_status = ValidationStatus::FAIL_NOT_SIMPLE;
    if (verbose)
      printf("  VALIDATION (simplicity) failed: T̄(0) simple_polygonal=%d, "
             "n_cells=%d (need F≥3), T(0) simplicial=%d (diagnostic).\n",
             stats_tbar_simple_polygonal, stats_tbar_n_cells,
             stats_t0_simplicial);
    return pos;
  }

  // Compute volume_norm = |V| / ⟨ℓ⟩³ for the well-formedness gate and
  // diagnostic stat.  (Reconstruction was performed up-front, before
  // the kappa-convergence gate, so pos is already populated.)
  double vol6 = 0;
  for (int f = 0; f < D.nf; f++) {
    if (D.f_he[f] < 0) continue;
    int ha = D.f_he[f];
    int hb = D.he_next[ha];
    int hc = D.he_next[hb];
    vol6 += pos[D.he_origin[ha]].dot(
              pos[D.he_origin[hb]].cross(pos[D.he_origin[hc]]));
  }
  double volume = std::abs(vol6) / 6.0;
  double sum_l = 0; int n_e = 0;
  for (int h = 0; h < D.nh; h += 2) {
    if (!D.alive(h)) continue;
    sum_l += D.he_length[h];
    n_e++;
  }
  double mean_l = (n_e > 0) ? sum_l / n_e : 1.0;
  stats_volume_norm = volume / (mean_l * mean_l * mean_l);

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
  bool well_formed_volume = std::isfinite(stats_volume_norm) &&
                              stats_volume_norm >= VOLUME_NORM_DEGENERATE;
  stats_polytope_no_self_intersect = !has_self_intersection(D, pos);
  // Volume gate first, then self-intersection gate.  Order is the order
  // of stats_status when multiple fail.
  if (!well_formed_volume) {
    stats_status = ValidationStatus::FAIL_VOLUME_DEGENERATE;
    if (verbose)
      printf("  VALIDATION (well-formedness/volume) failed: vol_norm=%.3e "
             "(threshold %.2e).\n",
             stats_volume_norm, VOLUME_NORM_DEGENERATE);
    return pos;
  }
  if (!stats_polytope_no_self_intersect) {
    stats_status = ValidationStatus::FAIL_SELF_INTERSECTING;
    if (verbose)
      printf("  VALIDATION (well-formedness/self-intersection) failed: "
             "two non-adjacent triangles cross in 3D.\n");
    return pos;
  }

  // ── CONVEXITY ───────────────────────────────────────────────────────
  //   Every non-face vertex on the inside of every face plane.
  //   Alexandrov's theorem requires the polytope to be convex; failure
  //   here indicates the reconstruction landed in a geometrically
  //   non-convex configuration.  Outward normal taken from the half-
  //   edge CCW convention; defensive precondition inside is_convex
  //   verifies signed volume is strictly positive.
  stats_polytope_convex = is_convex(D, pos);
  if (!stats_polytope_convex) {
    stats_status = ValidationStatus::FAIL_NOT_CONVEX;
    if (verbose)
      printf("  VALIDATION (convexity) failed: some vertex sticks out "
             "beyond a face plane.\n");
    return pos;
  }

  stats_status = ValidationStatus::OK;
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
    for (int v = 0; v < D.nv; v++) labels[v] = v;
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
    // For bigon edges (he_face[h] == he_face[h^1]), GCP::theta is not
    // meaningful and the edge is by definition not a flat-face diagonal
    // of any 2-face of P (it's part of a degenerate iji-bigon face).
    // Mark non-inessential.
    if (T.he_face[h] == T.he_face[h ^ 1]) continue;
    double theta = GCP::theta(T, r, h);
    if (std::isfinite(theta) && std::fabs(theta - M_PI) < eps) {
      tight[h]     = true;
      tight[h ^ 1] = true;
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
  // an isometric R³ embedding for a non-degenerate polytope.
  // Self-loops: any half-edge with origin == dest.
  // Multi-edges: two distinct edges with the same vertex pair.
  // Bigons: a face whose boundary has only 2 distinct edges (he_face[h] ==
  //   he_face[h^1] for some pair).
  std::map<std::pair<int,int>, int> pair_count;
  for (int h = 0; h < T.nh; h += 2) {
    if (!T.alive(h)) continue;
    int u = T.he_origin[h], v = T.dest(h);
    if (u == v) return false;                                 // self-loop
    if (T.he_face[h] == T.he_face[h ^ 1]) return false;       // bigon
    auto key = std::make_pair(std::min(u, v), std::max(u, v));
    if (++pair_count[key] > 1) return false;                  // multi-edge
  }
  return true;
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

  // Mean edge length, for relative tolerance.
  double sum_l = 0; int n_e = 0;
  for (int h = 0; h < T.nh; h += 2) {
    if (!T.alive(h)) continue;
    sum_l += T.he_length[h]; n_e++;
  }
  double mean_l = (n_e > 0) ? sum_l / n_e : 1.0;
  double abs_tol = tol * mean_l;

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
  double vol6 = 0;
  for (int f = 0; f < T.nf; f++) {
    if (T.f_he[f] < 0) continue;
    int ha = T.f_he[f], hb = T.he_next[ha], hc = T.he_next[hb];
    vol6 += pos[T.he_origin[ha]].dot(
              pos[T.he_origin[hb]].cross(pos[T.he_origin[hc]]));
  }
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
    int ha = T.f_he[f], hb = T.he_next[ha], hc = T.he_next[hb];
    int va = T.he_origin[ha], vb = T.he_origin[hb], vc = T.he_origin[hc];
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
    int ha = T.f_he[f], hb = T.he_next[ha], hc = T.he_next[hb];
    tris.push_back({T.he_origin[ha], T.he_origin[hb], T.he_origin[hc]});
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
