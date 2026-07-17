#include "fullerenes/wu-forcefield.hh"

#include <cmath>
#include <cstdio>
#include <deque>
#include <stdexcept>

using std::vector;

// ---------------------------------------------------------------------------
// Parameters (opt-standalone.f default_force_parameters, with SA_OptFF's unit
// conversions applied: angles/dihedrals deg -> rad, force constants N/m ->
// kJ/mol/A^2 by *6.02214129 = 1e-20 * N_A).
// ---------------------------------------------------------------------------

namespace {

constexpr double deg = M_PI / 180.0;
constexpr double unitconv = 6.02214129;  // 1e-20 * 6.02214129e20, as in SA_OptFF

// Wu (methods 1/2): force(1..9) = rp, rh, ap, ah, frp, frh, fap, fah, fco.
// The Wu field does not distinguish 1- from 2-pentagon bonds.
struct WuParams {
  double rp = 1.455, rh = 1.391;                       // A, from solid-state
  double ap = 108 * deg, ah = 120 * deg;
  double frp = 390.7 * unitconv, frh = 499.7 * unitconv;  // Ceulemans, Fowler
  double fap = 47.88 * 1.45 * 1.45 * unitconv;
  double fah = 80.86 * 1.45 * 1.37 * unitconv;
};

// Extended Wu (methods 3/4): force(1..19). Bond/dihedral classes are indexed
// below by the pentagon count np of the term's incident faces.
struct ExtWuParams {
  double r0[3] = {1.401, 1.458, 1.479};                // np = 0(hh), 1(hp), 2(pp)
  double kr[3] = {450 * unitconv, 390 * unitconv, 260 * unitconv};
  double a0[2] = {120 * deg, 108 * deg};               // hexagon, pentagon corner
  double ka[2] = {100 * unitconv, 100 * unitconv};
  double d0[4] = {0, 23.49 * deg, 29.20 * deg, 37.38 * deg};  // np = 0..3, "ideal_dihedral"
  double kd[4] = {270 * unitconv, 85 * unitconv, 65 * unitconv, 35 * unitconv};
};

}  // namespace

// ---------------------------------------------------------------------------
// Term compilation
// ---------------------------------------------------------------------------

WuForceField::WuForceField(const FullereneGraphView& g, int method) : N(g.N)
{
  if (method < 1 || method > 4)
    throw std::invalid_argument("WuForceField: method must be 1..4 (legacy iopt); "
                                "the hyperbolic variants 5/6 were never reachable "
                                "through this entry and are not ported");
  const bool extended = (method >= 3);
  const WuParams wu;
  const ExtWuParams xw;

  // legacy default fcoulomb = 0; methods 2/4 would activate it if nonzero.
  k_coulomb = 0;

  // --- Bonds: pentagon count of the edge's two faces (SA_get_edges) ---
  bonds.reserve(3 * size_t(N) / 2);
  for (const edge_t& e : g.undirected_edges()) {
    const int np = 12 - g.face_size(e.first, e.second) - g.face_size(e.second, e.first);
    if (extended)
      bonds.push_back({e.first, e.second, xw.r0[np], xw.kr[np]});
    else  // Wu: any pentagon-touching bond uses the pentagon parameters
      bonds.push_back({e.first, e.second, np ? wu.rp : wu.rh, np ? wu.frp : wu.frh});
  }

  // --- Angles: one per face corner, apex in the middle (SA_get_corners) ---
  angles.reserve(3 * size_t(N));
  for (const face_t& f : g.compute_faces(6)) {
    const int L = f.size(), pent = (L == 5);
    for (int j = 0; j < L; j++) {
      const node_t a = f[j], b = f[(j + 1) % L], c = f[(j + 2) % L];
      if (extended)
        angles.push_back({a, b, c, xw.a0[pent], xw.ka[pent]});
      else
        angles.push_back({a, b, c, pent ? wu.ap : wu.ah, pent ? wu.fap : wu.fah});
    }
  }

  // --- Dihedrals (extended only): one per vertex (SA_get_dihedrals) ---
  //
  //      t   B   s        With CCW neighbours (r,s,t) of u, face A lies
  //        \   /          between r and s, B between s and t, C between t
  //      C   u   A        and r. The term is the chain dihedral
  //          |            D(u, n[rot], n[rot+1], n[rot+2]) where rot keeps
  //          r            the legacy start: when exactly one face differs
  //                       in size from the other two, start at the
  //                       neighbour OPPOSITE that face; else start at r.
  if (extended) {
    dihedrals.reserve(N);
    for (node_t u = 0; u < N; u++) {
      const auto nb = g.nbrs(u);
      const node_t n3[3] = {nb[0], nb[1], nb[2]};
      // face_i is bounded by n3[i] and n3[(i+1)%3]; its size via the unique
      // shortest cycle through n3[(i+1)%3] -> u -> n3[i] (girth 5, faces <= 6),
      // exactly as the legacy SA_get_dihedrals probed through get_face.
      int l[3], np = 0;
      for (int i = 0; i < 3; i++) {
        l[i] = g.shortest_cycle({n3[(i + 1) % 3], u, n3[i]}, 6).size();
        np += (l[i] == 5);
      }
      int rot = 0;
      if (np == 1 || np == 2) {
        const int odd = (np == 1);  // the face size occurring once: 5 if np==1, 6 if np==2
        for (int i = 0; i < 3; i++)
          if ((l[i] == 5) == odd) rot = (i + 2) % 3;  // neighbour opposite the distinct face
      }
      dihedrals.push_back({u, n3[rot], n3[(rot + 1) % 3], n3[(rot + 2) % 3],
                           xw.d0[np], xw.kd[np]});
    }
  }
}

// ---------------------------------------------------------------------------
// Energy and gradient
// ---------------------------------------------------------------------------

double WuForceField::energy(std::span<const coord3d> x) const
{
  double E = 0;
  for (const Bond& t : bonds) {
    const double r = (x[t.a] - x[t.b]).norm();
    E += t.k * (r - t.r0) * (r - t.r0);
  }
  for (const Angle& t : angles) {
    const double a = coord3d::angle(x[t.a] - x[t.b], x[t.c] - x[t.b]);
    E += t.k * (a - t.a0) * (a - t.a0);
  }
  for (const Dihedral& t : dihedrals) {
    const double d = coord3d::dihedral(x[t.b] - x[t.a], x[t.c] - x[t.a], x[t.d] - x[t.a]);
    E += t.k * (d - t.d0) * (d - t.d0);
  }
  if (k_coulomb != 0)
    for (int u = 0; u < N; u++) E += k_coulomb / x[u].norm();
  return 0.5 * E;
}

void WuForceField::gradient(std::span<const coord3d> x, std::span<coord3d> grad) const
{
  for (int u = 0; u < N; u++) grad[u] = coord3d(0, 0, 0);

  for (const Bond& t : bonds) {
    const coord3d ab = x[t.a] - x[t.b];
    const double r = ab.norm();
    const coord3d dE = ab * ((t.k * (r - t.r0)) / r);
    grad[t.a] += dE;
    grad[t.b] -= dE;
  }
  for (const Angle& t : angles) {
    const coord3d a = x[t.a] - x[t.b], c = x[t.c] - x[t.b];
    const double ang = coord3d::angle(a, c);
    // At 0 or pi the derivative is singular; the legacy DANGLE zeroes it
    // there, and so do we (the term's energy is still counted above).
    if (ang < 1e-7 || ang > M_PI - 1e-7) continue;
    coord3d da, dc;
    coord3d::dangle(a, c, da, dc);
    const double s = t.k * (ang - t.a0);
    grad[t.a] += da * s;
    grad[t.c] += dc * s;
    grad[t.b] -= (da + dc) * s;
  }
  for (const Dihedral& t : dihedrals) {
    const coord3d b = x[t.b] - x[t.a], c = x[t.c] - x[t.a], d = x[t.d] - x[t.a];
    const double dih = coord3d::dihedral(b, c, d);
    coord3d db, dc, dd;
    coord3d::ddihedral(b, c, d, db, dc, dd);
    // Collinear triples make the dihedral derivative singular; treat like
    // the angle singularity above (energy counted, derivative dropped).
    if (!std::isfinite(db[0] + db[1] + db[2] + dc[0] + dc[1] + dc[2] +
                       dd[0] + dd[1] + dd[2]))
      continue;
    const double s = t.k * (dih - t.d0);
    grad[t.b] += db * s;
    grad[t.c] += dc * s;
    grad[t.d] += dd * s;
    grad[t.a] -= (db + dc + dd) * s;  // translation invariance
  }
  // dE/dx of 1/2 k sum 1/|x_u|. (The legacy dwu/dextwu added +k x/r^3 here,
  // sign-inconsistent with its own energy -- irrelevant in practice since
  // k_coulomb defaults to 0; we implement the energy-consistent gradient.)
  if (k_coulomb != 0)
    for (int u = 0; u < N; u++) {
      const double r2 = x[u].dot(x[u]);
      grad[u] -= x[u] * (0.5 * k_coulomb / (r2 * std::sqrt(r2)));
    }
}

// ---------------------------------------------------------------------------
// L-BFGS minimizer (self-contained; two-loop recursion, m=10, Armijo
// backtracking). Written fresh -- deliberately NOT a port of the legacy
// frprmn3d/brent3d, which are Numerical-Recipes-derived. Only the minimum
// matters for parity: the force field defines it, the minimizer just finds it.
// ---------------------------------------------------------------------------

namespace {

using field_t = vector<coord3d>;  // one 3-vector per atom

double dot(const field_t& a, const field_t& b)
{
  double s = 0;
  for (size_t i = 0; i < a.size(); i++) s += a[i].dot(b[i]);
  return s;
}

double max_abs(const field_t& a)
{
  double m = 0;
  for (const coord3d& v : a)
    for (int i = 0; i < 3; i++) m = std::max(m, std::abs(v[i]));
  return m;
}

}  // namespace

double WuForceField::optimize(std::span<coord3d> x, double ftol, int max_iter) const
{
  const int m = 10;        // L-BFGS history length
  const double c1 = 1e-4;  // Armijo sufficient-decrease constant

  field_t p(x.begin(), x.end()), g(N), p_new(N), g_new(N), dir(N);
  gradient(p, g);
  double E = energy(p);
  std::deque<std::pair<field_t, field_t>> hist;  // (s, y)

  int it = 0, n_converged = 0;
  for (; it < max_iter; it++) {
    if (max_abs(g) < 1e-12) break;  // stationary to machine precision

    // Two-loop recursion: dir = -H*g with H from the (s,y) history.
    dir = g;
    vector<double> alpha(hist.size());
    for (int i = int(hist.size()) - 1; i >= 0; i--) {
      const auto& [s, y] = hist[i];
      alpha[i] = dot(s, dir) / dot(y, s);
      for (int u = 0; u < N; u++) dir[u] -= y[u] * alpha[i];
    }
    if (!hist.empty()) {
      const auto& [s, y] = hist.back();
      const double gamma = dot(y, s) / dot(y, y);
      for (int u = 0; u < N; u++) dir[u] *= gamma;
    }
    for (int i = 0; i < int(hist.size()); i++) {
      const auto& [s, y] = hist[i];
      const double beta = dot(y, dir) / dot(y, s);
      for (int u = 0; u < N; u++) dir[u] += s[u] * (alpha[i] - beta);
    }
    for (int u = 0; u < N; u++) dir[u] = -dir[u];

    double gd = dot(g, dir);
    if (gd >= 0) {  // not a descent direction (degenerate history): restart
      hist.clear();
      const double gnorm = std::sqrt(dot(g, g));
      for (int u = 0; u < N; u++) dir[u] = g[u] * (-1.0 / gnorm);
      gd = dot(g, dir);
    }

    // Armijo backtracking from the natural L-BFGS step, with the trial step
    // capped so no atom moves more than max_move at once: descending e.g. a
    // grossly inflated start geometry must follow a continuous path, not
    // tunnel through a ridge into a tangled basin (the Fortran's Brent line
    // minimization had the same effect implicitly).
    const double max_move = 0.3;  // Angstroem
    double dir_max = 0;
    for (int u = 0; u < N; u++)
      dir_max = std::max(dir_max, std::sqrt(dir[u].dot(dir[u])));
    double step = hist.empty() ? std::min(1.0, 1.0 / std::sqrt(dot(g, g))) : 1.0;
    step = std::min(step, max_move / dir_max);
    double E_new = E;
    bool accepted = false;
    for (int bt = 0; bt < 50; bt++, step *= 0.5) {
      for (int u = 0; u < N; u++) p_new[u] = p[u] + dir[u] * step;
      E_new = energy(p_new);
      if (std::isfinite(E_new) && E_new <= E + c1 * step * gd) { accepted = true; break; }
    }
    if (!accepted) break;  // no decrease possible along dir: converged/stuck
    gradient(p_new, g_new);

    // Curvature update (skip degenerate pairs to keep H positive definite).
    field_t s(N), y(N);
    for (int u = 0; u < N; u++) {
      s[u] = p_new[u] - p[u];
      y[u] = g_new[u] - g[u];
    }
    if (dot(y, s) > 1e-12 * std::sqrt(dot(s, s) * dot(y, y))) {
      hist.emplace_back(std::move(s), std::move(y));
      if (int(hist.size()) > m) hist.pop_front();
    }

    p.swap(p_new);
    g.swap(g_new);
    const double E_old = E;
    E = E_new;

    // Legacy convergence criterion (relative energy change), required on two
    // consecutive iterations so a single slow step cannot stop the descent.
    if (2 * std::abs(E_old - E) <= ftol * (std::abs(E_old) + std::abs(E) + 1e-10)) {
      if (++n_converged >= 2) break;
    } else
      n_converged = 0;
  }
  if (it >= max_iter)
    fprintf(stderr, "WuForceField::optimize: iteration budget (%d) exhausted at "
            "E=%g, |grad|_inf=%g -- convergence, not the budget, is the problem\n",
            max_iter, E, max_abs(g));

  std::copy(p.begin(), p.end(), x.begin());
  return E;
}

int WuForceField::separate_coincident(std::span<coord3d> x) const
{
  const double too_close = 0.1, displacement = 0.5;  // optimize_other's constants
  int n_displaced = 0;
  for (int pass = 0; pass < 5; pass++) {  // re-scan: a displacement could re-collide
    int moved = 0;
    for (int u = 0; u < N; u++)
      for (int v = u + 1; v < N; v++)
        if ((x[u] - x[v]).norm() < too_close) {
          const coord3d d(displacement, displacement, displacement);
          x[u] += d;
          x[v] -= d;
          moved++;
        }
    n_displaced += moved;
    if (!moved) break;
  }
  return n_displaced;
}

vector<coord3d> wu_optimized_geometry(const FullereneGraphView& g,
                                      std::span<const coord3d> initial_geometry,
                                      int method, double ftol)
{
  WuForceField ff(g, method);
  vector<coord3d> x(initial_geometry.begin(), initial_geometry.end());
  ff.separate_coincident(x);
  ff.optimize(x, ftol);
  return x;
}
