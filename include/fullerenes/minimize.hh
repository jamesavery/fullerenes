#pragma once

// =====================================================================
// minimize -- dimension-agnostic smooth unconstrained minimization
// building blocks:
//
//   * strong-Wolfe line search (bracketing + zoom, Nocedal & Wright
//     Algorithms 3.5 / 3.6, bisection-safeguarded interpolation)
//   * L-BFGS (limited-memory BFGS, two-loop recursion, ring buffer)
//
// The oracle is any callable
//
//     double fg(std::span<const double> x, std::span<double> g)
//
// returning E(x) and filling g = grad E(x).  State is a flat double
// span (3N coordinates are viewed as a reinterpreted span of 3N
// doubles by geometric callers).
//
// Convergence tests (either may be disabled):
//   * relative energy decrease  2|f1 - f0| <= ftol_rel (|f1|+|f0|+eps)
//     with eps = 1e-9 -- the same criterion as the legacy Fortran
//     SA_frprmn3d, so drivers replacing it inherit identical stopping
//     semantics;
//   * gradient infinity norm    ||g||_inf <= gtol_inf.
//
// max_iters is a safeguard against runaway loops, not a tuning knob:
// it should never be reached, and Outcome.converged == false reports
// the trip so callers can treat it as the defect it is.
// =====================================================================

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <limits>
#include <span>
#include <vector>

namespace minimize {

struct Options {
    int    max_iters = 20000;   // safeguard (never-reached budget)
    int    memory    = 10;      // L-BFGS history pairs
    double gtol_inf  = 0.0;     // stop when ||g||_inf <= gtol_inf (0 = off)
    double gtol_2    = 0.0;     // stop when ||g||_2   <= gtol_2   (0 = off)
    double ftol_rel  = 1e-12;   // stop on relative energy stagnation (0 = off)
    double wolfe_c1  = 1e-4;    // Armijo (sufficient decrease)
    double wolfe_c2  = 0.9;     // curvature
    double max_step  = 0.0;     // ||step||_inf cap per iteration (0 = off).
                                // Bounds the line search like a trust
                                // radius: far from the basin the model is
                                // meaningless and unbounded steps jump into
                                // tangled minima; at the cap an Armijo-only
                                // point is accepted (curvature waived).
};

struct Outcome {
    double f         = 0;       // final value
    double gmax      = 0;       // final ||g||_inf
    int    iters     = 0;       // L-BFGS iterations taken
    int    n_evals   = 0;       // oracle calls (energy+gradient)
    bool   converged = false;   // a tolerance fired (vs safeguard trip)
};

namespace detail {

inline double dot(std::span<const double> a, std::span<const double> b) {
    double s = 0;
    for (size_t i = 0; i < a.size(); ++i) s += a[i] * b[i];
    return s;
}

inline double inf_norm(std::span<const double> a) {
    double m = 0;
    for (double v : a) m = std::max(m, std::fabs(v));
    return m;
}

inline void axpy(std::span<double> y, double a, std::span<const double> x) {
    for (size_t i = 0; i < y.size(); ++i) y[i] += a * x[i];
}

// Strong-Wolfe line search on phi(t) = f(x0 + t d).
//
// In:  x0 (start point), d (descent direction), f0 = phi(0),
//      dg0 = phi'(0) < 0, t first trial step, t_max step cap
//      (infinity = unbounded); at a capped step satisfying Armijo the
//      curvature condition is waived (bounded line search).
// Out (on success): x = x0 + t* d, g = grad f(x), f_out = f(x); true.
// On failure x, g, f_out hold the best Armijo-satisfying point seen if
// any (zoom's last iterate); returns false only when no acceptable
// point was found within the evaluation budget.
template <typename FG>
bool wolfe_search(FG&& fg,
                  std::span<const double> x0, std::span<const double> d,
                  double f0, double dg0, double t_init, double t_max,
                  std::span<double> x, std::span<double> g,
                  double& f_out, int& n_evals,
                  double c1, double c2)
{
    const size_t n = x0.size();
    const int    max_evals = 64;

    auto phi = [&](double t, double& dphi) -> double {
        for (size_t i = 0; i < n; ++i) x[i] = x0[i] + t * d[i];
        const double fv = fg(std::span<const double>(x.data(), n), g);
        ++n_evals;
        dphi = dot(std::span<const double>(g.data(), n), d);
        return fv;
    };

    // Zoom phase (Alg 3.6): maintain bracket [t_lo, t_hi] with
    // phi(t_lo) the lowest Armijo-satisfying value seen.
    auto zoom = [&](double t_lo, double f_lo, double dg_lo,
                    double t_hi, int evals_left) -> bool {
        for (int i = 0; i < evals_left; ++i) {
            // Quadratic interpolation using (t_lo, f_lo, dg_lo) and the
            // bisection midpoint as safeguard.
            double t = 0.5 * (t_lo + t_hi);
            double dphi = 0;
            const double fv = phi(t, dphi);
            if (!std::isfinite(fv) || fv > f0 + c1 * t * dg0 || fv >= f_lo) {
                t_hi = t;
                continue;
            }
            if (std::fabs(dphi) <= -c2 * dg0) { f_out = fv; return true; }
            if (dphi * (t_hi - t_lo) >= 0) t_hi = t_lo;
            t_lo = t; f_lo = fv; dg_lo = dphi;
        }
        // Budget exhausted: accept the best Armijo point if we have one.
        if (t_lo > 0) {
            double dphi = 0;
            f_out = phi(t_lo, dphi);
            return true;
        }
        return false;
    };

    double t_prev = 0, f_prev = f0, dg_prev = dg0;
    double t = std::min(t_init, t_max);
    for (int i = 0; i < max_evals; ++i) {
        double dphi = 0;
        const double fv = phi(t, dphi);
        if (!std::isfinite(fv)
            || fv > f0 + c1 * t * dg0
            || (i > 0 && fv >= f_prev))
            return zoom(t_prev, f_prev, dg_prev, t, max_evals - i - 1);
        if (std::fabs(dphi) <= -c2 * dg0) { f_out = fv; return true; }
        if (dphi >= 0)
            return zoom(t, fv, dphi, t_prev, max_evals - i - 1);
        if (t >= t_max) { f_out = fv; return true; }   // capped: Armijo holds
        t_prev = t; f_prev = fv; dg_prev = dphi;
        t = std::min(2 * t, t_max);
        if (t > 1e8) return false;
    }
    return false;
}

}  // namespace detail

// L-BFGS with strong-Wolfe line search.  Minimizes fg in place on x.
template <typename FG>
Outcome lbfgs(FG&& fg, std::span<double> x, const Options& opt = {})
{
    const size_t n = x.size();
    const int    M = std::max(1, opt.memory);
    const double EPS = 1e-9;    // Fortran SA_frprmn3d's zero-energy rectifier

    std::vector<double> g(n), gp(n), xp(n), d(n), q(n);
    std::vector<std::vector<double>> S(M, std::vector<double>(n));
    std::vector<std::vector<double>> Y(M, std::vector<double>(n));
    std::vector<double> rho(M), alpha(M);

    Outcome out;
    double f = fg(std::span<const double>(x.data(), n), std::span<double>(g));
    ++out.n_evals;

    int stored = 0, head = 0;   // ring buffer of (s, y) pairs
    for (int k = 0; k < opt.max_iters; ++k) {
        out.iters = k;
        out.gmax  = detail::inf_norm(g);
        if (opt.gtol_inf > 0 && out.gmax <= opt.gtol_inf) {
            out.converged = true;
            break;
        }
        if (opt.gtol_2 > 0 && std::sqrt(detail::dot(g, g)) <= opt.gtol_2) {
            out.converged = true;
            break;
        }

        // Descent direction d = -H g by the two-loop recursion.
        std::copy(g.begin(), g.end(), q.begin());
        for (int i = 0; i < stored; ++i) {                    // newest -> oldest
            const int j = ((head - 1 - i) % M + M) % M;
            alpha[j] = rho[j] * detail::dot(S[j], q);
            detail::axpy(q, -alpha[j], Y[j]);
        }
        if (stored > 0) {
            const int jn = ((head - 1) % M + M) % M;
            const double gamma = 1.0 / (rho[jn] * detail::dot(Y[jn], Y[jn]));
            for (double& v : q) v *= gamma;
        }
        for (int i = stored - 1; i >= 0; --i) {               // oldest -> newest
            const int j = ((head - 1 - i) % M + M) % M;
            const double beta = rho[j] * detail::dot(Y[j], q);
            detail::axpy(q, alpha[j] - beta, S[j]);
        }
        for (size_t i = 0; i < n; ++i) d[i] = -q[i];

        double dg = detail::dot(d, g);
        if (!(dg < 0)) {                 // not a descent direction: restart
            for (size_t i = 0; i < n; ++i) d[i] = -g[i];
            dg = -detail::dot(g, g);
            stored = 0;
            if (!(dg < 0)) { out.converged = true; break; }   // g == 0 exactly
        }

        // Line search from (xp, gp) = current point.
        std::copy(x.begin(), x.end(), xp.begin());
        std::copy(g.begin(), g.end(), gp.begin());
        const double t0 = (stored == 0)
                            ? std::min(1.0, 1.0 / std::max(out.gmax, 1e-12))
                            : 1.0;
        const double t_max = (opt.max_step > 0)
                               ? opt.max_step / std::max(detail::inf_norm(d), 1e-300)
                               : std::numeric_limits<double>::infinity();
        double f_new = f;
        const bool ok = detail::wolfe_search(fg, xp, d, f, dg,
                                             std::min(t0, t_max), t_max,
                                             x, std::span<double>(g), f_new,
                                             out.n_evals,
                                             opt.wolfe_c1, opt.wolfe_c2);
        if (!ok) {
            // Restore the pre-search point; a failed search from a
            // steepest-descent restart means we cannot make progress.
            std::copy(xp.begin(), xp.end(), x.begin());
            std::copy(gp.begin(), gp.end(), g.begin());
            if (stored == 0) break;
            stored = 0;
            continue;
        }

        // Fortran-compatible relative-energy stagnation test.
        if (opt.ftol_rel > 0
            && 2.0 * std::fabs(f_new - f)
                 <= opt.ftol_rel * (std::fabs(f_new) + std::fabs(f) + EPS)) {
            f = f_new;
            out.converged = true;
            ++out.iters;
            break;
        }

        // Store the (s, y) pair when it carries curvature information.
        auto& s = S[head];
        auto& y = Y[head];
        double sn2 = 0, yn2 = 0;
        for (size_t i = 0; i < n; ++i) {
            s[i] = x[i] - xp[i];
            y[i] = g[i] - gp[i];
            sn2 += s[i] * s[i];
            yn2 += y[i] * y[i];
        }
        const double sy = detail::dot(s, y);
        if (sy > 1e-10 * std::sqrt(sn2 * yn2)) {
            rho[head] = 1.0 / sy;
            head = (head + 1) % M;
            stored = std::min(stored + 1, M);
        }
        f = f_new;
    }

    out.f    = f;
    out.gmax = detail::inf_norm(g);
    return out;
}

}  // namespace minimize
