#pragma once

// =====================================================================
// optim::LineSearch -- the line-search globalization paradigm.  Owns
// the outer iteration loop; the Step policy (SteepestDescent, LBFGS,
// PR-CG) proposes descent directions, and one of two acceptance
// conditions shrinks the step along them:
//
//   STRONG_WOLFE -- a line-for-line port of minimize::detail::
//     wolfe_search (Nocedal-Wright Alg 3.5/3.6, bisection-safeguarded
//     zoom, 64-eval budget, capped step accepts Armijo-only).  The
//     oracle is the combined energy+gradient call, evaluated at every
//     trial point.  Bit-identical trajectories to minimize::lbfgs are
//     a gate requirement of the wu re-expression.
//
//   ARMIJO_QUAD -- Deltahedron's backtracking Armijo with quadratic
//     interpolation (alpha_q = -slope a^2 / 2(E_t - E - slope a),
//     clamped to [0.1a, 0.5a], <= 60 backtracks, energy-only trials);
//     the gradient is evaluated once at the accepted point.
//
// The loop reproduces minimize::lbfgs's control flow exactly when the
// extra Criteria (step_tol, domain predicates, work budget, stagnation
// window) are disabled: gtol tests at the top, descent-restart with
// LBFGS_RESET, failed-search restore/restart, and the Fortran-
// compatible relative-energy stagnation test required on two
// consecutive iterations.
// =====================================================================

#include "fullerenes/optim/core.hh"
#include "fullerenes/optim/preconditioner.hh"

#include <algorithm>
#include <cmath>
#include <limits>
#include <span>
#include <vector>

namespace optim {

enum class LineSearchCond { STRONG_WOLFE, ARMIJO_QUAD };

namespace detail {

// Strong-Wolfe line search on phi(t) = f(x0 + t d) -- the
// minimize::detail::wolfe_search algorithm, verbatim (bisection zoom,
// capped-step Armijo acceptance).  On success x, g, f_out hold the
// accepted point; on failure the best Armijo point seen, and false only
// when none was found within the budget.
// @anchor  optim-wolfe-search
// @pre     dg0 == dot(g(x0), d) < 0 (d is a descent direction)
// @pre     0 < t_init <= t_max
// @post    true  => Armijo holds at the returned point (and curvature
//          too, unless the step was capped at t_max or the zoom budget
//          ran out on a bracketed Armijo point)
// @post    false => x, g, f_out are unspecified; caller restores
// @variant 64-evaluation budget shared by bracket and zoom
template <typename FG>
bool wolfe_search(FG&& fg, std::span<const double> x0,
                  std::span<const double> d, double f0, double dg0,
                  double t_init, double t_max, std::span<double> x,
                  std::span<double> g, double& f_out,
                  double c1, double c2) {
  const std::size_t n = x0.size();
  const int max_evals = 64;

  auto phi = [&](double t, double& dphi) -> double {
    for (std::size_t i = 0; i < n; ++i) x[i] = x0[i] + t * d[i];
    const double fv = fg(std::span<const double>(x.data(), n), g);
    dphi = la::dot(std::span<const double>(g.data(), n), d);
    return fv;
  };

  auto zoom = [&](double t_lo, double f_lo, double t_hi,
                  int evals_left) -> bool {
    for (int i = 0; i < evals_left; ++i) {
      double t = 0.5 * (t_lo + t_hi);
      double dphi = 0;
      const double fv = phi(t, dphi);
      if (!std::isfinite(fv) || fv > f0 + c1 * t * dg0 || fv >= f_lo) {
        t_hi = t;
        continue;
      }
      if (std::fabs(dphi) <= -c2 * dg0) { f_out = fv; return true; }
      if (dphi * (t_hi - t_lo) >= 0) t_hi = t_lo;
      t_lo = t; f_lo = fv;
    }
    if (t_lo > 0) {            // budget exhausted: best Armijo point
      double dphi = 0;
      f_out = phi(t_lo, dphi);
      return true;
    }
    return false;
  };

  double t_prev = 0, f_prev = f0;
  double t = std::min(t_init, t_max);
  for (int i = 0; i < max_evals; ++i) {
    double dphi = 0;
    const double fv = phi(t, dphi);
    if (!std::isfinite(fv)
        || fv > f0 + c1 * t * dg0
        || (i > 0 && fv >= f_prev))
      return zoom(t_prev, f_prev, t, max_evals - i - 1);
    if (std::fabs(dphi) <= -c2 * dg0) { f_out = fv; return true; }
    if (dphi >= 0)
      return zoom(t, fv, t_prev, max_evals - i - 1);
    if (t >= t_max) { f_out = fv; return true; }   // capped: Armijo holds
    t_prev = t; f_prev = fv;
    t = std::min(2 * t, t_max);
    if (t > 1e8) return false;
  }
  return false;
}

// Deltahedron's backtracking Armijo with quadratic interpolation:
// energy-only trials, alpha clamped to [0.1a, 0.5a]; gradient evaluated
// once at the accepted point.  Like the parent, the final alpha is
// accepted even if the Armijo test never fired (the caller's stagnation
// logic then stops the descent).
// @anchor  optim-armijo-quad-search
// @pre     dg0 == dot(g(x0), d) < 0 (d is a descent direction)
// @post    x = x0 + alpha d and f_out = f(x), g = grad f(x) for the
//          final alpha; false iff f_out is non-finite
// @variant 60-backtrack budget
template <typename EOnly, typename FG>
bool armijo_quad_search(EOnly&& e_only, FG&& fg, std::span<const double> x0,
                        std::span<const double> d, double f0, double dg0,
                        double t_init, std::span<double> x,
                        std::span<double> g, double& f_out, double c1) {
  const std::size_t n = x0.size();
  double alpha = t_init;
  double f_trial = f0;
  for (int ls = 0; ls < 60; ++ls) {
    for (std::size_t i = 0; i < n; ++i) x[i] = x0[i] + alpha * d[i];
    f_trial = e_only(std::span<const double>(x.data(), n));
    if (f_trial <= f0 + c1 * alpha * dg0) break;
    const double denom = 2.0 * (f_trial - f0 - dg0 * alpha);
    if (denom > 1e-30) {
      const double alpha_q = -dg0 * alpha * alpha / denom;
      alpha = std::max(0.1 * alpha, std::min(0.5 * alpha, alpha_q));
    } else {
      alpha *= 0.5;
    }
  }
  for (std::size_t i = 0; i < n; ++i) x[i] = x0[i] + alpha * d[i];
  f_out = fg(std::span<const double>(x.data(), n), g);
  return std::isfinite(f_out);
}

}  // namespace detail

template <LineSearchStep Step, class Precond = Identity>
struct LineSearch {
  Step           step{};
  Precond        M{};
  LineSearchCond cond     = LineSearchCond::STRONG_WOLFE;
  double         c1       = 1e-4;   // Armijo (sufficient decrease)
  double         c2       = 0.9;    // curvature (strong Wolfe only)
  double         max_step = 0;      // ||step||_inf cap per iteration (0 = off)

  // Minimize prob.E in place on x.
  // @anchor  optim-linesearch-solve
  // @pre     x.size() is the model's DOF count (3N for geometric models)
  // @post    result.outcome names why iteration stopped; x holds the
  //          final iterate and result.f/gmax its value/gradient norm
  // @variant stop.max_iters outer iterations (a safeguard, never a knob)
  template <SmoothModel Model>
  Result solve(SeqCtx ctx, const Problem<Model>& prob, std::span<double> x,
               const Projection& project = {}) const {
    const std::size_t n = x.size();
    const double EPS = 1e-9;        // the Fortran zero-energy rectifier
    const Criteria& stop = prob.stop;

    Result out;
    auto fg = [&](std::span<const double> xv, std::span<double> gv) {
      ++out.n_grad;
      const double e = prob.E.gradient(ctx, xv, gv);
      prob.gauge.project(gv);
      return e;
    };
    auto e_only = [&](std::span<const double> xv) {
      ++out.n_energy;
      return prob.E.energy(ctx, xv);
    };

    std::vector<double> g(n), gp(n), xp(n), d(n), s(n), y(n);
    typename Step::State st(n, step);

    double f = fg(x, g);

    int n_stag = 0;        // consecutive ftol_rel events (minimize.hh)
    int stag_count = 0;    // iterations without meaningful decrease (Deltahedron)
    for (int k = 0; k < stop.max_iters; ++k) {
      out.iters = k;
      out.gmax  = la::inf_norm(g);
      if (auto stopped = converged_or_stop(stop, out.gmax, g, x, stag_count,
                                           out.work<Model>())) {
        out.outcome = *stopped;
        break;
      }

      // Descent direction from the step policy; restart on non-descent.
      step.direction(ctx, st, g, d);
      prob.gauge.project(d);
      double dg = la::dot(d, g);
      if (!(dg < 0)) {
        for (std::size_t i = 0; i < n; ++i) d[i] = -g[i];
        prob.gauge.project(d);
        dg = -la::dot(g, g);
        step.reset(st);
        out.flags |= Result::LBFGS_RESET;
        if (!(dg < 0)) {               // g == 0 exactly
          out.outcome = Outcome::CONVERGED;
          break;
        }
      }

      // Line search from (xp, gp) = current point.
      std::copy(x.begin(), x.end(), xp.begin());
      std::copy(g.begin(), g.end(), gp.begin());
      const double t0 = step.first_trial(st, out.gmax);
      const double t_max = (max_step > 0)
          ? max_step / std::max(la::inf_norm(d), 1e-300)
          : std::numeric_limits<double>::infinity();

      double f_new = f;
      bool ok;
      if (cond == LineSearchCond::STRONG_WOLFE)
        ok = detail::wolfe_search(fg, xp, d, f, dg, std::min(t0, t_max), t_max,
                                  x, std::span<double>(g), f_new, c1, c2);
      else
        ok = detail::armijo_quad_search(e_only, fg, xp, d, f, dg,
                                        std::min(t0, t_max), x,
                                        std::span<double>(g), f_new, c1);
      if (!ok) {
        // Restore the pre-search point; a failed search from a fresh
        // restart means no progress is possible.
        std::copy(xp.begin(), xp.end(), x.begin());
        std::copy(gp.begin(), gp.end(), g.begin());
        if (!step.has_history(st)) {
          out.outcome = Outcome::STEP_FAILED;
          break;
        }
        step.reset(st);
        continue;
      }

      // Fortran-compatible relative-energy stagnation test, required on
      // two consecutive iterations (minimize.hh semantics).
      if (stop.ftol_rel > 0
          && 2.0 * std::fabs(f_new - f)
               <= stop.ftol_rel * (std::fabs(f_new) + std::fabs(f) + EPS)) {
        if (++n_stag >= 2) {
          f = f_new;
          out.outcome = Outcome::CONVERGED;
          ++out.iters;
          break;
        }
      } else
        n_stag = 0;

      // Deltahedron-style stagnation window (off unless enabled).
      if (stop.stagnation_window > 0) {
        if (f - f_new > 1e-15 * std::max(1.0, std::fabs(f))) stag_count = 0;
        else ++stag_count;
      }

      for (std::size_t i = 0; i < n; ++i) {
        s[i] = x[i] - xp[i];
        y[i] = g[i] - gp[i];
      }
      step.accepted(ctx, st, s, y, g);
      f = f_new;

      if (project) {
        if (project(ctx, x) > 0) f = fg(x, g);   // re-anchor after projection
      }
      if (stop.step_tol > 0 && la::inf_norm(s) <= stop.step_tol) {
        out.outcome = Outcome::CONVERGED;
        ++out.iters;
        break;
      }
    }

    out.f    = f;
    out.gmax = la::inf_norm(g);
    return out;
  }
};

}  // namespace optim
