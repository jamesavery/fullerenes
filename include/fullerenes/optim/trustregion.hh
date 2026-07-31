#pragma once

// =====================================================================
// optim::TrustRegion -- the trust-region globalization paradigm, and
// TrustRadius, THE rho-ratio radius update (previously written twice:
// delaunay_alexandrov's TrustRegion::update with floor 1e-14, and
// Deltahedron's inlined variant with floor 1e-14 L -- structurally
// identical, same 0.1 / 0.75 / 0.5 / x2 / x0.25 constants, differing
// only in the floor, which is therefore a parameter).
//
// The paradigm owns the outer loop; the step is a subproblem-solver
// policy dispatched on its oracle:
//
//   * lsq_oracle_tag (LevenbergBisect) -- model provides residual R and
//     the SQUARE system matrix J; the loop is the generic form of
//     delaunay_alexandrov's Newton::polish (minus the feasibility clip
//     and post-accept re-flip, which stay with the Alexandrov polish
//     functor by the signed-off "monolithic for now" decision).
//     Convergence on Criteria::rtol_inf (max |R| < tol).
//
//   * hvp_oracle_tag (SteihaugCG) -- model provides gradient + HVP; the
//     loop is Deltahedron's STEIHAUG block (minus the phase machinery,
//     which stays with the delta driver): subproblem, predicted
//     reduction -(g'z + 1/2 z'Hz), rho-acceptance, optional acceptance
//     veto (the convex-constraint variant: reject steps a caller-
//     supplied predicate refuses, e.g. convex->concave transitions),
//     stagnation and work-budget criteria.
//
// An invalid pairing (TrustRegion<LBFGS>) is rejected by the
// TrustRegionStep constraint; a valid step over a model lacking its
// oracle face (TrustRegion<SteihaugCG> over a gradient-only model)
// fails the solve() requires-clause.
// =====================================================================

#include "fullerenes/optim/core.hh"
#include "fullerenes/optim/preconditioner.hh"
#include "fullerenes/optim/steps/levenberg.hh"
#include "fullerenes/optim/steps/steihaug.hh"

#include <algorithm>
#include <cmath>
#include <functional>
#include <span>
#include <vector>

namespace optim {

// The rho-ratio accept/grow/shrink policy.  rho/accept/reject are split
// so an acceptance veto can force the reject branch AFTER the rho test
// without the radius having already grown (the parent's exact order);
// update() is the plain Alexandrov TrustRegion::update shape.
// @anchor optim-trust-radius
// @inv    floor <= Delta <= Delta_max between calls
// @post   update: true iff rho > 0.1; on true Delta grew (x2, capped)
//         iff rho > 0.75 and the step was near the boundary; on false
//         Delta shrank (x0.25, floored)
struct TrustRadius {
  double Delta;
  double Delta_max;
  double floor = 1e-14;

  static double rho(double actual, double predicted) {
    return (predicted > 1e-30) ? actual / predicted : -1;
  }
  void accept(double rho_val, double dnorm) {
    if (rho_val > 0.75 && dnorm > 0.5 * Delta)
      Delta = std::min(2.0 * Delta, Delta_max);
  }
  void reject() { Delta = std::max(Delta * 0.25, floor); }

  bool update(double actual, double predicted, double dnorm) {
    const double r = rho(actual, predicted);
    if (r > 0.1) { accept(r, dnorm); return true; }
    reject();
    return false;
  }
};

// Acceptance veto: a predicate on the trial point the paradigm consults
// after the rho test (Deltahedron's convex-constraint gate).  Returning
// false rejects the step (CVX_REJECTED).
using AcceptVeto = std::function<bool(std::span<const double> x_trial)>;

template <TrustRegionStep Step, class Precond = Identity>
struct TrustRegion {
  Step       step{};
  Precond    M{};
  double     Delta0      = 0;      // required: initial radius
  double     Delta_max   = 0;      // required: radius cap
  double     delta_floor = 1e-14;  // reject-shrink floor (1e-14 L for AET)
  AcceptVeto accept_veto;          // optional (constrained variant)

  // --- LSQ route: Levenberg over residual + square system matrix -----
  // @anchor  optim-trustregion-solve-lsq
  // @pre     Delta0 > 0 and Delta_max >= Delta0 (caller-set radii)
  // @post    result.outcome: CONVERGED iff max|R| < stop.rtol_inf fired;
  //          STEP_FAILED after 20 consecutive rejections
  // @variant stop.max_iters outer iterations
  template <class Model>
    requires std::same_as<typename Step::oracle, lsq_oracle_tag> && HasLSQ<Model>
  Result solve(SeqCtx ctx, const Problem<Model>& prob, std::span<double> x,
               const Projection& project = {}) const {
    const Criteria& stop = prob.stop;
    const std::size_t m = prob.E.residual_size();
    Result out;

    TrustRadius tr{Delta0, Delta_max, delta_floor};
    la::V R(m), R_trial(m);
    Mat J(m, m);
    std::vector<double> x_trial(x.size());
    int rejects = 0;

    for (int iter = 0; iter < stop.max_iters; ++iter) {
      out.iters = iter;
      prob.E.residual(ctx, x, R);
      ++out.n_energy;
      if (stop.rtol_inf > 0 && la::max_abs(R) < stop.rtol_inf) {
        out.outcome = Outcome::CONVERGED;
        break;
      }
      if (rejects > 20) {
        out.outcome = Outcome::STEP_FAILED;
        break;
      }

      const double E = la::energy(R);
      prob.E.system_matrix(ctx, x, J);
      ++out.n_grad;

      la::V delta = Step::solve_subproblem(J, R, tr.Delta);
      const double pred = Step::predicted_reduction(J, R, delta);

      for (std::size_t i = 0; i < x.size(); ++i) x_trial[i] = x[i] + delta[i];
      prob.E.residual(ctx, x_trial, R_trial);
      ++out.n_energy;
      const double E_trial = la::energy(R_trial);

      if (tr.update(E - E_trial, pred, la::norm(delta))) {
        std::copy(x_trial.begin(), x_trial.end(), x.begin());
        rejects = 0;
        if (project) project(ctx, x);
      } else {
        ++rejects;
        out.flags |= Result::STEP_REJECTED;
      }
    }

    prob.E.residual(ctx, x, R);
    ++out.n_energy;
    out.f    = la::energy(R);
    out.gmax = la::max_abs(R);
    return out;
  }

  // --- HVP route: Steihaug-Toint over gradient + Hessian-vector ------
  // @anchor  optim-trustregion-solve-hvp
  // @pre     Delta0 > 0 and Delta_max >= Delta0 (caller-set radii)
  // @post    result.outcome names why iteration stopped; x holds the
  //          final iterate
  // @variant stop.max_iters outer iterations
  template <class Model>
    requires std::same_as<typename Step::oracle, hvp_oracle_tag> && HasHVP<Model>
  Result solve(SeqCtx ctx, const Problem<Model>& prob, std::span<double> x,
               const Projection& project = {}) const {
    const std::size_t n = x.size();
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

    std::vector<double> g(n), z(n), Hz(n), x_trial(n);
    TrustRadius tr{Delta0, Delta_max, delta_floor};

    double f = fg(x, g);
    int stag_count = 0;
    for (int iter = 0; iter < stop.max_iters; ++iter) {
      out.iters = iter;
      out.gmax  = la::inf_norm(g);
      if (auto stopped = converged_or_stop(stop, out.gmax, g, x, stag_count,
                                           out.work<Model>())) {
        out.outcome = *stopped;
        break;
      }

      Step::solve_subproblem(ctx, prob.E, prob.gauge, x, g, tr.Delta, z,
                             out.flags, out.n_hv);

      // Predicted reduction -(g'z + 1/2 z'Hz).
      ++out.n_hv;
      prob.E.hvp(ctx, x, z, Hz);
      const double pred = -(la::dot(g, z) + 0.5 * la::dot(z, Hz));

      for (std::size_t i = 0; i < n; ++i) x_trial[i] = x[i] + z[i];
      const double E_trial = e_only(x_trial);

      const double rho = TrustRadius::rho(f - E_trial, pred);
      bool accepted = rho > 0.1;
      if (accepted && accept_veto && !accept_veto(x_trial)) {
        accepted = false;
        out.flags |= Result::CVX_REJECTED;
      }

      if (accepted) {
        const double f_old = f;
        std::copy(x_trial.begin(), x_trial.end(), x.begin());
        f = fg(x, g);
        tr.accept(rho, la::nrm2(z));
        if (project) {
          if (project(ctx, x) > 0) f = fg(x, g);
        }
        if (f_old - f > 1e-15 * std::max(1.0, std::fabs(f_old))) stag_count = 0;
        else ++stag_count;
      } else {
        tr.reject();
        ++stag_count;                        // rejected step = no progress
        out.flags |= Result::STEP_REJECTED;
      }
    }

    out.f    = f;
    out.gmax = la::inf_norm(g);
    return out;
  }
};

}  // namespace optim
