#pragma once

// =====================================================================
// optim::delta_optimize -- Deltahedron::optimize() re-expressed on the
// framework: the same control flow (two-phase E_flat schedule, the
// gmax_L convergence test, the 50-iteration stagnation window, the
// unified 1:2:7 work budget, the convex-constrained acceptance gate,
// the post-optimization hull-projection cleanup) driving framework
// components over the AET model:
//
//   CG       -- Polak-Ribiere + Armijo-quad (detail::armijo_quad_search)
//   LBFGS    -- optim::LBFGS two-loop + Armijo-quad
//   STEIHAUG -- SteihaugCG subproblem + TrustRadius (floor 1e-14 L)
//
// Known, accepted arithmetic deviations from the parent (gate is
// quality-DISTRIBUTION parity, bench_delta_reexpr, not trajectory
// identity): the L-BFGS two-loop is the minimize.hh transcription
// (gamma as 1/(rho y'y), curvature skip 1e-10 sqrt(|s|^2 |y|^2)) where
// the parent Deltahedron carries a deque variant (gamma as ys/yy, skip
// 1e-10 s's); both line searches start at alpha = 1 as the parent does.
//
// Convexity restoration (reflect/hull loops BETWEEN optimize calls)
// stays with the caller exactly as in the parent pipeline; the
// post-optimization strict-convexity cleanup lives here, as it does in
// the parent optimize().
//
// In-loop reflection (cfg.reflect_threshold): at the top of every
// iteration, vertices with h < -reflect_threshold L are mirrored
// through their neighbour-centroid plane BEFORE the step (the parent's
// pre-2026-03-08 "periodic reflection", the DESIGN S3.3 between-steps
// Projection realized at this driver's altitude).  It is the
// fold-escape operator for cold starts: a descent from a fragile
// (Tutte-sphere) start folds a pentagon cap on the way down and, left
// alone, converges into a CONVEX wrong-basin minimum that no
// post-convergence reflect/hull pass can touch (C50 ang_max 94.13 vs
// the 75-degree guard; claude-projects/optimize/validation/
// C50-ANGMAX-INVESTIGATION.md).  Reflecting the fold as it forms keeps
// the descent in the equilateral basin.  0 disables.
// =====================================================================

#include "fullerenes/optim/core.hh"
#include "fullerenes/optim/linesearch.hh"
#include "fullerenes/optim/trustregion.hh"
#include "fullerenes/optim/steps/lbfgs.hh"
#include "fullerenes/optim/steps/steihaug.hh"
#include "fullerenes/optim/models/aet.hh"

#include "fullerenes/graphview.hh"
#include "fullerenes/stats.hh"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <span>
#include <vector>

namespace optim {

enum class DeltaMethod { CG, LBFGS, STEIHAUG };

struct DeltaConfig {
  DeltaMethod method = DeltaMethod::LBFGS;
  double target_L = 0;            // 0 = mean initial edge length
  double grad_tol = 1e-10;        // on gmax_L = L max_v |grad_v|_2 (free v)
  std::vector<bool> fixed;        // frozen vertices
  long long max_work = 0;         // 0 = 400 N^2 (n_e + 2 n_g + 7 n_hv)
  double angle_tol = 0;           // domain success: max angle rel-err, convex
  double k_flat = 2.0;            // phase-1 flatness (0 = single phase)
  double k_conv = 0;
  bool convex_constraint = false; // STEIHAUG: reject convex->concave steps
  bool skip_post_reflect = false; // caller handles convexity cleanup
  double reflect_threshold = 0.05; // in-loop reflection of h < -thr L (0 = off)
  FILE* log = nullptr;            // parent-format diagnostic lines (~20/run)
};

struct DeltaResult {
  Outcome outcome = Outcome::BUDGET_EXHAUSTED;
  double E = 0, gmax_L = 0, angle_relerr = 0;
  int n_concave = 0;
  int iters = 0;
  int n_energy = 0, n_grad = 0, n_hv = 0;
  int n_reflected = 0;            // vertices mirrored by in-loop reflection
  uint32_t flags = 0;
  // The unified budget in AET eval-cost units (the model's 1:2:7).
  double work() const {
    return eval_costs<AET>::energy * n_energy
         + eval_costs<AET>::gradient * n_grad
         + eval_costs<AET>::hvp * n_hv;
  }
};

// @anchor  optim-delta-optimize
// @pre     D.points holds the initial geometry (mutated in place)
// @post    result.outcome: CONVERGED (gmax_L or angle criterion),
//          STAGNATED (50 iters without meaningful decrease) or
//          BUDGET_EXHAUSTED; D.points holds the final geometry,
//          hull-cleaned unless cfg.skip_post_reflect
// @variant the work budget (default 400 N^2 in 1:2:7 eval units)
inline DeltaResult delta_optimize(DeltahedronView<double> D,
                                  const DeltaConfig& cfg = {}) {
  const int N = D.N;
  const bool has_fixed = !cfg.fixed.empty();
  DeltaResult out;

  // --- Model + configuration (the parent's fixed constants) ----------
  // One edge-list computation per call (the model's cache), as the
  // parent had; the default target L is its mean initial length.
  AET model(D, cfg.target_L);
  double L = cfg.target_L;
  if (L <= 0) {
    double sum = 0;
    for (const edge_t& e : model.edges)
      sum += coord3d::dist(D.points[e.first], D.points[e.second]);
    L = sum / model.edges.size();
    model.L = L;
  }
  model.k_flat = cfg.k_flat;
  model.k_conv = cfg.k_conv;
  model.fixed = cfg.fixed;

  long long max_work = cfg.max_work > 0 ? cfg.max_work : 400LL * N * N;
  const long long phase1_budget = max_work / 4;

  Gauge gauge;
  if (has_fixed) {
    gauge.free.assign(3 * N, 1);
    for (int v = 0; v < N; ++v)
      if (cfg.fixed[v]) gauge.free[3 * v] = gauge.free[3 * v + 1] =
                            gauge.free[3 * v + 2] = 0;
  }

  std::span<double> x = as_flat(D.points);
  const std::size_t n = x.size();

  auto fg = [&](std::span<const double> xv, std::span<double> gv) {
    ++out.n_grad;
    const double e = model.gradient({}, xv, gv);
    gauge.project(gv);
    return e;
  };
  auto e_only = [&](std::span<const double> xv) {
    ++out.n_energy;
    return model.energy({}, xv);
  };
  auto gmax_L = [&](std::span<const double> gv) {
    const auto gc = as_coords(gv);
    double m = 0;
    for (int v = 0; v < N; ++v) {
      if (has_fixed && cfg.fixed[v]) continue;
      m = std::max(m, gc[v].norm());
    }
    return m * L;
  };

  std::vector<double> g(n), gp(n), xp(n), d(n), h_current, h_trial;

  // Diagnostic logging, parent line formats (~20 lines over the budget).
  const char* method_name = (cfg.method == DeltaMethod::CG)      ? "CG"
                          : (cfg.method == DeltaMethod::LBFGS)   ? "LBFGS"
                                                                 : "ST";
  const int log_interval =
      cfg.log ? std::max(1, (int)(max_work / (20LL * N))) : 0;
  auto edge_cv = [&]() {              // log-path only (~20 calls per run)
    std::vector<double> lens;
    lens.reserve(model.edges.size());
    for (const edge_t& e : model.edges)
      lens.push_back(coord3d::dist(D.points[e.first], D.points[e.second]));
    return cv_twopass(lens);
  };

  int phase = (model.k_flat > 0) ? 1 : 2;
  double f = fg(x, g);
  const double phase1_g0 = la::nrm2(g);

  if (cfg.log)
    fprintf(cfg.log, "  %s start: E=%.6f |g|=%.4e L=%.4f cv=%.4f ph=%d tol=%.2e\n",
            method_name, f, phase1_g0, L, edge_cv(), phase, cfg.grad_tol);

  // Phase 1 ends when its work quarter is exhausted or |g| dropped 100x;
  // checked every 50th iteration, exactly as the parent.
  auto phase_transition = [&](int iter) {
    if (phase != 1 || iter % 50 != 49) return false;
    bool advance = out.work() >= phase1_budget;
    if (!advance && phase1_g0 > 0 && la::nrm2(g) < 0.01 * phase1_g0)
      advance = true;
    if (advance) {
      model.k_flat = 0;
      phase = 2;
      return true;
    }
    return false;
  };

  const int stag_window = 50;
  int stag_count = 0;
  bool converged = false;

  // The shared battery over delta's conventions: gtol_inf gated on the
  // L-scaled per-vertex norm (passed as gmax), the angle/convexity
  // success as a domain predicate, the parent's 50-iteration stagnation
  // window.  The work budget stays at the loop top (parent order:
  // budget before phase transition), so it is disabled here.
  Criteria stop_batt;
  stop_batt.gtol_inf = cfg.grad_tol;
  stop_batt.stagnation_window = stag_window;
  if (cfg.angle_tol > 0)
    stop_batt.domain.push_back([&](std::span<const double>) {
      return D.max_angle_relerr() < cfg.angle_tol && D.count_concave() == 0;
    });
  auto stopped = [&]() {
    auto s = converged_or_stop(stop_batt, gmax_L(g), g, x, stag_count, 0.0);
    if (s) converged = (*s == Outcome::CONVERGED);
    return s.has_value();
  };
  auto track_stagnation = [&](double f_old, double f_new) {
    if (f_old - f_new > 1e-15 * std::max(1.0, std::fabs(f_old))) stag_count = 0;
    else ++stag_count;
  };

  // In-loop reflection (header comment): mirror the vertices folded
  // deeper than reflect_threshold L, re-evaluate at the mirrored point.
  // D.points aliases x, so the mirrored geometry is the iterate.
  // Returns the number of vertices moved (0: nothing changed).
  auto reflect_step = [&]() {
    if (!(cfg.reflect_threshold > 0)) return 0;
    const int n_r = D.reflect_concave(D.points, cfg.reflect_threshold * L, cfg.fixed);
    if (n_r > 0) {
      out.n_reflected += n_r;
      f = fg(x, g);
    }
    return n_r;
  };

  // Polak-Ribiere direction update d <- -g + beta d from (g, gp),
  // shared by the CG method and the post-cleanup CG polish.
  auto pr_direction = [&]() {
    const double gg_old = la::dot(gp, gp);
    double beta = 0;
    if (gg_old > 1e-30) {
      double num = 0;
      for (std::size_t i = 0; i < n; ++i) num += g[i] * (g[i] - gp[i]);
      beta = std::max(0.0, num / gg_old);
    }
    for (std::size_t i = 0; i < n; ++i) d[i] = -g[i] + beta * d[i];
    gauge.project(d);
  };

  // Armijo-quad line search + move, shared by CG and the L-BFGS route.
  // d must be a descent direction; returns the new energy.
  auto search_and_move = [&](double f0) {
    std::copy(x.begin(), x.end(), xp.begin());
    std::copy(g.begin(), g.end(), gp.begin());
    const double dg0 = la::dot(d, g);
    double f_new = f0;
    detail::armijo_quad_search(e_only, fg, xp, d, f0, dg0, /*t_init=*/1.0,
                               x, std::span<double>(g), f_new, /*c1=*/1e-4);
    return f_new;
  };

  // ==================== CG (Polak-Ribiere) ====================
  if (cfg.method == DeltaMethod::CG) {
    for (std::size_t i = 0; i < n; ++i) d[i] = -g[i];
    gauge.project(d);

    for (int iter = 0;; ++iter) {
      out.iters = iter + 1;
      if (out.work() >= max_work) break;
      // Reflect before the step (parent order: reflect, phase, stop);
      // the direction is kept -- the descent check below resets it if
      // the mirrored gradient made it an ascent direction.
      const bool reflected = iter > 0 && reflect_step() > 0;
      if (phase_transition(iter)) {
        f = fg(x, g);
        for (std::size_t i = 0; i < n; ++i) d[i] = -g[i];
        gauge.project(d);
      }
      if (stopped()) break;

      if (la::dot(d, g) > 0) {                    // ensure descent
        for (std::size_t i = 0; i < n; ++i) d[i] = -g[i];
        gauge.project(d);
      }

      const double f_old = f;
      f = search_and_move(f);
      track_stagnation(f_old, f);
      if (log_interval > 0 && iter % log_interval == 0) {
        const double dd = la::dot(d, d);
        const double alpha = dd > 0 ? (la::dot(x, d) - la::dot(xp, d)) / dd : 0;
        fprintf(cfg.log, "  CG %4d: E=%.6f |g|=%.4e gmax*L=%.4e a=%.3e cv=%.4f ang=%.2e ph=%d%s\n",
                iter, f, la::nrm2(g), gmax_L(g), alpha, edge_cv(),
                D.max_angle_relerr(), phase, reflected ? " R" : "");
      }
      pr_direction();
    }
  }

  // ==================== L-BFGS ====================
  else if (cfg.method == DeltaMethod::LBFGS) {
    LBFGS policy{10};
    LBFGS::State st(n, policy);
    std::vector<double> s(n), y(n);

    for (int iter = 0;; ++iter) {
      out.iters = iter + 1;
      if (out.work() >= max_work) break;
      // Reflect before the step; the geometry changed outside the
      // quasi-Newton model, so the history is dropped (as the parent).
      const bool reflected = iter > 0 && reflect_step() > 0;
      if (reflected) policy.reset(st);
      if (phase_transition(iter)) {
        f = fg(x, g);
        policy.reset(st);
      }
      if (stopped()) break;

      policy.direction({}, st, g, d);
      gauge.project(d);
      if (la::dot(d, g) > 0) {       // parent predicate: reset only on ascent
        for (std::size_t i = 0; i < n; ++i) d[i] = -g[i];
        gauge.project(d);
        policy.reset(st);
        out.flags |= Result::LBFGS_RESET;
      }

      const double f_old = f;
      f = search_and_move(f);
      track_stagnation(f_old, f);

      for (std::size_t i = 0; i < n; ++i) {
        s[i] = x[i] - xp[i];
        y[i] = g[i] - gp[i];
      }
      policy.accepted({}, st, s, y, g);
      if (log_interval > 0 && iter % log_interval == 0) {
        const double dd = la::dot(d, d);
        const double alpha = dd > 0 ? la::dot(s, d) / dd : 0;
        std::vector<double> h_log;
        D.aet_h_values(D.points, h_log, cfg.fixed);
        double hmin = 1.0;
        for (int v = 0; v < N; ++v)
          if (!(has_fixed && cfg.fixed[v])) hmin = std::min(hmin, h_log[v]);
        fprintf(cfg.log, "  LB %4d: E=%.6f |g|=%.4e gmax*L=%.4e a=%.3e cv=%.4f ang=%.2e hmin=%+.4f ph=%d h=%d%s\n",
                iter, f, la::nrm2(g), gmax_L(g), alpha, edge_cv(),
                D.max_angle_relerr(), hmin, phase, st.stored, reflected ? " R" : "");
      }
    }
  }

  // ==================== Steihaug-Toint ====================
  else {
    TrustRadius tr{0.5 * L, L, 1e-14 * L};
    std::vector<double> z(n), Hz(n), x_trial(n);
    if (cfg.convex_constraint) D.aet_h_values(D.points, h_current, cfg.fixed);

    for (int iter = 0;; ++iter) {
      out.iters = iter + 1;
      if (out.work() >= max_work) break;
      // Reflect before the step; the constrained variant's h reference
      // follows the mirrored geometry.  The radius is kept (parent).
      const bool reflected = iter > 0 && reflect_step() > 0;
      if (reflected && cfg.convex_constraint) D.aet_h_values(D.points, h_current, cfg.fixed);
      if (phase_transition(iter)) {
        f = fg(x, g);
        tr.Delta = 0.5 * L;                        // reset radius on phase change
        if (cfg.convex_constraint) D.aet_h_values(D.points, h_current, cfg.fixed);
      }
      if (stopped()) break;

      const int hv_before = out.n_hv;
      SteihaugCG::solve_subproblem({}, model, gauge, x, g, tr.Delta, z,
                                   out.flags, out.n_hv);

      ++out.n_hv;
      model.hvp({}, x, z, Hz);
      const double pred = -(la::dot(g, z) + 0.5 * la::dot(z, Hz));

      for (std::size_t i = 0; i < n; ++i) x_trial[i] = x[i] + z[i];
      const double f_trial = e_only(x_trial);

      const double rho = TrustRadius::rho(f - f_trial, pred);
      bool accepted = rho > 0.1;
      if (accepted && cfg.convex_constraint) {
        D.aet_h_values(as_coords(std::span<const double>(x_trial)), h_trial,
                       cfg.fixed);
        for (int v = 0; v < N && accepted; ++v)
          if (h_current[v] > 0 && h_trial[v] < 0) {
            accepted = false;
            out.flags |= Result::CVX_REJECTED;
          }
      }

      if (accepted) {
        const double f_old = f;
        std::copy(x_trial.begin(), x_trial.end(), x.begin());
        f = fg(x, g);
        if (cfg.convex_constraint) D.aet_h_values(D.points, h_current, cfg.fixed);
        tr.accept(rho, la::nrm2(z));
        track_stagnation(f_old, f);
      } else {
        tr.reject();
        ++stag_count;
        out.flags |= Result::STEP_REJECTED;
      }
      if (log_interval > 0 && iter % log_interval == 0)
        fprintf(cfg.log, "  ST %4d: E=%.6f |g|=%.4e gmax*L=%.4e |z|=%.2e D=%.2e rho=%.2f ang=%.2e in=%d ph=%d %s%s\n",
                iter, f, la::nrm2(g), gmax_L(g), la::nrm2(z), tr.Delta, rho,
                D.max_angle_relerr(), out.n_hv - hv_before - 1, phase,
                accepted ? "acc" : "REJ", reflected ? " R" : "");
    }
  }

  out.gmax_L = gmax_L(g);
  const bool stagnated = !converged && stag_count >= stag_window;
  if (cfg.log)
    fprintf(cfg.log, "  %s done: %d iters, E=%.6f gmax*L=%.4e cv=%.4f refl=%d %s\n",
            method_name, out.iters, f, out.gmax_L, edge_cv(), out.n_reflected,
            converged ? "CONVERGED" : stagnated ? "STAGNATED" : "budget");

  // Post-optimization strict-convexity cleanup, exactly as the parent:
  // hull-project, then a brief CG polish if anything moved, then a final
  // projection.
  if (!cfg.skip_post_reflect) {
    const int projected = D.project_onto_convex_hull(D.points);
    if (projected > 0) {
      f = fg(x, g);
      for (std::size_t i = 0; i < n; ++i) d[i] = -g[i];
      gauge.project(d);
      for (int iter = 0; iter < 50; ++iter) {
        if (gmax_L(g) < cfg.grad_tol) break;
        f = search_and_move(f);
        pr_direction();
      }
      D.project_onto_convex_hull(D.points);
      if (cfg.log)
        fprintf(cfg.log, "  Post-project polish: projected=%d ang=%.4e\n",
                projected, D.max_angle_relerr());
    }
  }

  out.E = f;
  out.angle_relerr = D.max_angle_relerr();
  out.n_concave = D.count_concave();
  out.outcome = converged    ? Outcome::CONVERGED
                : stagnated  ? Outcome::STAGNATED
                             : Outcome::BUDGET_EXHAUSTED;
  return out;
}

}  // namespace optim
