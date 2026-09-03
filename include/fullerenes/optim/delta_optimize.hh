#pragma once

// =====================================================================
// optim::delta_optimize -- Deltahedron::optimize() as a framework
// COMPOSITION over the AET model:
//
//   Continuation{k_flat two-phase schedule}
//     ▸ LineSearch<LBFGS, ARMIJO QUAD_INTERP>     (method LBFGS)
//     ▸ LineSearch<ConjugateGradient, ARMIJO>     (method CG)
//     ▸ TrustRegion<SteihaugCG>                   (method STEIHAUG)
//   with Hooks { post_accept = in-loop reflection,
//                veto        = the convex-constraint gate (STEIHAUG),
//                observe     = the parent-format log lines }
//   + the post-optimization hull-projection cleanup epilogue
//     (LineSearch<ConjugateGradient> polish, 50 iterations).
//
// No hand-rolled loop remains: the gmax_L convention is
// Criteria::gnorm = BLOCK_INF{3, L}, the 50-iteration stagnation window
// and the unified 1:2:7 work budget are Criteria, the phase-1 exit
// (work quarter spent, or ||g||_2 dropped 100x) is the level's own
// Criteria (max_work, gtol_2_rel) and the schedule's advance().
//
// Known, accepted deviations from the pre-composition parent (the gate
// is quality parity -- bench_epopt A/B on the production pipeline,
// which runs single-phase LBFGS / STEIHAUG; deltahedron-test and the
// diag_c50_prep Tutte-start sweep on the standalone path):
//   * the phase-1 exit is tested every iteration, not every 50th;
//   * the stagnation counter restarts at the phase transition;
//   * CG after an in-loop reflection builds its next direction from the
//     reflected gradient (the parent reused the pre-reflection one);
//   * the L-BFGS two-loop is the minimize.hh transcription (gamma as
//     1/(rho y'y), curvature skip 1e-10 sqrt(|s|^2 |y|^2)) where the
//     parent Deltahedron carried a deque variant (gamma as ys/yy, skip
//     1e-10 s's); both line searches start at alpha = 1 as the parent
//     does (LBFGS::unit_first_trial);
//   * Result::iters counts completed iterations (the parent counted
//     loop entries: +1).
//
// In-loop reflection (cfg.reflect_threshold): after every accepted
// step, vertices with h < -reflect_threshold L are mirrored through
// their neighbour-centroid plane (the parent's pre-2026-03-08
// "periodic reflection", the DESIGN.md 3.3 between-steps projection as
// the post_accept hook).  It is the fold-escape operator for cold
// starts: a descent from a fragile (Tutte-sphere) start folds a
// pentagon cap on the way down and, left alone, converges into a
// CONVEX wrong-basin minimum that no post-convergence reflect/hull
// pass can touch (C50 ang_max 94.13 vs the 75-degree guard;
// claude-projects/optimize/validation/C50-ANGMAX-INVESTIGATION.md).
// Reflecting the fold as it forms keeps the descent in the equilateral
// basin.  0 disables.
//
// Convexity restoration BETWEEN optimize calls (the pipeline's
// reflect/hull loops) stays with the caller exactly as before.
// =====================================================================

#include "fullerenes/optim/core.hh"
#include "fullerenes/optim/driver.hh"
#include "fullerenes/optim/linesearch.hh"
#include "fullerenes/optim/trustregion.hh"
#include "fullerenes/optim/steps/directions.hh"
#include "fullerenes/optim/steps/lbfgs.hh"
#include "fullerenes/optim/steps/steihaug.hh"
#include "fullerenes/optim/models/aet.hh"

#include "fullerenes/graphview.hh"
#include "fullerenes/stats.hh"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <functional>
#include <limits>
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

// The two-phase E_flat schedule as a Continuation (driver.hh):
//   level 1  k_flat on; exits when its work quarter is spent
//            (max_work) or ||g||_2 dropped 100x (gtol_2_rel) -- or
//            when the whole solve is done (gmax_L / angle criterion /
//            stagnation), in which case there is no level 2;
//   level 2  k_flat off, the remaining budget.
// k_flat0 == 0 is the identity (one level).  Always terminal.
struct KFlatSchedule {
  double k_flat0;
  double phase1_budget;          // work ceiling of level 1
  double total_budget;
  double grad_tol;               // gmax_L tolerance (to tell which CONVERGED fired)
  std::function<bool()> domain_done;   // the angle criterion at the current geometry
  int    phase = 1;
  double work_used = 0;

  template <class M>
  void retarget(Problem<M>& p) const {
    if (phase == 1 && k_flat0 > 0) {
      p.E.k_flat = k_flat0;
      p.stop.gtol_2_rel = 0.01;
      p.stop.max_work = phase1_budget;
    } else {
      p.E.k_flat = 0;
      p.stop.gtol_2_rel = 0;
      p.stop.max_work = std::max(0.0, total_budget - work_used);
    }
  }
  template <class M>
  bool advance(Problem<M>&, std::span<const double>, const Result& r) {
    work_used += r.work<M>();
    if (phase != 1 || k_flat0 <= 0) return false;
    phase = 2;
    // Phase 1 ran out of its quarter, or its gradient dropped 100x
    // (CONVERGED, but neither the gmax_L tolerance nor the angle
    // criterion holds): continue with E_flat off.  Anything else ended
    // the whole solve, as the parent did.
    if (r.outcome == Outcome::BUDGET_EXHAUSTED) return true;
    if (r.outcome == Outcome::CONVERGED && r.gmax > grad_tol && !domain_done())
      return true;
    return false;
  }
  bool terminal() const { return true; }
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

  const long long max_work = cfg.max_work > 0 ? cfg.max_work : 400LL * N * N;

  Gauge gauge;
  if (has_fixed) {
    gauge.free.assign(3 * N, 1);
    for (int v = 0; v < N; ++v)
      if (cfg.fixed[v]) gauge.free[3 * v] = gauge.free[3 * v + 1] =
                            gauge.free[3 * v + 2] = 0;
  }

  std::span<double> x = as_flat(D.points);

  // --- Criteria: delta's conventions --------------------------------
  Criteria stop;
  stop.gnorm = GradNorm::BLOCK_INF;
  stop.gnorm_block = 3;
  stop.gnorm_scale = L;
  stop.gtol_inf = cfg.grad_tol;
  stop.stagnation_window = 50;
  stop.max_iters = std::numeric_limits<int>::max();   // the work budget is the guard
  auto domain_done = [&]() {
    return D.max_angle_relerr() < cfg.angle_tol && D.count_concave() == 0;
  };
  if (cfg.angle_tol > 0)
    stop.domain.push_back([&](std::span<const double>) { return domain_done(); });
  Problem<AET> prob{model, gauge, nullptr, stop};

  // --- Hooks ---------------------------------------------------------
  Hooks hooks;
  if (cfg.reflect_threshold > 0)
    hooks.post_accept = [&](SeqCtx, std::span<double>) {
      const int n_r = D.reflect_concave(D.points, cfg.reflect_threshold * L, cfg.fixed);
      out.n_reflected += n_r;
      return n_r;
    };
  std::vector<double> h_cur, h_trial;
  if (cfg.method == DeltaMethod::STEIHAUG && cfg.convex_constraint)
    hooks.veto = [&](std::span<const double> x_trial) {
      D.aet_h_values(D.points, h_cur, cfg.fixed);
      D.aet_h_values(as_coords(x_trial), h_trial, cfg.fixed);
      for (int v = 0; v < N; ++v)
        if (h_cur[v] > 0 && h_trial[v] < 0) return false;
      return true;
    };

  // Diagnostic logging, parent line formats (~20 lines over the budget).
  const char* method_name = (cfg.method == DeltaMethod::CG)      ? "CG"
                          : (cfg.method == DeltaMethod::LBFGS)   ? "LB"
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
  if (log_interval > 0)
    hooks.observe = [&](const IterationView& v) {
      if (v.iter % log_interval != 0) return;
      const int phase = model.k_flat > 0 ? 1 : 2;
      if (cfg.method == DeltaMethod::STEIHAUG)
        fprintf(cfg.log, "  ST %4d: E=%.6f |g|=%.4e gmax*L=%.4e |z|=%.2e D=%.2e rho=%.2f ang=%.2e in=%d ph=%d %s%s\n",
                v.iter, v.f, la::nrm2(v.g), v.gmax, v.step_len, v.radius, v.rho,
                D.max_angle_relerr(), v.n_inner, phase,
                v.accepted ? "acc" : "REJ", v.n_post > 0 ? " R" : "");
      else
        fprintf(cfg.log, "  %s %4d: E=%.6f |g|=%.4e gmax*L=%.4e a=%.3e cv=%.4f ang=%.2e ph=%d%s\n",
                method_name, v.iter, v.f, la::nrm2(v.g), v.gmax, v.step_len, edge_cv(),
                D.max_angle_relerr(), phase, v.n_post > 0 ? " R" : "");
    };

  if (cfg.log)
    fprintf(cfg.log, "  %s start: E=%.6f L=%.4f cv=%.4f kflat=%.1f tol=%.2e\n",
            method_name, model.energy({}, x), L, edge_cv(), model.k_flat, cfg.grad_tol);

  // --- The composition -----------------------------------------------
  KFlatSchedule schedule{cfg.k_flat, (double)(max_work / 4), (double)max_work,
                         cfg.grad_tol, domain_done};
  Result r;
  if (cfg.method == DeltaMethod::LBFGS) {
    LineSearch<LBFGS> method;
    method.step = LBFGS{10, /*unit_first_trial=*/true};
    method.cond = LineSearchCond::ARMIJO;
    r = minimize(SeqCtx{}, prob, method, x, {}, schedule, {}, hooks);
  } else if (cfg.method == DeltaMethod::CG) {
    LineSearch<ConjugateGradient> method;
    method.cond = LineSearchCond::ARMIJO;
    r = minimize(SeqCtx{}, prob, method, x, {}, schedule, {}, hooks);
  } else {
    TrustRegion<SteihaugCG> method;
    method.Delta0 = 0.5 * L;
    method.Delta_max = L;
    method.delta_floor = 1e-14 * L;
    r = minimize(SeqCtx{}, prob, method, x, {}, schedule, {}, hooks);
  }
  out.iters = r.iters;
  out.n_energy = r.n_energy;
  out.n_grad = r.n_grad;
  out.n_hv = r.n_hv;
  out.flags = r.flags;
  out.gmax_L = r.gmax;
  double f = r.f;

  if (cfg.log)
    fprintf(cfg.log, "  %s done: %d iters, E=%.6f gmax*L=%.4e cv=%.4f refl=%d %s\n",
            method_name, out.iters, f, out.gmax_L, edge_cv(), out.n_reflected,
            outcome_name(r.outcome));

  // --- Post-optimization strict-convexity cleanup, as the parent:
  // hull-project, then a brief CG polish if anything moved, then a
  // final projection.
  if (!cfg.skip_post_reflect) {
    const int projected = D.project_onto_convex_hull(D.points);
    if (projected > 0) {
      Criteria polish_stop = stop;
      polish_stop.domain.clear();
      polish_stop.stagnation_window = 0;
      polish_stop.max_work = 0;
      polish_stop.max_iters = 50;
      Problem<AET> polish_prob{model, gauge, nullptr, polish_stop};
      LineSearch<ConjugateGradient> polish;
      polish.cond = LineSearchCond::ARMIJO;
      const Result pr = polish.solve(SeqCtx{}, polish_prob, x);
      out.n_energy += pr.n_energy;
      out.n_grad += pr.n_grad;
      f = pr.f;
      D.project_onto_convex_hull(D.points);
      if (cfg.log)
        fprintf(cfg.log, "  Post-project polish: projected=%d ang=%.4e\n",
                projected, D.max_angle_relerr());
    }
  }

  out.E = f;
  out.angle_relerr = D.max_angle_relerr();
  out.n_concave = D.count_concave();
  out.outcome = r.outcome;
  return out;
}

}  // namespace optim
