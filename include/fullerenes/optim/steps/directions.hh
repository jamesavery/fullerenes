#pragma once

// =====================================================================
// Direction-generator step policies for the LineSearch paradigm.
//
// The line-search step contract (what LineSearch<Step> requires):
//
//   using paradigm = line_search_tag;
//   struct State { State(std::size_t n, const Step&); };
//       -- per-solve mutable history, sized for n DOF
//   void reset(State&) const;
//       -- the paradigm replaced the proposed direction by -g (restart)
//   void invalidate(State&) const;
//       -- the iterate moved outside the line search (post_accept hook):
//          drop whatever history no longer describes the new point
//   bool has_history(const State&) const;
//       -- false iff a failed search from here cannot be retried with
//          less history (the paradigm then reports STEP_FAILED)
//   void direction(SeqCtx, State&, span<const double> g,
//                  span<double> d) const;
//       -- propose the descent direction d from g and the history.
//          The paradigm verifies <d,g> < 0 and restarts otherwise.
//   void accepted(SeqCtx, State&, span<const double> s,
//                 span<const double> y, span<const double> g) const;
//       -- observe an accepted step: s = x1-x0, y = g1-g0, g = g1.
//   double first_trial(const State&, double gmax) const;
//       -- initial line-search step length t0 for this iteration.
//
// This file: SteepestDescent (the contract's simplest instance, used by
// the unit tests) and ConjugateGradient (Polak-Ribiere+, Deltahedron's
// CG method).  L-BFGS lives in steps/lbfgs.hh.
// =====================================================================

#include "fullerenes/optim/core.hh"

#include <algorithm>
#include <span>
#include <vector>

namespace optim {

struct SteepestDescent {
  using paradigm = line_search_tag;

  struct State {
    State(std::size_t, const SteepestDescent&) {}
  };

  void reset(State&) const {}
  void invalidate(State&) const {}
  bool has_history(const State&) const { return false; }

  void direction(SeqCtx, State&, std::span<const double> g,
                 std::span<double> d) const {
    for (std::size_t i = 0; i < g.size(); ++i) d[i] = -g[i];
  }

  void accepted(SeqCtx, State&, std::span<const double>,
                std::span<const double>, std::span<const double>) const {}

  double first_trial(const State&, double gmax) const {
    return std::min(1.0, 1.0 / std::max(gmax, 1e-12));
  }
};

// Polak-Ribiere+ conjugate gradient, the Deltahedron parent's formula
// in its floating-point order:
//   beta = max(0, sum_i g_i (g_i - gp_i) / <gp, gp>)   (0 if <gp,gp> <= 1e-30)
//   d    = -g + beta d_prev
// gp is a COPY of the gradient the previous direction was computed
// from (taken in direction(), not reconstructed as g - y, which differs
// in rounding).  On a paradigm restart (reset) the direction actually
// used was -gp, which the next direction() takes as d_prev -- exactly
// the parent, whose restart overwrote d in place and kept the PR
// momentum of the reset direction.  First trial is always 1 (the
// parent's Armijo-quad line search starts at alpha = 1).
struct ConjugateGradient {
  using paradigm = line_search_tag;

  struct State {
    std::vector<double> gp, d_prev;
    bool has_prev = false;   // gp/d_prev valid
    bool restarted = false;  // the paradigm replaced the last direction by -gp
    State(std::size_t n, const ConjugateGradient&) : gp(n), d_prev(n) {}
  };

  void reset(State& st) const { st.restarted = true; }
  // After an external move the parent kept its PR momentum (gp, d_prev)
  // and built the next direction from the new gradient; so does this.
  void invalidate(State&) const {}

  bool has_history(const State& st) const { return st.has_prev && !st.restarted; }

  void direction(SeqCtx, State& st, std::span<const double> g,
                 std::span<double> d) const {
    const std::size_t n = g.size();
    if (!st.has_prev) {
      for (std::size_t i = 0; i < n; ++i) d[i] = -g[i];
    } else {
      if (st.restarted) {                       // the step taken was -gp
        for (std::size_t i = 0; i < n; ++i) st.d_prev[i] = -st.gp[i];
        st.restarted = false;
      }
      const double gg_old = la::dot(st.gp, st.gp);
      double beta = 0;
      if (gg_old > 1e-30) {
        double num = 0;
        for (std::size_t i = 0; i < n; ++i) num += g[i] * (g[i] - st.gp[i]);
        beta = std::max(0.0, num / gg_old);
      }
      for (std::size_t i = 0; i < n; ++i) d[i] = -g[i] + beta * st.d_prev[i];
    }
    std::copy(g.begin(), g.end(), st.gp.begin());
    std::copy(d.begin(), d.end(), st.d_prev.begin());
    st.has_prev = true;
  }

  void accepted(SeqCtx, State&, std::span<const double>,
                std::span<const double>, std::span<const double>) const {}

  double first_trial(const State&, double) const { return 1.0; }
};

}  // namespace optim
