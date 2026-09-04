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
//       -- drop history (steepest-descent restart)
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
// the unit tests).  PR-CG joins it with the Deltahedron migration;
// L-BFGS lives in steps/lbfgs.hh.
// =====================================================================

#include "fullerenes/optim/core.hh"

#include <algorithm>
#include <span>

namespace optim {

struct SteepestDescent {
  using paradigm = line_search_tag;

  struct State {
    State(std::size_t, const SteepestDescent&) {}
  };

  void reset(State&) const {}

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

}  // namespace optim
