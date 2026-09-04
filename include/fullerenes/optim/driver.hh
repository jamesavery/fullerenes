#pragma once

// =====================================================================
// optim::minimize -- the framework driver (DESIGN.md 2, the spine):
//
//     continue( schedule ):                     outer homotopy (optional)
//         paradigm( step, preconditioner )      loop-owner + step policy
//             .solve( problem, x )
//
// Continuation is identity (a single pass) for every unconstrained
// path; the mu/eps/t schedules of the constrained embed arrive with
// migration step 5.  Prep is the ordered start-geometry pipeline
// (points -> points transforms) applied once before the first solve.
// Projection is the between-steps feasibility hook (DESIGN.md 3.3),
// forwarded to the paradigm.
//
// Workspace is accepted for signature stability (external storage is
// the parallel-primitives convention, needed by the step-6 batch
// lowering); the CPU paradigms currently allocate internally.
// =====================================================================

#include "fullerenes/optim/core.hh"

#include <functional>
#include <span>
#include <vector>

namespace optim {

using Prep = std::vector<std::function<void(std::span<double>)>>;

struct Workspace {};   // external-storage slot (step 6); unused on CPU

// Outer homotopy driver: `levels` inner solves, `retarget` called with
// the level index before each.  The default is the identity
// continuation -- one level, no retarget.
struct Continuation {
  int levels = 1;
  std::function<void(SeqCtx, int level)> retarget;
};

// @anchor optim-minimize
// @pre    x.size() is the model's DOF count; cont.levels >= 1
// @post   x holds the final iterate of the last continuation level;
//         result is that level's solve outcome
template <class Ctx, SmoothModel Model, class Paradigm>
Result minimize(const Ctx& ctx, Problem<Model> prob, const Paradigm& method,
                std::span<double> x, Workspace = {},
                const Continuation& cont = {}, const Prep& prep = {},
                const Projection& project = {}) {
  for (const auto& p : prep) p(x);

  Result out;
  for (int level = 0; level < cont.levels; ++level) {
    if (cont.retarget) cont.retarget(ctx, level);
    out = method.solve(ctx, prob, x, project);
  }
  return out;
}

}  // namespace optim
