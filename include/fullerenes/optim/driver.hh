#pragma once

// =====================================================================
// optim::minimize -- the framework driver (DESIGN.md 2, the spine):
//
//     continue( schedule ):                     outer continuation
//         paradigm( step, preconditioner )      loop-owner + step policy
//             .solve( problem, x, hooks )
//
// A Continuation is the model's parameter schedule with memory
// (DESIGN.md 10.1.2): three operations, no report object --
//
//   retarget(prob)            install the current level into the
//                             Problem: model targets (blend t, eps, mu,
//                             k_flat, ...) AND the level's stop criteria
//   advance(prob, x, result)  observe the solved level (the policy
//                             measures what it needs from prob.E and x
//                             itself; reads result.iters / n_accepted /
//                             outcome), move the schedule; false = done
//   terminal()                are the hard-invariant parameters at
//                             their terminal value?
//
// The driver never reports CONVERGED for a schedule that stopped short:
// !terminal() at exit is Outcome::SATURATED -- the t=1 invariant of the
// curvature-flow embed as a checked contract.  The identity
// continuation (one level, no retarget) degenerates the driver to a
// single solve.  A mesh swap (cascadic subdivision) is NOT a retarget --
// it changes the DOF count -- and stays a driver-of-drivers outside.
//
// Prep is the ordered start-geometry pipeline (points -> points
// transforms) applied once before the first solve.  Hooks are the
// paradigm extension points (core.hh), forwarded to every solve.
//
// Workspace is accepted for signature stability (external storage is
// the parallel-primitives convention, needed by the step-6 batch
// lowering); the CPU paradigms currently allocate internally.
// =====================================================================

#include "fullerenes/optim/core.hh"

#include <concepts>
#include <functional>
#include <span>
#include <vector>

namespace optim {

using Prep = std::vector<std::function<void(std::span<double>)>>;

struct Workspace {};   // external-storage slot (step 6); unused on CPU

template <class C, class Model>
concept Continuation =
    requires(C c, Problem<Model>& p, std::span<const double> x, const Result& r) {
      c.retarget(p);
      { c.advance(p, x, r) } -> std::convertible_to<bool>;
      { c.terminal() } -> std::convertible_to<bool>;
    };

// One level, no retarget, always terminal.
struct IdentityContinuation {
  template <class M> void retarget(Problem<M>&) const {}
  template <class M>
  bool advance(Problem<M>&, std::span<const double>, const Result&) const { return false; }
  bool terminal() const { return true; }
};

// Sum the per-level solves into one Result: counters add, flags OR,
// value/outcome/iteration-of-last-level are the last level's.
inline void accumulate(Result& acc, const Result& r) {
  acc.outcome    = r.outcome;
  acc.f          = r.f;
  acc.gmax       = r.gmax;
  acc.iters     += r.iters;
  acc.n_accepted += r.n_accepted;
  acc.n_energy  += r.n_energy;
  acc.n_grad    += r.n_grad;
  acc.n_hv      += r.n_hv;
  acc.flags     |= r.flags;
}

// @anchor optim-minimize
// @pre    x.size() is the model's DOF count
// @post   x holds the final iterate of the last continuation level;
//         result accumulates the levels' counters; result.outcome is
//         the last level's, or SATURATED if the schedule stopped with
//         !cont.terminal()
template <class Ctx, SmoothModel Model, class Paradigm,
          Continuation<Model> Cont = IdentityContinuation>
Result minimize(const Ctx& ctx, Problem<Model> prob, const Paradigm& method,
                std::span<double> x, Workspace = {}, Cont cont = {},
                const Prep& prep = {}, const Hooks& hooks = {}) {
  for (const auto& p : prep) p(x);

  Result acc;
  Result r;
  do {
    cont.retarget(prob);
    r = method.solve(ctx, prob, x, hooks);
    accumulate(acc, r);
  } while (cont.advance(prob, x, r));
  if (!cont.terminal()) acc.outcome = Outcome::SATURATED;
  return acc;
}

}  // namespace optim
