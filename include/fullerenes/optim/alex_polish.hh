#pragma once

// =====================================================================
// optim::alexandrov_polish -- the Alexandrov stage-4 trust-region
// Newton polish as a framework COMPOSITION: TrustRegion<LevenbergBisect>
// over the cell-resolved AlexKappaWD model, with two hooks --
//
//   clip         AlexandrovSolver::feasible_fraction (the F(T)-feasible
//                step rule; the paradigm scales the step by it and caps
//                the radius at the length actually taken)
//   post_accept  AlexKappaWD::commit (adopt the accepted trial's
//                weighted-Delaunay cell -- the parent's "adopt the copy
//                atomically on acceptance")
//
// and the GN-normal-equations subproblem + rho-ratio radius policy
// (floor 1e-14) from the paradigm.  This is a faithful re-expression of
// the library's internal Newton::polish (delaunay_alexandrov.cc, post
// the 2026-07-25 stall cure); the transcription-era functor that
// hand-rolled the same loop (the signed-off "monolithic for now"
// interim) is gone -- the hooks it was waiting for are the ones of
// DESIGN.md 10.2.
//
// kappa is only piecewise-smooth across flip boundaries, so the
// parent's two flip invariants live in the model: ENTRY (weighted-
// Delaunay restored before the first kappa/J evaluation -- the
// AlexKappaWD constructor) and TRIAL (every evaluation on the query
// point's own flipped copy; commit on acceptance).  Installed as
// AlexandrovSolver::polish_override, it runs inside the full production
// pipeline on the exact post-extrapolation state -- the head-to-head
// arrangement validate_alex_reexpr gates on (and solve() re-checks
// kappa < 1e-10 after the override, so a wrong verdict cannot leak).
// On this path AlexandrovSolver::stats_flips excludes the polish's own
// flips (documented at the seam).
// =====================================================================

#include "fullerenes/optim/models/alex_kappa.hh"
#include "fullerenes/optim/steps/levenberg.hh"
#include "fullerenes/optim/trustregion.hh"

#include "fullerenes/delaunay_alexandrov.hh"

#include <span>
#include <vector>

namespace optim {

// @anchor  optim-alexandrov-polish
// @pre     r in F(T) (every incident pyramid closes); r.size() == T.nv
// @post    true iff max|kappa(T, r)| < tol; T is weighted-Delaunay for
//          the final r (entry flip + per-trial flipped copies)
// @variant max_iter outer iterations; 20 consecutive rejections abort
inline bool alexandrov_polish(DelaunayTriangulation& T, std::vector<double>& r,
                              double tol = 1e-10, int max_iter = 50) {
  const std::size_t n = r.size();
  const double r_avg = la::dot(r, la::V(n, 1.0)) / n;

  AlexKappaWD model(T, r);                        // ENTRY invariant
  Criteria stop;
  stop.rtol_inf  = tol;
  stop.max_iters = max_iter;
  Problem<AlexKappaWD> prob{model, {}, nullptr, stop};

  TrustRegion<LevenbergBisect> method;
  method.Delta0      = 0.5 * r_avg;
  method.Delta_max   = 2.0 * r_avg;
  method.delta_floor = 1e-14;

  Hooks hooks;
  hooks.clip = [&](std::span<const double> x, std::span<const double> d) {
    return AlexandrovSolver::feasible_fraction(T, la::V(x.begin(), x.end()),
                                              la::V(d.begin(), d.end()));
  };
  hooks.post_accept = [&](SeqCtx, std::span<double> x) {
    model.commit(x);                              // the model's cell changed, x did not
    return 1;
  };

  const Result res = method.solve({}, prob, std::span<double>(r), hooks);
  return res.outcome == Outcome::CONVERGED;
}

}  // namespace optim
