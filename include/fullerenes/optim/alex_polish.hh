#pragma once

// =====================================================================
// optim::alexandrov_polish -- the Alexandrov stage-4 trust-region
// Newton polish re-expressed on the framework's components: the loop is
// a faithful transcription of the library's internal Newton::polish
// (delaunay_alexandrov.cc, post the 2026-07-25 stall cure), with
//
//   * the GN-normal-equations subproblem -> LevenbergBisect::solve_subproblem
//   * the rho-ratio radius policy       -> TrustRadius (floor 1e-14)
//   * the F(T)-feasibility clip         -> AlexandrovSolver::feasible_step
//   * the weighted-Delaunay flips       -> AlexandrovSolver::flip_to_weighted_delaunay
//
// the last two through the library's public statics (single source of
// truth; the seam exposed them for exactly this).
//
// kappa is only piecewise-smooth across flip boundaries, so the parent's
// two flip invariants are kept verbatim:
//   - ENTRY: restore weighted-Delaunay before the first kappa/J
//     evaluation (the endgame extrapolation moves r without flipping);
//   - TRIAL: evaluate each trial on its own flipped copy (T_trial,
//     r_trial), adopt the copy atomically on acceptance, discard on
//     rejection -- acceptance compares true energies across cells.
//
// Per the signed-off "monolithic for now" decision the functor owns the
// loop -- the clip feeds back into the radius (Delta capped at the
// clipped step length), which the generic TrustRegion paradigm has no
// hook for.  Installed as AlexandrovSolver::polish_override, it runs
// inside the full production pipeline on the exact post-extrapolation
// state -- the head-to-head arrangement validate_alex_reexpr gates on
// (and solve() re-checks kappa < 1e-10 after the override, so a wrong
// verdict cannot leak).  On this path AlexandrovSolver::stats_flips
// excludes the polish's own flips (documented at the seam).
// =====================================================================

#include "fullerenes/optim/steps/levenberg.hh"
#include "fullerenes/optim/trustregion.hh"

#include "fullerenes/delaunay_alexandrov.hh"

#include <algorithm>
#include <utility>
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
  TrustRadius tr{0.5 * r_avg, 2.0 * r_avg, 1e-14};
  int rejects = 0;

  AlexandrovSolver::flip_to_weighted_delaunay(T, r);   // ENTRY invariant

  for (int iter = 0; iter < max_iter; ++iter) {
    const auto kappa = AlexandrovSolver::kappa(T, r);
    if (la::max_abs(kappa) < tol) return true;
    if (rejects > 20) break;

    const double E = la::energy(kappa);
    const Mat J = AlexandrovSolver::jacobian(T, r);

    const la::V delta_raw = LevenbergBisect::solve_subproblem(J, kappa, tr.Delta);
    bool clipped;
    const auto delta = AlexandrovSolver::feasible_step(T, r, delta_raw, &clipped);
    const double pred = LevenbergBisect::predicted_reduction(J, kappa, delta);

    la::V r_trial(n);
    for (std::size_t i = 0; i < n; ++i) r_trial[i] = r[i] + delta[i];
    // TRIAL invariant: kappa(r_trial) on r_trial's own weighted-Delaunay
    // cell; adopt the copy atomically on acceptance.
    DelaunayTriangulation T_trial = T;
    AlexandrovSolver::flip_to_weighted_delaunay(T_trial, r_trial);
    const double E_trial = la::energy(AlexandrovSolver::kappa(T_trial, r_trial));

    const bool ok = tr.update(E - E_trial, pred, la::norm(delta));
    // A clipped step caps the radius at the length actually taken, so
    // the next subproblem does not repropose the infeasible direction.
    if (clipped) tr.Delta = std::min(tr.Delta, la::norm(delta));

    if (ok) {
      r = std::move(r_trial);
      T = std::move(T_trial);
      rejects = 0;
    } else
      ++rejects;
  }
  return false;
}

}  // namespace optim
