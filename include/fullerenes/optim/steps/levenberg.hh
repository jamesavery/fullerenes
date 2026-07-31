#pragma once

// =====================================================================
// optim::LevenbergBisect -- the Gauss-Newton Levenberg trust-region
// subproblem solver over a SQUARE SYMMETRIC system matrix, a verbatim
// transcription of delaunay_alexandrov's TrustRegion::solve /
// predicted_reduction (post the 2026-07-25 Newton-polish cure):
//
//   try the pure Newton root step delta = -J^{-1} R first (lambda = 0);
//   outside the radius (or with singular J), fall back to
//   Levenberg-Marquardt on the NORMAL equations
//       (J'J + lambda I) delta = -J'R,
//   bisecting on lambda >= 0 (10 bracketing probes x4, 20 bisections)
//   for ||delta|| ~ Delta.
//
// The normal-equations form is essential, not cosmetic (the parent's
// analysis, kept with the transcription): J is the Lorentzian Hessian
// of the B-I functional (one positive eigenvalue, the rest negative),
// so the former shifted-J system (J + lambda I) delta = -R was singular
// at every lambda in the negative spectrum and its large-lambda limit
// is an ASCENT direction for E = 1/2 ||R||^2 -- the hard-stall mode.
// J'J + lambda I is SPD for every lambda > 0 (no poles, monotone
// ||delta(lambda)||, well-posed bisection) and the damped limit is
// steepest descent on E.  The model contract requires J symmetric to
// roundoff (AlexKappa's J is symmetric BITWISE), so J' = J and
// J'J = J * J with no transpose machinery.
//
// Consumed by TrustRegion<LevenbergBisect> over a HasLSQ model, and
// standalone by the Alexandrov polish functor (which keeps the
// feasibility clip and the flip discipline in its own loop, per the
// signed-off "monolithic for now" decision).
// =====================================================================

#include "fullerenes/optim/core.hh"

namespace optim {

struct lsq_oracle_tag {};
struct hvp_oracle_tag {};

struct LevenbergBisect {
  using paradigm = trust_region_tag;
  using oracle   = lsq_oracle_tag;

  // Pure Newton if it fits; else GN normal equations
  // (J'J + lambda I) delta = -J'R with ||delta|| <= Delta.
  // @anchor  optim-levenberg-subproblem
  // @pre     J is square symmetric with J.n == R.size(); Delta > 0
  // @post    ||result|| <= Delta up to the bisection tolerance, and
  //          result solves (J'J + lambda I) delta = -J'R for some
  //          lambda >= 0 (or is the pure Newton root step when it fits)
  // @variant 10 bracketing probes + 20 bisections
  static la::V solve_subproblem(const Mat& J, const la::V& R, double Delta) {
    la::V negR(R.size());
    for (std::size_t i = 0; i < R.size(); ++i) negR[i] = -R[i];

    // Pure Newton (lambda = 0) -- identical to the pre-GN behaviour.
    auto delta = la::solve(J, negR);
    if (la::is_usable_step(delta) && la::norm(delta) <= Delta) return delta;

    // J symmetric => J'J = J * J (SPD), -J'R = -matvec(J, R) = -grad E.
    const Mat JtJ = J * J;
    la::V minus_Jtk = la::matvec(J, negR);

    double lo = 0, hi = la::max_abs(minus_Jtk) / Delta + 1.0;
    for (int probe = 0; probe < 10; ++probe) {
      delta = la::solve_shifted(JtJ, minus_Jtk, hi);
      if (la::is_usable_step(delta) && la::norm(delta) <= Delta) break;
      hi *= 4;
    }
    for (int bis = 0; bis < 20; ++bis) {
      const double mid = 0.5 * (lo + hi);
      delta = la::solve_shifted(JtJ, minus_Jtk, mid);
      if (!la::is_usable_step(delta) || la::norm(delta) > Delta) lo = mid;
      else hi = mid;
    }
    return la::solve_shifted(JtJ, minus_Jtk, hi);
  }

  // E(R) - E(R + J delta), E = half the squared norm.
  // @anchor optim-levenberg-predicted-reduction
  // @pre    as optim-levenberg-subproblem, delta.size() == R.size()
  // @post   the linear-model reduction; may be negative when the pure
  //         Newton root step overshoots (the paradigm's rho test
  //         rejects such steps)
  static double predicted_reduction(const Mat& J, const la::V& R,
                                    const la::V& delta) {
    la::V Jd = la::matvec(J, delta);
    la::V pred(R.size());
    for (std::size_t i = 0; i < R.size(); ++i) pred[i] = R[i] + Jd[i];
    return la::energy(R) - la::energy(pred);
  }
};

}  // namespace optim
