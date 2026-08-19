#pragma once

// =====================================================================
// optim::core -- the problem vocabulary of the unified optimizer
// framework (claude-projects/optimize/DESIGN.md, first landing).
//
// An optimizer is  globalize(step)  over an EnergyModel: a globalization
// PARADIGM (LineSearch / TrustRegion) owns the outer iteration loop, and
// the STEP (L-BFGS, PR-CG, Steihaug-CG, Levenberg-bisect) is a policy
// valid only within its paradigm.  This header defines everything both
// sides share:
//
//   * the model contract (concepts SmoothModel / HasHVP / HasHessian /
//     HasLSQ) -- a scheme requiring a face a model lacks fails to
//     COMPILE, not at runtime;
//   * the paradigm/step pairing tags (a Steihaug step in a line search
//     is a compile error);
//   * Criteria / Outcome / Result -- convergence tests and the named
//     exit status;
//   * Gauge -- the free-DOF mask (eliminated coordinates, not a
//     constraint);
//   * Problem -- model + gauge + (future) constraints + criteria.
//
// State is a flat std::span<double>; geometric models reinterpret it
// (3N doubles as N coord3d, or nv cone radii), exactly as minimize.hh
// does in the parent library.
//
// Convergence semantics (matching every optimizer being consolidated,
// see minimize.hh and Deltahedron::optimize): each ENABLED tolerance is
// a SUFFICIENT condition for CONVERGED -- a disjunction, with 0
// disabling a test.  Budget and stagnation trips exit with their own
// named outcomes instead.  Iteration/work budgets are safeguards, not
// tuning knobs: they should never be reached, and a trip is reported as
// the defect it is.
// =====================================================================

#include "fullerenes/optim/seq_ctx.hh"
#include "fullerenes/optim/linalg.hh"

#include "fullerenes/geometry.hh"

#include <algorithm>
#include <cmath>
#include <concepts>
#include <type_traits>
#include <cstdint>
#include <functional>
#include <optional>
#include <span>
#include <vector>

namespace optim {

// The flat-state convention: 3N doubles viewed as N coord3d, the same
// reinterpretation minimize.hh's geometric callers perform.  THE shared
// primitive -- models and drivers use these, never a bare cast.
static_assert(sizeof(coord3d) == 3 * sizeof(double) &&
              std::is_standard_layout_v<coord3d> &&
              alignof(coord3d) == alignof(double),
              "flat spans are reinterpreted as coord3d[] -- requires "
              "size, layout AND alignment compatibility");

inline std::span<const coord3d> as_coords(std::span<const double> x) {
  return { reinterpret_cast<const coord3d*>(x.data()), x.size() / 3 };
}
inline std::span<coord3d> as_coords(std::span<double> x) {
  return { reinterpret_cast<coord3d*>(x.data()), x.size() / 3 };
}
inline std::span<const double> as_flat(std::span<const coord3d> x) {
  return { reinterpret_cast<const double*>(x.data()), 3 * x.size() };
}
inline std::span<double> as_flat(std::span<coord3d> x) {
  return { reinterpret_cast<double*>(x.data()), 3 * x.size() };
}

// ---------------------------------------------------------------------
// Paradigm/step pairing.  Every step policy declares the paradigm it is
// valid in via  `using paradigm = line_search_tag;`  and the paradigm
// class templates constrain on these concepts, so an invalid pairing
// (LineSearch<SteihaugCG>, TrustRegion<LBFGS>) is rejected at the
// template constraint with a readable diagnostic.
// ---------------------------------------------------------------------
struct line_search_tag {};
struct trust_region_tag {};
struct sqp_tag {};              // ActiveSetSQP -- migration step 5, absent here

template <class S>
concept LineSearchStep = std::same_as<typename S::paradigm, line_search_tag>;
template <class S>
concept TrustRegionStep = std::same_as<typename S::paradigm, trust_region_tag>;

// ---------------------------------------------------------------------
// EnergyModel capability concepts (the "faces" of DESIGN.md 3.1).
//
// General smooth face (always):
//   double energy  (Ctx, span<const double> x) const;
//   double gradient(Ctx, span<const double> x, span<double> g) const;
//     -- returns E(x) and fills g = grad E(x) in one pass.
// Optional faces:
//   hvp(Ctx, x, v, Hv)          matrix-free Hessian-vector product
//   hessian(Ctx, x, Mat& H)     assembled dense Hessian
//   residual/system_matrix      the least-squares face: R(x) and the
//     SQUARE system matrix J = dR/dx consumed by LevenbergBisect --
//     the parent solvers' (J + lambda I) delta = -R scheme, NOT a
//     textbook JtJ Gauss-Newton (that step arrives with migration
//     step 5 for genuinely rectangular Jacobians).
// ---------------------------------------------------------------------
template <class M>
concept SmoothModel =
    requires(const M m, SeqCtx c, std::span<const double> x, std::span<double> g) {
      { m.energy(c, x) } -> std::convertible_to<double>;
      { m.gradient(c, x, g) } -> std::convertible_to<double>;
    };

template <class M>
concept HasHVP = SmoothModel<M> &&
    requires(const M m, SeqCtx c, std::span<const double> x,
             std::span<const double> v, std::span<double> Hv) {
      m.hvp(c, x, v, Hv);
    };

template <class M>
concept HasHessian = SmoothModel<M> &&
    requires(const M m, SeqCtx c, std::span<const double> x, Mat& H) {
      m.hessian(c, x, H);
    };

template <class M>
concept HasLSQ =
    requires(const M m, SeqCtx c, std::span<const double> x,
             std::span<double> R, Mat& J) {
      { m.residual_size() } -> std::convertible_to<std::size_t>;
      m.residual(c, x, R);
      m.system_matrix(c, x, J);
    };

// Per-model evaluation-cost weights for the work budget.  A model may
// override by defining  static constexpr double cost_energy / _gradient
// / _hvp;  the defaults are Deltahedron's calibrated 1 : 2 : 7 --
// explicitly NOT universal (DESIGN.md 4).
template <class M>
struct eval_costs {
  static constexpr double energy   = 1;
  static constexpr double gradient = 2;
  static constexpr double hvp      = 7;
};
template <class M>
  requires requires { M::cost_energy; M::cost_gradient; M::cost_hvp; }
struct eval_costs<M> {
  static constexpr double energy   = M::cost_energy;
  static constexpr double gradient = M::cost_gradient;
  static constexpr double hvp      = M::cost_hvp;
};

// ---------------------------------------------------------------------
// Outcome / Result -- the named exit status (style-failures.md).
// ---------------------------------------------------------------------
enum class Outcome {
  CONVERGED,          // an enabled tolerance fired
  STAGNATED,          // energy stopped decreasing (stagnation_window trips)
  BUDGET_EXHAUSTED,   // work/iteration safeguard tripped (a defect, not a result)
  STEP_FAILED,        // line search / subproblem could not make progress
  INFEASIBLE,         // iterate left the feasible region (constrained paths)
  CONSTRAINT_STUCK,   // feasibility restoration failed (reflect/hull, SQP)
  SATURATED           // a Continuation stopped before its terminal parameters
                      // (driver.hh: !cont.terminal() at exit; schedule-level
                      // exhaustion, not stagnation of an inner solve)
};

constexpr const char* outcome_name(Outcome o) {
  switch (o) {
    case Outcome::CONVERGED:        return "CONVERGED";
    case Outcome::STAGNATED:        return "STAGNATED";
    case Outcome::BUDGET_EXHAUSTED: return "BUDGET_EXHAUSTED";
    case Outcome::STEP_FAILED:      return "STEP_FAILED";
    case Outcome::INFEASIBLE:       return "INFEASIBLE";
    case Outcome::CONSTRAINT_STUCK: return "CONSTRAINT_STUCK";
    case Outcome::SATURATED:        return "SATURATED";
  }
  return "?";
}

struct Result {
  Outcome  outcome  = Outcome::BUDGET_EXHAUSTED;
  double   f        = 0;      // final value
  double   gmax     = 0;      // final gradient norm in Criteria::gnorm's convention
  int      iters    = 0;      // outer iterations taken
  int      n_accepted = 0;    // accepted steps (line search: = iters; trust region: iters - rejections)
  int      n_energy = 0;      // energy-only oracle calls
  int      n_grad   = 0;      // energy+gradient oracle calls
  int      n_hv     = 0;      // Hessian-vector products
  uint32_t flags    = 0;      // diagnostic events (Result flag bits below)

  // Optimizer-level diagnostic flags (a subset of the parent
  // PipelineDiag vocabulary, same semantics).
  static constexpr uint32_t NEG_CURVATURE = 1u << 0;  // Steihaug hit negative curvature
  static constexpr uint32_t TR_BOUNDARY   = 1u << 1;  // Steihaug truncated at the radius
  static constexpr uint32_t STEP_REJECTED = 1u << 2;  // trust region rejected a step
  static constexpr uint32_t CVX_REJECTED  = 1u << 3;  // convexity acceptance gate fired
  static constexpr uint32_t LBFGS_RESET   = 1u << 4;  // history dropped on non-descent

  bool converged() const { return outcome == Outcome::CONVERGED; }

  template <class M>
  double work() const {
    return eval_costs<M>::energy   * n_energy
         + eval_costs<M>::gradient * n_grad
         + eval_costs<M>::hvp      * n_hv;
  }
};

// ---------------------------------------------------------------------
// Criteria -- which convergence tests are enabled (0 = off), plus the
// safeguards.  Domain predicates are model-space success tests over the
// current iterate (e.g. max angle relative error), each sufficient for
// CONVERGED like the norm tests.
// ---------------------------------------------------------------------
using DomainPredicate = std::function<bool(std::span<const double> x)>;

// The gradient-norm CONVENTION the gtol_inf test (and Result::gmax) use:
//   FLAT_INF   max_i |g_i|                        (minimize.hh / wu)
//   FLAT_2     ||g||_2
//   BLOCK_INF  scale * max_b ||g_b||_2 over consecutive blocks of
//              `block` entries (Deltahedron's gmax_L: block = 3,
//              scale = L -- the per-vertex 2-norm, scaled to be
//              dimensionless in the edge length)
// Frozen (gauge) entries are zero in g, so a block norm over all
// blocks equals the parent's max over FREE vertices.
enum class GradNorm { FLAT_INF, FLAT_2, BLOCK_INF };

struct Criteria {
  GradNorm gnorm = GradNorm::FLAT_INF;   // convention for gtol_inf / Result::gmax
  int      gnorm_block = 1;              // BLOCK_INF: entries per block
  double   gnorm_scale = 1;              // BLOCK_INF: multiplier

  double gtol_inf = 0;      // gnorm(g) <= gtol_inf
  double gtol_2   = 0;      // ||g||_2   <= gtol_2
  double gtol_2_rel = 0;    // ||g||_2 <= gtol_2_rel * ||g_0||_2, g_0 = the solve's
                            // first gradient (Deltahedron's phase-1 exit: 0.01)
  double ftol_rel = 0;      // 2|df| <= ftol_rel (|f1|+|f0|+1e-9), required on
                            // TWO consecutive iterations (minimize.hh semantics)
  double step_tol = 0;      // ||step||_inf <= step_tol
  double rtol_inf = 0;      // LSQ face: max |R| <= rtol_inf (Alexandrov kappa)

  std::vector<DomainPredicate> domain;  // each sufficient for CONVERGED

  // Safeguards (never-reached budgets; a trip is reported, not tuned away).
  int       max_iters         = 20000;
  double    max_work          = 0;   // ceiling on eval_costs-weighted work (0 = off)
  int       stagnation_window = 0;   // iters without meaningful decrease -> STAGNATED (0 = off)
};

// The gradient norm in the Criteria's convention (header note above).
// @post  FLAT_INF: max|g_i|; FLAT_2: ||g||_2; BLOCK_INF: scale * max over
//        blocks of the block 2-norm (a trailing partial block counts)
inline double grad_norm(const Criteria& c, std::span<const double> g) {
  switch (c.gnorm) {
    case GradNorm::FLAT_INF: return la::inf_norm(g);
    case GradNorm::FLAT_2:   return la::nrm2(g);
    case GradNorm::BLOCK_INF: {
      const std::size_t b = (std::size_t)std::max(1, c.gnorm_block);
      double m = 0;
      for (std::size_t i = 0; i < g.size(); i += b) {
        double s = 0;
        for (std::size_t j = i; j < std::min(i + b, g.size()); ++j) s += g[j] * g[j];
        m = std::max(m, std::sqrt(s));
      }
      return m * c.gnorm_scale;
    }
  }
  return la::inf_norm(g);
}

// ---------------------------------------------------------------------
// Gauge -- eliminated degrees of freedom: rigid-motion gauge fixing and
// caller pins, as a free-DOF mask.  NOT a constraint: frozen entries
// simply never move (their gradient and step components are zeroed).
// Empty mask = all free.
// ---------------------------------------------------------------------
struct Gauge {
  std::vector<uint8_t> free;   // free[i] == 0 freezes coordinate i; empty = all free

  bool trivial() const { return free.empty(); }
  void project(std::span<double> v) const {   // zero frozen components
    if (free.empty()) return;
    for (std::size_t i = 0; i < v.size(); ++i)
      if (!free[i]) v[i] = 0;
  }
};

// ---------------------------------------------------------------------
// Problem -- WHAT to minimize, and where.  ConstraintSet is the
// SQP-only object of migration step 5; unconstrained problems leave it
// null and never pay for it.
// ---------------------------------------------------------------------
struct ConstraintSet;   // defined with the ActiveSetSQP paradigm (step 5)

template <SmoothModel Model>
struct Problem {
  Model&               E;
  Gauge                gauge{};
  const ConstraintSet* C = nullptr;
  Criteria             stop{};
};

// ---------------------------------------------------------------------
// Hooks -- the paradigm loops' extension points (DESIGN.md 10.2).  Each
// has at least two customers among the drivers that used to hand-roll
// their loops around the missing one; all default-empty, zero cost
// when unused.
//
//   clip         before the trial: the max admissible fraction s in
//                (0, 1] of the proposed step d from x (line search:
//                caps the first trial t; trust region: scales the
//                subproblem step and caps the radius at the length
//                taken).  Alexandrov's F(T)-feasible step, the embed's
//                fraction-to-boundary.
//   veto         on a trial point: false = the point is inadmissible --
//                a trust region rejects the step (CVX_REJECTED), a line
//                search treats it as infinite energy and keeps
//                backtracking.  Deltahedron's convex-constraint gate,
//                the embed's nonlinear feasibility check.
//   post_accept  after an accepted step: may move x (a points->points
//                projection: Deltahedron's reflect/hull, DESIGN.md 3.3)
//                or mutate the MODEL through its own capture
//                (Alexandrov's weighted-Delaunay re-flip); returns the
//                number of changes.  Nonzero => the paradigm re-anchors
//                f, g at x and drops the step's history.
//   observe      once per iteration, a read-only IterationView (all
//                diagnostic logging lives here, not in the loops).
// ---------------------------------------------------------------------
struct IterationView {
  int      iter      = 0;
  double   f         = 0;      // energy after this iteration
  double   gmax      = 0;      // gradient norm (Criteria::gnorm convention)
  std::span<const double> x, g;
  double   step_len  = 0;      // line search: t taken; trust region: ||step||_2
  double   radius    = 0;      // trust region: Delta after this iteration's update
  double   rho       = 0;      // trust region: reduction ratio (0 for line search)
  int      n_inner   = 0;      // trust region: subproblem HVPs this iteration
  int      n_post    = 0;      // changes reported by post_accept this iteration
  bool     accepted  = true;
  uint32_t flags     = 0;      // Result flag bits raised so far
};
using StepClip   = std::function<double(std::span<const double> x, std::span<const double> d)>;
using AcceptVeto = std::function<bool(std::span<const double> x_trial)>;
using PostAccept = std::function<int(SeqCtx, std::span<double> x)>;
using Observer   = std::function<void(const IterationView&)>;

struct Hooks {
  StepClip   clip;
  AcceptVeto veto;
  PostAccept post_accept;
  Observer   observe;
};

// THE convergence battery: the disjunction-of-enabled-tests semantics
// (header note above) written once, shared by every paradigm loop AND
// the delta driver.  Returns the outcome that fired, or nullopt to
// keep iterating.  gmax is the caller's gradient norm in ITS
// convention (grad_norm(stop, g)); g2_0 is ||g||_2 at the solve's
// start for the relative test.  The ftol_rel two-consecutive test
// stays with the accept path (it needs f0/f1), as does step_tol.
// @anchor optim-converged-or-stop
// @post   result == CONVERGED    iff an enabled tolerance/predicate fired
// @post   result == STAGNATED    iff the stagnation window tripped
// @post   result == BUDGET_EXHAUSTED iff the work ceiling tripped
inline std::optional<Outcome> converged_or_stop(
    const Criteria& stop, double gmax, std::span<const double> g,
    std::span<const double> x, int stag_count, double work_now,
    double g2_0 = 0) {
  if (stop.gtol_inf > 0 && gmax <= stop.gtol_inf) return Outcome::CONVERGED;
  if (stop.gtol_2 > 0 && std::sqrt(la::dot(g, g)) <= stop.gtol_2)
    return Outcome::CONVERGED;
  if (stop.gtol_2_rel > 0 && g2_0 > 0 && la::nrm2(g) <= stop.gtol_2_rel * g2_0)
    return Outcome::CONVERGED;
  if (!stop.domain.empty() &&
      std::any_of(stop.domain.begin(), stop.domain.end(),
                  [&](const DomainPredicate& pred) { return pred(x); }))
    return Outcome::CONVERGED;
  if (stop.stagnation_window > 0 && stag_count >= stop.stagnation_window)
    return Outcome::STAGNATED;
  if (stop.max_work > 0 && work_now >= stop.max_work)
    return Outcome::BUDGET_EXHAUSTED;
  return std::nullopt;
}

}  // namespace optim
