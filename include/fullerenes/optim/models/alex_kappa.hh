#pragma once

// =====================================================================
// optim::AlexKappa -- the Alexandrov B-I curvature residual kappa(r) as
// a framework EnergyModel with the least-squares face:
//
//   residual       R(r) = kappa(T, r)      (per-cone angle deficit)
//   system_matrix  J(r) = d kappa / d r    (dense symmetric, SQUARE)
//
// both delegated to the library's AlexandrovSolver statics (the single
// source of truth for the GCP geometry).  The smooth face is derived:
// E = 1/2 ||kappa||^2, grad E = J' kappa.  The state span is the radius
// vector r (one double per cone vertex), NOT 3D coordinates.
//
// The triangulation reference is mutable on purpose: moving r can make
// T non-weighted-Delaunay, and the re-flip
// (AlexandrovSolver::flip_to_weighted_delaunay) mutates T between
// accepted steps.  By the signed-off "monolithic for now" decision that
// re-flip and the feasibility clip live in the POLISH FUNCTOR
// (alex_polish.hh), not in model hooks -- this model is the pure
// residual/Jacobian face.
// =====================================================================

#include "fullerenes/optim/core.hh"

#include "fullerenes/delaunay_alexandrov.hh"

#include <span>
#include <vector>

namespace optim {

struct AlexKappa {
  DelaunayTriangulation& T;

  std::size_t residual_size() const { return (std::size_t)T.nv; }

  void residual(SeqCtx, std::span<const double> r, std::span<double> R) const {
    const auto k = AlexandrovSolver::kappa(T, la::V(r.begin(), r.end()));
    la::copy(k, R);
  }

  void system_matrix(SeqCtx, std::span<const double> r, Mat& J) const {
    J = AlexandrovSolver::jacobian(T, la::V(r.begin(), r.end()));
  }

  // Smooth face derived from the LSQ face.
  double energy(SeqCtx, std::span<const double> r) const {
    return la::energy(AlexandrovSolver::kappa(T, la::V(r.begin(), r.end())));
  }

  double gradient(SeqCtx, std::span<const double> r,
                  std::span<double> g) const {
    const la::V rv(r.begin(), r.end());
    const auto k = AlexandrovSolver::kappa(T, rv);
    const Mat J = AlexandrovSolver::jacobian(T, rv);
    // g = J' kappa = J kappa (J symmetric up to roundoff).
    la::copy(la::matvec(J, k), g);
    return la::energy(k);
  }
};

}  // namespace optim
