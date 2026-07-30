#pragma once

// =====================================================================
// optim::Identity -- the trivial preconditioner, and the lifecycle
// contract every preconditioner follows.
//
// A preconditioner is a factorization the STEP holds (DESIGN.md 3, not
// a free-standing axis):
//
//   void build(Ctx, const Problem&, span<const double> x);
//     -- factor once (or refresh per continuation level);
//   void apply(Ctx, span<const double> r, span<double> z) const;
//     -- z = M^{-1} r.
//
// Factored preconditioners (Jacobi diagonal, cotan-Laplacian IC0,
// barrier-augmented) arrive with the paths that need them; this landing
// uses only Identity, whose apply is a copy.
// =====================================================================

#include "fullerenes/optim/seq_ctx.hh"
#include "fullerenes/optim/linalg.hh"

#include <span>

namespace optim {

struct Identity {
  static constexpr bool refresh_per_level = false;

  template <class Problem>
  void build(SeqCtx, const Problem&, std::span<const double>) {}

  void apply(SeqCtx, std::span<const double> r, std::span<double> z) const {
    la::copy(r, z);
  }
};

}  // namespace optim
