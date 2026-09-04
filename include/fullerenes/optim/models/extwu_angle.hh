#pragma once

// =====================================================================
// optim::ExtWuAngle -- the Wu / extended-Wu molecular force field
// (angle-space) as a framework EnergyModel.
//
// Borrows the immutable classified term lists of a wu::ForceField
// (fullerenes/wu_forcefield.hh) and delegates energy / gradient to it;
// the model owns no state and is shareable across threads like the
// field itself.  The flat state span holds 3N doubles viewed as N
// coord3d, the same reinterpretation wu::optimize performs.
//
// All three faces delegate to the field (energy / energy_gradient /
// hvp -- the last promoted into wu::ForceField at graduation,
// implemented on the geometry_hessians.hh dual-number primitives and
// FD-verified by test_extwu_hvp); the assembled hessian face is the
// framework-side small-N diagnostic.
//
// NOTE the angle-space / cosine-space distinction (DESIGN.md 3.1): the
// SYCL cosine-space ExtWu has the same equilibria but different
// stiffness -- it is a SEPARATE model, not this one with different k.
// =====================================================================

#include "fullerenes/optim/core.hh"

#include "fullerenes/wu_forcefield.hh"

#include <span>

namespace optim {

struct ExtWuAngle {
  const wu::ForceField& FF;

  double energy(SeqCtx, std::span<const double> x) const {
    return FF.energy(as_coords(x));
  }

  double gradient(SeqCtx, std::span<const double> x,
                  std::span<double> g) const {
    return FF.energy_gradient(as_coords(x), as_coords(g));
  }

  // Matrix-free Hessian-vector product, delegated to the field's
  // promoted implementation (wu::ForceField::hvp -- graduated from this
  // model's incubation; exact to roundoff, FD-verified by
  // test_extwu_hvp).
  void hvp(SeqCtx, std::span<const double> x, std::span<const double> v,
           std::span<double> Hv) const {
    FF.hvp(as_coords(x), as_coords(v), as_coords(Hv));
  }

  // Assembled dense Hessian, column by column via the HVP.  Small-N
  // diagnostics and tests only -- O(3N) HVP passes.
  void hessian(SeqCtx ctx, std::span<const double> x, Mat& Hm) const {
    const std::size_t n = x.size();
    Hm = Mat(n, n);
    std::vector<double> e(n, 0.0), col(n);
    for (std::size_t j = 0; j < n; ++j) {
      e[j] = 1;
      hvp(ctx, x, e, col);
      e[j] = 0;
      for (std::size_t i = 0; i < n; ++i) Hm(i, j) = col[i];
    }
  }
};

}  // namespace optim
