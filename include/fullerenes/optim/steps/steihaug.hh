#pragma once

// =====================================================================
// optim::SteihaugCG -- the Steihaug-Toint truncated-CG trust-region
// subproblem solver, transcribed from Deltahedron's inner CG
// (deltahedron.cc, OptMethod::STEIHAUG):
//
//   approximately solve  min_z g'z + 1/2 z'Hz  s.t. ||z|| <= Delta
//   by CG on H z = -g, matrix-free through the model's HVP, stepping to
//   the trust-region boundary on negative curvature
//   (kappa <= 1e-15 rr) or when an iterate would leave the region, and
//   stopping the recurrence at ||r|| < 0.01 ||g|| or
//   max_inner = min(n, 200) iterations.
//
// Gauge-frozen coordinates are masked out of the residual and search
// direction exactly as the parent masks its fixed vertices.  The
// NEG_CURVATURE / TR_BOUNDARY events are reported through Result::flags.
// =====================================================================

#include "fullerenes/optim/core.hh"
#include "fullerenes/optim/steps/levenberg.hh"   // oracle tags

#include <algorithm>
#include <cmath>
#include <span>
#include <vector>

namespace optim {

struct SteihaugCG {
  using paradigm = trust_region_tag;
  using oracle   = hvp_oracle_tag;

  // Solve the subproblem; z is the returned step, flags accumulates
  // NEG_CURVATURE / TR_BOUNDARY, n_hv counts HVP evaluations.
  // @anchor  optim-steihaug-subproblem
  // @pre     g is the gradient at x, gauge-projected; Delta > 0
  // @post    ||z|| <= Delta; z approximately minimizes the quadratic
  //          model (exactly on the boundary when curvature is negative
  //          or the CG iterate leaves the region)
  // @variant min(n, 200) inner CG iterations
  template <class Model>
    requires HasHVP<Model>
  static void solve_subproblem(SeqCtx ctx, const Model& model,
                               const Gauge& gauge,
                               std::span<const double> x,
                               std::span<const double> g, double Delta,
                               std::span<double> z, uint32_t& flags,
                               int& n_hv) {
    const std::size_t n = x.size();
    const int max_inner = (int)std::min<std::size_t>(n, 200);

    std::vector<double> r(n), d(n), Hd(n), z_new(n);
    std::fill(z.begin(), z.end(), 0.0);
    for (std::size_t i = 0; i < n; ++i) r[i] = -g[i];
    gauge.project(r);
    std::copy(r.begin(), r.end(), d.begin());
    double rr = la::dot(r, r);
    const double gnorm = std::sqrt(rr);

    // ||z + tau d|| = Delta boundary step (positive root).
    auto boundary_tau = [&]() {
      const double zz = la::dot(z, z);
      const double zd = la::dot(z, d);
      const double dd = la::dot(d, d);
      const double disc = 4.0 * zd * zd - 4.0 * dd * (zz - Delta * Delta);
      return (-2.0 * zd + std::sqrt(std::max(0.0, disc))) / (2.0 * dd);
    };

    for (int j = 0; j < max_inner; ++j) {
      ++n_hv;
      model.hvp(ctx, x, d, Hd);

      const double kappa = la::dot(d, Hd);
      if (kappa <= 1e-15 * rr) {           // negative/zero curvature
        flags |= Result::NEG_CURVATURE;
        la::axpy(z, boundary_tau(), d);
        return;
      }

      const double alpha = rr / kappa;
      for (std::size_t i = 0; i < n; ++i) z_new[i] = z[i] + alpha * d[i];
      if (la::nrm2(z_new) >= Delta) {      // truncate to the boundary
        flags |= Result::TR_BOUNDARY;
        la::axpy(z, boundary_tau(), d);
        return;
      }
      std::copy(z_new.begin(), z_new.end(), z.begin());

      for (std::size_t i = 0; i < n; ++i) r[i] -= alpha * Hd[i];
      gauge.project(r);
      const double rr_new = la::dot(r, r);
      if (std::sqrt(rr_new) < 0.01 * gnorm) return;   // inner convergence

      const double beta = rr_new / rr;
      for (std::size_t i = 0; i < n; ++i) d[i] = r[i] + beta * d[i];
      gauge.project(d);
      rr = rr_new;
    }
  }
};

}  // namespace optim
