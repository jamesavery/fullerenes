#pragma once

// =====================================================================
// optim::AET -- the equilateral/Alexandrov-embedding-target energy of
// the Deltahedron optimizer (E_bond + E_angle + E_curv [+ E_flat]
// [+ E_conv]) as a framework EnergyModel.
//
// Delegates every evaluation to the library's AET seam
// (DeltahedronView<double>::aet_energy_gradient / aet_hv_product --
// thin wrappers over the internal term implementations, the single
// source of truth), with the flat double span viewed as coord3d[] and
// the undirected edge list cached at construction exactly as
// Deltahedron::optimize() caches it.
//
// Force constants default to optimize()'s fixed values
// (k_bond = k_angle = 1, k_curv = 2); k_flat / k_conv are the caller's
// phase configuration (k_flat > 0 only in the standalone phase 1,
// k_conv only in the constrained variants).  Eval-cost weights are the
// parent's calibrated 1 : 2 : 7.
//
// NOT shareable across threads (unlike ExtWuAngle): the const eval
// methods reuse mutable scratch vectors.  One AET instance per solve.
// =====================================================================

#include "fullerenes/optim/core.hh"

#include "fullerenes/graphview.hh"

#include <span>
#include <vector>

namespace optim {

struct AET {
  const DeltahedronView<double>& D;
  std::vector<edge_t> edges;          // cached undirected_edges()
  double L;                           // target edge length (> 0)
  double k_bond = 1.0, k_angle = 1.0, k_curv = 2.0;
  double k_flat = 0.0;
  double k_conv = 0.0, sigma_conv = 0.0;
  std::vector<bool> fixed;            // HVP gauge mask (parent convention)

  static constexpr double cost_energy = 1, cost_gradient = 2, cost_hvp = 7;

  AET(const DeltahedronView<double>& D_, double L_)
      : D(D_), edges(D_.undirected_edges()), L(L_),
        gscratch(D_.N), vscratch(D_.N), hscratch(D_.N) {}
  // The reference member + mutable scratch make copies alias D and
  // shear the one-instance-per-solve contract; a temporary view would
  // dangle.  Both are compile errors.
  AET(DeltahedronView<double>&&, double) = delete;
  AET(const AET&) = delete;
  AET& operator=(const AET&) = delete;

  double energy(SeqCtx, std::span<const double> x) const {
    return D.aet_energy_gradient(edges, as_coords(x), nullptr, L,
                                 k_bond, k_angle, k_curv, k_flat,
                                 k_conv, sigma_conv);
  }

  double gradient(SeqCtx, std::span<const double> x,
                  std::span<double> g) const {
    const double E = D.aet_energy_gradient(edges, as_coords(x), &gscratch, L,
                                           k_bond, k_angle, k_curv, k_flat,
                                           k_conv, sigma_conv);
    const double* p = reinterpret_cast<const double*>(gscratch.data());
    for (std::size_t i = 0; i < g.size(); ++i) g[i] = p[i];
    return E;
  }

  void hvp(SeqCtx, std::span<const double> x, std::span<const double> v,
           std::span<double> Hv) const {
    const auto vc = as_coords(v);
    std::copy(vc.begin(), vc.end(), vscratch.begin());
    D.aet_hv_product(edges, as_coords(x), vscratch, hscratch, L,
                     k_bond, k_angle, k_curv, k_flat, fixed,
                     k_conv, sigma_conv);
    const double* p = reinterpret_cast<const double*>(hscratch.data());
    for (std::size_t i = 0; i < Hv.size(); ++i) Hv[i] = p[i];
  }

 private:
  mutable std::vector<coord3d> gscratch, vscratch, hscratch;
};

}  // namespace optim
