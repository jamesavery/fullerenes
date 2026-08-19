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
// AlexKappa is the FIXED-CELL face: kappa on the triangulation as
// given, for any r.  kappa is only piecewise-smooth in r -- moving r
// can make T non-weighted-Delaunay -- so the polish needs the
// CELL-RESOLVED function AlexKappaWD below: every evaluation at a
// point r is taken on r's own weighted-Delaunay cell (a flipped scratch
// copy of T; the parent polish's "trial on its own flipped copy"
// invariant), and commit(r) adopts that cell as the base (the parent's
// "adopt the copy atomically on acceptance").  The two are distinct
// models of distinct functions, not one model with a flag.
// =====================================================================

#include "fullerenes/optim/core.hh"

#include "fullerenes/delaunay_alexandrov.hh"

#include <algorithm>
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

// Cell-resolved kappa: R(r) = kappa(WD(T, r), r), where WD(T, r) is the
// weighted-Delaunay triangulation reached from the BASE cell by
// AlexandrovSolver::flip_to_weighted_delaunay.  The base cell T is
// weighted-Delaunay for the committed point r_base (established by the
// entry flip in the constructor); an evaluation at r == r_base uses T
// directly, any other r a flipped scratch copy (cached on the last
// query point, so the accepted trial's cell is adopted without a
// second flip).  commit(r) is the post_accept hook's body.
// @inv  T is weighted-Delaunay for r_base
struct AlexKappaWD {
  DelaunayTriangulation& T;

  // @post T is weighted-Delaunay for r0 (the parent's ENTRY invariant)
  AlexKappaWD(DelaunayTriangulation& T_, std::span<const double> r0)
      : T(T_), r_base(r0.begin(), r0.end()) {
    AlexandrovSolver::flip_to_weighted_delaunay(T, r_base);
  }

  std::size_t residual_size() const { return (std::size_t)T.nv; }

  void residual(SeqCtx, std::span<const double> r, std::span<double> R) const {
    la::copy(AlexandrovSolver::kappa(cell(r), la::V(r.begin(), r.end())), R);
  }
  void system_matrix(SeqCtx, std::span<const double> r, Mat& J) const {
    J = AlexandrovSolver::jacobian(cell(r), la::V(r.begin(), r.end()));
  }
  double energy(SeqCtx, std::span<const double> r) const {
    return la::energy(AlexandrovSolver::kappa(cell(r), la::V(r.begin(), r.end())));
  }
  double gradient(SeqCtx, std::span<const double> r, std::span<double> g) const {
    const la::V rv(r.begin(), r.end());
    const DelaunayTriangulation& Tc = cell(r);
    const auto k = AlexandrovSolver::kappa(Tc, rv);
    la::copy(la::matvec(AlexandrovSolver::jacobian(Tc, rv), k), g);
    return la::energy(k);
  }

  // Adopt r's cell as the base (the accepted trial's flipped copy when
  // it is the cached query point; a fresh flip otherwise).
  // @post T is weighted-Delaunay for r; r_base == r
  void commit(std::span<const double> r) {
    if (same(r, r_scratch)) {
      T = std::move(T_scratch);
      r_scratch.clear();
    } else {
      AlexandrovSolver::flip_to_weighted_delaunay(T, la::V(r.begin(), r.end()));
    }
    r_base.assign(r.begin(), r.end());
  }

 private:
  la::V r_base;                               // the committed point
  mutable la::V r_scratch;                    // last non-base query point (empty: none)
  mutable DelaunayTriangulation T_scratch;    // its weighted-Delaunay cell

  static bool same(std::span<const double> a, const la::V& b) {
    return a.size() == b.size() && std::equal(a.begin(), a.end(), b.begin());
  }
  const DelaunayTriangulation& cell(std::span<const double> r) const {
    if (same(r, r_base)) return T;
    if (!same(r, r_scratch)) {
      T_scratch = T;
      r_scratch.assign(r.begin(), r.end());
      AlexandrovSolver::flip_to_weighted_delaunay(T_scratch, r_scratch);
    }
    return T_scratch;
  }
};

}  // namespace optim
