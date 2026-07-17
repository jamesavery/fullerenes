#pragma once

#include "fullerenes/graphview.hh"

#include <span>
#include <vector>

// The Wu / extended-Wu force field for fullerene 3D geometry, ported from
// the legacy Fortran (force.f + opt-standalone.f SA_OptFF) to replace the
// last C++ -> Fortran dependence, FullereneGraphView::optimized_geometry ->
// sa_optff_.
//
// References: Z. C. Wu, D. A. Jelski, T. F. George, Chem. Phys. Lett. 137,
// 291-295 (1987) (harmonic bonds + angles); the dihedral extension and the
// parameter set are the legacy program's (Fowler/Ceulemans force constants,
// see default parameters in the .cc).
//
// The topology is compiled ONCE at construction into flat term lists --
// every bond, angle and dihedral with its zero value and force constant
// already selected by pentagon content:
//
//   bonds     3N/2 terms, classified by the pentagon count of the edge's
//             two faces (hh / hp / pp),
//   angles    3N terms, one per face corner (60 pentagon + 3N-60 hexagon),
//   dihedrals N terms, one per vertex u: the chain dihedral D(u,b,c,d)
//             over u's CCW-ordered neighbours, rotated so the chain starts
//             at the neighbour opposite the "distinct" face when exactly
//             one incident face differs in size (the legacy convention the
//             zero values are calibrated to),
//
// with energy E = 1/2 sum_t k_t (c_t - c0_t)^2 (+ the 1/r-from-origin
// regularizer of methods 2/4, off by default).
//
// Units as in the legacy code: coordinates in Angstroem, angles in rad,
// force constants converted from N/m by *6.02214129 (1e-20 * N_A, i.e.
// kJ/mol/A^2). Geometry primitives are the library's coord3d::angle /
// dihedral, which are sign-exact equivalents of the Fortran's ANGLE /
// DIHEDRAL (verified: the Fortran's two normal-vector sign flips cancel in
// both atan2 arguments).
struct WuForceField {
  // Method numbering follows the legacy opt_method / iopt:
  //   1 = Wu (bonds + angles)         2 = Wu + Coulomb-from-origin
  //   3 = extended Wu (+ dihedrals)   4 = extended Wu + Coulomb
  // 3 is the default of the whole legacy pipeline. The Fortran's hyperbolic
  // variants 5/6 are not ported: default_force_parameters never filled
  // their parameters, so they were unreachable (uninitialized) through this
  // entry point.

  struct Bond     { node_t a, b;       double r0, k; };  // stretch, dist(a,b)
  struct Angle    { node_t a, b, c;    double a0, k; };  // bend at apex b
  struct Dihedral { node_t a, b, c, d; double d0, k; };  // chain dihedral abcd

  std::vector<Bond> bonds;
  std::vector<Angle> angles;
  std::vector<Dihedral> dihedrals;
  double k_coulomb = 0;  // methods 2/4; legacy default 0
  int N = 0;

  // Compile g's topology into terms. @pre g is an oriented fullerene graph.
  // @throws std::invalid_argument for method outside 1..4.
  WuForceField(const FullereneGraphView& g, int method = 3);

  double energy(std::span<const coord3d> x) const;
  // grad is overwritten with dE/dx.
  void gradient(std::span<const coord3d> x, std::span<coord3d> grad) const;

  // Minimize E from x in place (self-contained L-BFGS, m=10, Armijo
  // backtracking with the per-atom step capped at 0.3 A so descent from a
  // grossly inflated start follows a continuous path) and return the final
  // energy. Stops on the legacy relative-energy-change criterion
  // 2|dE| <= ftol (|E1|+|E2|+1e-10) on two consecutive iterations, or on a
  // vanishing gradient. max_iter is a safeguard, not a tuning knob; a trip
  // warns on stderr.
  double optimize(std::span<coord3d> x, double ftol = 1e-12,
                  int max_iter = 50000) const;

  // zero_order_geometry occasionally places two vertices at exactly the
  // same position, which makes angle/dihedral terms NaN and sent the legacy
  // sa_optff_ into a NaN spiral (observed: C40 isomer 8). Displace such
  // pairs apart (same remedy and constants as PolyhedronView::
  // optimize_other) and let the force field do the rest. Returns the number
  // of displaced pairs; deterministic.
  int separate_coincident(std::span<coord3d> x) const;
};

// Drop-in implementation of FullereneGraphView::optimized_geometry: compile
// the force field for `method` and minimize from initial_geometry.
std::vector<coord3d> wu_optimized_geometry(const FullereneGraphView& g,
                                           std::span<const coord3d> initial_geometry,
                                           int method = 3, double ftol = 1e-12);
