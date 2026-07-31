#pragma once

// =====================================================================
// wu -- the Wu / extended-Wu harmonic force field for fullerene 3D
// geometry, with its minimization driver.
//
//   Z. C. Wu, D. A. Jelski, T. F. George, "Vibrational Motions of
//   Buckminsterfullerene", Chem. Phys. Lett. 137, 291-295 (1987).
//
// This is the C++ replacement for the legacy Fortran optimizer
// (SA_OptFF in opt-standalone.f + extwu/dextwu in force.f), built as
// three layers:
//
//   1. Harmonic terms.  The energy is a sum of identical-shape terms
//         E = 1/2 sum_t k_t h(q_t(x) - q0_t)
//      where q_t is a geometric coordinate (bond length, interior
//      angle, signed dihedral), q0_t its zero value, k_t its force
//      constant, and h(d) = d^2 (hard) or d^2/sqrt(4+d^2) (soft,
//      variants 5/6).  Optionally + 1/2 k_c sum_u 1/|x_u| (Coulomb
//      origin repulsion, variants 2/4/6).
//
//   2. Term classification.  forcefield(G, variant) reads the face
//      structure of the fullerene graph G once and emits the typed
//      term lists: bonds keyed by adjacent-pentagon count (pp/hp/hh),
//      face corners keyed by face size (pentagon/hexagon), and
//      per-vertex dihedrals keyed by the sizes of the three incident
//      faces (555/556/566/666), each with its zero value and force
//      constant from the parameter table.
//
//   3. Minimization.  optimize() runs minimize::lbfgs on the summed
//      energy with analytic gradient, converging on the same
//      relative-energy criterion as the Fortran CG driver.
//
// The geometric coordinates and their derivatives are the library
// primitives coord3d::angle/dangle and coord3d::dihedral/ddihedral,
// which implement exactly the Fortran ANGLE/DANGLE/DIHEDRAL/DDIHEDRAL
// conventions (angle at the middle atom by the law of cosines,
// signed dihedral by atan2).  Energies therefore agree with the
// Fortran force field to floating-point roundoff at identical
// coordinates -- see delaunay-fillin/test_wu_forcefield.cc.
//
// Unlike the Fortran (module-global number_vertices), evaluation is
// thread-safe: a ForceField is immutable after construction and may
// be shared across threads; optimize() mutates only its span.
// =====================================================================

#include "fullerenes/geometry.hh"
#include "fullerenes/graphview.hh"
#include "fullerenes/minimize.hh"

#include <array>
#include <span>
#include <vector>

namespace wu {

// One harmonic term 1/2 k h(q(x) - q0) over K atom sites.  q is the
// bond length |x1 - x0| (K=2), the interior angle at atoms[1] (K=3),
// or the signed dihedral of the 4-tuple (K=4, Fortran quadruple
// convention: the angle between planes (a0,a1,a2) and (a1,a2,a3)).
template <int K>
struct Harmonic {
    std::array<node_t, K> atoms;
    double q0 = 0;    // zero value      [rad or Aangstroem]
    double k  = 0;    // force constant  [kJ/mol/rad^2 or kJ/mol/A^2]
};

using Bond     = Harmonic<2>;
using Corner   = Harmonic<3>;
using Dihedral = Harmonic<4>;

// Force-field parameter table in paper units: bond lengths in
// Aangstroem, angles/dihedrals in degrees, force constants in N/m --
// exactly the table of default_force_parameters (opt-standalone.f).
// forcefield() converts to internal units (rad, kJ/mol via the
// N/m -> kJ/mol/A^2 factor N_A * 1e-3 = 6.02214129).
//
// Bond classes are keyed by the number of pentagons adjacent to the
// edge; dihedral classes by the number of pentagons among the three
// faces at the vertex.
struct Parameters {
    double R55, R56, R66;             // bond zero values by pentagon count 2/1/0
    double A5, A6;                    // face interior angles
    double D555, D556, D566, D666;    // vertex dihedrals by pentagon count 3/2/1/0
    double fR55, fR56, fR66;
    double fA5, fA6;
    double fD555, fD556, fD566, fD666;
    double fCoulomb;

    // Extended Wu (Fortran iopt 3..6): 3 bond + 2 angle + 4 dihedral classes.
    static constexpr Parameters extwu() {
        return { 1.479, 1.458, 1.401,
                 108.0, 120.0,
                 37.38, 29.20, 23.49, 0.0,
                 260.0, 390.0, 450.0,
                 100.0, 100.0,
                 35.0, 65.0, 85.0, 270.0,
                 0.0 };
    }
    // Original Wu (Fortran iopt 1/2): a pentagon-adjacent bond class
    // (R5 = both pp and hp) and a hexagon-hexagon class; no dihedrals.
    static constexpr Parameters wu() {
        return { 1.455, 1.455, 1.391,
                 108.0, 120.0,
                 0, 0, 0, 0,
                 390.7, 390.7, 499.7,
                 47.88 * 1.45 * 1.45, 80.86 * 1.45 * 1.37,
                 0, 0, 0, 0,
                 0.0 };
    }
};

// The classified force field of one fullerene graph: immutable term
// lists + the two functional-form switches.
struct ForceField {
    std::vector<Bond>     bonds;
    std::vector<Corner>   corners;
    std::vector<Dihedral> dihedrals;
    double k_coulomb = 0;    // 1/2 k_c sum_u 1/|x_u| (0 = off)
    bool   soft      = false; // h(d) = d^2/sqrt(4+d^2) instead of d^2

    // E(x).
    // @anchor wu-energy
    // @pre  x.size() == number of graph vertices the terms index
    // @pure
    double energy(std::span<const coord3d> x) const;

    // E(x) and g = grad E(x) in one pass.
    // @anchor wu-energy-gradient
    // @pre  g.size() == x.size()
    // @post result == energy(x)
    double energy_gradient(std::span<const coord3d> x,
                           std::span<coord3d>       g) const;

    // Hv = (grad^2 E)(x) . v -- the matrix-free Hessian-vector product,
    // exact to roundoff for all four term classes (bond, corner,
    // dihedral, Coulomb).  The angle/dihedral second derivatives come
    // from a forward-mode dual pass through the gradient formulas
    // (fullerenes/optim/models/geometry_hessians.hh, pinned to
    // coord3d::dangle/ddihedral at 1e-12); FD-verified against
    // energy_gradient by claude-projects/optimize's test_extwu_hvp.
    // Feeds Newton/Steihaug steps of the optimizer framework.
    // @anchor wu-hvp
    // @pre  v.size() == x.size() && Hv.size() == x.size()
    // @post Hv is the exact directional derivative of energy_gradient's
    //       g along v (symmetric: v.hvp(x,w) == w.hvp(x,v))
    void hvp(std::span<const coord3d> x, std::span<const coord3d> v,
             std::span<coord3d> Hv) const;
};

// Build the classified force field of G under the legacy variant
// numbering (the iopt of the Fortran path, kept for API continuity):
//
//   1 Wu       2 Wu + Coulomb
//   3 ExtWu    4 ExtWu + Coulomb      (the library default is 3)
//   5 ExtWu soft   6 ExtWu soft + Coulomb
//
// @anchor wu-forcefield
// @pre  fullerene: G.is_a_fullerene() (cubic, faces of size 5 and 6 only)
// @pre  variant >= 1 && variant <= 6
// @post result.bonds.size() == 3*G.N/2 &&
//           result.corners.size() == 3*G.N &&
//           (variant <= 2 || result.dihedrals.size() == G.N)
// @throws std::invalid_argument when !(variant >= 1 && variant <= 6)
ForceField forcefield(const FullereneGraphView& G, int variant = 3,
                      const Parameters* params = nullptr);

// Minimize FF over x in place (L-BFGS, analytic gradient), stopping
// on the Fortran-compatible relative-energy criterion ftol.  For a
// Coulomb-carrying variant, the repulsion is switched off once
// ||g||_2 <= 10 and minimization continues without it -- the same
// two-phase schedule as the Fortran driver (SA_frprmn3d, iopt 4).
// @anchor wu-optimize
// @pre  x.size() == number of graph vertices FF's terms index
minimize::Outcome optimize(const ForceField& FF, std::span<coord3d> x,
                           double ftol = 1e-12);

// Start-geometry repair: zero_order_geometry occasionally places two
// vertices at exactly the same position (observed: C40 isomer 8), which
// makes bond/angle terms NaN and poisons any minimization from that
// start.  Displace such pairs apart (same remedy and constants as
// PolyhedronView::optimize_other) and let the force field do the rest.
// Deterministic; no-op on well-separated geometry.  Returns the number
// of displaced pairs.
// @anchor wu-separate-coincident
int separate_coincident(std::span<coord3d> x);

}  // namespace wu
