#include "fullerenes/wu_forcefield.hh"

#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace wu {

// =====================================================================
// Deviation profile h(d) and its derivative: the harmonic well d^2, or
// the hyperbolically softened d^2/sqrt(4+d^2) (Fortran variants 5/6,
// hyperbolic_par_2 = 4) whose walls grow linearly instead of
// quadratically far from the minimum.
// =====================================================================

namespace {

struct Profile { double h, dh; };

inline Profile profile(double d, bool soft) {
    if (!soft) return { d * d, 2 * d };
    const double s = 1.0 / std::sqrt(4 + d * d);
    return { d * d * s, d * (8 + d * d) * s * s * s };
}

}  // namespace

// =====================================================================
// Energy and gradient: plain sums over the typed term lists.  The
// geometric coordinates and their gradients are the library primitives
// coord3d::angle/dangle and coord3d::dihedral/ddihedral, which
// implement exactly the Fortran ANGLE/DANGLE/DIHEDRAL/DDIHEDRAL
// conventions (angle at the middle atom, signed-atan2 dihedral, both
// with the first atom translated to the origin).  The corner and
// dihedral terms therefore read identically: value + gradient of a
// geometric coordinate, its deviation profiled, scattered onto the
// term's atoms.
// =====================================================================

double ForceField::energy(std::span<const coord3d> x) const {
    double E = 0;

    for (const Bond& t : bonds) {
        const auto [a, b] = t.atoms;
        E += t.k * profile((x[a] - x[b]).norm() - t.q0, soft).h;
    }
    for (const Corner& t : corners) {
        const auto [a, b, c] = t.atoms;
        const double theta = coord3d::angle(x[a] - x[b], x[c] - x[b]);
        E += t.k * profile(theta - t.q0, soft).h;
    }
    for (const Dihedral& t : dihedrals) {
        const auto [a, b, c, d] = t.atoms;
        const double phi =
            coord3d::dihedral(x[b] - x[a], x[c] - x[a], x[d] - x[a]);
        E += t.k * profile(phi - t.q0, soft).h;
    }
    if (k_coulomb != 0)
        for (const coord3d& xu : x) E += k_coulomb / xu.norm();

    return 0.5 * E;
}

double ForceField::energy_gradient(std::span<const coord3d> x,
                                   std::span<coord3d>       g) const {
    std::fill(g.begin(), g.end(), coord3d());
    double E = 0;

    for (const Bond& t : bonds) {
        const auto [a, b] = t.atoms;
        const coord3d u = x[a] - x[b];
        const double  r = u.norm();
        const auto [h, dh] = profile(r - t.q0, soft);
        E += t.k * h;
        const coord3d f = u * (0.5 * t.k * dh / r);
        g[a] += f;
        g[b] -= f;
    }
    for (const Corner& t : corners) {
        const auto [a, b, c] = t.atoms;
        const coord3d A = x[a] - x[b], C = x[c] - x[b];
        const double theta = coord3d::angle(A, C);
        coord3d da, dc;
        coord3d::dangle(A, C, da, dc);
        const auto [h, dh] = profile(theta - t.q0, soft);
        E += t.k * h;
        const double w = 0.5 * t.k * dh;
        g[a] += da * w;
        g[c] += dc * w;
        g[b] -= (da + dc) * w;
    }
    for (const Dihedral& t : dihedrals) {
        const auto [a, b, c, d] = t.atoms;
        const coord3d B = x[b] - x[a], C = x[c] - x[a], D = x[d] - x[a];
        const double phi = coord3d::dihedral(B, C, D);
        coord3d db, dc, dd;
        coord3d::ddihedral(B, C, D, db, dc, dd);
        const auto [h, dh] = profile(phi - t.q0, soft);
        E += t.k * h;
        const double w = 0.5 * t.k * dh;
        g[b] += db * w;
        g[c] += dc * w;
        g[d] += dd * w;
        g[a] -= (db + dc + dd) * w;
    }
    if (k_coulomb != 0)
        for (size_t u = 0; u < x.size(); ++u) {
            const double rinv = 1.0 / x[u].norm();
            E += k_coulomb * rinv;
            g[u] -= x[u] * (0.5 * k_coulomb * rinv * rinv * rinv);
        }

    return 0.5 * E;
}

// =====================================================================
// Term classification.
// =====================================================================

ForceField forcefield(const FullereneGraphView& G, int variant,
                      const Parameters* params) {
    if (variant < 1 || variant > 6)
        throw std::invalid_argument("wu::forcefield: variant "
                                    + std::to_string(variant) + " not in 1..6");
    const bool extended = variant >= 3;
    const Parameters P =
        params ? *params : (extended ? Parameters::extwu() : Parameters::wu());

    // Internal units: N/m -> kJ/mol/A^2 (N_A * 1e-3), degrees -> rad.
    constexpr double kJ_mol = 6.02214129;
    const auto rad = [](double deg) { return deg * M_PI / 180.0; };

    // Class tables indexed by pentagon count.
    const double r0[3] = { P.R66, P.R56, P.R55 };
    const double kr[3] = { P.fR66 * kJ_mol, P.fR56 * kJ_mol, P.fR55 * kJ_mol };
    const double a0[2] = { rad(P.A6), rad(P.A5) };            // [face_size == 5]
    const double ka[2] = { P.fA6 * kJ_mol, P.fA5 * kJ_mol };
    const double d0[4] = { rad(P.D666), rad(P.D566), rad(P.D556), rad(P.D555) };
    const double kd[4] = { P.fD666 * kJ_mol, P.fD566 * kJ_mol,
                           P.fD556 * kJ_mol, P.fD555 * kJ_mol };

    ForceField FF;
    FF.soft      = (variant == 5 || variant == 6);
    FF.k_coulomb = (variant % 2 == 0) ? P.fCoulomb : 0;

    // Bonds: one term per undirected edge, keyed by the number of
    // pentagons among the edge's two adjacent faces.
    FF.bonds.reserve(3 * G.N / 2);
    for (node_t u = 0; u < G.N; ++u)
        for (node_t v : G.nbrs(u)) {
            if (v < u) continue;
            const int np = (G.face_size(u, v) == 5) + (G.face_size(v, u) == 5);
            FF.bonds.push_back({{u, v}, r0[np], kr[np]});
        }

    // Corners: each face of size m contributes its m interior angles
    // (the angle sits at the middle atom), keyed by face size.
    FF.corners.reserve(3 * G.N);
    for (const face_t& f : G.compute_faces(6)) {
        const int m = (int)f.size(), pent = (m == 5);
        for (int j = 0; j < m; ++j)
            FF.corners.push_back({{f[j], f[(j + 1) % m], f[(j + 2) % m]},
                                  a0[pent], ka[pent]});
    }

    // Dihedrals (extended field only): one per vertex.  With CCW
    // neighbours (r, s, t) the three incident faces are the wedges
    // A = s-u-r, B = t-u-s, C = r-u-t (each traced from its leading
    // vertex's arc into u).  The quadruple is (u, w1, w2, w3) with
    // (w1, w2, w3) the CCW rotation of (r, s, t) starting at the
    // LEADING neighbour of the odd-sized face (s for A, t for B, r
    // for C); the uniform classes 555/666 keep the stored order.
    //
    // Note this is the Fortran SA_get_dihedrals' EFFECTIVE behaviour,
    // which differs from its comment diagram: get_face(a,u,b) removes
    // the arcs a-u and u-b and walks the cycle leaving u by the third
    // neighbour, so each of its three wedge queries reports the NEXT
    // wedge cyclically.  The published zero values (D555/D556/D566)
    // are fitted to this quadruple choice, so we reproduce it --
    // verified by per-class energy agreement at Fortran-converged
    // geometries (delaunay-fillin/tools/diag_wu_fortran).
    if (extended) {
        FF.dihedrals.reserve(G.N);
        for (node_t u = 0; u < G.N; ++u) {
            const auto nu = G.nbrs(u);
            const node_t r = nu[0], s = nu[1], t = nu[2];
            const int lA = G.face_size(s, u);
            const int lB = G.face_size(t, u);
            const int lC = G.face_size(r, u);
            const int npent = (lA == 5) + (lB == 5) + (lC == 5);

            std::array<node_t, 4> quad;
            if      (lA == lB && lB == lC) quad = {u, r, s, t};   // 555 / 666
            else if (lB == lC)             quad = {u, s, t, r};   // A odd
            else if (lA == lC)             quad = {u, t, r, s};   // B odd
            else                           quad = {u, r, s, t};   // C odd
            FF.dihedrals.push_back({quad, d0[npent], kd[npent]});
        }
    }
    return FF;
}

// =====================================================================
// Minimization driver.
// =====================================================================

static_assert(sizeof(coord3d) == 3 * sizeof(double),
              "wu::optimize reinterprets coord3d[] as double[]");

minimize::Outcome optimize(const ForceField& FF, std::span<coord3d> x,
                           double ftol) {
    const size_t n = x.size();
    std::span<double> flat(reinterpret_cast<double*>(x.data()), 3 * n);

    // Minimize field F over the shared flat span with options o.
    const auto run = [&](const ForceField& F, const minimize::Options& o) {
        return minimize::lbfgs(
            [&](std::span<const double> xf, std::span<double> gf) {
                return F.energy_gradient(
                    { reinterpret_cast<const coord3d*>(xf.data()), n },
                    { reinterpret_cast<coord3d*>(gf.data()), n });
            }, flat, o);
    };

    minimize::Options opt;
    opt.ftol_rel = ftol;
    // Per-iteration step cap: no atom moves more than half a bond
    // length per L-BFGS step.  Inactive near the minimum (converged
    // steps are far smaller); from a pathological start (e.g. the
    // legacy sphere-wrapped zero-order geometry) it forces gradual
    // untangling instead of a jump into a tangled local minimum.
    opt.max_step = 0.7;

    if (FF.k_coulomb == 0) return run(FF, opt);

    // Coulomb-carrying variant: the Fortran schedule (SA_frprmn3d)
    // inflates the start geometry with the origin repulsion until
    // ||g||_2 <= 10, then switches it off and converges the pure field.
    minimize::Options pre = opt;
    pre.gtol_2 = 10.0;
    run(FF, pre);
    ForceField F0 = FF;
    F0.k_coulomb = 0;
    return run(F0, opt);
}

}  // namespace wu
