#pragma once

// =====================================================================
// dwu -- derived-geometry Wu force field for fullerene 3D structure.
//
// A successor to wu::forcefield (wu_forcefield.hh) differing from it in three
// ways, each of which was measured against a corpus of 507 GFN2-xTB reference
// geometries rather than assumed:
//
//   1. REST VALUES ARE DERIVED, NOT TABULATED.  extwu carries nine independent
//      rest values (three bond lengths, two corner angles, four dihedrals).
//      Only the bond lengths are independent: a corner angle is the interior
//      angle of the polygon a face's own edges make, and a vertex dihedral is
//      the fold three faces with those edges actually have.  Both follow from
//      the bond lengths by exact construction -- cyclic_polygon_angles and
//      coord3d::ideal_dihedral -- so this field fits TWO lengths and derives
//      the other seven quantities.
//
//      This is not a simplification for its own sake.  Fitting the rest values
//      freely reproduces the reference geometry marginally better while
//      producing a field whose parameters contradict each other: D566 = 0.18
//      deg against bond lengths implying 23.84, with the geometry rescued by
//      21 deg of cancelling strain in the dihedral term.  The derived field
//      cannot do that.
//
//   2. A FACE-FLATNESS TERM.  E = 1/2 k lambda_min(S_f), with S_f the centred
//      face scatter matrix, so lambda_min is exactly the summed squared
//      distance from the least-squares face plane.  This is the single most
//      valuable addition: it takes held-out RMSD from 0.078 to 0.028 A for two
//      parameters, a factor of 2.8.
//
//   3. NO 666 DIHEDRAL TERM.  extwu penalises the fold at a vertex of three hexagons
//      against a rest value of 0.  The rest value is right -- three regular hexagons tile
//      the plane, so 0 is what the construction gives -- but the TERM earns nothing, and
//      measurably so: its force constant sat on an arbitrary lower bound of 5 N/m in every
//      fit, and freeing that bound sent it to 0.005 N/m while held-out error IMPROVED by
//      4 %.  The face-flatness term already carries the bending physics the dihedral was
//      meant to, and on a better coordinate: lambda_min(S_f) is how far a face is from its
//      own plane, which is what bending a graphene sheet costs, while a vertex dihedral
//      measures one fold of three faces and is blind to how each face bends internally.
//
//      The Gaussian-curvature term, by contrast, STAYS -- and its status is honestly
//      provisional rather than settled.  On the IPR corpus it is neutral: +0.4 % without the
//      disclination field, -0.0 % with it, and nothing above 1.4 % in any size or
//      pentagon-distance band.  But an IPR corpus is the regime where the term has least to
//      say, because its target (pi/15)*n_pent takes only two values there -- 12 deg at a 566
//      vertex and 0 at a 666 one.  Non-IPR cages add 556 and 555 vertices at 24 and 36 deg,
//      four times the range, and that fit is where it will be decided.  Note also that fK
//      fits to 2-4 N/m rather than to zero: a useless term goes to zero, as the raw-angle
//      version did (0.011 N/m), so this one is EXCHANGEABLE with flatness rather than inert,
//      and which of the two is the better-conditioned description is not a question
//      geometry RMSD can answer.
//
//   4. A DISCLINATION FIELD ON THE 66 BOND.  A pentagon is a disclination, and
//      the strain it imposes decays with p, the topological distance to the
//      nearest pentagon (BFS from all pentagon vertices).  Measured, the 66
//      bond runs 1.370 A at p = 0 through 1.438 at p = 2 and decays to 1.430;
//      the rise over the first two shells is the classical bond alternation and
//      the decay beyond is the disclination's strain field.  So R66 takes free
//      values at p = 0 and 1 and q_inf + delta*exp(-p/xi) beyond, and the fitted
//      decay length reproduces to 1.5 % across independent training folds.
//
// So the field is 13 fitted constants: two bond rest lengths (one of them carrying the
// disclination field), five force constants, two flatness constants, and the curvature
// stiffness.  Everything else is derived or dropped.  Held-out RMSD 0.0254 A against the
// published field's 0.235 -- a factor of 9 -- on 507 IPR cages, C80-C180.
//
// PERFORMANCE.  The energy and gradient are wu::ForceField's, unchanged: this
// header only computes the rest values and adds the flatness term, so a dwu
// relaxation costs what a wu one costs plus one 3x3 eigen-decomposition per
// face per evaluation.
//
// STATUS.  The functional form is settled and validated; the CONSTANTS in
// Parameters::fitted() are provisional -- fitted on 507 IPR cages, C80-C180,
// against GFN2-xTB.  A larger fit (2,542 cages, C60-C300, including 1,320
// non-IPR) is in progress and will change the numbers, not the structure.  The
// R55 entry is a placeholder: an isolated-pentagon corpus contains no 55 bond,
// so it is unidentified, and non-IPR cages are what will fix it.
//
// Header-only deliberately, so it is usable against an already-built
// libfullerenes; it can move to src/c++ when the library is next rebuilt.
// =====================================================================

#include "fullerenes/wu_forcefield.hh"
#include "fullerenes/fullerenegraph.hh"
#include "fullerenes/geometry.hh"
#include "fullerenes/minimize.hh"

#include <algorithm>
#include <array>
#include <cmath>
#include <map>
#include <span>
#include <stdexcept>
#include <vector>

namespace dwu {

inline constexpr double kJ_mol = 6.02214129;   // N/m -> kJ/mol/A^2, wu_forcefield.cc:138

// --- the disclination field -------------------------------------------------
// A rest length as a function of p, the hop distance to the nearest pentagon:
// free values at the listed shells, an exponential decay to q_inf beyond.
//
// Why a free head rather than a decay everywhere.  The measured 66 bond RISES
// over the first shells -- 1.370 A at p = 0 through 1.394 to 1.438 at p = 2 --
// and only then decays.  p in {0, 1, 2} are the chemically distinct cases (both
// endpoints pentagon-adjacent, one, neither), i.e. the classical bond
// alternation, which no single monotone decay can represent.  So those shells
// get their own values and the decay describes the disclination's strain field
// beyond them.
//
// delta = 0 with no listed shells is a plain constant, which is how the classes
// with no p-dependence are expressed -- and for R55, R56 and A5 that is
// structural, not a modelling choice: every instance of those classes is
// pentagon-incident, so p takes ONE value across the class and a decay in it has
// nothing to vary over.
struct Field {
    std::vector<int>    shell;    // p values carrying their own value
    std::vector<double> value;    // parallel to shell
    double q_inf = 0, delta = 0, xi = 1;
    int    p0    = 0;             // the decay applies from this p upward

    // @pre shell.size() == value.size()
    // @pre the listed shells cover 0..p0-1 (see the clamp below)
    double at(int p) const {
        for (std::size_t i = 0; i < shell.size(); ++i)
            if (shell[i] == p) return value[i];
        // A p below p0 that no shell lists is a MISCONFIGURED field -- the head is
        // meant to cover every shell before the decay starts.  Clamping to the
        // first listed value keeps it continuous and bounded rather than
        // evaluating a decay outside its stated domain, where exp(-p/xi) > 1 would
        // overshoot q_inf + delta.  Stated because it is a fallback, not a feature:
        // with p0 = 2 and shells {0, 1}, as fitted(), it is unreachable.
        if (p < p0 && !value.empty()) return value.front();
        return q_inf + delta * std::exp(-double(p) / std::max(xi, 1e-9));
    }
    static Field constant(double q) { Field f; f.q_inf = q; return f; }
};

struct Parameters {
    Field  R55, R56, R66;                 // bond rest lengths, by pentagon count
    double fR55, fR56, fR66;              // bond force constants, N/m
    double fA5, fA6;                      // corner force constants, N/m
    double fD555, fD556, fD566, fD666;    // dihedral force constants, N/m
    double fFlat5, fFlat6;                // face flatness, N/m
    double fK;                            // vertex Gaussian curvature, N/m

    // PROVISIONAL constants: the 13-parameter fit of 2026-09-18, 507 IPR cages C80-C180
    // against GFN2-xTB, trained on 406 and sealed-tested on 101 -- held-out RMSD 0.02538 A
    // (train 0.02423), worst per-atom deviation 0.0542 A, nothing on a bound.  The published
    // extended-Wu field scores 0.235 A on the same references.
    //
    // Three groups are NOT fitted, for structural reasons rather than practical ones:
    //
    //   fR66      the GAUGE.  Scaling every force constant leaves the minimiser exactly
    //             where it was, so geometry determines only their ratios and one must be
    //             pinned.  Pinned at extwu's 450 N/m, which fixes the units of the rest.
    //             This is the one arbitrary number here, and it is arbitrary because
    //             geometry cannot see it.  A GFN2-xTB Hessian can: matching curvature rather
    //             than positions measures the scale.  Done for one cage, the answer spans
    //             3.9x depending on which modes are weighted (fR66 163 to 634 N/m), because
    //             this field's stiffness DISTRIBUTION is wrong vibrationally -- its highest
    //             mode is 2660 cm^-1 against xTB's 1521.  Fitting Hessians alongside
    //             geometry would settle the ratios and the scale together; until then 450
    //             stands as a convention, not a measurement.
    //   fD666     zero, and not a free parameter: the term is dropped (see 3 above).
    //   R55 fR55  an isolated-pentagon corpus has no two pentagons sharing an edge or a
    //   fD555     vertex, so these classes have no instances and keep their extwu values.
    //   fD556     A non-IPR fit moves them a long way -- fR55 260 -> 1063 N/m, fD555
    //             35 -> 555, fD556 65 -> 735 -- and moves the fitted ones with them
    //             (fR56 1752 -> 338), so these constants are NOT transferable to cages with
    //             adjacent pentagons.  The full-corpus fit supersedes them.
    //   R55 shell R55 and R56 are pentagon-incident by construction, so p takes one value
    //   R56 shell across each class and a decay in it is unidentifiable.
    //
    // Quoted to the precision at which they matter; claude-projects/forcefield-fit/tools/
    // dwu_check re-scores this header against the corpus the fit used, and should be run
    // after any change to these numbers.
    static Parameters fitted() {
        Parameters P;
        P.R55 = Field::constant(1.479);                  // extwu; no IPR instances
        P.R56 = Field::constant(1.42790);
        P.R66 = Field{{0, 1}, {1.38942, 1.40808}, 1.40138, 0.05703, 3.31101, 2};
        P.fR55 = 260.0;    P.fR56 = 1278.2;  P.fR66 = 450.0;   // fR66: the gauge
        P.fA5  = 365.48;   P.fA6  = 42.650;
        P.fD555 = 35.0;    P.fD556 = 65.0;                     // extwu; no IPR instances
        P.fD566 = 7.1019;  P.fD666 = 0.0;                      // fD666: term dropped
        P.fFlat5 = 168.44; P.fFlat6 = 80.369;
        P.fK     = 3.8594;                                     // centroid-fan curvature
        return P;
    }
};

// --- derived rest geometry --------------------------------------------------

// Interior angles (degrees) of the cyclic polygon with these sides: the ideal
// corner angles a face's own edge lengths imply.  Vertices on one circle of
// radius R, so side i subtends 2 asin(s_i/2R) and closure is
// sum_i asin(s_i/2R) = pi, monotone in R, hence unconditional bisection from
// R = max(s)/2.  Equal sides give the regular polygon exactly (108 deg at m=5,
// 120 at m=6), so this generalises the Wu convention rather than competing with
// it: it differs only at a face whose edges differ, i.e. a hexagon bordering a
// pentagon.  Result i is the corner between sides i-1 and i.
// @throws std::domain_error when the sides form no cyclic polygon -- checked,
//         because the bisection's bracket fails once the circumcentre leaves the
//         polygon and the clamped arcsines then return plausible angles with no
//         other signal.
inline std::vector<double> cyclic_polygon_angles(const std::vector<double>& sides) {
    const std::size_t m = sides.size();
    const double smax = *std::max_element(sides.begin(), sides.end());
    auto closure = [&](double R) {
        double t = 0;
        for (double s : sides) t += std::asin(std::min(1.0, s / (2 * R)));
        return t - M_PI;
    };
    double lo = smax / 2, hi = smax;
    while (closure(hi) > 0 && hi < 1e6 * smax) hi *= 2;
    for (int it = 0; it < 200; ++it) {
        const double mid = 0.5 * (lo + hi);
        (closure(mid) > 0 ? lo : hi) = mid;
    }
    const double R = 0.5 * (lo + hi);
    std::vector<double> theta(m), ang(m);
    for (std::size_t i = 0; i < m; ++i) theta[i] = 2 * std::asin(std::min(1.0, sides[i] / (2 * R)));
    for (std::size_t i = 0; i < m; ++i)
        ang[i] = (M_PI - 0.5 * (theta[(i + m - 1) % m] + theta[i])) * 180.0 / M_PI;
    double sum = 0;
    for (double a : ang) sum += a;
    if (std::abs(sum - 180.0 * (double(m) - 2)) > 1e-6)
        throw std::domain_error("dwu::cyclic_polygon_angles: sides form no cyclic polygon");
    return ang;
}

// Hop distance from every vertex to the nearest pentagon vertex.  Purely
// topological, so a rest value depending on it stays predictable from the graph
// alone -- no geometry and no reference data enter.
inline std::vector<int> pentagon_distance(const FullereneGraph& G,
                                          const std::vector<face_t>& faces) {
    std::vector<int> d(G.N, -1);
    std::vector<node_t> q;
    for (const face_t& f : faces)
        if (f.size() == 5)
            for (node_t v : f) if (d[v] != 0) { d[v] = 0; q.push_back(v); }
    for (std::size_t i = 0; i < q.size(); ++i)
        for (node_t w : G.nbrs(q[i])) if (d[w] < 0) { d[w] = d[q[i]] + 1; q.push_back(w); }
    return d;
}

// --- the curvature term -----------------------------------------------------
// E = 1/2 k (K_v - K_v*)^2 per vertex, with K_v = 2*pi - (angles around v) and the target
// K_v* = (pi/15)*n_pent(v): 36 deg at a vertex of three pentagons, 24 at two, 12 at one, and
// ZERO at a vertex of three hexagons.  The target is topological, not fitted -- each pentagon
// has five vertices, so it sums to 12*5*pi/15 = 4*pi, which is Gauss-Bonnet exactly -- so the
// term localises curvature onto the pentagons rather than setting its magnitude.
//
// The angles must come from a TRIANGULATION.  An angle defect is an intrinsic curvature only
// on a surface built of flat pieces, and a relaxed fullerene's faces are not flat: summing
// raw POLYGON angles charges the atoms for face non-planarity as though it were curvature,
// whose per-atom sum is 841.8 deg per cage on these references.  That is NOT a violation of
// Gauss-Bonnet: the theorem constrains a surface of FLAT faces, and a fullerene's faces are not
// planar, so the raw per-atom sum defines no surface and is not a curvature.  So
// every face is joined to its own centroid and the fan triangles supply the angles; then
// atoms plus centroids sum to 720.000 deg to 1e-12 on every cage measured.
//
// Triangulating through the face centroids DOES give a polyhedral surface, and there the theorem
// holds exactly.  The centroid share is not zero: the atoms carry 760.9
// deg and the centroids -40.9, so a target summing to 4*pi over the ATOMS presumes flat
// faces.  Holding the atoms to 4*pi therefore also pushes the faces flat, which is why this
// term and the flatness term are partly exchangeable (see the note at the top).
struct CurvatureField {
    std::vector<face_t> faces;      // oriented
    std::vector<double> target;     // per vertex, radians
    double k = 0;                   // internal kJ/mol/rad^2; 0 disables

    double energy_gradient(std::span<const coord3d> x, std::span<coord3d> g) const {
        if (k == 0) return 0;
        const std::size_t N = target.size();
        std::vector<coord3d> cen(faces.size(), coord3d(0, 0, 0));
        for (std::size_t f = 0; f < faces.size(); ++f) {
            for (node_t v : faces[f]) cen[f] += x[v];
            cen[f] /= double(faces[f].size());
        }
        std::vector<double> around(N, 0.0);
        for (std::size_t f = 0; f < faces.size(); ++f) {
            const face_t& fc = faces[f];
            const int m = int(fc.size());
            for (int j = 0; j < m; ++j) {
                const node_t p = fc[(j + m - 1) % m], v = fc[j], q = fc[(j + 1) % m];
                around[v] += coord3d::angle(x[p] - x[v], cen[f] - x[v])
                           + coord3d::angle(cen[f] - x[v], x[q] - x[v]);
            }
        }
        double E = 0;
        std::vector<double> w(N);              // dE/d(angle at v) = -k * residual
        for (std::size_t v = 0; v < N; ++v) {
            const double d = (2.0 * M_PI - around[v]) - target[v];
            E   += 0.5 * k * d * d;
            w[v] = -k * d;
        }
        // Every fan angle depends on the CENTROID as well as on its two rim atoms, and the
        // centroid on all m vertices of the face -- so a centroid derivative is spread over
        // the whole face with weight 1/m.
        for (std::size_t f = 0; f < faces.size(); ++f) {
            const face_t& fc = faces[f];
            const int m = int(fc.size());
            const double inv_m = 1.0 / m;
            for (int j = 0; j < m; ++j) {
                const node_t p = fc[(j + m - 1) % m], v = fc[j], q = fc[(j + 1) % m];
                coord3d d1, d2;
                coord3d::dangle(x[p] - x[v], cen[f] - x[v], d1, d2);
                g[p] += d1 * w[v];
                g[v] -= (d1 + d2) * w[v];
                for (node_t u : fc) g[u] += d2 * (w[v] * inv_m);
                coord3d::dangle(cen[f] - x[v], x[q] - x[v], d1, d2);
                for (node_t u : fc) g[u] += d1 * (w[v] * inv_m);
                g[q] += d2 * w[v];
                g[v] -= (d1 + d2) * w[v];
            }
        }
        return E;
    }
};

// --- the flatness term ------------------------------------------------------
// E = 1/2 k lambda_min(S_f); dE/dx_i = k ((x_i - c).n) n, with n the least
// eigenvector.  lambda_min carries A^2, so k converts like a bond constant and
// needs no angular treatment.
struct FlatnessField {
    std::vector<face_t> faces;
    std::vector<double> k;

    double energy_gradient(std::span<const coord3d> x, std::span<coord3d> g) const {
        double E = 0;
        for (std::size_t fi = 0; fi < faces.size(); ++fi) {
            if (k[fi] == 0) continue;
            const face_t& f = faces[fi];
            coord3d c(0, 0, 0);
            for (node_t u : f) c += x[u];
            c /= double(f.size());
            matrix3d A;
            for (node_t u : f) {
                const coord3d q = x[u] - c;
                for (int i = 0; i < 3; ++i)
                    for (int j = 0; j < 3; ++j) A(i, j) += q[i] * q[j];
            }
            // eigensystem() orders by |lambda| smallest first and returns an
            // orthonormal basis, so entry 0 is lambda_min with a usable
            // eigenvector even where it is degenerate -- the case eigenvector()
            // alone can only signal by returning the zero vector.
            const auto [lam, C] = A.eigensystem();
            const coord3d n(C(0, 0), C(0, 1), C(0, 2));
            E += 0.5 * k[fi] * std::max(0.0, lam[0]);
            for (node_t u : f) g[u] += n * (k[fi] * (x[u] - c).dot(n));
        }
        return E;
    }
};

// --- the field --------------------------------------------------------------

struct ForceField {
    wu::ForceField  wu;
    CurvatureField  curv;
    FlatnessField   flat;

    // @anchor dwu-energy-gradient
    // @pre g.size() == x.size()
    double energy_gradient(std::span<const coord3d> x, std::span<coord3d> g) const {
        // Sequenced deliberately.  wu::ForceField::energy_gradient ZERO-FILLS g before
        // accumulating, so it must run first -- and the operands of `+` are UNSEQUENCED
        // in C++, so writing the sum as one expression would leave an evaluation order
        // that erases the flatness gradient legal, and silently wrong.
        const double E_wu   = wu.energy_gradient(x, g);
        const double E_curv = curv.energy_gradient(x, g);
        const double E_flat = flat.energy_gradient(x, g);
        return E_wu + E_curv + E_flat;
    }
    double energy(std::span<const coord3d> x) const {
        std::vector<coord3d> g(x.size());
        return energy_gradient(x, {g.data(), g.size()});
    }
};

// Build the field of G under P.
//
// The term lists come from wu::forcefield, so the dihedral quadruple convention
// is the library's and not a re-derivation of it; only the rest values are
// replaced.  Each is computed from the bonds AT THE TERM'S OWN GEOMETRY: a
// face's angles from that face's own side lengths, a vertex's dihedral from its
// own three incident edges, each edge at its own p.
//
// The dihedral argument order is load-bearing.  ideal_dihedral's figure puts
// face A between edges ur and us, B between us and ut, C between ut and ur,
// while wu_forcefield.cc rotates its quadruple to start at the leading neighbour
// of the ODD face -- the two agree exactly when the odd face sits in slot C, so
// every mixed class is posed that way.  The quadruple is not symmetric under every
// rotation of the three faces: posing a 566 vertex as (5,6,6) reads 19.77 deg where
// the correct (6,6,5) gives 23.94.  (Permuting only the edge LENGTHS, faces held
// fixed, moves it by 0.04 deg -- what matters here is the posing of the faces.)  ideal_dihedral returns RADIANS, already the unit wu::forcefield converted
// its table into, so the value is assigned unscaled.
//
// Evaluated at extwu's OWN bond lengths -- each vertex class taking the edges it
// actually has, (R55,R55,R55) at a 555 vertex, (R56,R55,R56) at 556, (R56,R66,R56)
// at 566 -- the derivation reproduces three of the four published constants to the
// precision they are printed at: 37.3774 against 37.38, 29.2021 against 29.20, and
// 0 against 0.  The fourth comes out 23.9430 where the table prints 23.49, so the
// published D566 is a transposition of the digits of its own construction.
//
// @anchor dwu-forcefield
// @pre fullerene: G.is_a_fullerene()
// @post result.wu.bonds.size() == 3*G.N/2 && result.wu.corners.size() == 3*G.N
inline ForceField forcefield(const FullereneGraph& G, const Parameters& P) {
    ForceField M;
    // Seeded from the published table so a field this model does not set is
    // extwu's value rather than zero; every rest value is overwritten below.
    wu::Parameters k = wu::Parameters::extwu();
    k.fR55 = P.fR55; k.fR56 = P.fR56; k.fR66 = P.fR66;
    k.fA5  = P.fA5;  k.fA6  = P.fA6;
    k.fD555 = P.fD555; k.fD556 = P.fD556; k.fD566 = P.fD566; k.fD666 = P.fD666;
    k.fCoulomb = 0;
    M.wu = wu::forcefield(G, /*variant=*/3, &k);

    const std::vector<face_t> faces = G.compute_faces_oriented(6);
    const std::vector<int>    pd    = pentagon_distance(G, faces);
    const double D2R = M_PI / 180.0;
    const Field* bond[3] = {&P.R66, &P.R56, &P.R55};       // by pentagon count on the edge
    auto edge_rest = [&](node_t u, node_t v) {
        const int np = (G.face_size(u, v) == 5) + (G.face_size(v, u) == 5);
        return bond[np]->at(pd[u] + pd[v]);
    };

    for (wu::Bond& b : M.wu.bonds) b.q0 = edge_rest(b.atoms[0], b.atoms[1]);

    std::map<std::pair<node_t, node_t>, int> arcface;
    for (std::size_t fi = 0; fi < faces.size(); ++fi)
        for (std::size_t j = 0; j < faces[fi].size(); ++j)
            arcface[{faces[fi][j], faces[fi][(j + 1) % faces[fi].size()]}] = int(fi);

    std::vector<std::vector<double>> face_ang(faces.size());
    for (std::size_t fi = 0; fi < faces.size(); ++fi) {
        const face_t& f = faces[fi];
        std::vector<double> sides(f.size());
        for (std::size_t j = 0; j < f.size(); ++j)
            sides[j] = edge_rest(f[j], f[(j + 1) % f.size()]);
        face_ang[fi] = cyclic_polygon_angles(sides);
    }
    for (wu::Corner& c : M.wu.corners) {
        const auto [a, b, cc] = c.atoms;
        auto i1 = arcface.find({a, b}), i2 = arcface.find({b, cc});
        int fi = (i1 != arcface.end() && i2 != arcface.end() && i1->second == i2->second)
                 ? i1->second : (arcface.count({cc, b}) ? arcface[{cc, b}] : -1);
        int pos = -1;
        if (fi >= 0)
            for (std::size_t j = 0; j < faces[fi].size(); ++j)
                if (faces[fi][j] == b) { pos = int(j); break; }
        c.q0 = (fi >= 0 && pos >= 0) ? face_ang[fi][pos] * D2R
                                     : (G.face_size(a, b) == 5 ? 108.0 : 120.0) * D2R;
    }

    for (wu::Dihedral& d : M.wu.dihedrals) {
        const node_t u = d.atoms[0];
        const auto nu = G.nbrs(u);
        const node_t r = nu[0], s = nu[1], t = nu[2];
        const int lA = G.face_size(s, u), lB = G.face_size(t, u), lC = G.face_size(r, u);
        std::array<int, 3>    F{lA, lB, lC};
        std::array<node_t, 3> E{r, s, t};
        if (!(lA == lB && lB == lC)) {                 // rotate the odd face into slot C
            if      (lB == lC) { F = {lB, lC, lA}; E = {s, t, r}; }
            else if (lA == lC) { F = {lC, lA, lB}; E = {t, r, s}; }
        }
        d.q0 = coord3d::ideal_dihedral(F[0], F[1], F[2], edge_rest(u, E[0]),
                                       edge_rest(u, E[1]), edge_rest(u, E[2]));
    }

    M.flat.faces = faces;
    for (const face_t& f : faces)
        M.flat.k.push_back((f.size() == 5 ? P.fFlat5 : P.fFlat6) * kJ_mol);

    if (P.fK != 0) {
        M.curv.k     = P.fK * kJ_mol;
        M.curv.faces = faces;
        M.curv.target.assign(G.N, 0.0);
        for (node_t u = 0; u < G.N; ++u) {
            int npent = 0;
            for (node_t v : G.nbrs(u)) npent += (G.face_size(u, v) == 5);
            M.curv.target[u] = (M_PI / 15.0) * npent;
        }
    }
    return M;
}

// Minimise FF over x in place, stopping on a SCALE-FREE gradient criterion:
// ||g||_inf <= gtol_rel * ||g_0||_inf.  An absolute threshold is meaningless
// while force constants vary over orders of magnitude -- both sides here scale
// identically under any rescaling of energy or length, so the criterion, and the
// geometry it accepts, are invariant.
// @anchor dwu-optimize
// @pre x.size() == number of graph vertices FF's terms index
inline minimize::Outcome optimize(const ForceField& FF, std::span<coord3d> x,
                                  double gtol_rel = 1e-6) {
    const std::size_t n = x.size();
    std::vector<coord3d> g0(n);
    FF.energy_gradient(x, {g0.data(), n});
    double gs = 0;
    for (const coord3d& g : g0) for (int i = 0; i < 3; ++i) gs = std::max(gs, std::fabs(g[i]));

    minimize::Options opt;
    opt.max_step = 0.7;
    opt.ftol_rel = 1e-14;                  // backstop; the gradient decides
    opt.gtol_inf = gtol_rel * gs;
    std::span<double> flat(reinterpret_cast<double*>(x.data()), 3 * n);
    return minimize::lbfgs(
        [&](std::span<const double> xf, std::span<double> gf) {
            return FF.energy_gradient({reinterpret_cast<const coord3d*>(xf.data()), n},
                                      {reinterpret_cast<coord3d*>(gf.data()), n});
        }, flat, opt);
}

}  // namespace dwu
