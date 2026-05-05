#pragma once

// =====================================================================
// Eisenstein-lattice operations.
//
// Free functions and small helper types that extend the `Eisenstein`
// class with the number-theoretic primitives needed by lattice-based
// geometry algorithms (rasterisation, integer barycentric
// interpolation, D6-action alignment between norm-equal Eisensteins).
//
// Conventions:
//   - (a, b) in Z^2 represents a + b*omega where omega = exp(i*pi/3).
//   - Cartesian coords of (a, b): (a + b/2, b*sqrt(3)/2).
//   - Lattice directions d in {0..5} are the 6 unit Eisenstein integers
//     at angles d*60 degrees CCW, starting at the +a axis.
//   - Wedge w(u, v) = u.a*v.b - u.b*v.a, integer signed area in shear
//     coords; sign matches the Cartesian wedge.
// =====================================================================

#include "fullerenes/eisenstein.hh"
#include "fullerenes/geometry.hh"

#include <utility>
#include <vector>

// The 6 Eisenstein unit vectors, CCW starting at the +a axis.
//   d=0: ( 1,  0)  angle   0
//   d=1: ( 0,  1)  angle  60
//   d=2: (-1,  1)  angle 120
//   d=3: (-1,  0)  angle 180
//   d=4: ( 0, -1)  angle 240
//   d=5: ( 1, -1)  angle 300
inline Eisenstein unit_direction(int d) {
    static const Eisenstein U[6] = {
        Eisenstein( 1,  0), Eisenstein( 0,  1), Eisenstein(-1,  1),
        Eisenstein(-1,  0), Eisenstein( 0, -1), Eisenstein( 1, -1)
    };
    return U[((d % 6) + 6) % 6];
}

// Integer signed area in shear coords.  Sign matches the Cartesian
// wedge (Cartesian-wedge = (sqrt(3)/2) * integer-wedge), so the value
// is exact and trigonometry-free.
inline long wedge(Eisenstein u, Eisenstein v) {
    return (long)u.first * v.second - (long)u.second * v.first;
}

// i-reflection: complex conjugation in C, restricted to Z[ω]:
//   complex_conj(a + b*ω) = a + b*ω̄ = (a + b) - b*ω
// Companion to Eisenstein::eis_conj() (the ω-reflection in Z[ω]'s
// native (1, ω) basis).  Forwards to Eisenstein::complex_conj() so
// callers can use either the free-function or method form.
inline Eisenstein complex_conj(Eisenstein z) {
    return z.complex_conj();
}

// Cartesian (x, y) of Eisenstein (a, b) in the (1, ω) basis.
inline std::pair<double, double> to_cartesian(Eisenstein z) {
    return { z.first + 0.5 * z.second, z.second * 0.8660254037844386 };
}

// Inverse of to_cartesian, rounded to the nearest Eisenstein integer:
//   b = round(2 * y / sqrt(3));   a = round(x - b/2)
Eisenstein from_cartesian(double x, double y);

// Some Eisenstein integer (a, b) with a >= 0, b >= 0 and
// a^2 + a*b + b^2 == N.  Returns the first sector-0 representative
// found by scanning b = 0, 1, ...; aborts if no solution exists.
//
// Precondition: N >= 0 is a valid Eisenstein norm.
Eisenstein eisenstein_of_norm(int N);

// Enumerate ALL sector-0 Eisenstein reps (a >= 0, b >= 0) of norm N.
// Generic norms return 1 entry; split-prime norms return 2 entries in
// distinct rotation orbits.  Used to enumerate the candidate placements
// of a lattice triangle's edge endpoint when the edge length squares
// to a split-prime norm.
std::vector<Eisenstein> sector0_reps_of_norm(int N);

// =====================================================================
// D6 action between norm-equal Eisensteins.
// =====================================================================

// A D6 affine transform of Z[ω] onto itself:
//   T(z) = unit * z                  (if reflect == false)
//        = unit * complex_conj(z)    (if reflect == true)
// `unit` is one of the 6 Eisenstein units.
struct D6Affine {
    Eisenstein unit;
    bool       reflect;

    Eisenstein apply(Eisenstein z) const {
        return (reflect ? complex_conj(z) : z) * unit;
    }
};

// Find the D6Affine T with T(z_from) = z_to.  Both inputs must have
// the same norm.  Exactly one of the two branches (rotation,
// reflection) gives a unit by the D6 symmetry of Z[ω]; align returns
// that branch.  Aborts on norm mismatch or non-divisibility.
D6Affine align(Eisenstein z_from, Eisenstein z_to);

// =====================================================================
// Integer barycentric over a CCW lattice triangle.
// =====================================================================

// Integer barycentric weights (n0, n1, n2) of lattice point p inside
// triangle (P0, P1, P2):
//   denom = wedge(P1-P0, P2-P0) = 2 * signed area
// p is inside (or on boundary) iff all three are >= 0 (assuming CCW).
struct IntBary { long n0, n1, n2, denom; };

IntBary integer_barycentric(Eisenstein p,
                            Eisenstein P0,
                            Eisenstein P1,
                            Eisenstein P2);

// Integer barycentric weights divided by their common factor.  For an
// ON-EDGE lattice point (one weight = 0) reduce by gcd(non-zero pair)
// so that both adjacent cells produce bit-identical 3D output across
// the shared edge (the cell-wedge denominator drops out).  For
// interior points this is a no-op.
struct ReducedBary { long m0, m1, m2, denom; };

ReducedBary reduce_to_lowest_terms(IntBary b);

// Affine combination of three points with reduced rational weights:
//   v = (m0*C0 + m1*C1 + m2*C2) / denom
coord3d barycentric_combine(ReducedBary b,
                            const coord3d& C0,
                            const coord3d& C1,
                            const coord3d& C2);
