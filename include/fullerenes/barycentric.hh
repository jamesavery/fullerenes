#pragma once

// =====================================================================
// Integer barycentric coordinates over a CCW triangle in Z[w].
//
// For a lattice point p inside (or on the boundary of) a CCW triangle
// (P0, P1, P2):
//   IntBary { n0, n1, n2, denom } where denom = wedge(P1-P0, P2-P0)
//                                  and n_i are signed wedges.
// p is inside (or on boundary) iff n0, n1, n2 are all >= 0.
//
// For an ON-EDGE point (one of n_i is zero), reduce_to_lowest_terms
// divides through by gcd(non-zero pair) so that adjacent triangles
// produce bit-identical 3D output across the shared edge (the
// face-wedge denominator drops out of the FP arithmetic).
//
// barycentric_combine evaluates v = (m0*C0 + m1*C1 + m2*C2) / denom
// using doubles, with the rational weights kept in integer form until
// the final scalar divide.
//
// These primitives are generic to "barycentric over a CCW triangle in
// Z[w]"; they aren't fillin- or iDT-specific.
// =====================================================================

#include "fullerenes/eisenstein.hh"
#include "fullerenes/geometry.hh"

struct IntBary { long n0, n1, n2, denom; };

IntBary integer_barycentric(Eisenstein p,
                            Eisenstein P0,
                            Eisenstein P1,
                            Eisenstein P2);

struct ReducedBary { long m0, m1, m2, denom; };

ReducedBary reduce_to_lowest_terms(IntBary b);

coord3d barycentric_combine(ReducedBary b,
                            const coord3d& C0,
                            const coord3d& C1,
                            const coord3d& C2);

// The same affine combination for an already-normalized double triple
// (b0 + b1 + b2 == 1): v = b0*C0 + b1*C1 + b2*C2.  The FP-barycentric
// counterpart used by the (non-Loeschian) cubic transport path, so both
// number systems name the one "interpolate a triangle" operation.
coord3d barycentric_combine(const double b[3],
                            const coord3d& C0,
                            const coord3d& C1,
                            const coord3d& C2);
