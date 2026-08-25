#pragma once
// ============================================================================
// diamond_forms.hh -- DiamondForms<R>: the diamond's shared arithmetic
// skeleton over any ring with +, -.
//
// The diamond (delaunay_geometry.hh carries the picture):
//
//      B          upper triangle: sides (e, a, b), apex B
//     a b         a adjacent to u, b adjacent to v
//  u---e---v      e = the diagonal being tested / flipped
//     c d         lower triangle: sides (e, c, d), apex D
//      D          c adjacent to u, d adjacent to v
//
// Every exact diamond classifier -- DiamondSq (Eisenstein squared lengths,
// delaunay_geometry.hh), cyclotomic::Diamond (the flattened kis surface's
// wedge carry, cyclotomic.hh), cyclotomic::Diamond60 (the conductor-60
// cross-check, cyclotomic_ambient.hh), and the layer-2 CyclotomicMetric
// policy to come -- shares this skeleton: the five side fields, the four
// law-of-cosines numerators, and the reversal involution.  The verdict
// composition and the area accessors (tau vs the Heron form vs the carried
// wedge) genuinely differ per ring and stay with the derived classifiers.
//
// sigma, the REVERSAL INVOLUTION: the diamond read from the diagonal's
// other endpoint v, sigma(e,a,b,c,d) = (e,b,a,d,c); sigma(sigma(D)) = D.
// Whole-diamond quantities (the Delaunay form, the faces and their areas)
// are sigma-invariant; per-endpoint quantities (P, Q, endpoint convexity)
// transform by sigma -- the second endpoint's test IS the first one's on
// the reversed diamond.  Stated once, here.
//
// Fields are the SQUARED side lengths in each regime's own convention
// (DiamondSq: Loeschian integers; the cyclotomic diamonds: Z[gamma] at the
// x kLsqScale scaling).  The double-precision banded Diamond
// (delaunay_geometry.hh) carries real lengths and band constants instead
// and deliberately does not derive from this skeleton.
// ============================================================================

template <class R>
struct DiamondForms {
  R e, a, b, c, d;

  // The law-of-cosines numerators (squared coordinates): s_* are the
  // cotangent numerators at the two apexes, P/Q the endpoint-convexity
  // numerators at the origin endpoint u.
  R s_upper() const { return a + b - e; }
  R s_lower() const { return c + d - e; }
  R P() const { return e + a - b; }
  R Q() const { return e + c - d; }

  DiamondForms reversed() const { return {e, b, a, d, c}; }
};
