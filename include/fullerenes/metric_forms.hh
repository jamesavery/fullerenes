#pragma once
// ============================================================================
// metric_forms.hh -- the metric identities shared across coordinate rings.
//
// heron_product_sq: the Heron product in SQUARED-side coordinates,
//   H(x, y, z) = 2(xy + yz + zx) - (x^2 + y^2 + z^2) = 16 * Area^2,
// equivalently 4xy - (x + y - z)^2; negative iff the triangle inequality
// fails.  ONE spelling for every ring that has +, -, and *: long long
// (the Eisenstein exact regime, via eisenstein.hh's global overload, which
// delegates here), cyclotomic::Real30 (the flattened kis surface,
// cyclotomic.hh), cyclotomic::Zeta60 (the verification tier,
// cyclotomic_ambient.hh).
//
// The template lives in its own namespace so that INTEGER arguments always
// bind to eisenstein.hh's long long overload (promotion) instead of
// deducing R = int and silently narrowing -- the 2026-08-25 review
// exhibited exactly that hijack for a bare global template.  Ring callers
// use the qualified name.
//
// (delaunay_detail::heron_product is the SAME form in real-length
// coordinates, factored as (a+b+c)(-a+b+c)(a-b+c)(a+b-c) and clamped at
// the triangle-inequality boundary for floating-point conditioning -- a
// deliberately separate evaluation, not a second spelling of this one.)
// ============================================================================

namespace metric_forms {

template <class R>
constexpr R heron_product_sq(const R& x, const R& y, const R& z) {
  return 2 * (x * y + y * z + z * x) - (x * x + y * y + z * z);
}

}  // namespace metric_forms
