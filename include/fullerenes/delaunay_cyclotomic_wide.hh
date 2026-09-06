#pragma once
// ============================================================================
// delaunay_cyclotomic_wide.hh -- the 64/128 tier of the kis metric policy:
// 128-bit arithmetic (cyclotomic_wide.hh) over the 64-bit storage of the
// carry, for single huge isomers beyond the 64-bit envelope (about
// C15,000; through N of order 10^9).  HOST ONLY -- the library's own
// single-isomer owner path (DelaunayTriangulation::
// remove_and_complete_cyclotomic_kis) runs on this tier, since one isomer
// at a time costs nothing extra and needs no width decision.
// ============================================================================

#include "cyclotomic_wide.hh"
#include "delaunay_cyclotomic.hh"

namespace cyclotomic {

using CyclotomicMetric64 = CyclotomicMetricT<Real30, Real30Wide>;
using CarryStore64 = CarryStoreT<Real30, Real30Wide>;
using CyclotomicKisCarry64 = CyclotomicKisCarryT<Real30, Real30Wide>;

inline CyclotomicKisCarry64 derive_cyclotomic_kis_carry_wide(
    const DelaunayView& V, int n_centres, std::span<const int> centre_size,
    const char* op = "derive_cyclotomic_kis_carry") {
  return derive_cyclotomic_kis_carry<Real30, Real30Wide>(V, n_centres, centre_size, op);
}

}  // namespace cyclotomic
