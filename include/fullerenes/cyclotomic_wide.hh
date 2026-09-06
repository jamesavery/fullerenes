#pragma once
// ============================================================================
// cyclotomic_wide.hh -- the 64/128 tier of the cyclotomic ring: 128-bit
// arithmetic coefficients over a 64-bit storage of the carry, for single
// isomers beyond the 64-bit envelope (about C15,000; cyclotomic.hh's
// WIDTHS paragraph).  HOST ONLY: this is the one header of the module
// that spells a 128-bit integer type; the batch tier's headers never do.
//
// Everything is cyclotomic.hh's templates instantiated at Coef =
// __int128 -- the ring, its guard (a 128-bit overflow-checked magnitude
// product), the exact division (a 128-bit numerator's residue is its two
// 64-bit words folded), the point ring and the diamond -- plus the sign
// oracle's parameters at this width: the Liouville bound needs
// p >= 4 * 127 + 17 = 525 bits, the module takes 576 (18 words), with the
// fixed-point gamma derived by the same exact bisection as the 64-bit
// one and verified through psi at compile time.
//
// The envelope at this width, carry_coeff_max<__int128>(), is
// 5.4 * 10^17: the (5,0) nanotube's largest coordinate, 0.56 N^2, stays
// inside it through N of order 10^9.  A stored coordinate is 64-bit, so
// the storage binds first, at N ~ 4 * 10^9 -- neither is reached.
// ============================================================================

#include "cyclotomic.hh"

namespace cyclotomic {

using WideCoef = __int128;
using Real30Wide = Real30T<WideCoef>;
using Zeta30Wide = Zeta30T<Real30Wide>;
using DiamondWide = DiamondT<Real30Wide>;

inline constexpr WideCoef kCarryCoeffMaxWide = carry_coeff_max<WideCoef>();
static_assert(kCarryCoeffMaxWide > (WideCoef)500'000'000 * 1'000'000'000,
              "the 128-bit envelope exceeds 5 * 10^17");
static_assert(envelope_holds<WideCoef>(), "the 128-bit envelope's binding terms");

namespace detail {

template <>
struct SignParams<WideCoef> {
  static constexpr int kBits = 576;                  // p
  static constexpr int kConstLimbs = 55;             // m^k 2^{p(3-k)} < 2^{3p+3} = 2^1731
  static constexpr int kAccLimbs = 59;               // + 128-bit coefficients, four terms: < 2^1861
  static constexpr int kVerifyLimbs = 73;            // psi at scale 2^{4p}: < 2^2309
  // floor(gamma 2^576), little-endian 64-bit words
  // (claude-projects/delaunay/tools/derive_gamma_fixed.py 576).
  static constexpr std::array<uint64_t, 10> kGammaWords = {
      0x90995d93a2f56f52ULL, 0x958bef41838f42deULL, 0x3d71909d6767f5a6ULL,
      0xce1cc5d3311c3213ULL, 0x45effbeef3b55e15ULL, 0x10c1517b2ab04fe0ULL,
      0x5af6b5d8ea29e11eULL, 0xf5be43c6e4340270ULL, 0xf4cfc327a007f8a9ULL,
      0x0000000000000001ULL};
};

static_assert(sign_params_sound<WideCoef>(), "sign oracle, 128-bit: parameters");
static_assert(gamma_fixed_isolates<WideCoef>(),
              "sign oracle, 128-bit: psi(m / 2^p) < 0 < psi((m + 1) / 2^p) must hold");

}  // namespace detail

}  // namespace cyclotomic
