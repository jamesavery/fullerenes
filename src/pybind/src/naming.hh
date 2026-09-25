#pragma once
//
// naming.hh — the canonical fullerene name as the Python bindings write it.

#include <string>

#include "fullerenes/triangulation.hh"   // FullereneDualView
#include "fullerenes/spiral.hh"

namespace pyf {

// Returns the canonical name of the fullerene whose dual triangulation is
// `dual`, prefixed with its carbon count 2*Nv - 4: "C<N>-" followed by
// FullereneDualView::name().to_string(), e.g. "C60-[1,7,9,...]-fullerene".
// The spiral is the canonical general spiral (pentagon starts); the name
// carries a search tag only when that spiral has jumps (spiral.hh, the name
// grammar).
inline std::string canonical_name(const FullereneDualView& dual) {
    return "C" + std::to_string(2 * (int)dual.N - 4) + "-"
         + dual.name(/*rarest_start=*/true).to_string();
}

}  // namespace pyf
