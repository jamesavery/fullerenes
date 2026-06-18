#pragma once
//
// naming.hh — the project-wide canonical-name convention in one place.

#include <string>

#include "fullerenes/triangulation.hh"   // FullereneDualView
#include "fullerenes/spiral.hh"

namespace pyf {

// Canonical generalized-spiral name of a fullerene dual triangulation, e.g.
// "C60-[GS:1,7,9,...]-fullerene". The single source of the C-prefix + GS-spiral
// convention (CLAUDE.md mandates it everywhere). Carbon count = 2*Nv - 4.
inline std::string canonical_name(const FullereneDualView& dual) {
    return "C" + std::to_string(2 * (int)dual.N - 4) + "-"
         + dual.name(/*rarest_start=*/true).to_string();
}

}  // namespace pyf
