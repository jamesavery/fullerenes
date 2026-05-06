#pragma once

// =====================================================================
// TikZ helpers for visualising Eisenstein-lattice geometry.
//
// Header-only.  Each consumer TU gets its own inlined copy; the code
// is small and only called from cold (diagnostic) paths, so the
// multi-instantiation cost is negligible compared to the benefit of
// having one source of truth.
// =====================================================================

#include "fullerenes/eisenstein.hh"

#include <algorithm>
#include <cmath>
#include <ostream>
#include <utility>

namespace tikz {

constexpr double SQ3_2 = 0.8660254037844386;     // sin(60deg) = sqrt(3)/2

// Eisenstein -> Cartesian.  (a + b w) maps to (a + b/2, b*sqrt(3)/2)
// with w = exp(i*pi/3).
inline std::pair<double, double> cart(Eisenstein z) {
    return { z.first + 0.5 * z.second, z.second * SQ3_2 };
}

// Cartesian bounding box accumulated by bump(); span() is the larger
// of width / height (used to scale figure margins).
struct BBox {
    double xmin = 1e9, xmax = -1e9;
    double ymin = 1e9, ymax = -1e9;

    void bump(Eisenstein p) {
        const auto [x, y] = cart(p);
        if (x < xmin) xmin = x;  if (x > xmax) xmax = x;
        if (y < ymin) ymin = y;  if (y > ymax) ymax = y;
    }
    double span() const {
        return std::max(xmax - xmin, ymax - ymin);
    }
};

// Faint Eisenstein-lattice background filling the box.  For every
// lattice point in (and just outside) the box, draws the 3 outgoing
// unit edges (east, NE, NW) so every lattice edge whose midpoint is
// inside the box appears exactly once.
inline void emit_grid(std::ostream& os, BBox b) {
    const int j_lo = (int)std::floor(b.ymin / SQ3_2) - 1;
    const int j_hi = (int)std::ceil (b.ymax / SQ3_2) + 1;
    static constexpr int DA[3] = {  1, 0, -1 };
    static constexpr int DB[3] = {  0, 1,  1 };
    for (int j = j_lo; j <= j_hi; ++j) {
        const double y = j * SQ3_2;
        const int a_lo = (int)std::floor(b.xmin - 0.5 * j) - 1;
        const int a_hi = (int)std::ceil (b.xmax - 0.5 * j) + 1;
        for (int a = a_lo; a <= a_hi; ++a) {
            const double x = a + 0.5 * j;
            for (int d = 0; d < 3; ++d) {
                const double xn = x  + DA[d] + 0.5 * DB[d];
                const double yn = y  + DB[d] * SQ3_2;
                const double mx = 0.5 * (x + xn);
                const double my = 0.5 * (y + yn);
                if (mx < b.xmin || mx > b.xmax ||
                    my < b.ymin || my > b.ymax) continue;
                os << "  \\draw[black!30,line width=0.35pt] ("
                   << x  << "," << y  << ") -- ("
                   << xn << "," << yn << ");\n";
            }
        }
    }
}

// Directed outline of the lattice triangle P0 -> P1 -> P2 -> P0 in red.
inline void emit_lattice_triangle(std::ostream& os,
                                   Eisenstein P0, Eisenstein P1, Eisenstein P2)
{
    const auto [x0, y0] = cart(P0);
    const auto [x1, y1] = cart(P1);
    const auto [x2, y2] = cart(P2);
    os << "  \\draw[->,red!80!black,line width=1pt] ("
       << x0 << "," << y0 << ") -- (" << x1 << "," << y1 << ");\n"
       << "  \\draw[->,red!80!black,line width=1pt] ("
       << x1 << "," << y1 << ") -- (" << x2 << "," << y2 << ");\n"
       << "  \\draw[->,red!80!black,line width=1pt] ("
       << x2 << "," << y2 << ") -- (" << x0 << "," << y0 << ");\n";
}

}  // namespace tikz
