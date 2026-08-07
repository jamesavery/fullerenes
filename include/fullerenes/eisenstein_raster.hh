#pragma once

// =====================================================================
// Rasterisation primitives over the Eisenstein lattice Z[ω].
//
// walk_line(T, ...)     -- raster a directed integer line through a
//                          CCW-oriented equilateral triangulation,
//                          unfolding each crossed face into the
//                          line's lattice frame as it goes.
// scan_triangle(P0..P2) -- per-scanline a-range of the integer
//                          lattice points in a CCW triangle.
//
// Both depend only on the Eisenstein primitives in eisenstein.hh
// and on Triangulation's CCW neighbour-list contract.  No Alexandrov,
// no fillin-specific machinery.
// =====================================================================

#include "fullerenes/eisenstein.hh"
#include "fullerenes/triangulation.hh"

#include <algorithm>
#include <array>
#include <limits>
#include <span>
#include <utility>
#include <vector>

// =====================================================================
// Line walk through a CCW-oriented equilateral triangulation.
// =====================================================================

// The walk_line step budget, a safeguard: a primitive walk to
// displacement s = |a|+|b| crosses <= 3s faces, so valid inputs never
// approach it; exceeding it is reported through final_vertex = -1 (a
// modelled outcome, not an abort).  THE one spelling of this bound --
// the walker loop and the scratch-capacity formulas derived from it
// (eisenstein_paint.cc) read this name.
inline constexpr int walk_max_steps = 4096;

struct WalkResult {
    // Lattice points on the line, in order, with their T vertex IDs.
    // Always includes the starting vertex u at (0, 0).  For a line with
    // gcd = g, contains g+1 entries: the start, g-1 intermediate hex
    // vertices, and the end.
    std::vector<std::pair<Eisenstein, int>> walk;

    // T face triples (a, b, c) CCW traversed by the walk, in order.
    // Always includes the initial face (u, v, w).
    std::vector<std::array<int, 3>> face_path;

    // Lattice positions of each face's three vertices in the global
    // walk frame (start vertex u at (0, 0)).  Parallel to face_path.
    std::vector<std::array<Eisenstein, 3>> face_pos;

    // Vertex ID at the endpoint, or -1 if the walk failed (line does
    // not land at a face vertex within walk_max_steps, or an invalid exit
    // is encountered).
    int final_vertex = -1;

    // Final face's vertex IDs and lattice positions in the global
    // walk frame.  Used by callers that need to continue a walk past
    // the endpoint (e.g. through an intermediate hex vertex on a
    // gcd > 1 line).
    std::array<int, 3>        final_face_verts = {-1, -1, -1};
    std::array<Eisenstein, 3> final_face_pos   = { Eisenstein(0, 0),
                                                   Eisenstein(0, 0),
                                                   Eisenstein(0, 0) };
};

// Walk the integer line from (0, 0) to `endpoint` through T's faces,
// starting from face (u, v, w) with:
//   pos(u) = (0, 0)
//   pos(v) = unit_direction(dir_uv)
//   pos(w) = unit_direction(dir_uv + 1)
// The face (u, v, w) must be CCW-oriented in T (i.e.,
// w == T.next_on_face(u, v)).  Caller is responsible for this.
//
// Supports both gcd == 1 and gcd > 1 endpoints.  For gcd = g > 1, the
// line passes through g-1 intermediate lattice points (hex vertices on
// the geodesic); the walk returns g+1 entries in `walk`.  Internally
// the walk is g concatenated primitive (gcd == 1) sub-walks, with
// frame propagation at each intermediate hex vertex.
//
// @pre  For a correctly-developed (CCW, orientation-preserving) strip corridor,
//       `endpoint` should lie in the sector0 wedge (endpoint.first >= 0 and
//       endpoint.second >= 0).  The walk itself is orientation-preserving (it
//       develops every crossed face CCW), so the frame is never reflected; but
//       an endpoint outside sector0 -- especially with a negative component --
//       may make the marcher exit the wrong edge and reach the target via a
//       DIFFERENT same-length geodesic.  Callers that cannot guarantee the
//       sector must therefore VALIDATE the returned development, not trust it:
//       find_chains (eisenstein_paint.cc) probes the raw cell-frame displacement
//       and its sector0 rotation-representative over all start combinations and
//       KEEPS EVERY walk that lands the target exactly AND folds no non-corner
//       cone into the face; embed_cell then selects, per cell, the triple of
//       edge developments that agree wherever their strips share a lattice
//       position (the multi-edge / obtuse-cell cross-edge disambiguator).
WalkResult walk_line(const TriangulationView& T,
                     int u, int v, int w,
                     int dir_uv,
                     Eisenstein endpoint);

// =====================================================================
// Triangle scan-conversion over Z[ω].
// =====================================================================

struct ScanLine {
    int a_left  = 1;        // inclusive; valid range iff a_left <= a_right
    int a_right = 0;
    bool empty() const { return a_left > a_right; }
};

struct ScanLines {
    int b_min = 0;
    int b_max = -1;             // empty iff b_max < b_min
    std::vector<ScanLine> lines;   // size = b_max - b_min + 1
};

namespace raster_detail {
// Signed floor/ceiling division (round toward -inf / +inf) — the scan's
// exact chain arithmetic, ONE spelling for the owner and span bodies.
inline long floor_div(long n, long d) {
    long q = n / d, r = n % d;
    if (r != 0 && ((r < 0) != (d < 0))) --q;
    return q;
}
inline long ceil_div(long n, long d) {
    long q = n / d, r = n % d;
    if (r != 0 && ((r < 0) == (d < 0))) ++q;
    return q;
}
}  // namespace raster_detail

// The scan body over CALLER storage (device-legal: no allocation, no
// throw) — Abrash 39.2-style polygon scan of the CCW triangle
// (P0, P1, P2), filling lines[0 .. b_max-b_min] with the inclusive
// a-range per integer scanline (skip-scanlines keep a_left > a_right).
// Returns the scanline count; -1 for a degenerate/CW triangle (the
// owner's throw), -2 when `lines` is too small.  The owning
// scan_triangle below wraps this — ONE body.
inline long scan_triangle_into(Eisenstein P0, Eisenstein P1, Eisenstein P2,
                               std::span<ScanLine> lines,
                               int& b_min_out, int& b_max_out) {
    using raster_detail::ceil_div;
    using raster_detail::floor_div;
    if (wedge(P1 - P0, P2 - P0) <= 0) return -1;

    const int b_min = std::min({P0.second, P1.second, P2.second});
    const int b_max = std::max({P0.second, P1.second, P2.second});
    const long n = (long)b_max - b_min + 1;
    if ((long)lines.size() < n) return -2;
    for (long i = 0; i < n; ++i)
        lines[i] = ScanLine{std::numeric_limits<int>::max(),
                            std::numeric_limits<int>::min()};

    // For a CCW triangle, each non-flat edge in the cycle P0->P1->P2->P0
    // is on the LEFT chain iff its second endpoint has STRICTLY LOWER b
    // than its first endpoint (b decreasing in CCW direction).  The
    // other edges are on the RIGHT chain.
    const Eisenstein verts[3] = { P0, P1, P2 };
    for (int i = 0; i < 3; ++i) {
        const Eisenstein A = verts[i];
        const Eisenstein B = verts[(i + 1) % 3];
        if (A.second == B.second) continue;     // flat edge: skip
        const bool left_chain = (B.second < A.second);
        const long Aa = A.first, Ab = A.second;
        const long Ba = B.first, Bb = B.second;
        const long db = Bb - Ab;                // non-zero by the check above
        const long da = Ba - Aa;
        const int b_lo = (int)std::min(Ab, Bb);
        const int b_hi = (int)std::max(Ab, Bb);
        for (int b = b_lo; b <= b_hi; ++b) {
            // a_real = Aa + (b - Ab) * da / db
            const long num     = (long)(b - Ab) * da;
            const long a_floor = Aa + floor_div(num, db);
            const long a_ceil  = Aa + ceil_div(num, db);
            ScanLine& sl = lines[b - b_min];
            if (left_chain) {
                if ((int)a_ceil  < sl.a_left ) sl.a_left  = (int)a_ceil;
            } else {
                if ((int)a_floor > sl.a_right) sl.a_right = (int)a_floor;
            }
        }
    }
    b_min_out = b_min;
    b_max_out = b_max;
    return n;
}

// Run Abrash 39.2-style polygon scan over the CCW triangle (P0, P1, P2).
// Returns the inclusive a-range per integer scanline b in [b_min, b_max].
// Skip-scanlines (where the triangle's a-range at b spans less than one
// lattice unit) are marked by a_left > a_right (ScanLine::empty()).
//
// Preconditions: P0, P1, P2 non-collinear AND CCW
// (wedge(P1 - P0, P2 - P0) > 0).  Aborts otherwise.
// Owning wrapper of scan_triangle_into — one body.
ScanLines scan_triangle(Eisenstein P0, Eisenstein P1, Eisenstein P2);
