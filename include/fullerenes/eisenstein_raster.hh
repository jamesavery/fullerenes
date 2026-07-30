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

#include <array>
#include <utility>
#include <vector>

// =====================================================================
// Line walk through a CCW-oriented equilateral triangulation.
// =====================================================================

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
    // not land at a face vertex within MAX_STEPS, or an invalid exit
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

// Run Abrash 39.2-style polygon scan over the CCW triangle (P0, P1, P2).
// Returns the inclusive a-range per integer scanline b in [b_min, b_max].
// Skip-scanlines (where the triangle's a-range at b spans less than one
// lattice unit) are marked by a_left > a_right (ScanLine::empty()).
//
// Preconditions: P0, P1, P2 non-collinear AND CCW
// (wedge(P1 - P0, P2 - P0) > 0).  Aborts otherwise.
ScanLines scan_triangle(Eisenstein P0, Eisenstein P1, Eisenstein P2);
