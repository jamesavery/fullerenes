#pragma once

#include "fullerenes/eisenstein.hh"
#include "fullerenes/triangulation.hh"
#include "fullerenes/delaunay.hh"

#include <iosfwd>
#include <vector>

// =====================================================================
// Per directed iDT arc: unfold the triangle strip the geodesic
// traverses into a canonical Eisenstein-lattice frame.
//
// For each live iDT arc h (half-edge), `unfold_arc_strip` computes the
// strip of T_sorted triangles the geodesic from u = he_origin[h] to
// v = dest(h) passes through, in the arc's canonical Eisenstein frame:
// u at (0, 0), v at v_position (sector-0 rep of |he_length[h]|^2).
//
// Arc indexing (rather than undirected-edge indexing) lets a face-
// placement consumer take each strip in the orientation of the cell's
// CCW boundary directly: for a cell's arc c_i -> c_j, the canonical
// frame already has c_i at (0, 0) and c_j at v_position, so the
// cell-frame affine is always a pure rotation (never a reflection).
// =====================================================================

// One vertex of the triangle strip in the arc's canonical frame.
struct StripVertex {
    Eisenstein position;        // canonical-frame lattice position
    int        vertex_id;       // T_sorted vertex
    int        frame_offset;    // k0 in the canonical frame, 0..5
};

// Per-arc strip data.
//
// `left` lists strip vertices on the LEFT of directed arc u -> v
// (i.e., on the side of the iDT face he_face[h]), plus the two cone
// endpoints, plus any intermediate on-line hex (gcd > 1 case).
// `right` is symmetric for the RIGHT of u -> v (the twin face's side).
// Endpoints and intermediate on-line hex appear in both lists.
//
// Invariant: left.size() == right.size() == v_position.first +
// v_position.second + 1 (sector-0 v_position with both components >= 0).
struct Strip {
    int        u = -1;             // = he_origin[h]
    int        v = -1;             // = dest(h)
    Eisenstein v_position{0, 0};   // canonical position of v (with u at origin)
    std::vector<StripVertex> left;
    std::vector<StripVertex> right;
    bool       ok = false;
};

// Unfold the strip for one iDT arc (half-edge).
//
// Preconditions:
//   - 0 <= arc_h < D.nh
//   - The half-edge is live (D.alive(arc_h))
//
// Aborts on alignment failure (contract violation between iDT and the
// walking-line model -- should not fire on valid iDT input).
Strip unfold_arc_strip(const DelaunayTriangulation& D,
                       const Triangulation& T_sorted,
                       int arc_h);

// Unfold strips for every live iDT arc in parallel.
// Returns vector<Strip> of size D.nh; dead slots have ok == false.
std::vector<Strip> unfold_all_arc_strips(const DelaunayTriangulation& D,
                                         const Triangulation& T_sorted);

// Diagnostic flag: when set, unfold_arc_strip prints a per-edge trace
// line to stderr.  Off by default.  Set via test driver CLI flags
// rather than environment variables.
void set_unfold_arc_strip_debug(bool on);
bool unfold_arc_strip_debug();

// Self-contained TikZ document showing the unfolded strip in the
// canonical Eisenstein frame: T_sorted edges between strip vertices,
// the directed iDT geodesic u -> v, and the strip vertices coloured
// by side (left / right / on-line) and shaped by type (square =
// endpoint cone, circle = interior).
void dump_strip_tikz(const Strip& S,
                     const Triangulation& T_sorted,
                     std::ostream& os,
                     double scale = 0.7);
