#pragma once
// seam_step.hh -- fold Step 3 primitives:
//
//   1. EisFrame -- the affine lattice isometry a track or holonomy walk carries,
//   2. SeamAtlas -- the outline's chart-changing structure: mate pairing,
//      glue frames, corner material sectors, exit-crossing detection,
//   3. cone_holonomy -- the CCW neighbour-ring walk around a cone, chaining
//      through the outline glues.
//
// The Z[w] direction algebra (dot2, same_dir, primitive_dir, first_ccw,
// first_cw) and the half-open reflex-safe Sector membership this builds on
// live in eisenstein.hh.  Algorithm reference: doc/ALGORITHM-seam-step.md.
// Everything is exact integer arithmetic; failures are diagnosed by throw,
// never patched over.

#include "fullerenes/triangulation.hh"
#include "fullerenes/eisenstein.hh"
#include <map>
#include <numeric>
#include <stdexcept>
#include <vector>


// ---------------------------------------------------------------------------
// EisFrame: the affine lattice isometry  p |-> rot*p + t  (rot a unit) carried by
// a walk; composed with a glue at every seam crossing.
// ---------------------------------------------------------------------------
struct EisFrame {
  Eisenstein rot{1, 0}, t{0, 0};

  Eisenstein operator()(Eisenstein p) const { return rot * p + t; }
  Eisenstein inv(Eisenstein q) const { return rot.complex_conj() * (q - t); }  // valid since |rot| = 1
  EisFrame      after(const EisFrame& f) const { return {rot * f.rot, (*this)(f.t)}; }  // this composed after f
  bool       is_isometry() const { return rot.norm2() == 1; }
};

// ---------------------------------------------------------------------------
// SeamAtlas: the outline as a chart-changing structure.
//
// The outline is the CW (or CCW -- decided from the signed area, never
// assumed) polygon of cone-labelled corners; each directed seam segment
// appears exactly once with its label-reversed mate elsewhere (the seam is a
// tree).  The atlas derives, per DESIGN.md section 3, the mate pairing and
// the glue isometries, plus the two geometric predicates every walk needs.
// ---------------------------------------------------------------------------
class SeamAtlas {
public:
  // Throws if any segment lacks a unique label-reversed mate (a non-tree seam).
  explicit SeamAtlas(const std::vector<std::pair<Eisenstein, node_t>>& outline);

  int        n() const { return (int)O.size(); }
  Eisenstein pos(int i) const { return O[wrap(i)].first; }
  node_t     lab(int i) const { return O[wrap(i)].second; }

  // The Eisenstein::turn sign of a direction pointing into the material.
  int interior_turn() const { return int_turn; }

  // Mate segment index of edge i, and the glue mapping edge i's lip frame
  // onto the mate's: pos(i) |-> pos(mate+1), pos(i+1) |-> pos(mate).
  int   mate(int i) const;
  EisFrame glue(int i) const;

  // Is this development position an outline corner (i.e. a cone copy)?
  bool is_corner(Eisenstein p) const { return corner_ix.count(p) != 0; }

  // Unit direction u points into the material at corner m.  Half-open, so
  // each boundary lattice cell is owned by exactly one corner copy.  The
  // ownership side is the traversal convention: CCW walks own (incoming,
  // outgoing] (the default), CW walks the mirrored [incoming, outgoing) --
  // each excludes the far copy of a just-crossed seam arc, which for a
  // LEAF cone (single copy, e.g. every non-source cone of a star unfolding)
  // sits exactly on the other boundary ray of the same corner.
  bool material(int m, Eisenstein u, bool ccw = true) const;

  // The first outline edge the segment G -> Gn crosses from inside-or-on to
  // strictly outside, at exact rational parameter num/den STRICTLY beyond
  // (min_num/min_den) -- the caller walks a segment through successive charts
  // by gluing at each exit and passing the last crossing's parameter back in
  // (parameters are preserved by the glues).  min_num < 0 marks the segment's
  // true start and admits parameter-0 exits (a mid-edge boundary departure).
  //
  // Corners: at the true start, departure from an outline corner is governed
  // by that corner's material sector -- material means an interior departure
  // (the incident edges are not crossing candidates), non-material means no
  // unique glue exists and kind == corner is reported (e.g. a track trying to
  // continue straight through a cone).  A segment passing exactly through any
  // other corner also reports kind == corner: the caller must stop, never
  // guess a glue.
  struct Exit {
    enum Kind { none, edge, corner } kind;
    int index;                                  // the crossed edge / the corner's outline index
    long num = 0, den = 1;                      // crossing parameter along G -> Gn (kind == edge)
  };
  Exit first_exit(Eisenstein G, Eisenstein Gn, long min_num = -1, long min_den = 1) const;

private:
  int wrap(int i) const { const int m = (int)O.size(); return ((i % m) + m) % m; }

  std::vector<std::pair<Eisenstein, node_t>> O;
  int int_turn;
  std::map<std::pair<node_t, node_t>, int> mate_ix;
  std::map<Eisenstein, int> corner_ix;          // corner position -> outline index
};

// ---------------------------------------------------------------------------
// cone_holonomy: assemble the CCW neighbour ring of the cone at outline
// corner m0, starting from the anchor direction (which must be material
// there).  Rotating CCW exits through the OUTGOING edge of the current copy;
// each exit records the wedge-boundary ray and composes the glue.  The walk
// runs until it closes back onto (m0, anchor); non-closure, an undeveloped
// ring cell, or a runaway glue chain throw.
// ---------------------------------------------------------------------------
struct ConeRing {
  // The unit arc cone -> neighbour may itself cross NON-incident seams (a thin
  // neck or spike of the outline within unit distance of the cone) -- it is
  // then one of the seam-crossing arcs of DESIGN.md section 4.  The direction
  // is material at the copy, but the cell is resolved by following the unit
  // segment's glue chain: F is that chain's composed isometry (identity when
  // the arc is developed in-lip) and cell = F(pos(copy) + dir).
  struct Entry {
    int        copy;   // outline corner index owning the arc's direction
    Eisenstein dir;    // global unit direction from that copy to the neighbour
    EisFrame      F;      // glue chain of the unit segment (identity if in-lip)
    Eisenstein cell;   // resolved development cell = F(pos(copy) + dir)
    node_t     id;     // the neighbour (grid lookup at cell)
  };
  // A wedge-boundary seam ray (for the figures): sits CCW after ring entry
  // [after]; dir is the outgoing-edge direction in the frame current at the
  // event, w_acc the rotation accumulated since entry [after] was emitted.
  struct BoundaryRay {
    int        after;
    Eisenstein w_acc;
    Eisenstein dir;
  };
  std::vector<Entry>       ring;
  std::vector<BoundaryRay> rays;
};

ConeRing cone_holonomy(const SeamAtlas& atlas, const std::map<Eisenstein, node_t>& grid,
                       int m0, Eisenstein anchor);

// The chart walk of one straight UNIT arc from the development position P0
// (a cone corner, an on-seam copy, or an interior cell) in direction `dir`:
// glue at each exit crossing under the strictly-increasing parameter
// protocol, land, look up.  A unit arc contains no interior lattice point,
// hence no interior cone, so the walk stays within flat charts regardless of
// any cones near the seams -- this is THE primitive from which fold completes
// every arc that is not developed in-sheet.  F is the composed glue chain
// (identity if in-sheet), cell = F(P0 + dir).  Throws on an undeveloped
// landing or a corner-degenerate exit; both are impossible for real arcs (a
// unit edge's interior contains no vertex), so they are invariants.
struct ArcLanding { EisFrame F; Eisenstein cell; node_t id; };
ArcLanding resolve_arc(const SeamAtlas& atlas, const std::map<Eisenstein, node_t>& grid,
                       Eisenstein P0, Eisenstein dir, const std::string& who);

// The six unit arcs of the flat vertex developed at p, in p's own chart:
// result[j] = the neighbour in direction unit[j].  For an on-seam copy the
// on-line directions walk the boundary and the far side glues across, so the
// row is coherent in that one copy's frame.
std::array<node_t, 6> unit_ring(const SeamAtlas& atlas, const std::map<Eisenstein, node_t>& grid,
                                Eisenstein p);

