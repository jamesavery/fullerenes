#pragma once

#include "fullerenes/planargraph.hh"

// Free functions for 2D layout operations on planar graphs.
// These were previously member functions of PlanarGraph that depended on a stored layout2d member.
// Now they take the layout as an explicit parameter, separating combinatorial graph structure from 2D embedding.

namespace layout2d {

  // ---------------------------------------------------------------------------
  // BOUNDARY-ONLY.  The two functions below WRITE the rotation system: they
  // reorder every neighbour row.  Orientation is an invariant held from
  // construction to destruction (parent CLAUDE.md), so a graph that is already
  // oriented has nothing for them to establish and one that lost its
  // orientation internally is a BUG -- repairing it there hides the bug and
  // silently swaps in the mirror embedding, which no combinatorial test can
  // tell from the original (Steinitz/Whitney: a 3-connected planar graph has
  // exactly two genus-0 rotation systems).
  //
  // Internal code must therefore NEVER call these.  They are reachable from
  // exactly four places, each one an input format that carries NO orientation
  // by nature, each documented as a boundary at its definition:
  //
  //   1. Polyhedron::from_mol2            (polyhedron-io.cc)  -- bare edge list
  //   2. oriented_graph_from_adjacency    (graph_fortran.cc)  -- adjacency matrix
  //   3. orient_polyhedron_neighbours     (polyhedron.cc, file-static) -- the
  //      Polyhedron-from-coordinates construction, whose only callers are the
  //      Polyhedron(graph, points) constructors
  //   4. set_layout2d_                    (graph_fortran.cc)  -- the Fortran C ABI
  //      re-import of an externally computed embedding (VERIFIED and rolled back
  //      on failure; see the note there -- it is the one caller whose input may
  //      already be oriented)
  //
  // Regression coverage for all four is tests/orientation-test.cc.
  // ---------------------------------------------------------------------------

  // Establish a rotation system on a graph that has none: Tutte-embed G and sort
  // each neighbour row CCW by angle in that embedding.  BOUNDARY-ONLY (above).
  //
  // The result is a planar embedding of unspecified handedness -- CW vs CCW
  // needs coordinates, which this signature does not have.  A caller that has
  // them must fix the handedness itself (orient_polyhedron_neighbours does, by
  // signed volume).
  //
  // @post   result == G.is_consistently_oriented()  -- CHECK IT.  A false means
  //         G's rows were rewritten and are still not a planar embedding, so the
  //         graph is left worse than it arrived; every boundary above throws.
  bool planar_orient(GraphView& G);

  // Sort each neighbour row of G CCW by angle around the vertex in `layout`.
  // BOUNDARY-ONLY (above).  Unchecked: the caller must verify the result with
  // G.oriented_surface() / require_oriented_surface().
  //
  // The rotation system produced is a planar embedding iff `layout` is a
  // crossing-free straight-line drawing of G.  It is NOT a layout setter for an
  // arbitrary layout: given a drawing with crossings it writes rows that embed G
  // in a higher-genus surface, silently.
  //
  // @pre    laid_out: layout.size() == G.N
  // @pre    planar_drawing: `layout` draws G without edge crossings
  void orient_neighbours(GraphView& G, const vector<coord2d>& layout);
  face_t find_outer_face(const PlanarGraphView& G, const vector<coord2d>& layout);
  bool layout_is_crossingfree(const PlanarGraphView& G, const vector<coord2d>& layout);

  vector<coord2d> spherical_projection(const PlanarGraphView& G, const vector<coord2d>& layout);
  bool optimize_layout(PlanarGraphView& G, vector<coord2d>& layout,
                       double zv_dist=0.2, double k_dist=10.0, double k_angle=10.0, double k_area=10.0);

  vector<double> edge_lengths(const PlanarGraphView& G, const vector<coord2d>& layout);
  coord2d width_height(const vector<coord2d>& layout);
  void scale(vector<coord2d>& layout, const coord2d& x);
  void move(vector<coord2d>& layout, const coord2d& x);

  string to_latex(const PlanarGraphView& G, const vector<coord2d>& layout,
                  double w_cm = 10, double h_cm = 10, bool show_dual = false, bool number_vertices = false, bool include_latex_header = false,
                  int edge_colour = 0x6a5acd, int path_colour = 0x6a5acd, int vertex_colour = 0x8b2500,
                  double edge_width = 0.1, double path_width = 0.1, double vertex_diameter = 2.0,
                  int Npath = 0, int *path = 0);

  string to_povray(const PlanarGraphView& G, const vector<coord2d>& layout,
                   double w_cm = 10, double h_cm = 10,
                   int line_colour = 0x6a5acd, int vertex_colour = 0x8b2500,
                   double line_width = 0.1, double vertex_diameter = 2.0);
}
