#pragma once

#include "fullerenes/planargraph.hh"

// Free functions for 2D layout operations on planar graphs.
// These were previously member functions of PlanarGraph that depended on a stored layout2d member.
// Now they take the layout as an explicit parameter, separating combinatorial graph structure from 2D embedding.

namespace layout2d {

  // Orient a planar graph using the Boyer-Myrvold planarity algorithm.
  // This works on any planar graph regardless of current neighbour ordering.
  // Returns true if the graph is planar, false otherwise.
  bool planar_orient(Graph& G);

  void orient_neighbours(Graph& G, const vector<coord2d>& layout);
  face_t find_outer_face(const PlanarGraph& G, const vector<coord2d>& layout);
  bool layout_is_crossingfree(const PlanarGraph& G, const vector<coord2d>& layout);

  vector<coord2d> spherical_projection(const PlanarGraph& G, const vector<coord2d>& layout);
  bool optimize_layout(PlanarGraph& G, vector<coord2d>& layout,
                       double zv_dist=0.2, double k_dist=10.0, double k_angle=10.0, double k_area=10.0);

  vector<double> edge_lengths(const PlanarGraph& G, const vector<coord2d>& layout);
  coord2d width_height(const vector<coord2d>& layout);
  void scale(vector<coord2d>& layout, const coord2d& x);
  void move(vector<coord2d>& layout, const coord2d& x);

  string to_latex(const PlanarGraph& G, const vector<coord2d>& layout,
                  double w_cm = 10, double h_cm = 10, bool show_dual = false, bool number_vertices = false, bool include_latex_header = false,
                  int edge_colour = 0x6a5acd, int path_colour = 0x6a5acd, int vertex_colour = 0x8b2500,
                  double edge_width = 0.1, double path_width = 0.1, double vertex_diameter = 2.0,
                  int Npath = 0, int *path = 0);

  string to_povray(const PlanarGraph& G, const vector<coord2d>& layout,
                   double w_cm = 10, double h_cm = 10,
                   int line_colour = 0x6a5acd, int vertex_colour = 0x8b2500,
                   double line_width = 0.1, double vertex_diameter = 2.0);
}
