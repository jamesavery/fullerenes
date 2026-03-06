#ifndef DELAUNAY_HH
#define DELAUNAY_HH

#include "triangulation.hh"

// Intrinsic Delaunay triangulation of a fullerene dual.
//
// Given a fullerene dual triangulation (all equilateral triangles with unit edge
// length), computes the intrinsic Delaunay triangulation whose vertices are only
// the cone points (degree != 6 vertices). This is the Delaunay triangulation on
// the piecewise-flat surface manifold, using only the combinatorial structure and
// the implicit equilateral metric -- no 3D embedding required.
//
// Algorithm: incrementally remove flat (degree-6) vertices from the back,
// re-triangulating each hole via fan triangulation + local Delaunay flipping.
// A final global Delaunay sweep ensures correctness.
//
// References:
//   Fisher, Springborn, Schroder, Bobenko. "An Algorithm for the Construction
//   of Intrinsic Delaunay Triangulations." SIGGRAPH 2006.
//   Bobenko, Springborn. "A discrete Laplace-Beltrami operator for simplicial
//   surfaces." 2005.

class FulleroidDelaunay: public Triangulation {
public:
  matrix<double> edge_lengths; // N_alloc x N_alloc; edge_lengths(u,v) = 0 if no edge u-v

  // Construct from a fullerene dual triangulation.
  // Sorts vertices so non-degree-6 come first, degree-6 last.
  // Initializes all edge lengths to 1.
  FulleroidDelaunay(const Triangulation& T);

  // --- Edge length access ---
  double get_length(node_t u, node_t v) const { return edge_lengths(u,v); }
  void   set_length(node_t u, node_t v, double len) {
    edge_lengths(u,v) = len;
    edge_lengths(v,u) = len;
  }

  // --- Intrinsic geometry from edge lengths ---

  // Cotangent of the angle at vertex vk in the triangle (vi, vj, vk),
  // opposite edge vi-vj. Uses law of cosines: cot = cos/sin.
  double cot_opposite(node_t vi, node_t vj, node_t vk) const;

  // --- Delaunay operations ---

  // Check if edge (u,v) satisfies the local Delaunay criterion:
  // cot(angle_B) + cot(angle_D) >= 0 where B,D are the vertices
  // opposite (u,v) in its two adjacent triangles.
  bool is_delaunay_edge(node_t u, node_t v) const;

  // Compute the length of the other diagonal in the diamond around edge (u,v).
  // Lays out the 4 diamond vertices in 2D and computes the Euclidean distance.
  double flipped_edge_length(node_t u, node_t v) const;

  // Flip edge (u,v): remove it and add the other diagonal (B,D).
  // Returns true if flip was performed.
  // Returns false if flip would create a self-loop, multi-edge, or invalid length.
  bool flip_edge(node_t u, node_t v);

  // Flip all non-Delaunay edges until the triangulation is Delaunay.
  // Uses Lawson's algorithm (stack-based edge flipping).
  // Returns the total number of flips performed.
  int flip_to_delaunay();

  // Check if ALL edges satisfy the Delaunay criterion.
  bool is_delaunay() const;

  // --- Vertex removal ---

  // Lay out the fan of triangles around vertex v in 2D coordinates.
  // Places v at the origin; neighbours[v][0] along the positive x-axis.
  // Returns 2D positions for each neighbor of v (in CCW order).
  vector<coord2d> layout_fan(node_t v) const;

  // Remove flat vertex v: remove it and all its edges, fan-triangulate
  // the resulting hole, then locally Delaunayify.
  void remove_flat_vertex(node_t v);

  // Remove all flat (originally degree-6) vertices.
  // Main entry point for computing the reduced intrinsic Delaunay triangulation.
  void remove_flat_vertices();

  // --- Validation ---
  bool edge_lengths_are_symmetric() const;
};

#endif
