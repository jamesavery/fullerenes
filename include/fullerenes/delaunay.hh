#ifndef DELAUNAY_HH
#define DELAUNAY_HH

#include "triangulation.hh"
#include "geometry.hh"

// Diamond: the local geometry around an edge in a metrized triangulation.
//
//      B
//     / \          upper triangle: sides (e, a, b)
//    a   b         a = side adjacent to u, b = side adjacent to v
//   /     \  .
//  u---e---v       e = diagonal being tested/flipped
//   \     /
//    c   d         lower triangle: sides (e, c, d)
//     \ /          c = side adjacent to u, d = side adjacent to v
//      D
//
// All geometric predicates (Delaunay, convexity, flipped length) depend only
// on these five edge lengths.  No vertex IDs or topology needed.
struct Diamond {
  double e, a, b, c, d;

  bool   is_delaunay()    const;  // cot(angle_B) + cot(angle_D) >= 0
  bool   is_convex()      const;  // angle sum < pi at both u and v
  double flipped_length() const;  // length of BD (the other diagonal)
};

// Intrinsic Delaunay triangulation of an equilateral triangulation.
//
// Given any equilateral triangulation (all edges unit length) on a closed
// orientable surface, computes the intrinsic Delaunay triangulation whose
// vertices are only the cone points (degree != 6 vertices).
//
// Algorithm: incrementally remove flat (degree-6) vertices by reducing their
// degree to 3 via local edge flips, then trivially removing them.  Each edge
// flip computes the new edge length from the diamond geometry.
// Lawson flipping after each removal restores the Delaunay property.
//
// References:
//   Fisher, Springborn, Schroder, Bobenko. "An Algorithm for the Construction
//   of Intrinsic Delaunay Triangulations." SIGGRAPH 2006.
//   Bobenko, Springborn. "A discrete Laplace-Beltrami operator for simplicial
//   surfaces." 2005.

class FulleroidDelaunay: public Triangulation {
public:
  matrix<double> edge_lengths; // edge_lengths(u,v) = length if edge exists, 0 otherwise

  // Construct from any equilateral triangulation.
  // Sorts vertices so cone points (degree != 6) come first, degree-6 last.
  // Initializes all edge lengths to 1.
  FulleroidDelaunay(const Triangulation& T);

  // --- Edge length access ---
  double get_length(node_t u, node_t v) const { return edge_lengths(u,v); }
  void   set_length(node_t u, node_t v, double len) {
    edge_lengths(u,v) = len;
    edge_lengths(v,u) = len;
  }

  // --- Intrinsic geometry ---

  // Extract the diamond geometry around edge (u,v).
  Diamond diamond(node_t u, node_t v) const;

  // --- Delaunay operations ---

  bool is_delaunay_edge(node_t u, node_t v) const;
  bool flip_edge(node_t u, node_t v, bool verbose = false);
  int  flip_to_delaunay();           // Full: Phase 1 + Phase 2 (blocker resolution)
  int  flip_to_delaunay_phase1();    // Phase 1 only: standard Lawson, no blocker resolution
  bool is_delaunay() const;

  // --- Vertex removal ---

  void remove_flat_vertex(node_t v);
  void remove_flat_vertices();

  // --- Validation ---
  bool edge_lengths_are_symmetric() const;

  // --- 3D Embedding ---

  // All-pairs shortest-path distances on the reduced graph (Floyd-Warshall).
  matrix<double> all_pairs_distances() const;

  // Embed the reduced triangulation in 3D to match edge lengths.
  // Uses classical MDS for initial guess, then stress refinement.
  // Returns 3D coordinates for each vertex.
  vector<coord3d> embed_3d() const;
};

#endif
