#ifndef DELAUNAY_HH
#define DELAUNAY_HH

#include "triangulation.hh"
#include "geometry.hh"
#include <memory>

class Symmetry;

// Paired permutations + O(3) matrices for symmetry-constrained embedding.
// perms[i] is a permutation on 0..nv-1, matrices[i] is the corresponding
// orthogonal matrix.  Invariant: matrices[i]*matrices[j] == matrices[k]
// whenever perms[i]*perms[j] == perms[k].
struct SymmetryConstraint {
  vector<vector<int>> perms;
  vector<matrix3d> matrices;
  bool empty() const { return perms.empty(); }
};

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

// ============================================================================
// DCEL-based intrinsic Delaunay triangulation (delta-complex).
//
// Half-edge (DCEL) representation that correctly handles multi-edges and
// self-loops.  Every edge is identified by a half-edge index, so flip_edge(h)
// is unambiguous even when multiple edges connect the same vertex pair.
//
// Twin convention: half-edges 2k and 2k+1 are always twins.
// Face orientation: he_next traverses each face CCW.
// Vertex circulation: cw(h) rotates CW around origin(h).
// ============================================================================

struct DelaunayTriangulation {
  // --- Counts ---
  int nv = 0;   // live vertices
  int nh = 0;   // allocated half-edges (including dead slots)
  int nf = 0;   // allocated faces (including dead slots)

  // --- Half-edge topology (indexed 0..nh-1) ---
  vector<int>    he_next;    // next half-edge CCW in same face
  vector<int>    he_origin;  // origin vertex (-1 = dead)
  vector<int>    he_face;    // face to the left

  // --- Metric (indexed 0..nh-1) ---
  vector<double> he_length;  // edge length (same for h and twin(h))
  vector<double> he_angle;   // angle at origin(h) in face(h)

  // --- Per-vertex (indexed 0..nv-1) ---
  vector<int>    v_out;          // one outgoing half-edge (-1 = dead vertex)
  vector<double> v_cone_angle;   // original_degree * pi/3
  vector<int>    v_orig_degree;  // degree in original equilateral triangulation

  // --- Per-face (indexed 0..nf-1) ---
  vector<int>         f_he;      // one boundary half-edge (-1 = dead face)
  vector<vector<int>> f_origin;  // original face IDs covered by this iDT face

  // --- Free lists ---
  vector<int> free_edges;  // recycled edge slots (half-edge id / 2)
  vector<int> free_faces;  // recycled face slots

  // --- Exact face-origin tracking ---
  // flip_edge() and remove_flat_vertex() use exact Eisenstein arithmetic
  // (turn predicate on the Z[omega] grid) to repartition f_origin.
  struct OriginTracker;
  std::shared_ptr<const OriginTracker> origin_tracker;

  // --- Clean accessors ---
  int  twin(int h)  const { return h ^ 1; }
  int  edge(int h)  const { return h >> 1; }
  int  prev(int h)  const { return he_next[he_next[h]]; }  // only for triangulations
  int  dest(int h)  const { return he_origin[h ^ 1]; }
  bool alive(int h) const { return he_origin[h] >= 0; }

  // CW rotation around origin(h): next outgoing half-edge clockwise.
  int cw(int h) const { return he_next[h ^ 1]; }

  // CCW rotation around origin(h): next outgoing half-edge counterclockwise.
  int ccw(int h) const { return (he_next[he_next[h]]) ^ 1; }

  int vertex_degree(int v) const;  // count outgoing half-edges from v

  // --- Construction ---
  static DelaunayTriangulation from_triangulation(const Triangulation& T);

  // --- Geometry ---
  Diamond diamond(int h) const;
  void recompute_face_angles(int f);
  void recompute_all_angles();

  // --- Delaunay operations ---
  bool is_delaunay_edge(int h) const;
  bool flip_edge(int h);
  int  lawson_sweep();
  int  count_non_delaunay() const;
  int  flip_to_delaunay();
  bool is_delaunay() const;

  // --- Vertex removal ---
  void remove_flat_vertex(int v);
  void remove_flat_vertices();

  // --- Edge/face allocation ---
  int  alloc_edge();         // returns half-edge id of first half-edge
  int  alloc_face();         // returns face id
  void dealloc_edge(int h);  // mark edge as dead, add to free list
  void dealloc_face(int f);  // mark face as dead, add to free list

  // --- Full algorithm ---
  // Computes the intrinsic Delaunay triangulation.
  // track_origins: when true, f_origin is computed exactly using Eisenstein
  // arithmetic during every flip and vertex removal.  When false, f_origin
  // is left empty (no origin tracking).
  static DelaunayTriangulation compute(const Triangulation& T,
                                       bool track_origins = false);

  // --- 3D Embedding ---
  // Embed the reduced triangulation in 3D to match intrinsic edge lengths.
  // When symmetry is provided (permutations + O(3) matrices from
  // representation_3d()), optimization runs in the reduced orbit-rep
  // parameter space: one coord3d per orbit, projected into the stabilizer's
  // fixed-point subspace.  All orbit members are generated by the group
  // action, so the embedding is exactly symmetric by construction.
  vector<coord3d> embed_3d(const SymmetryConstraint& sym = {}) const;

  static DelaunayTriangulation compute_symmetric(const Triangulation& T,
                                                  const Symmetry& S);
  int check_edge_symmetry(const vector<vector<int>>& cone_perms) const;

  // Bisect all multi-edges by inserting midpoint vertices.
  // Multi-edges (multiple geodesics between the same cone-point pair) can't be
  // realized as distinct straight edges in R³.  Bisecting each with a midpoint
  // makes them geometrically distinct.  Returns the number of vertices added.
  int bisect_multi_edges();

  // --- Validation ---
  bool check_consistency() const;
};

// Restrict full-triangulation automorphisms to iDT cone-point indices.
// G: automorphisms of the full triangulation (from Symmetry::G).
// T: the original triangulation (before sort_flat_last / iDT computation).
// Returns permutations on 0..11 (the iDT vertex ordering: cone points sorted
// by original index). Empty input → empty output.
vector<vector<int>> restrict_to_cone_points(
    const vector<vector<int>>& G, const Triangulation& T);

// Restrict full-triangulation symmetry (permutations + O(3) matrices) to
// iDT cone-point indices, maintaining the perm<->matrix pairing.
// G: automorphisms (Symmetry::G), R: 3D matrices (Representation3D::R).
SymmetryConstraint restrict_symmetry_to_cone_points(
    const vector<vector<int>>& G, const vector<matrix3d>& R,
    const Triangulation& T);

// Compute vertex orbits from a group of permutations (union-find).
vector<vector<int>> compute_orbits(int n, const vector<vector<int>>& G);

#endif
