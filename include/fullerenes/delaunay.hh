#ifndef DELAUNAY_HH
#define DELAUNAY_HH

#include "triangulation.hh"
#include "geometry.hh"
#include <memory>

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

struct IDTAudit;  // forward declaration

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
  IDTAudit* audit = nullptr;   // null = no checking; non-null = full invariant checking

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
  int  lawson_sweep();               // Standard Lawson: flip all flippable non-Delaunay edges
  int  count_non_delaunay() const;   // Count remaining non-Delaunay edges
  int  delaunay_resolve();           // Search-based escape from Lawson local minima
  int  flip_to_delaunay();           // Full: lawson_sweep + delaunay_resolve
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

// Invariant checker for iDT operations.  Attach to a FulleroidDelaunay via
// its `audit` pointer to enable comprehensive postcondition checking after
// every mutation (flip, vertex removal, Lawson sweep).  Null pointer = no
// checking = zero cost.
//
// Usage:
//   FulleroidDelaunay D(T);
//   IDTAudit audit(D);
//   D.audit = &audit;
//   D.remove_flat_vertices();
//   assert(audit.n_failures == 0);
struct IDTAudit {
  // --- Captured at construction ---
  vector<int> original_degrees;   // original degree of every vertex
  int    original_faces;          // face count of initial triangulation
  double expected_area;           // original_faces * sqrt(3)/4

  // --- Results ---
  int n_checks   = 0;
  int n_failures = 0;

  // --- Options ---
  bool stop_on_failure = true;    // abort on first failure (good for debugging)

  explicit IDTAudit(const FulleroidDelaunay& D);

  // --- Operation hooks (called from FulleroidDelaunay methods) ---
  void after_flip(const FulleroidDelaunay& D, node_t u, node_t v);
  void after_removal(const FulleroidDelaunay& D, node_t removed);
  void after_sweep(const FulleroidDelaunay& D, int n_flips);

private:
  void check_all(const FulleroidDelaunay& D, const char* context);

  bool verify_euler(const FulleroidDelaunay& D, const char* ctx);
  bool verify_orientation(const FulleroidDelaunay& D, const char* ctx);
  bool verify_edge_consistency(const FulleroidDelaunay& D, const char* ctx);
  bool verify_triangle_inequality(const FulleroidDelaunay& D, const char* ctx);
  bool verify_positive_area(const FulleroidDelaunay& D, const char* ctx);
  bool verify_cone_angles(const FulleroidDelaunay& D, const char* ctx);
  bool verify_total_area(const FulleroidDelaunay& D, const char* ctx);
  bool verify_loeschian(const FulleroidDelaunay& D, const char* ctx);
  bool verify_no_multi_edges(const FulleroidDelaunay& D, const char* ctx);

  void fail(const char* invariant, const char* ctx, const string& detail);
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

  // --- Optional exact face-origin tracking ---
  // When non-null, flip_edge() and remove_flat_vertex() use exact Eisenstein
  // arithmetic (turn predicate on the Z[omega] grid) to repartition f_origin.
  // When null, conservative merge is used (both faces get the union).
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
  int  delaunay_resolve();
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
  // exact_origins: when true, f_origin is computed exactly using Eisenstein
  // arithmetic during every flip and vertex removal.  When false, f_origin
  // uses a conservative merge (each face gets the union of its neighbors).
  static DelaunayTriangulation compute(const Triangulation& T,
                                       bool exact_origins = false);

  // --- 3D Embedding ---
  vector<coord3d> embed_3d() const;

  // --- Validation ---
  bool check_consistency() const;
};

#endif
