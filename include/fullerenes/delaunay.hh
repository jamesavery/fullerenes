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
  vector<int> f_he;        // one boundary half-edge (-1 = dead face)

  // --- Free lists ---
  vector<int> free_edges;  // recycled edge slots (half-edge id / 2)
  vector<int> free_faces;  // recycled face slots

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
  // Like flip_edge, but does not reject the B == D case.  A successful
  // flip with B == D produces a self-loop edge at B.  The result is a
  // well-formed delta-complex triangulation; the caller takes
  // responsibility for handling the self-loop.
  bool flip_edge_allow_self_loop(int h);
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

  // Allocate an edge and set its endpoints and length.
  // Returns the half-edge h with origin(h) = u, origin(twin h) = v,
  // length = L on both sides.  Faces remain unassigned.
  int  alloc_directed_edge(int u, int v, double L);

  // Wire three half-edges into a CCW triangle face and compute its angles
  // from the stored edge lengths.  Returns the new face id.
  // Preconditions: h0, h1, h2 already have their origin and length set;
  // their endpoints form a triangle with origin(h0)=u, dest(h0)=v=origin(h1),
  // dest(h1)=w=origin(h2), dest(h2)=u.
  int  wire_triangle(int h0, int h1, int h2);

  // Ensure v_out[v] points to a live outgoing half-edge from v.
  // Walks the CW ring at v if the current pointer is stale.  Leaves v_out[v]
  // at -1 iff v has no incident live edge (i.e. v is dead).
  void ensure_v_out(int v);

  // --- Full algorithm ---
  // Computes the unique intrinsic Delaunay triangulation of the input
  // surface (Bobenko-Springborn 2007).  The output is generally a
  // delta-complex and may contain multi-edges, self-loops, and bigons
  // around cone vertices (deg-2 cones) -- all valid iDT features.
  static DelaunayTriangulation compute(const Triangulation& T);

  // Smallest degree among live (non-removed) vertices, or INT_MAX if
  // none.  A value below 3 indicates the iDT has at least one bigon
  // (a non-simplicial feature, valid in delta-complex form).
  int min_live_degree() const;

  // Bisect all multi-edges by inserting midpoint vertices.
  // Multi-edges (multiple geodesics between the same cone-point pair) can't be
  // realized as distinct straight edges in R³.  Bisecting each with a midpoint
  // makes them geometrically distinct.  Returns the number of vertices added.
  int bisect_multi_edges();

  // --- Validation ---
  bool check_consistency() const;
};


#endif
