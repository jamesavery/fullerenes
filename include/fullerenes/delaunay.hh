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

  // Cocircular ("tight") test: cot(angle_B) + cot(angle_D) == 0 exactly,
  // i.e. the four points u, v, B, D are concyclic on the surface.  In this
  // case both triangulations of the diamond are equally-valid Delaunay
  // refinements of the same cell.  Uses exact integer arithmetic on
  // length-squared (valid when all five lengths square to non-negative
  // integers, e.g. equilateral triangulations and their flips).
  // See CANONICAL-TESSELATION.md for the derivation.
  bool is_cocircular() const;
};

// ============================================================================
// Canonical Delaunay tesselation
//
// On a piecewise-flat surface, the iDT triangulation is not unique: when
// k >= 4 cone points are concyclic, the cocircular k-gon admits multiple
// equally-valid Delaunay triangulations, and which one the algorithm emits
// depends on the input vertex labelling.  The Bobenko-Springborn (2007)
// theorem guarantees that the underlying TESSELATION (cell decomposition
// where every cocircular polygon is left intact) IS unique.
//
// CanonicalTesselation is exactly that invariant: a sorted multi-set of
// CCW-oriented polygons, each in user-supplied vertex labels with
// integer length-squared annotation per edge.  Two iDT outputs that
// differ only by trivial flips inside cocircular cells (and/or by a
// label-equivariant relabelling, once mapped through the same external
// label-map) compare equal here, even though their raw triangulations
// differ.
// ============================================================================

struct CanonicalTesselation {
  // Polygon = CCW boundary of one Delaunay cell.  Each entry is
  // (vertex_label, length²-of-outgoing-edge-to-next-vertex).  Polygons
  // are normalized to lex-min cyclic rotation.
  using Polygon = std::vector<std::pair<int, long long>>;

  // Sorted multi-set of cells (a deterministic linear order is enough for
  // equality / ordering).
  std::vector<Polygon> cells;

  bool operator==(const CanonicalTesselation& o) const { return cells == o.cells; }
  bool operator!=(const CanonicalTesselation& o) const { return cells != o.cells; }
  bool operator< (const CanonicalTesselation& o) const { return cells <  o.cells; }

  // Counts (for quick queries; empty iDT yields zeros).
  int n_cells() const { return (int)cells.size(); }

  // Stable 64-bit fingerprint for hash maps and at-a-glance comparison.
  size_t fingerprint() const;

  // Compact human-readable form; one cell per line.
  std::string to_string() const;
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
  // Flip the diagonal of the diamond around edge h.  Accepts any
  // non-Delaunay edge with a convex diamond, including the B == D case
  // (which produces a self-loop edge at B; see Lemma 1 of
  // claude-projects/delaunay/CORRECTNESS-PROOF.md).
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
  // none.  A value below 3 is one (but not the only) non-simplicial
  // signature -- use is_simplicial() for a complete check.
  int min_live_degree() const;

  // True iff the iDT's 1-skeleton is a simple graph: no self-loops,
  // no multi-edges.  Equivalently: the map h |-> (origin(h), dest(h))
  // is injective on live half-edges.  Self-loops fail because both
  // twins encode the arc (v,v); multi-edges fail because two non-twin
  // half-edges encode the same arc.
  // Non-simplicial outputs are valid iDT delta-complexes (Bobenko-
  // Springborn 2007), not algorithm failures.
  // Cost: O(E log E).
  bool is_simplicial() const;

  // True iff the DCEL is structurally well-formed: every live half-edge
  // belongs to exactly one face cycle traversed via he_next, and every
  // such cycle has length 3.  This is the DCEL counterpart of
  // Graph::is_consistently_oriented (which walks faces via next(v,u)).
  // A well-formed iDT can still be non-simplicial (multi-edges, self-
  // loops); a not-well-formed iDT signals a bug in the algorithm.
  // Cost: O(E).
  bool is_well_formed() const;

  // Bisect all multi-edges by inserting midpoint vertices.
  // Multi-edges (multiple geodesics between the same cone-point pair) can't be
  // realized as distinct straight edges in R³.  Bisecting each with a midpoint
  // makes them geometrically distinct.  Returns the number of vertices added.
  int bisect_multi_edges();

  // --- Tesselation invariant ---
  // True iff edge h is cocircular ("tight"): the diamond's four cone points
  // share a common circumcircle, so flipping h yields an equally-valid
  // Delaunay triangulation.  Both half-edges of an edge return the same value.
  bool is_cocircular_edge(int h) const;

  // Per-half-edge cocircular mask: tight[h] == tight[h^1]; dead half-edges
  // are false.  O(num_edges) integer-arithmetic predicates.
  std::vector<bool> cocircular_edges() const;

  // Canonical Delaunay tesselation in `vertex_labels` coordinates.
  // `vertex_labels[k]` is the external label assigned to D's live vertex k;
  // typically the vertex's position in the original Triangulation it came
  // from (as recovered through `Triangulation::sort_flat_last()`).
  // Cocircular interior edges are merged so each cell becomes one polygon;
  // polygons are normalized (lex-min rotation), the multi-set is sorted.
  CanonicalTesselation canonical_tesselation(
      const std::vector<int>& vertex_labels) const;

  // General tesselation extraction with a caller-supplied "interior" mask.
  // `tight[h]` true ⇒ edge h is interior to a cell (it is collapsed away);
  // `tight[h]` false ⇒ edge h is on a cell boundary.  Both half-edges of an
  // edge must agree (`tight[h] == tight[h^1]`).
  //
  // The cell-walk is identical to canonical_tesselation; only the source of
  // the interior-edge mask changes.  This lets the Alexandrov solver pass
  // a Bobenko-Izmestiev "inessential" mask (edges with θ_e = π in the
  // weighted-Delaunay tesselation at κ=0) to obtain the polytope's
  // 2-skeleton T̄ (the polygonal flat 2-faces of the convex polytope),
  // distinct from but compatible with the cocircular tesselation.
  //
  // Both notions agree exactly on flat-2-face diagonals whose four cone
  // points are also cocircular in the surface metric; they may differ
  // when the flat 2-face is non-cyclic, or when 4 surface-cocircular
  // points are spread across multiple polytope faces.  Cross-validation
  // between the two is a useful sanity check on the GCP geometry at κ=0.
  CanonicalTesselation canonical_tesselation(
      const std::vector<int>& vertex_labels,
      const std::vector<bool>& tight) const;

  // --- Validation ---
  bool check_consistency() const;
};


#endif
