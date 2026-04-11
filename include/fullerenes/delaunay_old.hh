#ifndef DELAUNAY_OLD_HH
#define DELAUNAY_OLD_HH

// Old adjacency-list iDT implementation.
// Superseded by the DCEL-based DelaunayTriangulation in delaunay.hh.
// Kept for reference and for the old test suite.

#include "triangulation.hh"
#include "geometry.hh"

struct Diamond;  // defined in delaunay.hh

struct IDTAudit;  // forward declaration

class FulleroidDelaunay: public Triangulation {
public:
  matrix<double> edge_lengths; // edge_lengths(u,v) = length if edge exists, 0 otherwise
  IDTAudit* audit = nullptr;   // null = no checking; non-null = full invariant checking

  FulleroidDelaunay(const Triangulation& T);

  // --- Edge length access ---
  double get_length(node_t u, node_t v) const { return edge_lengths(u,v); }
  void   set_length(node_t u, node_t v, double len) {
    edge_lengths(u,v) = len;
    edge_lengths(v,u) = len;
  }

  // --- Intrinsic geometry ---
  Diamond diamond(node_t u, node_t v) const;

  // --- Delaunay operations ---
  bool is_delaunay_edge(node_t u, node_t v) const;
  bool flip_edge(node_t u, node_t v, bool verbose = false);
  int  lawson_sweep();
  int  count_non_delaunay() const;
  int  delaunay_resolve();
  int  flip_to_delaunay();
  bool is_delaunay() const;

  // --- Vertex removal ---
  void remove_flat_vertex(node_t v);
  void remove_flat_vertices();

  // --- Validation ---
  bool edge_lengths_are_symmetric() const;

  // --- 3D Embedding ---
  matrix<double> all_pairs_distances() const;
  // embed_3d() not available — use embed_delaunay_3d() instead.
};

struct IDTAudit {
  // --- Captured at construction ---
  vector<int> original_degrees;
  int    original_faces;
  double expected_area;

  // --- Results ---
  int n_checks   = 0;
  int n_failures = 0;

  // --- Options ---
  bool stop_on_failure = true;

  explicit IDTAudit(const FulleroidDelaunay& D);

  // --- Operation hooks ---
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

#endif
