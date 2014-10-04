#ifndef DELAUYNAY_HH
# define DELAUYNAY_HH

#include "fullerenegraph.hh"
#include "triangulation.hh"
#include "polyhedron.hh"

class FulleroidDelaunay: public Triangulation {
public:
  matrix<double> distances;
  static constexpr double epsilon = 1e-5;

  struct Quad {
    node_t v[4];

    Quad(node_t A, node_t B, node_t C, node_t D) : v{A,B,C,D} {}
    Quad flipped() const { return Quad(v[1],v[2],v[3],v[0]); }
  };

  FulleroidDelaunay(const Triangulation& T) : Triangulation(T.sort_nodes()),
					      distances(surface_distances()) {}

  double angle(node_t A, node_t B, node_t C) const;
  double angle(const Quad& Q, int i, int subangle = 0) const;

  bool is_delaunay(const Quad& Q) const { return fabs(angle(Q,1) + angle(Q,3)) - epsilon <= M_PI; }
  bool is_consistent(const Quad& Q,int i) const;
  bool is_consistent(const Quad& Q) const;
  


  vector<dedge_t> triangulate_hole(const vector<node_t>& hole);
  vector<dedge_t> delaunayify_hole(const vector<dedge_t>& edges);

  void remove_last_vertex();
  void remove_flat_vertices();
};

#endif
