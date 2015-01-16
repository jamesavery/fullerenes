#ifndef DELAUYNAY_HH
#define DELAUYNAY_HH

#include "fullerenegraph.hh"
#include "triangulation.hh"
#include "polyhedron.hh"

class FulleroidDelaunay: public Triangulation {
public:
  matrix<double> distances;
  matrix<double> edge_lengths; // zero if !edge
  static constexpr double epsilon = 1e-5;

  struct Quad {
    node_t v[4];

    Quad(node_t A, node_t B, node_t C, node_t D) : v{A,B,C,D} {}
    Quad flipped() const { return Quad(v[1],v[2],v[3],v[0]); }
  };

  FulleroidDelaunay(const Triangulation& T) : Triangulation(T.sort_nodes()),
					      distances(surface_distances()), edge_lengths(N,N,0) {
    for(int i=0; i<N; i++){
      for(int j=i; j<N; j++){
        if(edge_exists(dedge_t(i,j))){
          edge_lengths(i,j) = distances(i,j);
          edge_lengths(j,i) = distances(j,i);
        }
      }
    }
  }

  double angle(node_t A, node_t B, node_t C) const;
  double angle_d6y(node_t A, node_t B, node_t C) const;
  double angle(const Quad& Q, int i, int subangle = 0) const;

  bool is_delaunay(const Quad& Q) const { return fabs(angle(Q,1) + angle(Q,3)) - epsilon <= M_PI; }
  bool is_consistent(const Quad& Q,int i) const;
  bool is_consistent(const Quad& Q) const;
  
  void remove_edge_d6y(node_t a, node_t c){
    remove_edge(edge_t(a,c));
    edge_lengths(a,c)=0;
    edge_lengths(c,a)=0;
  }

  void insert_edge_d6y(edge_t e, node_t b, node_t d, double length){
    insert_edge(e, b, d);
    edge_lengths(e.first,e.second)=length;
    edge_lengths(e.second,e.first)=length;
  }


  vector<dedge_t> triangulate_hole(const vector<node_t>& hole);
  vector<dedge_t> delaunayify_hole(const vector<dedge_t>& edges);
  void delaunayify_hole_2(const vector<dedge_t>& edges);

  void remove_flat_vertices();
  void remove_flat_vertex(node_t v);

};

#endif

