#ifndef DELAUYNAY_HH
#define DELAUYNAY_HH

#include "fullerenegraph.hh"
#include "triangulation.hh"
#include "polyhedron.hh"

class FulleroidDelaunay: public Triangulation {
public:
  matrix<double> distances;
  matrix<double> edge_lengths_d6y; // zero if !edge
  static constexpr double epsilon = 1e-5;

  struct Quad {
    node_t v[4];

    Quad(node_t A, node_t B, node_t C, node_t D) : v{A,B,C,D} {}
    Quad flipped() const { return Quad(v[1],v[2],v[3],v[0]); }

    friend ostream& operator<<(ostream& s, const Quad& q);
  };

  FulleroidDelaunay(const Triangulation& T) : Triangulation(T.sort_nodes()),
					      distances(surface_distances()), edge_lengths_d6y(N,N,0) {
    for(int i=0; i<N; i++){
      for(int j=i; j<N; j++){
        if(edge_exists(dedge_t(i,j))){
          edge_lengths_d6y(i,j) = distances(i,j);
          edge_lengths_d6y(j,i) = distances(j,i);
        }
      }
    }
  }

  double angle(node_t A, node_t B, node_t C) const;
  double angle_d6y(node_t A, node_t B, node_t C) const;
  double angle(const Quad& Q, int i, int subangle = 0) const;
  double angle_d6y(const Quad& Q, int i, int subangle = 0) const;

  bool is_delaunay(const Quad& Q) const { return fabs(angle(Q,1) + angle(Q,3)) - epsilon <= M_PI; }
  bool is_delaunay_d6y(const Quad& Q) const { return fabs(angle_d6y(Q,1) + angle_d6y(Q,3)) - epsilon <= M_PI; }
  bool is_consistent(const Quad& Q,int i) const;
  bool is_consistent(const Quad& Q) const;
  
  void remove_edge_d6y(const dedge_t e){
    cout << "remove " << e << endl;
    remove_edge(e);
    edge_lengths_d6y(e.first,e.second)=0;
    edge_lengths_d6y(e.second,e.first)=0;
  }

  // in the neighbour list of e.first, insert e.second in the position of b
  // in the neighbour list of e.second, insert e.first in the position of d
  void insert_edge_d6y(const edge_t e, const node_t b, const node_t d, const double length){
    cout << "--insert_edge_d6y--" << endl;
    cout << "  edge;length: " << e << "/" << length<< endl;
    insert_edge(e, b, d);
    edge_lengths_d6y(e.first,e.second)=length;
    edge_lengths_d6y(e.second,e.first)=length;
    cout << "--insert_edge_d6y--" << endl;
  }


  vector<dedge_t> triangulate_hole(const vector<node_t>& hole);
  vector<dedge_t> delaunayify_hole(const vector<dedge_t>& edges);
  void delaunayify_hole_2(const vector<edge_t>& edges);

  void remove_flat_vertices();
  void remove_flat_vertex(node_t v);

  bool edge_lengths_d6y_are_symmetric(){
    for(int i=0; i<N; i++){
      for(int j=i; j<N; j++){
        if(edge_lengths_d6y(i,j) != edge_lengths_d6y(j,i)) return false;
      }
    }
    return true;
  }
};

#endif

