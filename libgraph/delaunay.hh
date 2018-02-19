#ifndef DELAUNAY_HH
#define DELAUNAY_HH

#include "fullerenegraph.hh"
#include "triangulation.hh"
#include "polyhedron.hh"

class FulleroidDelaunay: public Triangulation {
public:
  matrix<double> edge_lengths; // zero if !edge

  FulleroidDelaunay(const Triangulation& T) : Triangulation(T.sort_nodes()), edge_lengths(N,N,0)
  {
    for(node_t u=0;u<N;u++)
      for(node_t v: neighbours[u])
        edge_lengths(u,v) = 1;
  }

  struct Quad {
    /*   v1
     *  /  \
     *v0---v2
     *  \  /
     *   v3
     */
    node_t v[4];

    Quad(node_t A, node_t B, node_t C, node_t D) : v{A,B,C,D} {}
    Quad flipped() const { return Quad(v[1],v[2],v[3],v[0]); }


    vector<int> to_vector() const { return vector<int>(&v[0],&v[4]); }

    friend ostream& operator<<(ostream& s, const Quad& q);
  };


  double tan_halfangle(node_t vi, node_t vj, node_t  vk) const;
  double cot_angle(node_t  vi, node_t  vj, node_t  vk) const;
  double tan_adh(const Quad& Q) const;
  double cos_ad(const Quad& Q) const;
  static double add_tan(double t_a, double t_b);

  /*  length of side f, where:
   *
   *  b /  f
   *   /   |
   *  v0-e--
   *   \   |
   *  c \  f
   */
  double flipped_length(const Quad& Q) const;


  bool is_delaunay(const Quad& Q) const;
  void flip(const Quad& Q);

  void remove_edge(const edge_t e){
    Graph::remove_edge(e);
    edge_lengths(e.first,e.second)=0;
    edge_lengths(e.second,e.first)=0;
  }

  void           align_hole(vector<node_t>& hole) const;
  vector<double>  new_distances(const node_t& v, const vector<node_t>& hole) const;
  vector<edge_t> triangulate_hole(const vector<node_t>& hole, const vector<double>& distances);

  void insert_edge(const dedge_t& uv, const node_t p, const node_t q, const double length){
    node_t u=uv.first, v = uv.second;
    /* Insert v in neighbours[u] just before p, and insert u in neighbours[v] just before q:
     *
     *     p
     *    /
     *   /
     *  u ---- v
     *        /
     *       /
     *      q
     */
    Graph::insert_edge(uv, p, q);
    assert(edge_lengths(u,v)==0);
    assert(edge_lengths(v,u)==0);
    edge_lengths(u,v)=length;
    edge_lengths(v,u)=length;
  }

  bool is_consistent(const Quad& Q,int i) const;
  bool is_consistent(const Quad& Q) const;


  vector<dedge_t> delaunayify_hole(const vector<edge_t>& edges);

  void remove_flat_vertices();
  void remove_flat_vertex(node_t v);

  bool edge_lengths_are_symmetric(){
    for(int i=0; i<N; i++)
      for(int j=i; j<N; j++)
        if(edge_lengths(i,j) != edge_lengths(j,i)) return false;

    return true;
  }

};

#endif

