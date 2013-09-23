#ifndef TRIANGULATION_HH
# define TRIANGULATION_HH

#include "planargraph.hh"
#include "unfold.hh"

class Triangulation : public PlanarGraph {
public:
  using Graph::N;
  using Graph::neighbours;

  map<dedge_t,node_t> nextCW;

  // Operations:
  //  1. Orient triangulation
  //  2. Unfold (assert deg(v) <= 6 for all v)
  //  3. GC     (assert deg(v) <= 6 for all v)
  //  4. Spirals (constructor + all_spirals + canonical_spiral)
  //  5. Embed in 2D
  //  6. Embed in 3D 
  Triangulation(const Graph& g = Graph()) : PlanarGraph(g) { orient_neighbours(); }
  Triangulation(const neighbours_t& neighbours) : PlanarGraph(Graph(neighbours)) { orient_neighbours();  }

  Triangulation(const vector<int>& spiral_string, const jumplist_t& jumps = jumplist_t());


  PlanarGraph& dual_graph() const;
  
  pair<node_t,node_t> adjacent_tris(const edge_t &e) const;

  vector<tri_t> compute_faces() const;          // Returns oriented triangulation
  void          orient_neighbours();		// Ensures that neighbours are ordered consistently
  vector<tri_t> compute_faces_oriented() const; // If orient_neighbours() has been called, compute faces more efficiently
  
  
  Unfolding unfold() const;
  Triangulation GCtransform(int k, int l) const;

  void get_spiral(const node_t f1, const node_t f2, const node_t f3, vector<int>& v, jumplist_t& j, bool general=true) const;
  void get_canonical_spiral(vector<int>& v, jumplist_t& j, bool general=true) const;

};

class FullereneDual : public Triangulation {
public:
  // 1. Fullerene get_dual()
  // 2. Construct with buckygen
  // 3. Spiral+gen. spiral special case
  // 4. Embed-in-3D special case

  void get_canonical_fullerene_rspi(vector<int>& r, jumplist_t& j, bool general=true) const;

};

#endif
