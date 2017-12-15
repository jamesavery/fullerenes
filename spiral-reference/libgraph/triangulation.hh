#pragma once

#include "spiral.hh"
#include "planargraph.hh"

class Triangulation : public PlanarGraph {
public:
  vector<tri_t> triangles;	// Faces

  // Operations:
  //  1. Orient triangulation
  //  2. Unfold (assert deg(v) <= 6 for all v)
  //  3. GC     (assert deg(v) <= 6 for all v)
  //  4. Spirals (constructor + all_spirals + canonical_spiral)
  //  5. Embed in 2D
  //  6. Embed in 3D
  Triangulation(int N) : PlanarGraph(Graph(N,true)) {}
  Triangulation(const Graph& g = Graph(), bool already_oriented = false) : PlanarGraph(g) { update(already_oriented); }
  Triangulation(const Graph& g, const vector<tri_t>& tris) : PlanarGraph(g), triangles(tris) { 
    orient_triangulation(triangles);
    orient_neighbours();
  }
  Triangulation(const neighbours_t& neighbours, bool already_oriented = false) : PlanarGraph(Graph(neighbours)) { update(already_oriented); }

  Triangulation(const vector<int>& spiral_string, const jumplist_t& jumps = jumplist_t(), const bool best_effort=false); // and the opposite of 'best-effort' is 'fast and robust'
  Triangulation(const spiral_nomenclature &fsn): Triangulation(fsn.spiral_code, fsn.jumps, true){} // best_effort = true

  PlanarGraph dual_graph() const;
  vector<face_t> dual_faces() const;

  // takes a triangulation, and returns the dual of the inverse leapfrog
  // this is cheap because we just remove a set of vertices
  PlanarGraph inverse_leapfrog_dual() const;
  
  pair<node_t,node_t> adjacent_tris(const dedge_t &e) const;

  vector<tri_t> compute_faces() const;          // Returns non-oriented triangles
  vector<tri_t> compute_faces_oriented() const; // Compute oriented triangles given oriented neighbours
  void          orient_neighbours();		// Ensures that neighbours are ordered consistently
  
  //  Unfolding unfold() const;
  Triangulation halma_transform(int m) const;

  // spiral stuff
  bool get_spiral_implementation(const node_t f1, const node_t f2, const node_t f3, vector<int>& v, jumplist_t& j, 
				 vector<node_t>& permutation, const bool general=true,
				 const vector<int>& S0=vector<int>(), const jumplist_t& J0=jumplist_t()) const;
  // the one defined by three nodes
  bool get_spiral(const node_t f1, const node_t f2, const node_t f3, vector<int>& v, jumplist_t& j, vector<node_t>& permutation, const bool general=true) const;


  // Get canonical general spiral and permutation of nodes compared to current triangulation
  bool get_spiral(vector<int>& v, jumplist_t& j, vector<vector<node_t>> &permutations, const bool only_rarest_special=true, const bool general=true) const;  
  // Get canonical general spiral
  bool get_spiral(vector<int>& v, jumplist_t& j, const bool rarest_start=true, const bool general=true) const;
  general_spiral get_general_spiral(const bool rarest_start=true) const;

  void get_all_spirals(vector< vector<int> >& spirals, vector<jumplist_t>& jumps,
		       vector<vector<node_t>>& permutations,
		       const bool only_special=false, const bool general=false) const;

  void symmetry_information(int N_generators, Graph& coxeter_diagram, vector<int>& coxeter_labels) const;

  vector<node_t> vertex_numbers(vector<vector<node_t>> &perms, const vector<node_t> &loc) const;
  
  void update(bool already_oriented) {
    //    renumber(); // TODO!
    if(count_edges() > 0){
      if(already_oriented){
        triangles = compute_faces_oriented();
      }
      else {
        triangles = compute_faces();
        orient_neighbours();
      }
    }
  }


  Triangulation sort_nodes() const;
};

class FullereneDual : public Triangulation {
public:
  // 1. Fullerene get_dual()
  // 2. Construct with buckygen
  // 3. Spiral+gen. spiral special case
  // 4. Embed-in-3D special case
  FullereneDual(const Triangulation& g = Triangulation()) : Triangulation(g) {}
  FullereneDual(const int N, const vector<int>& rspi, const jumplist_t& jumps = jumplist_t()) {
    vector<int> spiral(N/2+2,6);
    for(int i: rspi) spiral[i] = 5;
    *this = Triangulation(spiral,jumps);
  }
  
  bool get_rspi(const node_t f1, const node_t f2, const node_t f3, vector<int>& r, jumplist_t& j, const bool general=true) const;
  bool get_rspi(vector<int>& r, jumplist_t& j, const bool general=true, const bool pentagon_start=true) const;

};


