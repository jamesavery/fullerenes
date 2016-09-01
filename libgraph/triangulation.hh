#ifndef TRIANGULATION_HH
#define TRIANGULATION_HH

#include "matrix.hh"
#include "planargraph.hh"
#include "unfold.hh"

class Triangulation : public PlanarGraph {
public:
  vector<tri_t> triangles;	// Faces

  struct general_spiral { 
    vector<int> spiral;
    jumplist_t  jumps;
    
    bool operator< (const general_spiral& s){ 
      return (jumps.size() < s.jumps.size()) 
	|| (jumps.size() == s.jumps.size() && jumps < s.jumps)
	|| (jumps == s.jumps && spiral < s.spiral);
    }
    bool operator==(const general_spiral& s){ return jumps == s.jumps && spiral == s.spiral;   }
  };
  // Operations:
  //  1. Orient triangulation
  //  2. Unfold (assert deg(v) <= 6 for all v)
  //  3. GC     (assert deg(v) <= 6 for all v)
  //  4. Spirals (constructor + all_spirals + canonical_spiral)
  //  5. Embed in 2D
  //  6. Embed in 3D 
  Triangulation(const Graph& g = Graph(), bool already_oriented = false) : PlanarGraph(g) { update(already_oriented); }
  Triangulation(const Graph& g, const vector<tri_t>& tris) : PlanarGraph(g), triangles(tris) { 
    orient_triangulation(triangles);
    orient_neighbours();
  }
  Triangulation(const neighbours_t& neighbours, bool already_oriented = false) : PlanarGraph(Graph(neighbours)) { update(already_oriented); }

  Triangulation(const vector<int>& spiral_string, const jumplist_t& jumps = jumplist_t());

  PlanarGraph dual_graph() const;
  vector<face_t> dual_faces() const;
  
  pair<node_t,node_t> adjacent_tris(const edge_t &e) const;

  
  node_t nextCW(const dedge_t& uv)   const { return nextCW(uv.first,uv.second);  } // TODO: Remove.
  node_t nextCCW(const dedge_t& uv)  const { return nextCCW(uv.first,uv.second); }
  node_t nextCW(node_t u, node_t v) const;
  node_t nextCCW(node_t u, node_t v) const;


  vector<tri_t> compute_faces() const;          // Returns non-oriented triangles
  vector<tri_t> compute_faces_oriented() const; // Compute oriented triangles given oriented neighbours
  void          orient_neighbours();		// Ensures that neighbours are ordered consistently
  
  //  Unfolding unfold() const;
  Triangulation GCtransform(const unsigned k=1, const unsigned l=0) const;

  // spiral stuff
  bool get_spiral_implementation(const node_t f1, const node_t f2, const node_t f3, vector<int>& v, jumplist_t& j, 
				 vector<node_t>& permutation, const bool general=true, const vector<int>& S0=vector<int>()) const;
  // the one defined by three nodes
  bool get_spiral(const node_t f1, const node_t f2, const node_t f3, vector<int>& v, jumplist_t& j, vector<node_t>& permutation, const bool general=true) const;
  // the canonical one
  bool get_spiral(vector<int>& v, jumplist_t& j, const bool canonical=true, const bool only_special=false, const bool general=true) const;
  void get_all_spirals(vector< vector<int> >& spirals, vector<jumplist_t>& jumps, // TODO: Should only need to supply jumps when general=true
		       vector< vector<int> >& permutations,
		       const bool only_special=false, const bool general=false) const;

  void symmetry_information(int N_generators, Graph& coxeter_diagram, vector<int>& coxeter_labels) const;

  void update(bool already_oriented) {
    //    renumber(); // TODO!
    if(N>0){
      if(already_oriented){
        triangles = compute_faces_oriented();
      }
      else {
        triangles = compute_faces();
        orient_neighbours();
      }
    }
  }

  matrix<double> surface_distances() const;
  matrix<int>    convex_square_surface_distances() const;
  node_t         end_of_the_line(node_t u0, int i, int a, int b) const;

  Triangulation sort_nodes() const;
};

class FullereneDual : public Triangulation {
public:
  // 1. Fullerene get_dual()
  // 2. Construct with buckygen
  // 3. Spiral+gen. spiral special case
  // 4. Embed-in-3D special case
  FullereneDual(int N, const vector<int>& rspi) {
    vector<int> spiral(N/2+2,6);
    for(int i=0;i<12;i++) spiral[rspi[i]-1] = 5;
    *this = Triangulation(spiral);
  }
  FullereneDual(const Triangulation& g = Triangulation()) : Triangulation(g) {}

  vector<coord3d> zero_order_geometry(double scalerad=4) const;
  vector<coord3d> optimized_geometry(const vector<coord3d>& initial_geometry) const;

  bool get_rspi(const node_t f1, const node_t f2, const node_t f3, vector<int>& r, jumplist_t& j, bool general=false) const;
  bool get_rspi(vector<int>& r, jumplist_t& j, bool canonical=true, bool only_special = false, bool general=false) const;
  bool get_rspi_standard(const node_t f1, const node_t f2, const node_t f3, vector<int>& rspi) const;
  bool get_rspi_standard(vector<int>& rspi, bool canonical, bool only_special) const;
};

#endif
