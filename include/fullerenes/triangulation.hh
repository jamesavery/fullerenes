#pragma once

#include <functional>

#include "fullerenes/matrix.hh"
#include "fullerenes/spiral.hh"
#include "fullerenes/planargraph.hh"
#include "fullerenes/planargraph.hh"

// TODO: Easy correspondence between cubic and dual 
//  1. Triangle numbers correspond to cubic nodes
//  2. 
class Triangulation : public PlanarGraph {
public:
  typedef function<bool(Triangulation)> predicate_t;

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
  Triangulation(const neighbours_t& neighbours, bool already_oriented = true) : PlanarGraph(Graph(neighbours)) { update(already_oriented); }

  Triangulation(const vector<int>& spiral_string, const jumplist_t& jumps = jumplist_t(), const bool best_effort=false); // and the opposite of 'best-effort' is 'fast and robust'
  Triangulation(const spiral_nomenclature &fsn): Triangulation(fsn.spiral_code, fsn.jumps, true){} // best_effort = true
  
  PlanarGraph dual_graph() const;
  vector<face_t> cubic_faces() const;
  unordered_map<arc_t,arc_t> arc_translation() const;
  
  size_t max_degree() const {
    size_t max_degree = 0;
    for(auto &nu: neighbours) max_degree = ::max(max_degree, nu.size());
    return max_degree;
  }

  vector<uint8_t> n_degrees() const {
    vector<uint8_t> n_degrees(max_degree(),0);
    for(auto &nu: neighbours) n_degrees[nu.size()-1]++;
    return n_degrees;
  }
  
  // takes a triangulation, and returns a dual of the inverse leapfrog
  // this is cheap because we just remove a set of faces
  PlanarGraph inverse_leapfrog_dual() const;
  
  pair<node_t,node_t> adjacent_tris(const arc_t &e) const;

  vector<tri_t> compute_faces() const;          // Returns non-oriented triangles
  vector<tri_t> compute_faces_oriented() const; // Compute oriented triangles given oriented neighbours
  void          orient_neighbours();		// Ensures that neighbours are ordered consistently
  
  //  Unfolding unfold() const;
  Triangulation GCtransform(const unsigned k=1, const unsigned l=0) const;
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
  
  void update(bool already_oriented=true) {
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

  struct simple_geodesic {
    Eisenstein g;
    int axis;

    simple_geodesic(int a, int b=0, int axis=0) : g(a,b), axis(axis) {}
  };

  struct geodesic {
    vector<Eisenstein> g;
    double d;
    int axis;
  };
  
  matrix<int>    pentagon_distance_mtx() const;


  matrix<int>              simple_square_surface_distances(vector<node_t> only_nodes = {},bool calculate_self_geodesics=false) const;
  matrix<double>           surface_distances(vector<node_t> only_nodes = {},bool calculate_self_geodesics=false) const;
  
  matrix<geodesic>         surface_geodesics(vector<node_t> only_nodes = {},bool calculate_self_geodesics=false) const;    
  matrix<simple_geodesic>  simple_geodesics(vector<node_t> only_nodes = {},bool calculate_self_geodesics=false) const;
  
  node_t         end_of_the_line(node_t u0, int i, int a, int b) const;
  vector<vector<node_t>> quads_of_the_line(node_t u0, int i, int a, int b) const;  

  Triangulation sort_nodes() const;
};


class FullereneDual : public Triangulation {
public:
  // 1. Fullerene get_dual()
  // 2. Construct with buckygen
  // 3. Spiral+gen. spiral special case
  // 4. Embed-in-3D special case
  FullereneDual(const Triangulation& g = Triangulation()) : Triangulation(g) {}
  FullereneDual(const int N, const general_spiral& rspi) : FullereneDual(N,rspi.spiral,rspi.jumps) {}
  FullereneDual(const int N, const vector<int>& rspi, const jumplist_t& jumps = jumplist_t()) {
    vector<int> spiral(N/2+2,6);
    for(int i: rspi) spiral[i] = 5;
    *this = Triangulation(spiral,jumps);
  }
  
  bool get_rspi(const node_t f1, const node_t f2, const node_t f3, vector<int>& r, jumplist_t& j, const bool general=true) const;
  bool get_rspi(vector<int>& r, jumplist_t& j, const bool general=true, const bool pentagon_start=true) const;
  general_spiral get_rspi(const bool rarest_start=true) const; // TODO: Replace above by this simplified API

  static vector<general_spiral> isomer_search(const Triangulation::predicate_t& predicate, size_t N, size_t print_step=0,
					      bool IPR=false, bool only_nontrivial_symmetry=false, size_t N_chunks=1, size_t chunk_index=0);

  spiral_nomenclature name(bool rarest_start=true) const;  
};


class CubicPair {
  Triangulation T;
  PlanarGraph   G;
  IDCounter<tri_t> triangle_id;  
  vector<vector<arc_t>> CtoD, DtoC;
  
  int face_start(const face_t &f){
    node_t i_m = 0;
    for(int i=0, m=INT_MAX; i<int(f.size()); i++) if(f[i] < m){ i_m = i; m = f[i]; }
    return i_m;
  }
    
  CubicPair(const Triangulation &T) : G(T.dual_graph()), CtoD(G.N,vector<arc_t>(3)), DtoC(T.N)
  {
    for(const auto &t: T.triangles) triangle_id.insert(t.sorted());
  
    for(node_t u=0;u<T.N;u++){
      const auto& nu = T.neighbours[u];
      DtoC[u].resize(nu.size()); 

      // For each directed edge v->u
      for(size_t i=0;i<nu.size();i++){

	node_t v = nu[i];
	node_t s = nu[(i+1)%nu.size()];           // u->v->s is triangle associated with u->v
	node_t t = nu[(i+nu.size()-1)%nu.size()]; // v->u->t is triangle associated with v->u

	tri_t t1 = {u,v,s}, t2 = {v,u,t};
	
	// arcs get a unique number: id(u,i) = row_offset[u]+i
	node_t U   = triangle_id(t1.sorted()), V = triangle_id(t2.sorted());
	node_t i_V = G.arc_ix(U,V), i_U = G.arc_ix(V,U);
	node_t i_v = T.arc_ix(u,v), i_u = T.arc_ix(v,u);
	
	CtoD[U][i_V] = {u,i_v};
	CtoD[V][i_U] = {v,i_u};
	DtoC[u][i_v] = {U,i_V};
	DtoC[v][i_u] = {V,i_U};
      }
    }
  }
};
