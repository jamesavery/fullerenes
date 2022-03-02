#pragma once

// TODO FOR SIMPLIFIED & MORE EFFICIENT IMPLEMENTATION:
//  * neighbours_t for CubicGraph is new sparse_matrix datatype whose data is a flat N*d_max neighbours array + Nx1 degrees array
//    This allows dynamic insertion and deletion of edges/arcs without reallocation. Maybe add "compactify()" for when no more dynamic operations are needed.
//  * memory can be allocated at construction + freed at destruction ( bool owns_memory == true )
//    or can be a view of memory allocated elsewhere by passing a pointer at construction ( owns_memory == false)
//  * clean up interface (not everything is needed)
//  * layout2d() should no longer be a part of PlanarGraph - separate out
//  * dedge_t -> arc_t
//  * add more efficient interface that exploits arc representations of the form (u,i) instead of (u,v) (with v = neighbours[u][i])
//  * remove dependency on hashes (Use N*3 tables instead). arc2cubic (tri), arc2dual, carc2darc, etc.
//  * simplify computation of dual
//  * simplify Halma + GC
//  * make simple dynamic array class to replace std::vector (as std::vector doesn't allow initialization from existing memory -> no zero-copy init).
//    => this will let us put a Graph/FullereneGraph/CubicGraph/etc. view on existing flat-memory representation with no overhead.

#include <stdio.h>
#include <map>
#include <vector>
#include <list>
#include <set>
#include <algorithm>
#include <string>
#include <iostream>
#include <limits.h>
#include <math.h>
#include <assert.h>

using namespace std;
#include "geometry.hh"
#include "auxiliary.hh"
#include "matrix.hh"

struct Graph {
  int N;
  neighbours_t neighbours;
  bool is_oriented;
  
  Graph(size_t N=0, bool is_oriented=false) : N(N), neighbours(N), is_oriented(is_oriented) {}
  Graph(const set<edge_t>& edge_set) : is_oriented(false) {
    update_from_edgeset(edge_set);
  }
  Graph(const neighbours_t& neighbours, bool is_oriented=false) : N(neighbours.size()), neighbours(neighbours), is_oriented(is_oriented) { }
  Graph(const unsigned int N, const vector<int>& adjacency, bool is_oriented=false) : N(N), neighbours(N), is_oriented(is_oriented) {
    assert(adjacency.size() == N*N);
    set<edge_t> edge_set;

    for(unsigned int i=0;i<N;i++)
      for(unsigned int j=i+1;j<N;j++){
	if(adjacency[i*N+j]) edge_set.insert(edge_t(i,j));
      }

    update_from_edgeset(edge_set);
  }


  bool insert_edge(const dedge_t& e, const node_t suc_uv=-1, const node_t suc_vu=-1);
  bool remove_edge(const edge_t& e);
  bool edge_exists(const edge_t& e) const;
  void remove_isolated_vertices();
  void remove_vertices(set<int> &sv);

  int  dedge_ix(node_t u, node_t v) const;  
  node_t next(node_t u, node_t v) const;
  node_t prev(node_t u, node_t v) const;
  node_t next_on_face(node_t u, node_t v) const;
  node_t prev_on_face(node_t u, node_t v) const;

  bool is_consistently_oriented() const;
  bool adjacency_is_symmetric() const;
  bool has_separating_triangles() const;

  bool is_connected(const set<node_t>& subgraph = set<node_t>()) const;
  //  vector<int> shortest_path(node_t source, node_t dest, const vector<int>& D0) const;
  void      single_source_shortest_paths(node_t source, int *distances, size_t max_depth = INT_MAX) const;
  matrix<int> all_pairs_shortest_paths(const vector<node_t>& V,
				       const unsigned int max_depth = INT_MAX) const;
  matrix<int> all_pairs_shortest_paths(const unsigned int max_depth = INT_MAX) const;
  vector<vector<node_t> > connected_components() const;
  
  // Find shortest cycle of the form s->...->s
  vector<node_t> shortest_cycle(node_t s, const int max_depth) const;
  // Find shortest cycle of the form s->t->r->...->s
  vector<node_t> shortest_cycle(const vector<node_t> &prefix, const int max_depth) const;
  vector<int> multiple_source_shortest_paths(const vector<node_t>& sources, const unsigned int max_depth=INT_MAX) const;

  // vector<int> multiple_source_shortest_paths(const vector<node_t>& sources, const vector<bool>& used_edges, 
  //						      const vector<bool>& used_nodes, const unsigned int max_depth=INT_MAX) const;

  int hamiltonian_count() const;
  int hamiltonian_count(node_t current_node, vector<bool>& used_edges, vector<bool>& used_nodes, vector<node_t>& path, const vector<int>& distances) const;

  coord2d centre2d(const vector<coord2d>& layout) const; 
  coord3d centre3d(const vector<coord3d>& layout) const; 
  void orient_neighbours(const vector<coord2d>& layout);

  int degree(node_t u) const;
  int max_degree() const; 

  void update_from_edgeset(const set<edge_t>& edge_set); 
  vector<edge_t>  undirected_edges() const;
  vector<dedge_t> directed_edges()   const;

  size_t count_edges() const;
  
  friend ostream& operator<<(ostream& s, const Graph& g);
};


