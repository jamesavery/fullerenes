#ifndef GRAPH_HH
#define GRAPH_HH

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

  node_t next(const node_t& u, const node_t& v) const;
  node_t prev(const node_t& u, const node_t& v) const;
  node_t next_on_face(const node_t &u, const node_t &v) const;
  node_t prev_on_face(const node_t &u, const node_t &v) const;

  bool is_consistently_oriented() const;
  bool adjacency_is_symmetric() const;
  bool has_separating_triangles() const;

  bool is_connected(const set<node_t>& subgraph = set<node_t>()) const;
  vector<int> shortest_path(const node_t& source, const node_t& dest, const vector<int>& D0) const;
  vector<int> shortest_paths(const node_t& source, const vector<bool>& used_edges, 
				      const vector<bool>& used_nodes, const unsigned int max_depth = INT_MAX) const;
  vector<int> shortest_paths(const node_t& source, const unsigned int max_depth = INT_MAX) const;
  vector<int> all_pairs_shortest_paths(const unsigned int max_depth = INT_MAX) const;
  list< list<node_t> > connected_components() const;
  
  // Find shortest cycle of the form s->...->s
  vector<node_t> shortest_cycle(const node_t& s, const int max_depth) const;
  // Find shortest cycle of the form s->t->...->s
  vector<node_t> shortest_cycle(const node_t& s, const node_t& t, const int max_depth) const;
  // Find shortest cycle of the form s->t->r->...->s
  vector<node_t> shortest_cycle(const node_t& s, const node_t& t, const node_t& r, const int max_depth) const;
  vector<int> multiple_source_shortest_paths(const vector<node_t>& sources, const unsigned int max_depth=INT_MAX) const;

  vector<int> multiple_source_shortest_paths(const vector<node_t>& sources, const vector<bool>& used_edges, 
						      const vector<bool>& used_nodes, const unsigned int max_depth=INT_MAX) const;

  int hamiltonian_count() const;
  int hamiltonian_count(const node_t& current_node, vector<bool>& used_edges, vector<bool>& used_nodes, vector<node_t>& path, const vector<int>& distances) const;

  coord2d centre2d(const vector<coord2d>& layout) const; 
  coord3d centre3d(const vector<coord3d>& layout) const; 
  void orient_neighbours(const vector<coord2d>& layout);

  int degree(const node_t& u) const;
  int max_degree() const; 

  void update_from_edgeset(const set<edge_t>& edge_set); 
  set<edge_t>  undirected_edges() const;
  set<dedge_t> directed_edges()   const;

  size_t count_edges() const;
  
  friend ostream& operator<<(ostream& s, const Graph& g);
};

#endif
