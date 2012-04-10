#ifndef GRAPH_HH
# define GRAPH_HH

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


// TODO: Separate planar graph stuff from general graph stuff.
// TODO: Don't constantly pass layout2d around -- it's a member!
struct Graph {
  int N;
  neighbours_t neighbours;
  edges_t edges;
  set<edge_t> edge_set;
  string name;

  Graph(const set<edge_t>& edge_set = set<edge_t>()) : edge_set(edge_set) {
    update_auxiliaries();
  }
  Graph(const unsigned int N, const vector<int>& adjacency) : N(N), neighbours(N), edges(N*(N-1)/2) {
    assert(adjacency.size() == N*N);

    for(unsigned int i=0;i<N;i++)
      for(unsigned int j=i+1;j<N;j++){
	if(adjacency[i*N+j]) edge_set.insert(edge_t(i,j));
      }

    update_auxiliaries();
  }

  vector<unsigned int> shortest_paths(const node_t& source, const vector<bool>& used_edges, 
				      const vector<bool>& used_nodes, const unsigned int max_depth = INT_MAX) const;
  vector<unsigned int> shortest_paths(const node_t& source, const unsigned int max_depth = INT_MAX) const;
  vector<unsigned int> all_pairs_shortest_paths(const unsigned int max_depth = INT_MAX) const;
  
  // Find shortest cycle of the form s->t->...->s
  vector<node_t> shortest_cycle(const node_t& s, const node_t& t, const int max_depth=INT_MAX) const;
  // Find shortest cycle of the form s->t->r->...->s
  vector<node_t> shortest_cycle(const node_t& s, const node_t& t, const node_t& r, const int max_depth) const;
  vector<unsigned int> multiple_source_shortest_paths(const vector<node_t>& sources, const unsigned int max_depth=INT_MAX) const;

  vector<unsigned int> multiple_source_shortest_paths(const vector<node_t>& sources, const vector<bool>& used_edges, 
						      const vector<bool>& used_nodes, const unsigned int max_depth=INT_MAX) const;


  int hamiltonian_count() const;
  int hamiltonian_count(const node_t& current_node, vector<bool>& used_edges, vector<bool>& used_nodes, vector<node_t>& path, const vector<unsigned int>& distances) const;

  coord2d centre2d(const vector<coord2d>& layout) const; 
  coord3d centre3d(const vector<coord3d>& layout) const; 
  void orient_neighbours(const vector<coord2d>& layout);

  void update_auxiliaries(); 
  void update_from_neighbours(); 

  friend ostream& operator<<(ostream& s, const Graph& g);

};

#endif
