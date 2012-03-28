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
  unsigned int N;
  neighbours_t neighbours;
  edges_t edges;
  set<edge_t> edge_set;
  string name;

  Graph(const unsigned int N=0) : N(N) {}
  Graph(const unsigned int N, const set<edge_t>& edge_set) : N(N), edge_set(edge_set) {
    update_auxiliaries();
  }
  Graph(const unsigned int N, const vector<int>& adjacency) : N(N), neighbours(N), edges(N*(N-1)/2) {
    assert(adjacency.size() == N*N);

  //  printf(" New Graph(%d,adjacency)\n",N);
    for(unsigned int i=0;i<N;i++)
      for(unsigned int j=i+1;j<N;j++){
	if(adjacency[i*N+j]) edge_set.insert(edge_t(i,j));
      }

    printf(" %d edges\n",int(edge_set.size()));
    update_auxiliaries();
  }

  vector<coord2d> layout2d; 	// If graph is planar, we can associate a 2D layout

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

  typedef map<unsigned int,set<face_t> > facemap_t;
  facemap_t compute_faces(unsigned int Nmax=INT_MAX, const vector<coord2d> layout = vector<coord2d>()) const;
  facemap_t compute_faces_oriented(const vector<coord2d>& layout) const;
  vector<face_t> compute_faces_flat(unsigned int Nmax=INT_MAX, const vector<coord2d> layout = vector<coord2d>()) const;

  Graph dual_graph(unsigned int Fmax=INT_MAX, const vector<coord2d> layout = vector<coord2d>()) const;

  void    orient_neighbours(const vector<coord2d>& layout);
  coord2d centre2d(const vector<coord2d>& layout) const; // TODO: Move to static member of geometry.hh::coord2d
  coord3d centre3d(const vector<coord3d>& layout) const; // TODO: Move to geometry.hh::coord3d

  int hamiltonian_count() const;
  int hamiltonian_count(const node_t& current_node, vector<bool>& used_edges, vector<bool>& used_nodes, vector<node_t>& path, const vector<unsigned int>& distances) const;

  void update_auxiliaries(); 
  void update_from_neighbours(); 

  friend ostream& operator<<(ostream& s, const Graph& g);

  string to_latex(const vector<coord2d>& layout2d, double w_cm = 10, double h_cm = 10, bool show_dual = false, bool number_vertices = false, bool include_latex_header = false) const;
  string to_latex(const vector<coord3d>& layout3d, bool show_dual = false, bool number_vertices = false, bool include_latex_header = false) const;

  // TODO: Move.
  vector<coord2d> tutte_layout(const node_t s=0, const node_t t=0, const node_t r=0) const;
  vector<face_t> Graph::triangulation(int face_max = INT_MAX) const;
  vector<face_t> Graph::triangulation(const vector<face_t>& faces) const;

};

#endif
