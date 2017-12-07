#pragma once
#include <list>

#include "spiral.hh"
#include "cubicgraph.hh"
#include "geometry.hh"

class FullereneGraph : public CubicGraph {
public:
  FullereneGraph(const Graph& g, const vector<coord2d>& layout = vector<coord2d>()) : CubicGraph(g,layout) { if(N>0) fullerene_check();  }
  FullereneGraph(const PlanarGraph& g) : CubicGraph(g,g.layout2d) { if(N>0) fullerene_check(); }
  FullereneGraph(const set<edge_t>& edges=set<edge_t>(), const vector<coord2d>& layout = vector<coord2d>()) 
    : CubicGraph(Graph(edges),layout) { if(N>0) fullerene_check(); }
  FullereneGraph(FILE *file) : CubicGraph(file) { if(N>0) fullerene_check(); }
  FullereneGraph(const unsigned int index, FILE *file) : CubicGraph(index, file) { if(N>0) fullerene_check(); }
  FullereneGraph(const int N, const vector<int>& spiral_indices, const jumplist_t& jumps = jumplist_t()); 

  void fullerene_check() const
  {
    if(!is_a_fullerene()){
      fprintf(stderr,"Fullerene graph constructor called for non-fullerene graph.\n");
      abort();
    }
  }

  // Creates the m-point halma-fullerene from the current fullerene C_n with n(1+m)^2 vertices. 
  // (I.e. 4,9,16,25,... for n=1,2,3,4,...)
  FullereneGraph halma_fullerene(const int n, const bool do_layout=false) const;

  // Creates the next leapfrog fullerene C_{3n} from the current fullerene C_n
  FullereneGraph leapfrog_fullerene(const bool do_layout=false) const;

  // spiral from graph, with or without starting point
  bool get_rspi_from_fg(const node_t f1, const node_t f2, const node_t f3, vector<int> &rspi, jumplist_t &jumps, const bool general=true) const;
  bool get_rspi_from_fg(vector<int> &rspi, jumplist_t &jumps, const bool general=true, const bool pentagon_start=true) const;

  // create a matrix that holds the topological distances between all pentagons
  vector<int> pentagon_distance_mtx() const;

  vector<coord3d> zero_order_geometry(double scalerad=4) const;
};


