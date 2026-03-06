#pragma once

#include <list>

#include "fullerenes/spiral.hh"
#include "fullerenes/cubicgraph.hh"
#include "fullerenes/geometry.hh"
#include "fullerenes/triangulation.hh"

class FullereneGraph : public CubicGraph {
public:
  FullereneGraph(const Graph& g) : CubicGraph(g) { if(N>0) fullerene_check();  }
  FullereneGraph(const PlanarGraph& g) : CubicGraph(g) { if(N>0) fullerene_check(); }

  FullereneGraph(const int N, const vector<int>& spiral_indices, const jumplist_t& jumps = jumplist_t()); 
  FullereneGraph(const spiral_nomenclature &fsn){
    *this =  Triangulation(fsn).dual_graph();
  } 
  
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

  // Creates the (k,l)-Goldberg-Coxeter construction C_{(k^2+kl+l^2)n} of the current fullerene C_n
  FullereneGraph GCtransform(unsigned k=1, unsigned l=0) const;

  // spiral from graph, with or without starting point
  bool get_rspi_from_fg(const node_t f1, const node_t f2, const node_t f3, vector<int> &rspi, jumplist_t &jumps, const bool general=true) const;
  bool get_rspi_from_fg(vector<int> &rspi, jumplist_t &jumps, const bool general=true, const bool pentagon_start=true) const;

  // create a matrix that holds the topological distances between all pentagons
  matrix<int> pentagon_distance_mtx() const;

  vector<coord3d> zero_order_geometry(double scalerad=4) const;
  vector<coord3d> optimized_geometry(const vector<coord3d>& initial_geometry, int opt_method = 3, double ftol = 1e-12) const;

  static FullereneGraph C20() {
    // CW-oriented neighbour lists for dodecahedral C20, obtained from buckygen
    return FullereneGraph(Graph(neighbours_t{
      {1,4,7},   {2,0,9},    {3,1,10},   {4,2,13},
      {5,0,3},   {6,4,14},   {7,5,16},   {8,0,6},
      {9,7,17},  {12,1,8},   {11,2,12},  {13,10,19},
      {10,9,18}, {14,3,11},  {15,5,13},  {16,14,19},
      {17,6,15}, {18,8,16},  {19,12,17}, {15,11,18}
    }));
  }
};


