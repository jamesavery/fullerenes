#pragma once

#include <list>

#include "fullerenes/spiral.hh"
#include "fullerenes/cubicgraph.hh"
#include "fullerenes/geometry.hh"
#include "fullerenes/triangulation.hh"

// FullereneGraph: owned fullerene graph (3-regular, 12 pentagons, rest hex).
// Inherits algorithm methods from FullereneGraphView via Owned<FullereneGraphView>.
// Adds validation constructor and static C20().
// Fullerenes C_N exist for every even N >= 20 except N = 22.
inline bool is_fullerene_size(int N) { return N >= 20 && N != 22 && !(N & 1); }

class FullereneGraph : public Owned<FullereneGraphView> {
  using base_t = Owned<FullereneGraphView>;
public:
  FullereneGraph() = default;
  explicit FullereneGraph(int N) : base_t(N) {}

  FullereneGraph(const GraphView& g) : base_t(g) {
    if(N > 0 && dmax != 3) restride_inplace(3);
    if(N>0) fullerene_check();
  }

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

  static FullereneGraph C20() {
    // CW-oriented neighbour lists for dodecahedral C20, obtained from buckygen
    return FullereneGraph(Graph{
      {1,4,7},   {2,0,9},    {3,1,10},   {4,2,13},
      {5,0,3},   {6,4,14},   {7,5,16},   {8,0,6},
      {9,7,17},  {12,1,8},   {11,2,12},  {13,10,19},
      {10,9,18}, {14,3,11},  {15,5,13},  {16,14,19},
      {17,6,15}, {18,8,16},  {19,12,17}, {15,11,18}
    });
  }
};
