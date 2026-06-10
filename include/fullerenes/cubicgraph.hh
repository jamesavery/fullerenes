#pragma once

#include <iostream>
#include "fullerenes/owned.hh"
#include "fullerenes/spiral.hh"

class Triangulation;

// CubicGraph: owned cubic planar graph.
// Inherits algorithm methods from CubicGraphView via Owned<CubicGraphView>.
// Adds validation constructors that check degree-3 and restride.
struct CubicGraph : public Owned<CubicGraphView> {
  using base_t = Owned<CubicGraphView>;

  CubicGraph() {}

  CubicGraph(const GraphView& g) : base_t(g) {
    for(node_t u=0;u<N;u++)
      if((*this)[u].size() != 3){
        fprintf(stderr,"Graph not cubic: deg(%d) = %d\n",u,int((*this)[u].size()));
        abort();
      }
    if(N > 0 && dmax != 3) restride_inplace(3);
  }

  CubicGraph(const int N, const vector<int>& spiral_string, const jumplist_t& jumps = jumplist_t());
  CubicGraph(const spiral_nomenclature &fsn);
};
