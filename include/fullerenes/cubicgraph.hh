#pragma once

#include <iostream>
#include "fullerenes/planargraph.hh"
#include "fullerenes/spiral.hh"

class Triangulation;

// Cubic planar graphs (cubic polyhedral graphs)
// TODO: Make this more useful, merge cubic-specific stuff from planargraph, make triangulation dual produce cubicgraph
struct CubicGraph : public PlanarGraph {

  //  CubicGraph() {} // TODO: Does it make sense to have a default constructor for a span-view?
  CubicGraph(const PlanarGraph& g) : PlanarGraph(g) {
    for(node_t u=0;u<N;u++)
      if(neighbours[u].size() != 3){
        fprintf(stderr,"Graph not cubic: deg(%d) = %d\n",u,int(neighbours[u].size()));
        abort();
      }
  }

  CubicGraph(const Graph& g, const vector<coord2d>& layout) : PlanarGraph(g,layout) {}
  //TODO: Sane constructors with new graph format. But also make CubicGraph a useful class.
  //  CubicGraph(const int N, const vector<int>& spiral_string, const jumplist_t& jumps = jumplist_t());
  //  CubicGraph(const spiral_nomenclature &fsn);

  bool get_spiral_from_cg(const node_t f1, const node_t f2, const node_t f3, vector<int> &spiral, jumplist_t &jumps, const bool general=true) const;
  bool get_spiral_from_cg(vector<int> &spiral, jumplist_t &jumps, const bool canonical=true, const bool general=true, const bool pentagon_start=true) const;

  vector<node_t> vertex_numbers(const Triangulation &T, const vector<vector<node_t>> &perm, const vector<node_t>& loc) const;
};

