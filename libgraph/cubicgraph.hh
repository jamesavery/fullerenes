#ifndef GRAPH3_HH
# define GRAPH3_HH

#include "planargraph.hh"
#include <iostream>
// TODO: Assumes planarity. Should perhaps split into cubic class and planar class?
struct CubicGraph : public PlanarGraph {

  CubicGraph() {}
  CubicGraph(const PlanarGraph& g) : PlanarGraph(g) {
    for(node_t u=0;u<N;u++)
      if(neighbours[u].size() != 3){
	fprintf(stderr,"Graph not cubic: deg(%d) = %d\n",u,int(neighbours[u].size()));
	abort();
      }
  }

  CubicGraph(const Graph& g, const vector<coord2d>& layout) : PlanarGraph(g,layout) {}

  CubicGraph(FILE *file);

  CubicGraph(const unsigned int *index, FILE *file=stdin);
};

#endif
