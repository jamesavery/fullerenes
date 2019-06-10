#ifndef GRAPH3_HH
#define GRAPH3_HH

#include "planargraph.hh"
#include <iostream>

typedef list<pair<int,int> > jumplist_t;

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
  CubicGraph(const unsigned int index, FILE *file=stdin);
  CubicGraph(const int N, const vector<int>& spiral_string, const jumplist_t& jumps = jumplist_t());

  bool get_spiral_from_cg(const node_t f1, const node_t f2, const node_t f3, vector<int> &spiral, jumplist_t &jumps, const bool general=true) const;
  bool get_spiral_from_cg(vector<int> &spiral, jumplist_t &jumps, const bool canonical=true, const bool general=true) const;

  // creates the (k,l)-Goldberg-Coxeter construction C_{(k^2+kl+l^2)n} of the current C_n
  CubicGraph GCtransform(const unsigned k=1, const unsigned l=0, const bool do_layout=false) const;

};

#endif
