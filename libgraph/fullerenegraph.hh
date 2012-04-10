#ifndef FULLERENE_GRAPH_HH
# define FULLERENE_GRAPH_HH

#include "cubicgraph.hh"

class FullereneGraph : public CubicGraph {
public:
  FullereneGraph(const Graph& g, const vector<coord2d>& layout = vector<coord2d>()) : CubicGraph(g,layout) { fullerene_check();  }
  FullereneGraph(const set<edge_t>& edges=set<edge_t>(), const vector<coord2d>& layout = vector<coord2d>()) 
    : CubicGraph(Graph(edges),layout) { fullerene_check(); }
  FullereneGraph(FILE *file) : CubicGraph(file) { fullerene_check(); }

  void fullerene_check() const
  {
    if(!this_is_a_fullerene()){
      fprintf(stderr,"Fullerene graph constructor called for non-fullerene graph.\n");
      abort();
    }
  }

  // Creates the m-point halma-fullerene from the current fullerene C_n with n(1+m)^2 vertices. 
  // (I.e. 4,9,16,25,... for n=1,2,3,4,...)
  FullereneGraph halma_fullerene(const unsigned int n, const bool do_layout=false) const;

  // Creates the next leapfrog fullerene C_{3n} from the current fullerene C_n
  FullereneGraph leapfrog_fullerene(const bool do_layout=false) const;

  // Creates the (i,j)-Goldberg-Coxeter construction C_{(i^2+ij+k^2)n} of the current fullerene C_n
  FullereneGraph coxeter_fullerene(const unsigned int i, const unsigned int j, const bool do_layout=false) const;

  // Compute sets <P,H> of pentagonal and hexagonal faces.
  pair<set< face_t>, set<face_t> > compute_faces56() const;



  // static FullereneGraph C20() {
  //   Graph g;
  //   g.neighbours = vector<node_t>(C20_neighbours,C20_neighbours+3*20);
  //   g.update_from_neighbours();
  //   return FullereneGraph(g);
  // }
private:
  static node_t C20_neighbours[20*3];
};

#endif
