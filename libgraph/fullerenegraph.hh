#ifndef FULLERENE_GRAPH_HH
# define FULLERENE_GRAPH_HH

#include "cubicgraph.hh"

class FullereneGraph : public CubicGraph {
public:
  FullereneGraph(const Graph& g) : CubicGraph(g) { fullerene_check();  }
  FullereneGraph(const unsigned int N, const vector<node_t>& neighbours) : CubicGraph(N,neighbours) { fullerene_check(); }
  FullereneGraph(FILE *file = stdin) : CubicGraph(file) { fullerene_check(); }

  bool this_is_a_fullerene() const;
  void fullerene_check() const
  {
    if(!this_is_a_fullerene()){
      fprintf(stderr,"Fullerene graph constructor called for non-fullerene graph.\n");
      abort();
    }
  }

  // Creates the m-point halma-fullerene from the current fullerene C_n with n(1+m)^2 vertices. 
  // (I.e. 4,9,16,25,... for n=1,2,3,4,...)
  FullereneGraph halma_fullerene(const unsigned int n, const bool do_layout=true) const;

  // Creates the next leapfrog fullerene C_{3n} from the current fullerene C_n
  FullereneGraph leapfrog_fullerene(const bool do_layout=true) const;

  // Compute sets <P,H> of pentagonal and hexagonal faces.
  pair<set< face_t>, set<face_t> > compute_faces56() const;



  static FullereneGraph C20() {
    return CubicGraph(20,vector<node_t>(C20_neighbours,C20_neighbours+3*20));
  }
private:
  static node_t C20_neighbours[20*3];
};

#endif
