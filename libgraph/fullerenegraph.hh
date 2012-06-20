#ifndef FULLERENE_GRAPH_HH
# define FULLERENE_GRAPH_HH

#include "cubicgraph.hh"

class FullereneGraph : public CubicGraph {
public:
  FullereneGraph(const Graph& g, const vector<coord2d>& layout = vector<coord2d>()) : CubicGraph(g,layout) { if(N>0) fullerene_check();  }
  FullereneGraph(const set<edge_t>& edges=set<edge_t>(), const vector<coord2d>& layout = vector<coord2d>()) 
    : CubicGraph(Graph(edges),layout) { if(N>0) fullerene_check(); }
  FullereneGraph(FILE *file) : CubicGraph(file) { if(N>0) fullerene_check(); }
  FullereneGraph(int N, const vector<int>& spiral_indices, bool IPR); 

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

  // Creates the (i,j)-Goldberg-Coxeter construction C_{(i^2+ij+k^2)n} of the current fullerene C_n
  FullereneGraph coxeter_fullerene(const unsigned int i, const unsigned int j, const bool do_layout=false) const;

  // Compute sets <P,H> of pentagonal and hexagonal faces.
  pair<set< face_t>, set<face_t> > compute_faces56() const;


  static FullereneGraph C20() {
    PlanarGraph g;
    g.layout2d.resize(20);
    for(node_t u=0;u<20;u++)
      g.layout2d[u] = coord2d(C20_layout2d[u][0],C20_layout2d[u][1]);

    for(int i=0;i<30;i++)
      g.edge_set.insert(edge_t(C20_edges[i][0],C20_edges[i][1]));
    
    g.update_auxiliaries();
    return FullereneGraph(g,g.layout2d);
  }
private:
  static node_t C20_edges[30][2];
  static double C20_layout2d[20][2];
};

#endif
