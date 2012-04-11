#ifndef PLANARGRAPH_HH
# define PLANARGRAPH_HH

#include "graph.hh"

class PlanarGraph : public Graph {
public:
  using Graph::N;
  using Graph::neighbours;
  using Graph::edges;
  using Graph::edge_set;
  using Graph::name;

  mutable face_t outer_face;
  vector<coord2d> layout2d; 	// If graph is planar, we can associate a 2D layout
  vector<coord2d> spherical_layout;

  PlanarGraph() {}
  PlanarGraph(const PlanarGraph& g) : Graph(g), layout2d(g.layout2d), spherical_layout(g.spherical_layout) {  }
  PlanarGraph(const Graph& g, const node_t s=-1, const node_t t=0, const node_t r=0) : Graph(g)
  {
    if(s!=-1){ // Compute planar layout
      layout2d = tutte_layout(s,t,r);
      outer_face = find_outer_face();
    } 
  }

  PlanarGraph(const Graph& g, const vector<coord2d>& layout) : Graph(g), layout2d(layout) {  }


  bool this_is_a_fullerene() const;

  facemap_t compute_faces(unsigned int Nmax=INT_MAX, bool planar_layout=false) const;
  facemap_t compute_faces_oriented() const;
  vector<face_t> compute_faces_flat(unsigned int Nmax=INT_MAX, bool planar_layout=false) const;
  face_t find_outer_face() const; 

  PlanarGraph dual_graph(unsigned int Fmax=INT_MAX) const;



  vector<face_t>  triangulation(int face_max = INT_MAX) const;
  vector<face_t>  triangulation(const vector<face_t>& faces) const;

  vector<coord2d> tutte_layout(node_t s=0, node_t t=-1, node_t r=-1) const;
  vector<coord2d> spherical_projection() const;

  vector<double> edge_lengths() const;
  coord2d width_height() const;
  void scale(const coord2d& x);
  void move (const coord2d& x);

  string to_latex(double w_cm = 10, double h_cm = 10, bool show_dual = false, bool number_vertices = false, bool include_latex_header = false) const;

  friend ostream& operator<<(ostream& s, const PlanarGraph& g);
};

#endif
