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

  face_t outer_face;
  vector<coord2d> layout2d; 	// If graph is planar, we can associate a 2D layout
  vector<coord2d> spherical_layout;

  PlanarGraph() {}
  PlanarGraph(const Graph& g, const node_t s=-1, const node_t t=0, const node_t r=0) : Graph(g)
  {
    if(s!=-1){ // Compute planar layout
      layout2d = tutte_layout(s,t,r);
      outer_face = find_outer_face();
    } 
  }

  PlanarGraph(const Graph& g, const vector<coord2d>& layout) : Graph(g) { 
    if(layout.size() != N)
      layout2d = tutte_layout();
    else 
      layout2d = layout;

    outer_face = find_outer_face();
  }


  bool this_is_a_fullerene() const;

  facemap_t compute_faces(unsigned int Nmax=INT_MAX) const;
  facemap_t compute_faces_oriented() const;
  vector<face_t> compute_faces_flat(unsigned int Nmax=INT_MAX) const;
  face_t find_outer_face() const; 

  PlanarGraph dual_graph(unsigned int Fmax=INT_MAX) const;

  vector<double> edge_lengths() const;

  vector<face_t>  triangulation(int face_max = INT_MAX) const;
  vector<face_t>  triangulation(const vector<face_t>& faces) const;

  vector<coord2d> tutte_layout(const node_t s=0, const node_t t=0, const node_t r=0) const;
  vector<coord2d> spherical_projection(const vector< coord2d >& layout2d) const;

  void scale(const double s);
  void move (const coord2d s);

  string to_latex(double w_cm = 10, double h_cm = 10, bool show_dual = false, bool number_vertices = false, bool include_latex_header = false) const;
};

#endif
