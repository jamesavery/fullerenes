#ifndef PLANARGRAPH_HH
#define PLANARGRAPH_HH

#include "graph.hh"

// TODO: Separate planar cubic graph stuff away from general planar graph into CubicGraph class.
//       Exploit duality between triangulation and cubic planar graph.
class PlanarGraph : public Graph {
public:
  // TODO: Change jumplist_t to vector<pair<int,int>> - faster and easier access!
  typedef list<pair<int,int> > jumplist_t;

  mutable face_t outer_face;
  vector<coord2d> layout2d; 	// If graph is planar, we can associate a 2D layout

  // TODO: Get rid of layout_is_spherical
  PlanarGraph() {}
  PlanarGraph(const PlanarGraph& g) : Graph(g), layout2d(g.layout2d) {  }
  PlanarGraph(const Graph& g, const node_t s=-1, const node_t t=0, const node_t r=0) : Graph(g)
  {
    if(s!=-1){ // Compute planar layout
      layout2d = tutte_layout(s,t,r);
      outer_face = find_outer_face();
    } 
  }

  // Assumes file is at position of a graph start
  static Graph read_hog_planarcode(FILE *planarcode_file);
  static vector<Graph> read_hog_planarcodes(FILE *planarcode_file);
  
  PlanarGraph(const Graph& g, const vector<coord2d>& layout) : Graph(g), layout2d(layout) {  }

  
  bool is_a_fullerene(bool verbose=false) const; // TODO: Do something better with output
  bool is_cubic() const;
  bool is_triangulation() const;

  bool layout_is_crossingfree() const;

  bool is_cut_vertex(const node_t v) const;

  // This group of functions should be used whenever the graphs are oriented. Fmax is not necessary, just an extra back-stop.
  vector<face_t> compute_faces_oriented(int Fmax=INT_MAX) const; // TODO: This should replace the old layout-based method  
  face_t get_face_oriented(const dedge_t &e, int Fmax=INT_MAX) const; 
  dedge_t get_face_representation(dedge_t e, int Fmax=INT_MAX) const; 
  vector<dedge_t> compute_face_representations(int Fmax=INT_MAX) const; // Unique representation of face in oriented planar graph

  // This should all be phased out. For non-oriented graphs, better to compute planar embedding and orient once and for all,
  // than to use planar layout for orientation everywhere (or better yet, make sure graph is oriented in the first place).
  vector<face_t> compute_faces(unsigned int Nmax=INT_MAX, bool planar_layout=false) const;
  vector<face_t> compute_faces_layout_oriented() const;
  face_t get_face_layout_oriented(int u, int v) const; // TODO: This uses the layout, should work towards retirement


  void orient_neighbours();   // Ensures that neighbours are ordered CCW

  face_t find_outer_face() const;


  PlanarGraph dual_graph(unsigned int Fmax=INT_MAX, bool planar_layout=true) const;
  
  size_t count_perfect_matchings() const;

  Graph leapfrog_dual() const;

  vector<int> vertex_numbers(vector<int> &perm) const;
  
  vector<tri_t>  triangulation(int face_max=INT_MAX) const;
  vector<tri_t>  triangulation(const vector<face_t>& faces) const;
  vector<tri_t>  centroid_triangulation(const vector<face_t>& faces) const ;
  vector<tri_t>&  orient_triangulation(vector<tri_t>& tris) const;


  vector<coord2d> tutte_layout(node_t s=0, node_t t=-1, node_t r=-1, unsigned int face_max=INT_MAX) const;
  vector<coord2d> tutte_layout(const face_t& outer_face) const;
  vector<coord2d> tutte_layout_iterative(const face_t& outer_face, const vector<coord2d>& outer_coords) const;
  vector<coord2d> tutte_layout_direct(const face_t& outer_face, const vector<coord2d>& outer_coords) const;
  vector<coord2d> spherical_projection() const;
  bool optimize_layout(const double zv_dist=0.2, const double k_dist=10.0, const double k_angle=10.0, const double k_area=10.0);

  vector<double> edge_lengths() const;
  coord2d width_height() const;
  void scale(const coord2d& x);
  void move (const coord2d& x);

  

  vector<coord3d> zero_order_geometry(double scalerad=4) const;

  string to_latex(double w_cm = 10, double h_cm = 10, bool show_dual = false, bool number_vertices = false, bool include_latex_header = false,
		  int edge_colour = 0x6a5acd, int path_colour = 0x6a5acd, int vertex_colour = 0x8b2500,
		  double edge_width = 0.1,double path_width = 0.1, double vertex_diameter = 2.0,
		  int Npath = 0, int *path = 0
		  ) const;

  string to_povray(double w_cm = 10, double h_cm = 10, 
		   int line_colour = 0x6a5acd, int vertex_colour = 0x8b2500,
		   double line_width = 0.1, double vertex_diameter = 2.0
		   ) const;

  friend ostream& operator<<(ostream& s, const PlanarGraph& g);
};

#endif
