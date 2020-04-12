#pragma once
#include "fullerenes/triangulation.hh"
#include "fullerenes/eisenstein.hh"

class Unfolding {
public:
  typedef pair<Eisenstein,Eisenstein> dedgecoord_t;

  Triangulation graph;

  map<dedge_t,dedgecoord_t>      edgecoords; // If initialized by
					     // outline, only outline
					     // dedges are defined



  vector< pair<Eisenstein,node_t> > outline; // Polygon outline in the Eisenstein plane. This is always initialized.
  vector<int> degrees;

  // There are two ways to create an "unfolding":
  // 
  // 1. Provide a triangulation of the sphere, i.e., the dual of a planar cubic graph.
  Unfolding(const Triangulation& graph, bool planar_layout = false) : graph(graph), degrees(graph.N) 
  {
    // Store degrees of each node
    vector<tri_t> triangles = graph.compute_faces_oriented();
    edgecoords              = unfold(triangles);
    outline                 = get_outline(edgecoords);
    for(int u=0;u<graph.N;u++) degrees[u] = graph.neighbours[u].size();
  }

  // 2. Provide a simple polygon in the Eisenstein plane, with polygon vertices annotated with
  //    the corresponding node numbers.
  Unfolding(const vector< pair<Eisenstein,node_t> > &outline) : outline(outline) {
    // Calculate degrees of each node and store directed edge coordinates.
    Eisenstein x, y;
    node_t u,v;
    
    polygon P(get_keys(outline));
    for(int i=0;i<outline.size();i++){
      tie(x,u) = outline[i];
      tie(y,v) = outline[(i+1)%outline.size()];


      Eisenstein unit(1,0);
      for(int j=0;j<6;j++,unit = unit.nextCW())
	if(P.point_inside(x+unit)) degrees[u]++;
	
      edgecoords[dedge_t(u,v)] = {x,y};
    }
  }

  // Transform unfolding into a polygon with only straight lines between nodes of
  // degree != 6; i.e., degree 6 nodes can only be interior or lie on a straight line.
  Unfolding straighten_lines() const;

  // This function unfolds a triangulation and lays it out on an equilateral
  // triangular grid, such that if one were to cut along the outline
  // and glue together the nodes with the same labels, one would obtain
  // again the original fullerene dual.
  //
  // Preconditions: Triangles are oriented consistently, i.e. CW or CCW.
  static map<dedge_t,dedgecoord_t> unfold(const vector<tri_t> &triangulation);

  // Compute outline in CW order of the map returned from unfold().
  static vector< pair<Eisenstein,node_t> > get_outline(const map<dedge_t,dedgecoord_t>& edgecoords);

  // Simple transformations in the Eisenstein plane.
  Unfolding& operator *= (const Eisenstein& y){ for(int i=0;i<outline.size();i++) outline[i].first *= y; return *this;  }
  Unfolding operator*(const Eisenstein& y) const { Unfolding U(outline); return (U *= y); }

  Unfolding& operator += (const Eisenstein& y){ for(int i=0;i<outline.size();i++) outline[i].first += y; return *this;  }
  Unfolding operator+(const Eisenstein& y) const { Unfolding U(outline); return (U += y); }

  // Unfolding& operator /= (const Eisenstein& y){ for(int i=0;i<outline.size();i++) outline[i].first /= y; return *this;  }
  // Unfolding operator/(const Eisenstein& y) const { Unfolding U(outline); return (U /= y); }

  static void transform_line(const dedgecoord_t& l1, const dedgecoord_t& l2, Eisenstein& x0, Eisenstein& x0p, Eisenstein& w);



  // Output
  string to_latex(int K=1, int L=0,int label_vertices=1, bool draw_equilaterally=true, bool include_headers=false) const;
  string to_mathematica() const;
};




class Folding {
public:
  typedef Unfolding::dedgecoord_t dedgecoord_t;

  const vector< pair<Eisenstein,node_t> > outline;
  const polygon P;

  IDCounter<Eisenstein> grid;

  // Debug flags:
  enum {WRITE_FILE=1, DONT_ROTATE=2, DONT_CONNECT_POLYGON = 4, DONT_CONNECT_ACROSS = 8, DONT_IDENTIFY_NODES = 16, DO_NONPLANAR_LAYOUT = 32};
  int debug_flags;
  ostream &debug_file;

  Folding(const Unfolding& U, int debug_flags=0, ostream& debug_file = std::cerr) : 
    outline(U.outline), P(get_keys(outline)), debug_flags(debug_flags), debug_file(debug_file)
  {
    const vector<Eisenstein> allpoints(P.allpoints());
    for(auto &x: allpoints) grid.insert(x);
  }

  vector<edge_t> connect();    

  vector<edge_t> connect_polygon(const Eisenstein& w);   // Connect scanlines of polygon rotated by w.
  vector<edge_t> connect_cross(const Eisenstein& w);	 // Connect horizontal lines across edges or polygon rotated by w.
  vector<edge_t> connect(const Eisenstein& w);		 // Both of the above operations.

  // Collect nodes in unfolding that will correspond to the same nodes in the folded graph.
  vector<node_t> identify_nodes() const;
  vector<node_t> outline_nodes() const;

  PlanarGraph fold();
};
