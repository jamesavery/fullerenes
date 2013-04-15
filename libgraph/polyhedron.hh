#include "planargraph.hh"
#include "fullerenegraph.hh"
#include <fstream>
#include <sstream>
#include <list>

using namespace std;

struct Polyhedron : public PlanarGraph {
  int face_max;
  vector<coord3d> points;
  coord3d centre;
  vector<face_t> faces;

  //---- Constructors ----//
  // Default constructor
  Polyhedron(const int face_max = INT_MAX) : face_max(face_max) {  }
  Polyhedron(const PlanarGraph& G, const vector<coord3d>& points_ = vector<coord3d>(), const int face_max = INT_MAX,
	     const vector<face_t> faces = vector<face_t>());

  // Create polyhedron from skeleton graph and 3D vertex coordinates 
  // Read polyhedron from a .pol file
  Polyhedron(const string& path) {
    ifstream file(path.c_str());
    file >> *this;
    file.close();
  }
  

  double surface_area() const;

  double volume() const { return volume_divergence(); }
  double volume_tetra() const;
  double volume_divergence() const;

  Polyhedron convex_hull() const { return incremental_convex_hull(); }
  Polyhedron incremental_convex_hull() const;

  void scale(const coord3d& x) {
    for(node_t u=0;u<N;u++) points[u] *= x;
  }

  void move(const coord3d& x) {
    for(node_t u=0;u<N;u++) points[u] += x;
  }

  static Polyhedron C20() {
    vector<coord3d> points(20);
    for(node_t u=0;u<20;u++) points[u] = coord3d(C20_points[u][0],C20_points[u][1],C20_points[u][2]);
      
    return Polyhedron(FullereneGraph::C20(),points);
  }


  static vector<coord3d> polar_mapping(const vector<coord2d>& angles) {
    vector<coord3d> surface(angles.size());
    for(node_t u=0;u<surface.size();u++){
      const double &theta = angles[u].first, &phi = angles[u].second;
      surface[u] = coord3d(cos(theta)*sin(phi), sin(theta)*sin(phi), cos(phi));
    }
    return surface;
  }

  vector<coord2d> polar_angles() const {
    vector<coord2d> angles(N);
    for(node_t u=0;u<N;u++) angles[u] = points[u].polar_angle();
    return angles;
  }

  friend ostream& operator<<(ostream& s, const Polyhedron& P){
    vector<node_t> reachable_points;
    for(node_t u=0;u<P.N;u++) if(P.neighbours[u].size()!=0) reachable_points.push_back(u);
    s << "{{";
    for(unsigned int i=0;i<reachable_points.size();i++) s << reachable_points[i] << (i+1<reachable_points.size()?", ":"},{");
    for(unsigned int i=0;i<P.points.size();i++) s << P.points[i] << (i+1<P.points.size()?", ":"},");
    s << "{";
    for(int i=0;i<P.faces.size();i++){
	s << "{";
	for(int j=0; j < P.faces[i].size();j++) 
	      s << P.faces[i][j] << (j+1<P.faces[i].size()?", ":"}");
      s << (i+1<P.faces.size()?", ":"},");
    }
    s << static_cast<Graph>(P) << "}";
    return s;
  }

  friend istream& operator>>(istream& f, Polyhedron& P){
    string s;
    node_t u=0,v;
    coord3d x;

    while(getline(f,s)){
      stringstream l(s);
      l >> x;
      if(l.fail()) continue; // Invalid line
      l >> v;
      if(l.fail()) continue; // Invalid line
      while(!l.fail()){
	P.edge_set.insert(edge_t(u,v-1)); // File format numbers nodes from 1 
	l >> v;
      }
      P.points.push_back(x);
      u++;
    }
    P.update_auxiliaries();
    P.layout2d = P.tutte_layout();
    P.faces = P.compute_faces_flat(P.face_max);
    P.centre = P.centre3d(P.points);
    return f;
  }

  double  diameter() const;
  coord3d width_height_depth() const;

  string to_latex(bool show_dual = false, bool number_vertices = false, bool include_latex_header = false) const;
  string to_povray(double w_cm = 10, double h_cm = 10, 
		   int line_colour = 0x6a5acd, int vertex_colour = 0xb03500, int face_colour = 0x667744,
		   double line_width = 0.7, double vertex_diameter = 2.0, double face_opacity = 0.5) const;
private:
  static double C20_points[20][3];
};

