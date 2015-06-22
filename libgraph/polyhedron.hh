#ifndef POLYHEDRON_HH
#define POLYHEDRON_HH

#include "planargraph.hh"
#include "fullerenegraph.hh"
#include <fstream>
#include <sstream>
#include <list>

using namespace std;

struct Polyhedron : public PlanarGraph {
  int face_max;
  vector<coord3d> points;
  vector<face_t> faces;

  //---- Constructors ----//
  // Default constructor
  Polyhedron(const int face_max = INT_MAX) : face_max(face_max) {  }
  Polyhedron(const PlanarGraph& G, const vector<coord3d>& points_ = vector<coord3d>(), const int face_max = INT_MAX,
	     const vector<face_t> faces = vector<face_t>());

  // Create polyhedron from skeleton graph and 3D vertex coordinates 
  // Read polyhedron from a .xyz or .mol2 file.
  Polyhedron(const string& path);
  // Create polyhedron from point collection, assuming shortest distance is approximate bond length
  Polyhedron(const vector<coord3d>& xs, double tolerance = 1.2);

  double surface_area() const;

  double volume() const { return volume_divergence(); }
  double volume_tetra() const;
  double volume_divergence() const;

  Polyhedron convex_hull() const { return incremental_convex_hull(); }
  Polyhedron incremental_convex_hull() const;

  matrix3d inertia_matrix() const;
  matrix3d inertial_frame() const;

  void scale(const coord3d& x) {
    for(node_t u=0;u<N;u++) points[u] *= x;
  }

  void move(const coord3d& x) {
    for(node_t u=0;u<N;u++) points[u] += x;
  }

  void move_to_origin() {
    coord3d x0(centre3d(points));
    move(-x0);
  }
  void align_with_axes(){
    matrix3d If(inertial_frame());
    points = If*points;
  }

  void orient_neighbours();

  bool optimize(int opt_method = 3, double ftol = 1e-10);
  bool optimize_other(bool optimize_angles = true, map<edge_t, double> zero_values_dist=map<edge_t, double>());

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
    s << "{" << (reachable_points+1) << ", " << P.points << ", " << (vector<vector<int> >(P.faces.begin(),P.faces.end())+1) << ", " << static_cast<Graph>(P) << "}";
    return s;
  }

  Polyhedron dual(int Fmax=INT_MAX, bool planar_layout=true) const;
  bool is_triangulation() const;

  double  diameter() const;
  coord3d width_height_depth() const;

  string to_latex(bool show_dual = false, bool number_vertices = false, bool include_latex_header = false) const;
  string to_povray(double w_cm = -1, double h_cm = 10, 
		   int line_colour = 0x888888, int vertex_colour = 0x667744, int face_colour = 0xc03500,
		   double line_width = 0.7, double vertex_diameter = 2.0, double face_opacity = 0.4) const;


  string to_turbomole() const;
  string to_xyz() const;
  string to_mol2() const;
  string to_cc1() const;

  static Polyhedron from_xyz(const string& path);
  static Polyhedron from_mol2(const string& path);
  //  static Polyhedron from_rspi(const string& path); 


private:
  static double C20_points[20][3];
};

#endif
