#pragma once

#include <fstream>
#include <sstream>
#include <list>

#include "fullerenes/owned.hh"
#include "fullerenes/fullerenegraph.hh"

using namespace std;

// Polyhedron: owned planar graph with 3D vertex coordinates.
// Inherits geometry methods from PolyhedronView via Owned<PolyhedronView>.
struct Polyhedron : public Owned<PolyhedronView> {
  using base_t = Owned<PolyhedronView>;
  int face_max = INT_MAX;

  //---- Constructors ----//
  Polyhedron() = default;
  Polyhedron(const int face_max) : face_max(face_max) {}

  // Owning: copies points into owned storage
  Polyhedron(const PlanarGraphView& G, const vector<coord3d>& points_ = vector<coord3d>(), const int face_max = INT_MAX);

  // View: uses external coordinate memory (caller manages lifetime)
  Polyhedron(const PlanarGraphView& G, std::span<coord3d> points_, const int face_max = INT_MAX);

  // Create polyhedron from point collection, assuming shortest distance is approximate bond length
  Polyhedron(const vector<coord3d>& xs, double tolerance = 1.2);

  // Replace owned coordinate storage and repoint the span.
  void set_points(std::vector<coord3d> pts) {
    owned_points = std::move(pts);
    repoint();
  }

  Polyhedron convex_hull() const { return incremental_convex_hull(); }

  bool is_triangulation() const;

  static Polyhedron fullerene_polyhedron(FullereneGraph G);

  static Polyhedron C20() {
    constexpr double bond_length = 1.45;
    vector<coord3d> pts(20);
    for(node_t u=0;u<20;u++)
      pts[u] = coord3d(dodecahedron_points[u][0],dodecahedron_points[u][1],dodecahedron_points[u][2]) * bond_length;
    return Polyhedron(FullereneGraph::C20(),pts);
  }

  friend ostream& operator<<(ostream& s, const Polyhedron& P){
    vector<node_t> reachable_points;
    for(node_t u=0;u<P.N;u++) if(P.degree(u)!=0) reachable_points.push_back(u);
    auto fs = P.faces(P.face_max);
    s << "{" << (reachable_points+1) << ", " << P.points << ", " << (vector<vector<int> >(fs.begin(),fs.end())+1) << "}";
    return s;
  }

  // Graph I/O
  static vector<string> formats,format_alias, input_formats, output_formats;
  enum {ASCII,PLANARCODE,XYZ,MOL2,MATHEMATICA,LATEX,CC1,TURBOMOLE,GAUSSIAN,WAVEFRONT_OBJ,SPIRAL} formats_t;
  static int format_id(string id);

  static Polyhedron from_file(string path);
  static Polyhedron from_file(FILE *file, string format);
  static Polyhedron from_xyz(FILE *file);
  static Polyhedron from_mol2(FILE *file);

  static bool to_file(const Polyhedron &G, string path);
  static bool to_file(const Polyhedron &G, FILE *file, string format);
  static bool to_ascii(const Polyhedron &G, FILE *file);
  static bool to_wavefront_obj(const Polyhedron &G, FILE *file);
  static bool to_turbomole(const Polyhedron &G, FILE *file);
  static bool to_gaussian(const Polyhedron &P, FILE *file, string header="");
  static bool to_xyz(const Polyhedron &G, FILE *file);
  static bool to_mol2(const Polyhedron &G, FILE *file);
  static bool to_cc1(const Polyhedron &G, FILE *file);

  string to_latex(bool show_dual = false, bool number_vertices = false, bool include_latex_header = false) const;
  string to_povray(double w_cm = -1, double h_cm = 10,
                   int line_colour = 0x888888, int vertex_colour = 0x667744, int face_colour = 0xc03500,
                   double line_width = 0.7, double vertex_diameter = 2.0, double face_opacity = 0.4) const;

  static double dodecahedron_points[20][3];
};
