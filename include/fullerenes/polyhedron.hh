#pragma once

#include <fstream>
#include <sstream>
#include <list>
#include <stdexcept>

#include "fullerenes/owned.hh"
#include "fullerenes/fullerenegraph.hh"
#include "fullerenes/geo-format.hh"
#include "fullerenes/mesh-io-error.hh"

using namespace std;

// Polyhedron: owned planar graph with 3D vertex coordinates.
// Inherits geometry methods from PolyhedronView via Owned<PolyhedronView<double>>.
struct Polyhedron : public Owned<PolyhedronView<double>> {
  using base_t = Owned<PolyhedronView<double>>;
  int face_max = INT_MAX;

  //---- Constructors ----//
  Polyhedron() = default;
  Polyhedron(const int face_max) : face_max(face_max) {}

  // Owning: copies points into owned storage
  Polyhedron(const PlanarGraphView& G, const vector<coord3d>& points_ = vector<coord3d>(), const int face_max = INT_MAX);

  // Owning as well: the coordinates are copied in from the span (an owner never borrows)
  Polyhedron(const PlanarGraphView& G, std::span<coord3d> points_, const int face_max = INT_MAX);

  // Create polyhedron from point collection, assuming shortest distance is approximate bond length
  Polyhedron(const vector<coord3d>& xs, double tolerance = 1.2);

  // Write the coordinates: one per vertex, copied into the owner's storage.
  // @pre  size: pts.size() == size_t(N)
  // @throws std::invalid_argument on a size mismatch
  void set_points(std::span<const coord3d> pts) {
    if (pts.size() != size_t(N))
      throw std::invalid_argument("Polyhedron::set_points: " + std::to_string(pts.size())
                                  + " coordinates for " + std::to_string(N) + " vertices");
    std::copy(pts.begin(), pts.end(), points.begin());
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
  enum {ASCII,PLANARCODE,XYZ,MOL2,MATHEMATICA,LATEX,CC1,TURBOMOLE,GAUSSIAN,WAVEFRONT_OBJ,SPIRAL,PLY,GEO} formats_t;
  static int format_id(string id);

  static Polyhedron from_file(string path);
  static Polyhedron from_file(FILE *file, string format);
  static Polyhedron from_xyz(FILE *file);
  static Polyhedron from_mol2(FILE *file);
  // PLY (ascii or binary_little_endian) -> oriented Polyhedron, normalised to
  // outward-facing (CCW-on-outside) by a signed-volume check + flip if inward.
  // @post   result.is_consistently_oriented() && signed enclosed volume >= 0
  // @throws mesh_io_error  on a null/unsupported/malformed file or invalid topology
  static Polyhedron from_ply(FILE *file);

  static bool to_file(const Polyhedron &G, string path);
  static bool to_file(const Polyhedron &G, FILE *file, string format);
  static bool to_ascii(const Polyhedron &G, FILE *file);
  static bool to_wavefront_obj(const Polyhedron &G, FILE *file);
  static bool to_turbomole(const Polyhedron &G, FILE *file);
  static bool to_gaussian(const Polyhedron &P, FILE *file, string header="");
  static bool to_xyz(const Polyhedron &G, FILE *file);
  static bool to_mol2(const Polyhedron &G, FILE *file);
  // PLY writer (binary_little_endian by default, else ascii). n-gon faces and
  // triangles serialise through the same per-face (count, indices) record.
  // @post   result == (no write error occurred)
  // @throws mesh_io_error  on a null file (NullFile) or a face with > 255 vertices (FaceTooLarge)
  static bool to_ply(const Polyhedron &G, FILE *file, bool binary=true);
  static bool to_cc1(const Polyhedron &G, FILE *file);

  // Compact binary geometry (.geo, GEO-FORMAT.md): fixed-size records of fixed-point
  // or float coordinates and optional connectivity, one seek per record.
  //
  // from_geo: record `index`, with its own graph.
  //   @throws mesh_io_error as geo::read_record; UnsupportedFormat for a record
  //           without coordinates or graph; NonSimplicial for a self-loop or
  //           parallel edge (read those as DelaunayTriangulation)
  // from_geo(file, G, index): record `index` of a file without graph, on G.
  //   @throws as above; std::invalid_argument when the file has a graph or the
  //           record's vertex count differs from G.N
  // to_geo: P as one record, a fresh single-record file, or appended (GEO-FORMAT.md sec. 10).
  //   Rows are written as stored; they must be the counter-clockwise rotation system.
  //   @post   result == no stdio write failed
  //   @throws mesh_io_error NonSimplicial, InvalidTopology, CapacityExceeded,
  //           ValueOutOfRange, HeaderMismatch; std::invalid_argument
  // to_geo(Ps, file, opt): a fresh file with every polyhedron of Ps, the options'
  //   derive fields taken tightly from all of them.
  static geo_header read_geo_header(FILE *file);
  static bool       verify_geo(FILE *file);
  static Polyhedron from_geo(FILE *file, uint64_t index = 0);
  static Polyhedron from_geo(FILE *file, const PlanarGraphView &G, uint64_t index = 0);
  static bool       to_geo(const Polyhedron &P, FILE *file, bool append = false, const geo_options &opt = {});
  static bool       to_geo(std::span<const Polyhedron> Ps, FILE *file, const geo_options &opt = {});

  string to_latex(bool show_dual = false, bool number_vertices = false, bool include_latex_header = false) const;
  string to_povray(double w_cm = -1, double h_cm = 10,
                   int line_colour = 0x888888, int vertex_colour = 0x667744, int face_colour = 0xc03500,
                   double line_width = 0.7, double vertex_diameter = 2.0, double face_opacity = 0.4) const;

  static double dodecahedron_points[20][3];
};
