#pragma once

#include "fullerenes/graph.hh"
#include "fullerenes/owned.hh"
#include "fullerenes/spiral.hh"

// PlanarGraph: owned planar graph with oriented embedding.
// Inherits algorithm methods from PlanarGraphView via Owned<PlanarGraphView>.
class PlanarGraph : public Owned<PlanarGraphView> {
public:
  using base_t = Owned<PlanarGraphView>;

  PlanarGraph() = default;
  explicit PlanarGraph(size_t N, int dmax=6) : base_t(int(N), uint8_t(dmax)) {}

  // Deep copy from any GraphView (Graph, PlanarGraph, Triangulation, etc.)
  PlanarGraph(const GraphView& g) : base_t(g) {}

  PlanarGraph(const spiral_nomenclature &fsn);

  // PlanarGraph I/O.
  static vector<string> formats,input_formats, output_formats;
  enum {ASCII,PLANARCODE,XYZ,MOL2,MATHEMATICA,LATEX,SPIRAL} formats_t;
  static int format_id(string id);

  static PlanarGraph from_file(string path, int index=0);
  static PlanarGraph from_file(FILE *file, string format, int index=0);
  static PlanarGraph from_spiral(FILE *file, size_t index=0);
  static PlanarGraph from_ascii(FILE *file);
  static PlanarGraph from_planarcode(FILE *file, const size_t index=0);
  static PlanarGraph from_xyz(FILE *file);
  static PlanarGraph from_mol2(FILE *file);

  static bool to_file(const PlanarGraph &G, string path);
  static bool to_file(const PlanarGraph &G, FILE *file, string format);
  static bool to_spiral(const PlanarGraph &G, FILE *file);
  static bool to_ascii(const PlanarGraph &G, FILE *file);
  static bool to_mathematica(const PlanarGraph &G, FILE *file);
  static bool to_planarcode(const PlanarGraph &G, FILE *file);

  static PlanarGraph read_hog_planarcode(FILE *planarcode_file);
  static vector<PlanarGraph> read_hog_planarcodes(FILE *planarcode_file);
  static bool read_hog_metadata(FILE *file, size_t &graph_count, size_t &graph_size);

  friend ostream& operator<<(ostream& s, const PlanarGraph& g);
};
