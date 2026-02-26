#pragma once
#include "fullerenes/triangulation.hh"
#include "fullerenes/geometry.hh"

class Polyhedron;  // forward declaration

class Deltahedron : public Triangulation {
public:
  vector<coord3d> points;

  // Constructors
  Deltahedron() = default;
  Deltahedron(const Triangulation& T, const vector<coord3d>& points);
  Deltahedron(const Polyhedron& P);  // must be a triangulation

  // Faces from triangles on the fly (no storage)
  vector<face_t> compute_dual_faces() const;

  // Laplacian smoothing
  void smooth(double q);

  // GC transform with 3D coordinates (l=0 for now)
  Deltahedron GCtransform(unsigned k, unsigned l=0) const;
};
