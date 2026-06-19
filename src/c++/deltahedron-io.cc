// Deltahedron mesh I/O. A Deltahedron is the triangulation case of a Polyhedron,
// so both readers/writers route through Polyhedron::{from,to}_ply (the raw PLY
// serialisation + the graph-from-faces builder live in polyhedron-io.cc). This
// TU mirrors the polyhedron.cc / polyhedron-io.cc split: geometry I/O is kept
// out of the optimizer-heavy deltahedron.cc.

#include "fullerenes/deltahedron.hh"
#include "fullerenes/polyhedron.hh"

#include <vector>

Deltahedron Deltahedron::from_ply(FILE *file)
{
  // Build the oriented (outward-normalised) Polyhedron, then require it to be a
  // triangulation. The Deltahedron(const Polyhedron&) ctor asserts the same, so
  // checking here turns a malformed input into a catchable mesh_io_error rather
  // than an abort.
  Polyhedron P = Polyhedron::from_ply(file);
  if (!P.is_triangulation())
    throw mesh_io_error(mesh_io_error::Code::NotATriangulation, "Deltahedron::from_ply: PLY mesh is not a triangulation");
  return Deltahedron(P);
}

bool Deltahedron::to_ply(const Deltahedron &D, FILE *file, bool binary)
{
  // A Deltahedron is a Polyhedron whose faces happen to be triangles; wrap the
  // view as a Polyhedron (copies graph + points) and reuse its writer.
  Polyhedron P(static_cast<const PlanarGraphView&>(D),
               std::vector<coord3d>(D.points.begin(), D.points.end()));
  return Polyhedron::to_ply(P, file, binary);
}
