#include "fullerenes/planargraph.hh"
#include "fullerenes/polyhedron.hh"
#include "fullerenes/geometry.hh"

using namespace std;

// General force-field optimization of non-fullerene polyhedra (arbitrary
// cubic graphs and triangulations).  This was historically a GSL
// Polak-Ribiere conjugate-gradient optimizer over a uniform harmonic
// field; the GSL path was never actually built (the GSL_INCLUDE_DIR
// preprocessor macro was never defined by the build, so this stub was
// what compiled), and GSL has been removed as a dependency -- all dense
// linear algebra goes through the BLAS-free dense_linalg module.
//
// Fullerene graphs do NOT reach this method: PolyhedronView::optimize
// routes them to FullereneGraphView::optimized_geometry (the C++
// Wu/ExtWu force field, wu_forcefield.hh).  Only non-fullerene cubic
// graphs / triangulations / the leapfrog-dual fallback land here.
//
// TODO(optimizer-framework): re-provide this as the framework's general
// harmonic EnergyModel (uniform bond/angle/dihedral constants) driven by
// the shared minimizer, so non-fullerene polyhedra optimize through the
// same path as everything else.  Until then this reports "unavailable"
// rather than silently returning an unoptimized geometry.
template<>
bool PolyhedronView<double>::optimize_other(bool, map<edge_t, double>)
{
  cerr << "PolyhedronView::optimize_other: general (non-fullerene) polyhedron "
          "optimization is not currently available (awaiting the unified "
          "optimizer framework)." << endl;
  return false;
}
