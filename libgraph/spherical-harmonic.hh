#ifndef SPHERICAL_HARMONIC
# define SPHERICAL_HARMONIC

#include "geometry.hh"
#include "polyhedron.hh"

namespace RealSphericalHarmonic {
  
  struct Ylm_coefficient { int l, m; double coefficient; };
  
  double Y3 (int l, int m, const coord3d& u);

  vector<Ylm_coefficient> decompose_polyhedron(const Polyhedron& P, int Lmax);
}
#endif
