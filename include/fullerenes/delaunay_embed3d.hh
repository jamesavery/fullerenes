#ifndef DELAUNAY_EMBED3D_HH
#define DELAUNAY_EMBED3D_HH

#include "delaunay.hh"

// Paired permutations + O(3) matrices for symmetry-constrained embedding.
struct SymmetryConstraint {
  vector<vector<int>> perms;
  vector<matrix3d> matrices;
  bool empty() const { return perms.empty(); }
};

// Embed the reduced triangulation in 3D to match intrinsic edge lengths.
// MDS initialization + trust-region Steihaug CG optimization of
// E_edge (distance matching) + E_cone (cone angle matching).
// When symmetry is provided, uses Reynolds-projected MDS for initialization.
vector<coord3d> embed_delaunay_3d(const DelaunayTriangulation& D,
                                   const SymmetryConstraint& sym = {});

// Restrict full-triangulation symmetry (permutations + O(3) matrices) to
// iDT cone-point indices, maintaining the perm<->matrix pairing.
// G: automorphisms (Symmetry::G), R: 3D matrices (Representation3D::R).
SymmetryConstraint restrict_symmetry_to_cone_points(
    const vector<vector<int>>& G, const vector<matrix3d>& R,
    const Triangulation& T);

// Compute vertex orbits from a group of permutations (union-find).
vector<vector<int>> compute_orbits(int n, const vector<vector<int>>& G);

#endif
