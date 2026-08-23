#pragma once

#include <vector>
#include <iostream>

#include "fullerenes/spiral.hh"
#include "fullerenes/triangulation.hh"
#include "fullerenes/geometry.hh"
#include "fullerenes/permutation.hh"

class PointGroup {
public:
  typedef enum { UNKNOWN, C, D, T, S, O, I }      symmetry_type;
  typedef enum { NONE, REF_V, REF_H, REF_D, REF_I, REF_S } symmetry_reflection;

  symmetry_type sym_type;
  unsigned int  n;
  symmetry_reflection sym_reflection;

  PointGroup(symmetry_type t = UNKNOWN, symmetry_reflection r = NONE) : sym_type(t), n(0), sym_reflection(r) {}
  PointGroup(symmetry_type t, unsigned int n, symmetry_reflection r = NONE) :
    sym_type(t), n(n), sym_reflection(r) {}
  PointGroup(const string& name);

  static PointGroup FullereneSymmetries[28];

  string to_string() const;
  // One of the 28 point groups a fullerene can have (Fowler & Manolopoulos).
  bool is_fullerene_group() const;

  friend std::ostream& operator<<(ostream& S, const PointGroup& G){
    S << G.to_string();
    return S;
  }

  bool operator==(const PointGroup& G) const {
    return sym_type == G.sym_type && n == G.n && sym_reflection == G.sym_reflection;
  }
};


// 3D representation of a point group: one orthogonal matrix per group element.
// R[i] is the 3D matrix for the i-th group element (same indexing as Symmetry::G).
// Invariant: R[i]*R[j] == R[k] whenever G[i]*G[j] == G[k].
// det(R[i]) == +1 for orientation-preserving, -1 for orientation-reversing.
struct Representation3D {
  vector<matrix3d> R;
};

class Symmetry : public Triangulation {
public:
  vector<int> S0;
  jumplist_t  J0;
  vector< Permutation > G, Gedge, Garc, Gtri;
  IDCounter<edge_t>   edge_id;
  IDCounter<arc_t> arc_id;

  vector<Permutation> permutation_representation() const;
  vector<Permutation> tri_permutation(const vector<Permutation>& Gf)  const;
  vector<Permutation> edge_permutation(const vector<Permutation>& Gf) const;
  vector<Permutation> arc_permutation(const vector<Permutation>& Gf) const;

  // Returns the involutions *except* from the identity
  vector<int>           involutions() const;
  vector<int>           fixpoints(const Permutation& pi) const;
  vector<int>           group_fixpoints(const vector<Permutation>& pi) const;
  vector<int>           site_symmetry_counts(const vector<Permutation>& pi) const;
  vector< vector<int> > multiplication_table() const ;
  bool                  reverses_orientation(const Permutation& pi) const;

  vector< pair<int,int> > NMR_pattern() const;

  PointGroup point_group() const;

  // @anchor symmetry-equivalence-classes
  // @pre  nonempty: !G.empty() -- always holds (G contains the identity for any
  //                 Symmetry); stated because the body reads G[0].
  // @post result == orbits((int)G[0].size(), G) -- the partition of the vertices
  //                 into orbits under the permutation group G.
  vector<vector<node_t>> equivalence_classes(const vector<Permutation>& G) const;

  // Compute the 3D rotation/reflection matrices for each element of G.
  // Uses the point group type (from point_group()) to generate standard
  // matrices, then matches them to permutations via multiplication table
  // isomorphism.  Coordinate convention: principal Cn axis along z,
  // one C2' axis along x for dihedral groups.
  Representation3D representation_3d() const;
  
  void initialize(){
    for(node_t u=0;u<N;u++){
      auto nu = (*this)[u];
      for(int i=0;i<nu.size();i++){
	const node_t v = nu[i];
	arc_id.insert({u,v});
	edge_id.insert(edge_t{u,v});
      }
    }

    G = permutation_representation();
    Gedge  = edge_permutation(G);
    Garc = arc_permutation(G);
    Gtri  = tri_permutation(G);
  }
  
  Symmetry(const vector<int>& spiral, const jumplist_t& jumps) : Triangulation(spiral,jumps), S0(spiral), J0(jumps)
  {
    initialize();
  }

  Symmetry(const Triangulation& g) : Triangulation(g)
  {
    g.get_spiral(S0,J0);

    initialize();
  }  
  
};


