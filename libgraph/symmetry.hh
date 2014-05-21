#ifndef FULLERENE_SYMMETRY_HH
# define FULLERENE_SYMMETRY_HH

#include <vector>
#include <iostream> 
#include "triangulation.hh"

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

  friend std::ostream& operator<<(ostream& S, const PointGroup& G){
    S << G.to_string();
    return S;
  }

  bool operator==(const PointGroup& G){
    return sym_type == G.sym_type && n == G.n && sym_reflection == G.sym_reflection;
  }
};


struct Permutation : public vector<int> {
  Permutation(const vector<int>& p) : vector<int>(p){}
  Permutation(int N) : vector<int>(N){}

  static Permutation identity(int N);

  Permutation inverse() const;
  int order() const;

  // Permutation composition
  Permutation operator*(const Permutation& q) const;
  bool operator==(const Permutation& q) const;

};


class Symmetry : public Triangulation {
public:


  vector<int> S0;
  vector< Permutation > G, Gedge, Gtri;

  vector<Permutation> permutation_representation() const;
  vector<Permutation> tri_permutation(const vector<Permutation>& Gf)  const;
  vector<Permutation> edge_permutation(const vector<Permutation>& Gf) const;

  // Returns the involutions *except* from the identity
  vector<int>           involutions() const;
  vector<int>           fixpoints(const Permutation& pi) const;
  vector<int>           group_fixpoints(const vector<Permutation>& pi) const;
  vector<int>           site_symmetry_counts(const vector<Permutation>& pi) const;
  vector< vector<int> > multiplication_table() const ;
  bool                  reverses_orientation(const Permutation& pi) const;

  vector< pair<int,int> > NMR_pattern() const;
  
  PointGroup point_group() const;

  Symmetry(const vector<int>& spiral) : Triangulation(spiral), S0(spiral)
  {
    G = permutation_representation();
    Gedge = edge_permutation(G);
    Gtri  = tri_permutation(G);
  }
  
};

#endif
