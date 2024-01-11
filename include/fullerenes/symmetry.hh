#pragma once

#include <vector>
#include <iostream> 

#include "fullerenes/spiral.hh"
#include "fullerenes/triangulation.hh"

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
  Permutation(int N=0) : vector<int>(N){}

  static Permutation identity(int N);

  Permutation inverse() const;
  int order() const;

  // Permutation composition
  Permutation operator*(const Permutation& q) const;
  //  bool operator==(const Permutation& q) const;

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
  vector<vector<node_t>> equivalence_classes(const vector<Permutation>& G) const;
  
  void initialize(){
    for(node_t u=0;u<N;u++){
      const vector<node_t> &nu = neighbours[u];
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


