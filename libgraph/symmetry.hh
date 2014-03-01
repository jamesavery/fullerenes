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

  PointGroup(symmetry_type t = UNKNOWN, symmetry_reflection r = NONE) 
    : PointGroup(t,0,r) {}
  PointGroup(symmetry_type t, unsigned int n, symmetry_reflection r = NONE) :
    sym_type(t), n(n), sym_reflection(r) {}

  PointGroup(const Triangulation& g);

  static PointGroup FullereneSymmetries[28];

  string to_string() const;

  friend std::ostream& operator<<(ostream& S, const PointGroup& G){
    S << G.to_string();
    return S;
  }
};
