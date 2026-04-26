#pragma once

#include <vector>
#include <functional>

using std::vector;

// A permutation pi of {0, ..., N-1}.  Element-access semantics:
//   pi[u_old] == u_new
// I.e. pi maps an "old" label to its "new" label.  pi.inverse()[u_new]
// gives the old label that ended up at u_new.
//
// Composition: (pi * q)[u] == q[pi[u]] (apply pi first, then q).
//
// Defined in its own header so non-symmetry consumers (e.g. Graph
// primitives like argsort/apply_permutation) can use it without
// pulling in the symmetry/triangulation headers.
struct Permutation : public vector<int> {
  Permutation(const vector<int>& p) : vector<int>(p){}
  Permutation(int N=0) : vector<int>(N){}

  static Permutation identity(int N);

  Permutation inverse() const;
  int order() const;

  // Permutation composition
  Permutation operator*(const Permutation& q) const;
};

namespace std {
  template<> struct hash<Permutation> {
    size_t operator()(const Permutation &v) const {
      return std::hash<vector<int>>()(v);
    }
  };
}
