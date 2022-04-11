#pragma once
#include <utility> 
#include <vector>
#include <string>
#include <iostream>

#include "auxiliary.hh"

class PlanarGraph;		// Circular dependence to planargraph.hh

typedef pair<int,int>  jump_t;
typedef vector<jump_t> jumplist_t;

struct general_spiral {
  jumplist_t  jumps;
  vector<int> spiral;

  bool operator<(const general_spiral &s) const
  {
    return jumps.size() < s.jumps.size() ||
    (jumps.size() == s.jumps.size() && jumps < s.jumps) ||
      (jumps == s.jumps && spiral < s.spiral);
      // The following gives spiral strings precedence over jump content (but still prefers shorter jump lists)
      //    (jumps.size() == s.jumps.size() && spiral < s.spiral) ||
      //      (jumps.size() == s.jumps.size() && spiral == s.spiral && jumps < s.jumps);
  }

  bool operator==(const general_spiral &s) const
  {
    return jumps == s.jumps && spiral == s.spiral;
  }
  
  friend ostream &operator<<(ostream &s, const general_spiral &GS)
  {
    return s << make_pair(GS.jumps,GS.spiral); 
  }
  
};

// Make general spirals hashable
namespace std {
  template <> struct hash<general_spiral> {
    size_t operator()(const general_spiral& S) const {
      size_t seed(0);
      for(const auto &j: S.jumps){
	hash_combine(seed,j.first);
	hash_combine(seed,j.second);
      }
      for(const auto &d: S.spiral)
	hash_combine(seed,d);

      return seed;
    }
  };
}



// Specification:
// (PG-)?[SCHEMES JUMPS SPIRAL]-FORMULA-cage
// (PG-)?[SCHEMES JUMPS INDICES]-FORMULA-fullerene
// (PG-)?[SCHEMES JUMPS INDICES_LIST]-(FACEDEGREES)BASEFACE-FORMULA-fulleroid
//
//
// PG                  ::= String "-"  (* Point group *)
// SCHEMES             ::= EMPTY | (List(',') of SEARCH_SCHEME | CONSTRUCTION_SCHEME) ":"
// SEARCH_SCHEME       ::= "GS" | "CS" | nil      (* "General Spiral" | "Compatibility Spiral"    *)
// CONSTRUCTION_SCHEME ::= "C" | "T" | "LF" | nil (* "Cubic graph" | "Triangulation" | "Leapfrog" *)
// JUMPS               ::= EMPTY | (List of int ) ";"
// SPIRAL              ::= List(',') of int
// INDICES             ::= List(',') of int
// INDICES_LIST        ::= List(';') of INDICES
//
// 0. PG is redundant information. If given, verify after construction.
// 1. If SEARCH_SCHEME is nil, assume "GS" if first index is non-base-face, "CS" otherwise.
// 2. If CONSTRUCTION_SCHEME is nil, assume "C".
// 3. FORMULA is redundant information. If given, verify N is correct after construction.
//
//
// I believe this is not a context-free grammar, and has to be parsed by hand. The following
// is a looser but much simpler grammar (that accepts all correct spiral names plus some more).
//
// PG? '[' SPIRAL_SPEC ']' CAGE_TYPE
// SPIRAL_SPEC ::= (SCHEMES ':'?) List(';') of (List(',') of int)
// CAGE_TYPE ::= List('-') of ('(' (List(',') of int) ')') | FORMULA | ("fullerene" | "fulleroid" | "cage");


struct spiral_nomenclature {
  typedef enum { SS_UNSPECIFIED, CANONICAL_GENERALIZED_SPIRAL, COMPATIBILITY_CANONICAL_SPIRAL } search_scheme_t;
  typedef enum { CS_NONE, CUBIC, TRIANGULATION, LEAPFROG } construction_scheme_t; // -> naming_scheme?
  typedef enum { GT_NONE, FULLERENE, FULLEROID, CAGE } naming_scheme_t; // -> naming_scheme

  naming_scheme_t       naming_scheme;
  search_scheme_t       search_scheme;
  construction_scheme_t construction_scheme;

  string point_group, chemical_formula;
  int base_face_degree;
  vector<int> face_degrees;	// Non-base-face degrees

  // TODO: Change to general_spiral everywhere?
  jumplist_t  jumps;
  vector<int> spiral_code;

  void fulleroid_constructor(const vector<vector<int>> &spiral_numbers, vector<int> face_degrees = {3,4,5},
			     int base_face_degree=6);

  void cage_constructor(const vector<vector<int>> &spiral_numbers);

  spiral_nomenclature(const string &str);
  spiral_nomenclature(const PlanarGraph &G, const naming_scheme_t name_type=CAGE,
		      const construction_scheme_t construction_scheme=CS_NONE,		      
		      bool rarest_special_start = true);

  static string search_scheme_txt[4], construction_scheme_txt[4], naming_scheme_txt[4];

  string to_string(bool unpacked=false) const;
  
  friend ostream& operator<<(ostream& s, const spiral_nomenclature &sn)
  {
    s << sn.to_string();
    return s;
  }

  
};


