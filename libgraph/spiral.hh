#pragma once
//TODO: The plan is to gather the spiral stuff here, so it's not scattered around the place.

#include <utility> //required for pair
#include <vector>
#include <string>
#include <iostream>

#include "auxiliary.hh"
//#include "triangulation.hh"

typedef pair<int,int>  jump_t;
typedef vector<jump_t> jumplist_t;

#if 1
// TODO: Gather spiral stuff in spiral.hh

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
  
  friend ostream &operator<<(ostream &s, const general_spiral &GS)
  {
    return s << make_pair(GS.jumps,GS.spiral); 
  }
};
#endif




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


struct full_spiral_name {
  typedef enum { SS_UNSPECIFIED, CANONICAL_GENERALIZED_SPIRAL, COMPATIBILITY_CANONICAL_SPIRAL } search_scheme_t;
  typedef enum { CS_NONE, CUBIC, TRIANGULATION, LEAPFROG } construction_scheme_t;
  typedef enum { GT_NONE, FULLERENE, FULLEROID, CAGE } graph_type_t;

  graph_type_t          graph_type;
  search_scheme_t       search_scheme;
  construction_scheme_t construction_scheme;

  string point_group, chemical_formula;
  int base_face_degree;
  vector<int> face_degrees;	// Non-base-face degrees

  jumplist_t  jumps;
  vector<int> spiral_code;

  void fulleroid_constructor(const vector<vector<int>> &spiral_numbers, vector<int> face_degrees = {3,4,5},
			     int base_face_degree=6);

  void cage_constructor(const vector<vector<int>> &spiral_numbers);

  full_spiral_name(const string &str);

  static string search_scheme_txt[4], construction_scheme_txt[4], graph_type_txt[4];

  friend ostream& operator<<(ostream& s, const full_spiral_name &sn)
  {
    s << "{\n\t"
      << "graph_type: \""<<full_spiral_name::graph_type_txt[sn.graph_type]<<"\",\n\t"
      << "search_scheme: \""<<full_spiral_name::search_scheme_txt[sn.search_scheme]<<"\",\n\t"
      << "construction_scheme: \""<<full_spiral_name::construction_scheme_txt[sn.construction_scheme]<<"\",\n\t"
      << "point_group: \""<<(sn.point_group.empty()? "UNSPECIFIED" : sn.point_group) <<"\",\n\t"
      << "chemical_formula: \""<<sn.chemical_formula<<"\",\n\t"
      << "base_face_degree: "<<sn.base_face_degree<<",\n\t"
      << "face_degrees: " << sn.face_degrees << ",\n\t"
      << "jumps: " << sn.jumps << ",\n\t" // indices start counting at 0
      << "spiral_code: " << sn.spiral_code << ",\n\t"
      << "}";
    return s;
  }

};



#if 0

// This general spiral type can represent any polyhedral graph.
struct general_spiral {
  typedef enum { CUBIC, TRIANGULATION, GENERAL } spiral_type;

  vector<int> spiral;  
  jumplist_t  jumps;
  spiral_type type;

  general_spiral(const vector<int>& spiral,
		 const jumplist_t& jumps = jumplist_t(),
		 const polyhedron_type& type) : spiral(spiral),
						jumps(jumps),
						type(type) {}

  general_spiral(const PlanarGraph &g, bool compatibility=false);
  general_spiral(const string& name);
  
  bool operator<(const general_spiral &s) const;
  friend ostream &operator<<(ostream &s, const general_spiral &GS);
};


struct spiral_nomenclature {
  bool is_a_fullerene;
  string atom;
  string point_group;
  
  PlanarGraph   graph;
  Triangulation triangulation;
  bool compatibility;

  general_spiral GS;
  Permutation permutation;

  spiral_nomenclature(const PlanarGraph &g, const string& atom="", bool compatibility=false);
  spiral_nomenclature(const string& name);
  
  friend ostream& operator<<(ostream &s, const spiral_nomenclature &n);

};
#endif


