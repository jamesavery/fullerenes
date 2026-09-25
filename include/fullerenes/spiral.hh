#pragma once
#include <utility> 
#include <vector>
#include <string>
#include <iostream>
#include <sstream>
#include <stdexcept>

#include "auxiliary.hh"

class PlanarGraph;		// Circular dependence to planargraph.hh
struct TriangulationView;	// ... and to graphview.hh

typedef pair<int,int>  jump_t;
typedef vector<jump_t> jumplist_t;

struct general_spiral {
  jumplist_t  jumps;
  vector<int> spiral_code;

  bool operator<(const general_spiral &s) const
  {
    return jumps.size() < s.jumps.size() ||
    (jumps.size() == s.jumps.size() && jumps < s.jumps) ||
      (jumps == s.jumps && spiral_code < s.spiral_code);
      // The following gives spiral strings precedence over jump content (but still prefers shorter jump lists)
      //    (jumps.size() == s.jumps.size() && spiral < s.spiral) ||
      //      (jumps.size() == s.jumps.size() && spiral == s.spiral && jumps < s.jumps);
  }

  bool operator==(const general_spiral &s) const
  {
    return jumps == s.jumps && spiral_code == s.spiral_code;
  }
  
  friend ostream &operator<<(ostream &s, const general_spiral &GS)
  {
    return s << make_pair(GS.jumps,GS.spiral_code); 
  }

};

// Thrown by the spiral search when it finds no spiral or cannot search.
//   reason   -- NO_SPIRAL: the general-spiral search found no spiral.  A
//               general spiral (one that may jump past cut vertices) closes
//               from every starting triple of every oriented triangulation of
//               the sphere, so this never reports a legitimate outcome: it
//               signals a malformed input (not a sphere triangulation, not
//               oriented) or a bug in the search.
//               DEGREE_LIMIT: the triangulation has a vertex of degree above
//               degree_limit, the width of the search's per-vertex bit mask
//               of remaining neighbours; the search refuses it before starting.
//   N        -- vertex count of the triangulation searched
//   start    -- NO_SPIRAL: the starting triple whose general spiral did not
//               close, or {-1,-1,-1} when the search had no starting triple
//               to try; DEGREE_LIMIT: {-1,-1,-1}
//   only_rarest_special, CW_only -- NO_SPIRAL: the search parameters;
//               DEGREE_LIMIT: false
//   max_degree -- DEGREE_LIMIT: the largest vertex degree found; NO_SPIRAL: -1
struct SpiralSearchFailed : std::runtime_error {
  enum reason_t { NO_SPIRAL, DEGREE_LIMIT };
  static constexpr int degree_limit = 16;

  reason_t reason;
  int N;
  node_t start[3];
  bool only_rarest_special, CW_only;
  int max_degree;

  SpiralSearchFailed(int N, node_t f1, node_t f2, node_t f3,
                     bool only_rarest_special, bool CW_only)
    : std::runtime_error(message(N, f1, f2, f3, only_rarest_special, CW_only)),
      reason(NO_SPIRAL), N(N), start{f1, f2, f3},
      only_rarest_special(only_rarest_special), CW_only(CW_only), max_degree(-1) {}

  // The DEGREE_LIMIT refusal of an N-vertex triangulation of maximal degree max_degree.
  static SpiralSearchFailed degree_exceeded(int N, int max_degree)
  {
    return SpiralSearchFailed(N, max_degree);
  }

private:
  SpiralSearchFailed(int N, int max_degree)
    : std::runtime_error("Spiral search refused a " + std::to_string(N) + "-vertex triangulation: "
                         "a vertex has degree " + std::to_string(max_degree) + ", above the limit "
                         + std::to_string(degree_limit)),
      reason(DEGREE_LIMIT), N(N), start{-1, -1, -1},
      only_rarest_special(false), CW_only(false), max_degree(max_degree) {}

  static string message(int N, node_t f1, node_t f2, node_t f3,
                        bool only_rarest_special, bool CW_only)
  {
    ostringstream s;
    s << "General spiral search failed on a " << N << "-vertex triangulation ("
      << (only_rarest_special ? "rarest-special starts" : "all starts")
      << (CW_only ? ", clockwise only" : ", both orientations") << "): ";
    if(f1 < 0) s << "no starting triple to try";
    else       s << "no closing spiral from start (" << f1 << "," << f2 << "," << f3 << ")";
    return s.str();
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
      for(const auto &d: S.spiral_code)
	hash_combine(seed,d);

      return seed;
    }
  };
}



// Thrown by the name parser, spiral_nomenclature(const string&), when the
// string is not a name: it breaks the grammar below, or its numbers cannot be
// the spiral of a graph of the named type.  what() quotes the input and states
// the defect; `defect` classifies it.
struct SpiralNameMalformed : std::invalid_argument {
  enum defect_t {
    BRACKETS,       // not exactly one "[" followed by one "]"
    SCHEME_TAG,     // an unknown or repeated scheme tag, or a ":" with no tag
    NUMBER,         // a list entry that is not a decimal non-negative int
    SEGMENT_COUNT,  // the number of ";"-separated lists does not fit the graph type
    JUMPS,          // an odd jump list, a position below 3 or not increasing,
                    // a rotation count below 1, a position past a cage spiral
    INDEX_RANGE,    // a face index below 1, or too large for the spiral length
    INDEX_ORDER,    // face indices not strictly increasing, or one index under two degrees
    INDEX_COUNT,    // not 12 pentagon indices, or fulleroid faces failing Euler's relation
    DEGREES,        // a degree below 3, cage degrees failing Euler's relation, or a
                    // fulleroid degree list that is empty, repeated, or holds the base degree
    SUFFIX          // a graph type other than cage/fullerene/fulleroid, an empty or
                    // extra "-" segment, a malformed face-degree group or base degree
  };
  string name;
  defect_t defect;

  SpiralNameMalformed(const string& name, defect_t defect, const string& reason)
    : std::invalid_argument("spiral_nomenclature: malformed name \"" + name + "\": " + reason),
      name(name), defect(defect) {}
};

// The name grammar.
//
//   NAME     ::= (PG "-")? "[" (TAGS ":")? (JUMPS ";")? NUMBERS "]" SUFFIX
//   SUFFIX   ::= EMPTY                                   (a cage)
//              | "-"? (QUALIFIER "-")* TYPE
//   QUALIFIER::= FORMULA | FACES   (at most one FORMULA; one FACES exactly for a
//                                   fulleroid, none otherwise; either order)
//   TYPE     ::= "cage" | "fullerene" | "fulleroid"
//   FACES    ::= "(" List(',') of int ")" ("_"? "6")?  (the special face degrees; base 6)
//   TAGS     ::= List(',') of TAG, at most one construction and one search tag
//   TAG      ::= "C" | "T" | "LF"                      (construction: cubic, triangulation, leapfrog)
//              | "GS" | "CS"                           (search: general, compatibility spiral)
//   JUMPS    ::= List(',') of int: pairs (position, rotations), 1-based positions
//                >= 3 and strictly increasing, rotations >= 1
//   NUMBERS  ::= cage:      List(',') of int, the degree sequence of the spiral
//              | fullerene: List(',') of int, the 12 pentagon positions
//              | fulleroid: List(';') of (List(',') of int), the positions of the
//                           faces of each degree in FACES, in that order
//
// Every int is a decimal numeral without sign; face positions are 1-based and
// strictly increasing.  "-" may also be written as an en dash.
//
// 0. PG (point group) and FORMULA are redundant information, recorded but not
//    checked against the spiral.
// 1. The spiral is always a vertex spiral of a triangulation, the ENVELOPING
//    triangulation of the named graph G (PlanarGraphView::enveloping_triangulation):
//    G itself when G is a triangulation ("T"), G's dual when G is cubic (no tag,
//    or "C", which is never written), G's leapfrog dual otherwise ("LF").  The
//    construction scheme therefore says WHICH graph the name describes:
//    "[1,7,...]-fullerene" names the cubic fullerene (its face spiral),
//    "[T:1,7,...]-fullerene" names its dual triangulation (the same spiral,
//    read as a vertex spiral).  The two are different graphs with the same
//    spiral numbers.
// 2. The search scheme says which canonical spiral is written.
//    CANONICAL_GENERALIZED_SPIRAL, the minimal spiral over the starts at
//    vertices of degree != 6, is the default and is never written; a name
//    without a search tag parses to it.  COMPATIBILITY_CANONICAL_SPIRAL, the
//    minimal spiral over all starts, is always written, "CS".  Both rules hold
//    with and without jumps; jumps are recognised by the form
//    "<jumps>; <numbers>" alone.  A name thus records exactly what a reader
//    needs to rebuild the graph and to know which spiral it holds.
// 3. Legacy forms.  The parser also accepts "GS" (the default, redundant) on
//    any spiral, "[GS:1,7,...]-fullerene" naming the same graph as
//    "[1,7,...]-fullerene", and "[T,GS:1,7,...]-fullerene", which parses as
//    the "T" name it spells: a record in that form that was meant as the
//    fullerene name is read as the name of its dual triangulation, and
//    FullereneGraph(sn) recovers the fullerene.  Records from earlier writers
//    also hold untagged jump-free spirals from the all-starts search; for a
//    fullerene these coincide with the default search's spiral (the minimal
//    all-starts spiral starts at a pentagon whenever a pentagon start closes
//    without jumps), for a general cage they need not.
// 4. The parser refuses, with SpiralNameMalformed naming the defect, every
//    string outside this grammar and every spiral that cannot belong to a
//    triangulated sphere of the named type by counting alone: a fullerene
//    lists exactly 12 pentagons; a cage's degrees d >= 3 and fulleroid faces
//    around hexagons satisfy Euler's relation sum(6 - d) = 12.  Whether the
//    spiral winds up (closes) is decided only by building the graph.


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
  general_spiral spiral;
  
  void fulleroid_constructor(const vector<vector<int>> &spiral_numbers, vector<int> face_degrees = {3,4,5},
			     int base_face_degree=6);

  void cage_constructor(const vector<vector<int>> &spiral_numbers);

  // Parse a name (the grammar above, legacy forms included).
  // @throws SpiralNameMalformed when the string is not a name (4. above).
  // Nothing is constructed: build the graph with Triangulation(sn) (always the
  // enveloping triangulation, i.e. the DUAL of a named cubic graph and the named
  // graph itself for a "T" name), FullereneGraph(sn) (its dual) or
  // PlanarGraph(sn) (the named graph).
  spiral_nomenclature(const string &str);

  // Names the graph G: construction_scheme says how G
  // relates to the triangulation whose vertex spiral is written (1. above).
  // construction_scheme == CS_NONE lets G decide: TRIANGULATION if G is a
  // triangulation, CUBIC if G is cubic, LEAPFROG otherwise; the member records
  // that decision.  Any other value is recorded as given (the caller asserts
  // it; enveloping_triangulation does not check it).
  // @throws SpiralSearchFailed if the general-spiral search on G's
  //         enveloping triangulation finds no spiral (see its definition).
  spiral_nomenclature(const PlanarGraph &G, const naming_scheme_t name_type=CAGE,
		      const construction_scheme_t construction_scheme=CS_NONE,
		      bool rarest_special_start = true);

  // Names the FULLERENE (the cubic graph) whose dual triangulation is
  // `dual`.  The dual's canonical vertex spiral is the fullerene's canonical face
  // spiral, so the spiral is searched on `dual` directly (no dualization) and
  // the name carries no "T".  Ih-C60, for example, is named
  // "[1,7,9,11,13,15,18,20,22,24,26,32]-fullerene"; a spiral that needs jumps
  // is written "[<jumps>; ...]-fullerene" (2. above).
  // This is the library's canonical fullerene name.
  // @pre  dual is oriented and dual.is_fullerene_dual() (degrees 5 and 6,
  //       twelve 5s); not checked here (the name of a non-fullerene triangulation would carry
  //       the "-fullerene" suffix regardless).
  // @throws SpiralSearchFailed as the constructor above.
  static spiral_nomenclature fullerene_from_dual(const TriangulationView &dual,
                                                 bool rarest_special_start = true);

  static string search_scheme_txt[4], construction_scheme_txt[4], naming_scheme_txt[4];

  // Writes the name string: the construction scheme unless CUBIC (1. above),
  // the search scheme unless CANONICAL_GENERALIZED_SPIRAL (2. above).
  string to_string(bool unpacked=false) const;

private:
  spiral_nomenclature(naming_scheme_t naming_scheme, construction_scheme_t construction_scheme,
                      bool rarest_special_start);
  // Sets spiral to T's canonical general spiral and face_degrees to its
  // non-base degrees.
  void name_triangulation(const TriangulationView &T, bool rarest_special_start);
public:
  
  friend ostream& operator<<(ostream& s, const spiral_nomenclature &sn)
  {
    s << sn.to_string();
    return s;
  }

  
};


