#include "libgraph/triangulation.hh"
#include <string.h>

string trim(const string& str, const string& wschars)
{
    size_t first = str.find_first_not_of(wschars);
    if(first == string::npos)
    {
        return str;
    }
    size_t last = str.find_last_not_of(wschars);

    return str.substr(first, (last - first + 1));
}

template <typename T> vector<T> split(const string& parse_str, const string& delimiters, const string wschars=" \t\r\n")
{
  vector<string> string_result = split<string>(parse_str,delimiters,wschars);
  vector<T> result;

  for(string s: string_result) result.push_back(from_string<T>(s));

  return result;
}

template <> vector<string> split(const string& parse_str, const string& delimiters, const string wschars)
{
  vector<string> result;
  const char *del_pt = delimiters.c_str();
  char s[parse_str.size()+1];
  strncpy(s,parse_str.c_str(),parse_str.size()+1);
  char *pt = 0;

  //  cerr << "Split \"" << s << "\" on \"" << del_pt << "\" ";
  const char *next = strtok_r(s, del_pt, &pt);

  while(next != 0){
    result.push_back(trim(next,wschars));
    next = strtok_r(0,del_pt,&pt);
  }
  //  cerr << "yields " << result << ".\n";
  return result;
}

vector<string> find_parenthetical(const string& parse_str, const char parentheses[2])
{
  size_t start = parse_str.find(parentheses[0], 0);
  if (start == string::npos) return vector<string>();
  size_t end   = parse_str.find(parentheses[1], start);
  if (end == string::npos) return vector<string>();
  
  return vector<string>{{parse_str.substr(0, start), parse_str.substr(start+1, end-1), parse_str.substr(end+1,parse_str.size()-1)}};
}

// #include "libgraph/spiral.hh"
// #include "libgraph/spiral.cc"

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
  typedef enum { CANONICAL_GENERALIZED_SPIRAL, COMPATIBILITY_CANONICAL_SPIRAL } search_scheme_t;
  typedef enum { CUBIC, TRIANGULATION, LEAPFROG } construction_scheme_t;
  typedef enum { FULLERENE, FULLEROID, CAGE } graph_type_t;

  graph_type_t          graph_type;
  search_scheme_t       search_scheme;
  construction_scheme_t construction_scheme;  

  string point_group;
  int base_face_degree;
  vector<int> face_degrees;	// Non-base-face degrees

  vector<int>      jumps;
  vector<int>      spiral_code;

  void fulleroid_constructor(const vector<vector<int>> &spiral_numbers, vector<int> face_degrees = {3,4,5},
			     int base_face_degree=6);  
  
  void cage_constructor(const vector<vector<int>> &spiral_numbers);

  full_spiral_name(const string &str);

  friend ostream& operator<<(ostream& s, const full_spiral_name &sn)
  {
    s << "{\n\t"
      << "graph_type: \""<<(sn.graph_type==FULLERENE?"FULLERENE":
			  (sn.graph_type==FULLEROID?"FULLEROID":
			   (sn.graph_type==CAGE?"CAGE" : "null")))<<"\",\n\t"
      << "search_scheme: \""<<(sn.search_scheme==CANONICAL_GENERALIZED_SPIRAL?"CANONICAL_GENERALIZED_SPIRAL":
			     (sn.search_scheme==COMPATIBILITY_CANONICAL_SPIRAL?
			      "COMPATIBILITY_CANONICAL_SPIRAL":"null"))<<"\",\n\t"
      << "construction_scheme: \""<<(sn.construction_scheme==CUBIC?"CUBIC":
				   (sn.construction_scheme==TRIANGULATION?"TRIANGULATION":
				    (sn.construction_scheme==LEAPFROG?"LEAPFROG" : "null")))<<"\",\n\t"
      << "point_group: \""<<sn.point_group<<"\",\n\t"
      << "base_face_degree: "<<sn.base_face_degree<<",\n\t"
      << "face_degrees: " << sn.face_degrees << ",\n\t"
      << "jumps: " << sn.jumps << ",\n\t"
      << "spiral_code: " << sn.spiral_code << ",\n\t"
      << "}";
    return s;
  }
};

full_spiral_name::full_spiral_name(const string &str) : graph_type(CAGE), search_scheme(CANONICAL_GENERALIZED_SPIRAL),
							construction_scheme(CUBIC), point_group("null"),
							base_face_degree(6), face_degrees({5})
{
  vector<string> segments = find_parenthetical(str,"[]");
  cerr << "segments = " << segments << "\n";
    
  if(segments.size() != 3){
    cerr << "Error in spiral string \"" << str << "\": Number of \"[]\"-delimited segments is " << segments.size()
	 << ", not 3.\n";
    abort();
  }
    
  const string &prefix_string = segments[0], &spiral_spec = segments[1], &suffix_string = segments[2];

  //------------------------------ Parse Prefix ------------------------------
  // Is the point group specified?
  if(!segments[0].empty()) point_group = trim(prefix_string,"-– \t\r\n");
  // That's all we have stuck into the prefix.
    
  // ---------------------- Parse Spiral Specification -----------------------
  // TODO: Make runlength-coding parser (a,b,c)n -> a,b,c,a,b,c,...,a,b,c, and 
  //       parse number segments in case runlength-coding is used
  vector<string> schemes_and_spiral = split<string>(spiral_spec,":");
  assert(schemes_and_spiral.size() <= 2);
    
  bool schemes_present = schemes_and_spiral.size() == 2;

  // cerr << "schemes=\"" << (schemes_present? schemes_and_spiral[0] : "") <<"\";\n";
  // cerr << "spiral =\"" << schemes_and_spiral[schemes_present] <<"\";\n";    
  vector<string> spiral_segments  = split<string>(schemes_and_spiral[schemes_present],";");

  // Construction- and canonical-spiral-search-schemes are stuck inside the spiral spec
  // Graph type is in suffix, parsed later
  if(schemes_present){
    vector<string> schemes = split<string>(schemes_and_spiral[0],",");

    for(auto s: schemes){
      // Default construction scheme is CUBIC
      if(s == "LF") construction_scheme = LEAPFROG;
      if(s == "T")  construction_scheme = TRIANGULATION;

      // Default search scheme is GS
      if(s == "CS") search_scheme = COMPATIBILITY_CANONICAL_SPIRAL;
    }
  }  
  
  // Get spiral numbers: Optional jumps followed by either spiral face-degree string (spiral_code) or non-base-face indices
  vector<vector<int>> spiral_numbers(spiral_segments.size());    
  for(int i=0; i<spiral_segments.size(); i++){
    spiral_numbers[i] = split<int>(spiral_segments[i],",");
  }

  //------------------------------ Parse SUFFIX ------------------------------
  vector<string> suffix_segments = split<string>(suffix_string,"-–");
  string suffix = suffix_segments.empty()? "cage" : suffix_segments.back();

  // General cage
  if (suffix == "cage"){
    graph_type = CAGE;
    cage_constructor(spiral_numbers);
    return;
  } 

  if(!(suffix == "fullerene" || suffix == "fulleroid")){
    cerr << "Graph type is \""<<suffix<<"\", must be one of \"cage\", \"fullerene\", or \"fulleroid\".\n";
    abort();
  }
    
  if (suffix == "fullerene"){
    graph_type = FULLERENE;
    base_face_degree = 6;
    face_degrees     = vector<int>{{5}};
  }

  if(suffix == "fulleroid"){
    graph_type = FULLEROID;
    // TODO: Parse fulleroid specification
  }

  fulleroid_constructor(spiral_numbers,face_degrees,base_face_degree);
}

void full_spiral_name::cage_constructor(const vector<vector<int>> &spiral_numbers)
{
  assert(spiral_numbers.size() == 1 || spiral_numbers.size() == 2);
  bool has_jumps = (spiral_numbers.size() == 2);

  if(has_jumps)
    jumps = spiral_numbers[0];
  spiral_code = spiral_numbers[has_jumps];
}

void full_spiral_name::fulleroid_constructor(const vector<vector<int>> &spiral_numbers, vector<int> face_degrees, int base_face_degree)
{
  assert(spiral_numbers.size() ==  face_degrees.size() || spiral_numbers.size() == face_degrees.size()+1);

  int max_index = 0;
  for(const auto &v: spiral_numbers)
    for(const auto &ix: v)
      max_index = max(ix,max_index);

  spiral_code = vector<int>(2*max_index,6);

  bool has_jumps = (spiral_numbers.size() == face_degrees.size()+1);

  if(has_jumps)
    jumps = spiral_numbers[0];
  
  for(int i=0;i<face_degrees.size();i++)
    for(auto ix: spiral_numbers[i+has_jumps])
      spiral_code[ix-1] = face_degrees[i]; 
}




int main(int ac, char **av)
{
  full_spiral_name
    tutte_clean_name("[CS:11, 1, 17, 1; 5, 10, 5, 5, 5, 9, 5, 4, 5, 4, 4, 5, 4, 10, 5, 5, 5, 5, 5, 10, 4, 5, 5, 4, 5]"),
    tutte_name("C3–[CS: 11, 1, 17, 1; 5, 10, 5, 5, 5, 9, 5, 4, 5, 4, 4, 5, 4, 10, 5, 5, 5, 5, 5, 10, 4, 5, 5, 4, 5]-C46-cage");

  cout << "tutte_name = " << tutte_name << ";\n";
  cout << "tutte_clean_name = " << tutte_clean_name << ";\n";
  return 0;
}
