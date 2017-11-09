#include "libgraph/triangulation.hh"
#include <string.h>

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


vector<string> find_parenthetical(const string& parse_str, const char parentheses[2])
{
  size_t start = parse_str.find(parentheses[0], 0);
  if (start == string::npos) return vector<string>();
  size_t end   = parse_str.find(parentheses[1], start);
  if (end == string::npos) return vector<string>();
  
  return vector<string>{{parse_str.substr(0, start),
                         parse_str.substr(start+1, end-start-1),
                         parse_str.substr(end+1,parse_str.size()-end-1)}};
}

string trim(const string& str, const string& wschars)
{
    size_t first = str.find_first_not_of(wschars);
    if(first == string::npos)
      return "";
    size_t last = str.find_last_not_of(wschars);

    return str.substr(first, (last - first + 1));
}

template <> vector<string> split(const string& parse_str, const string& delimiters, const string wschars)
{
  vector<string> result;
  const char *del_pt = delimiters.c_str();
  char s[parse_str.size()+1];
  strncpy(s,parse_str.c_str(),parse_str.size()+1);
  char *pt = 0;

  const char *next = strtok_r(s, del_pt, &pt);

  while(next != 0){
    result.push_back(trim(next,wschars));
    next = strtok_r(0,del_pt,&pt);
  }
  return result;
}



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

  vector<int>      jumps;
  vector<int>      spiral_code;

  void fulleroid_constructor(const vector<vector<int>> &spiral_numbers, vector<int> face_degrees = {3,4,5},
			     int base_face_degree=6);  
  
  void cage_constructor(const vector<vector<int>> &spiral_numbers);

  full_spiral_name(const string &str);

  static string search_scheme_txt[4], construction_scheme_txt[4], graph_type_txt[4];
};

string full_spiral_name::search_scheme_txt[4]       = {"UNSPECIFIED","CANONICAL_GENERALIZED_SPIRAL","COMPATIBILITY_CANONICAL_SPIRAL"};
string full_spiral_name::construction_scheme_txt[4] = {"null","CUBIC","TRIANGULATION", "LEAPFROG"};
string full_spiral_name::graph_type_txt[4]          = {"null","FULLERENE", "FULLEROID", "CAGE"};

ostream& operator<<(ostream& s, const full_spiral_name &sn)
{
  s << "{\n\t"
    << "graph_type: \""<<full_spiral_name::graph_type_txt[sn.graph_type]<<"\",\n\t"
    << "search_scheme: \""<<full_spiral_name::search_scheme_txt[sn.search_scheme]<<"\",\n\t"
    << "construction_scheme: \""<<full_spiral_name::construction_scheme_txt[sn.construction_scheme]<<"\",\n\t"
    << "point_group: \""<<sn.point_group<<"\",\n\t"
    << "chemical_formula: \""<<sn.chemical_formula<<"\",\n\t"
    << "base_face_degree: "<<sn.base_face_degree<<",\n\t"
    << "face_degrees: " << sn.face_degrees << ",\n\t"
    << "jumps: " << sn.jumps << ",\n\t"
    << "spiral_code: " << sn.spiral_code << ",\n\t"
    << "}";
  return s;
}


full_spiral_name::full_spiral_name(const string &str) : graph_type(CAGE), search_scheme(SS_UNSPECIFIED),
							construction_scheme(CUBIC), point_group("null"),
							base_face_degree(6), face_degrees({5})
{
  vector<string> segments = find_parenthetical(str,"[]");
  // cerr << "segments = " << segments[0] << "; " << segments[1] << "; " << segments[2] << "\n";
    
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

      // Default search scheme is UNSPECIFIED
      if(s == "GS") search_scheme = CANONICAL_GENERALIZED_SPIRAL;
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
    if(suffix_segments.size()>=2)
      chemical_formula = suffix_segments[suffix_segments.size()-2];
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
    vector<string> fulleroid_face_spec = find_parenthetical(suffix_segments[0],"()");
    base_face_degree = from_string<int>(fulleroid_face_spec[2]);
    if(base_face_degree < 3) base_face_degree = 6;
    
    face_degrees = split<int>(fulleroid_face_spec[1],",");
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


map<string,string> spiral_paper_examples{{
    {"Tetrahedron", "[3,3,3,3]"},
    {"Truncated_tetrahedron","[CS: 6,3,6,3,6,3,6,3]"},
    {"Truncated_tetrahedron GS","[6,1; 3,6,6,6,3,3,6,3]"},
    {"Triakis_tetrahedron","[T,CS: 6,3,6,3,6,3,6,3]"},
    {"Omnitruncated_octahedron_spiral","[4, 6,6,6, 6,4,6,4,6,4,6,4,6,4]"},
    {"Omnitruncated_octahedron_full",  "Oh-[4, 6,6,6, 6,4,6,4,6,4,6,4,6,4]-24-cage"},
    {"Tutte_graph","[CS:11, 1, 17, 1; 5, 10, 5, 5, 5, 9, 5, 4, 5, 4, 4, 5, 4, 10, 5, 5, 5, 5, 5, 10, 4, 5, 5, 4, 5]"},
    {"Tutte_molecule","C3-[CS: 11,1,17,1; 5,10,5,5,5,9,5,4,5,4,4,5,4,10,5,5,5,5,5,10,4,5,5,4,5]-C46-cage"},
      // Different ways of writing the Td-C100
    {"Td-C100_shortest","[2,8,9,23,24,28,29,37,41,45,46,52]-fullerene"},
    {"Td-C100_full_GS","Td-[GS: 43,2; 1,4,5,26,27,31,32,40,43,47,48,52]-C100-fullerene"},
    {"Td-C100_roid_a","[GS: 43,2; 1,4,5,26,27,31,32,40,43,47,48,52]-(5)-C100-fulleroid"},
    {"Td-C100_roid_b","[GS: 43,2; 1,4,5,26,27,31,32,40,43,47,48,52]-(5)6-fulleroid"},
    {"Td-C100_roid_c","[GS: 43,2; 1,4,5,26,27,31,32,40,43,47,48,52]-(5)_6-fulleroid"},
    {"Td-C100_roid_d","[GS: 43,2; 1,4,5,26,27,31,32,40,43,47,48,52]-(4,5)-fulleroid"},
    {"Td-C100_roid_e","[GS: 43,2; 1,4,5,26,27,31,32,40,43,47,48,52]-(4,5,8)-fulleroid"},            
      // Non-spiralable fullerene examples
    {"NS-T-C380","[GS: 162,2,186,3; 1,4,5,126,127,136,137,155,162,171,172,186]-C380-fullerene"},
    {"NS-D3-C384","D3-[GS:171,8,178,9; 1,3,4,5,168,169,170,178,190,191,192,194]-C384-fullerene"},
    {"NS-D3-C440","D3-[CS:62,1;39,40,41,62,170,171,197,198,218,219,220,222]-C440-fullerene"},
    {"NS-D3-C672","D3-[GS:297,8,310,9; 1,10,12,14,260,262,264,310,324,326,328,338]-C672-fullerene"},
      // Non-cubic polyhedra nanoparticle-examples
    {"M12L24","Oh-[LF,GS: 9,1,12,1,15,1,24,3; "
	   "4, 8,8,8,8, "
	   "3, 8, 4, 3, 8, 4, 3, 8, 4, 3, 8, 4,"
 	   "8, 3,8, 3,8, 3,"
  	   "3, 8, 4]-12-cage"},
    {"M24L48","Oh-[LF,GS: 22,1,27,1,42,3; "
	   "3, 8,8,8, "
	   "4,8, 4,8, 4,8, 4,8, 4,8, 4,8, "
	   "3,8,4,8,4, 3,8,4,8,4, 3,8,4,8,4, "
	   "8,4,3, "
	   "8,4, 8,4, 8,4, 8,4, "
	   "8,3"
	"]-24-cage"},
    {"M30L60","O-[LF,GS: 54,1,57,1,59,1; "
	   "3,8,8,8, "
	   "4,8, 4,8, 4,8, 4,8, 4,8, 4,8, 4,8, "
	   "3,8,4,8,4,8, 3,8,4,8,4,8, 3,8,4,8,4,8, "
	   "4,8, "
	   "3,8,4,8,4,8, 3,8,4,8,4,8, "
	   "8,3, "
	   "8,4,4, 8,4,4, 8,4,4, "
	   "8,3"
	"]-30-cage"}
  }};


int main(int ac, char **av)
{
  for(auto example: spiral_paper_examples){
    cout << example.first << "_name = \"" << example.second << "\";\n";
    cout << example.first << " = " << full_spiral_name(example.second) << ";\n\n";
  }
  
  return 0;
}
