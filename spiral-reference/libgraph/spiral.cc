
#include <utility>
#include <vector>
#include <string>


#include "spiral.hh"
#include "triangulation.hh"

using namespace std;

string trim(const string& str, const string& wschars)
{
    size_t first = str.find_first_not_of(wschars);
    if(first == string::npos)
      return "";
    size_t last = str.find_last_not_of(wschars);

    return str.substr(first, (last - first + 1));
}

vector<string> find_parenthetical(const string& parse_str, const char parentheses[2])
{
  size_t start = parse_str.find(parentheses[0], 0);
  if (start == string::npos) return vector<string>();
  size_t end   = parse_str.find(parentheses[1], start);
  if (end == string::npos) return vector<string>();

  return vector<string>{{trim(parse_str.substr(0, start)," \t\r\n"),
	                 trim(parse_str.substr(start+1, end-start-1)," \t\r\n"),
	                 trim(parse_str.substr(end+1,parse_str.size()-end-1)," \t\r\n")}};
}

template <typename T> vector<T> split(const string& parse_str, const string& delimiters, const string wschars=" \t\r\n")
{
  vector<string> string_result = split<string>(parse_str,delimiters,wschars);
  vector<T> result;

  for(string s: string_result) result.push_back(from_string<T>(s));

  return result;
}

// Version of split that handles empty strings ("a;;b;c;" -> {"a","","b","c",""} in stead of {"a","b"}.
template <> vector<string> split(const string& parse_str, const string& delimiters, const string wschars)
{
  // Unlike Python's split(), "" splits to [] instead of [""].
  if(parse_str.empty()) return vector<string>();
  
  vector<string> result;
  size_t
    start = 0,
    end   = parse_str.find_first_of(delimiters);
  
  while(end != string::npos){
    result.push_back(trim(parse_str.substr(start,end-start),wschars));
    start = end+1;
    end   = parse_str.find_first_of(delimiters,start);
  }
  result.push_back(trim(parse_str.substr(start,string::npos),wschars));  

  return result;
}


string full_spiral_name::search_scheme_txt[4]       = {"UNSPECIFIED","CANONICAL_GENERALIZED_SPIRAL","COMPATIBILITY_CANONICAL_SPIRAL"};
string full_spiral_name::construction_scheme_txt[4] = {"UNSPECIFIED","CUBIC","TRIANGULATION", "LEAPFROG"};
string full_spiral_name::graph_type_txt[4]          = {"null","FULLERENE", "FULLEROID", "CAGE"};

full_spiral_name::full_spiral_name(const string &str) : graph_type(CAGE), search_scheme(SS_UNSPECIFIED),
							construction_scheme(CUBIC), 
							base_face_degree(6), face_degrees({5})
{
  vector<string> segments = find_parenthetical(str,"[]");
  //  cerr << "segments = " << segments[0] << "; " << segments[1] << "; " << segments[2] << "\n";

  if(segments.size() != 3){
    cerr << "Error in spiral string \"" << str << "\": Number of \"[]\"-delimited segments is " << segments.size()
	 << ", not 3.\n";
    abort();
  }

  string &prefix_string = segments[0], &spiral_spec = segments[1], &suffix_string = segments[2];

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
  //  cerr << "suffix_segments = " << suffix_segments << ";\n";

  string suffix = suffix_segments.empty()? "cage" : suffix_segments.back();

  int suffix_start = !suffix_segments.empty() && suffix_segments[0].empty()? 1 : 0; // first '-' doesn't separate, it's just there to look nice
  
  // I think this should hold in general, whether it's 'cage', 'fullerene', or 'fulleroid',
  // except in one particlar case: [...]-(f1,..,fp)-fulleroid, handled below.
  if(suffix_segments.size()>=2)
    chemical_formula = suffix_segments[suffix_segments.size()-2];
  
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
    vector<string> fulleroid_face_spec = find_parenthetical(suffix_segments[suffix_start],"()");
    if(fulleroid_face_spec[2].size()==0)
      base_face_degree = 6;
    else{
      base_face_degree = from_string<int>(trim(fulleroid_face_spec[2],"_"));
    }
    if(base_face_degree < 3) return;

    face_degrees = split<int>(fulleroid_face_spec[1],",");
    if(suffix_segments.size() == 3) chemical_formula = "";
  }

  // cerr << "spiral_numbers   = " << spiral_numbers << ";\n"
  //      << "face_degrees     = " << face_degrees << ";\n"
  //      << "base_face_degree = " << base_face_degree << ";\n\n";
  
  fulleroid_constructor(spiral_numbers,face_degrees,base_face_degree);
}

void full_spiral_name::cage_constructor(const vector<vector<int>> &spiral_numbers)
{
  assert(spiral_numbers.size() == 1 || spiral_numbers.size() == 2);
  bool has_jumps = (spiral_numbers.size() == 2);

  if(has_jumps){
    for(int i=0; i<spiral_numbers[0].size()/2; i++){
      jumps.push_back(make_pair(spiral_numbers[0][2*i]-1,spiral_numbers[0][2*i+1])); // -1 because indices start counting at 0
    }
  }
  spiral_code = spiral_numbers[has_jumps];

  // Set face_degrees to something sensible
  set<int> face_degree_set;
  for(int f: spiral_code)
    if(f!=base_face_degree) face_degree_set.insert(f);

  face_degrees = vector<int>(face_degree_set.begin(), face_degree_set.end());
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

  if(has_jumps){
    for(int i=0; i<spiral_numbers[0].size()/2; i++){
      jumps.push_back(make_pair(spiral_numbers[0][2*i]-1,spiral_numbers[0][2*i+1])); // -1 because indices start counting at 0
    }
  }

  for(int i=0;i<face_degrees.size();i++){
    for(auto ix: spiral_numbers[i+has_jumps]){
      spiral_code[ix-1] = face_degrees[i];
    }
  }
}

// TODO: Should it be possible to specify base_face_degree?
full_spiral_name::full_spiral_name(const PlanarGraph &G, const graph_type_t graph_type, bool rarest_special_start) : graph_type(graph_type), search_scheme(rarest_special_start? CANONICAL_GENERALIZED_SPIRAL : COMPATIBILITY_CANONICAL_SPIRAL), base_face_degree(6)
{
  Triangulation T(G.enveloping_triangulation(construction_scheme));
  general_spiral spiral = T.get_general_spiral(rarest_special_start);

  // Which face degrees appear?
  set<int> face_degree_set;
  for(int d: spiral.spiral) if(d != base_face_degree) face_degree_set.insert(d);
  face_degrees = vector<int>(face_degree_set.begin(), face_degree_set.end());

  spiral_code = spiral.spiral;
  jumps       = spiral.jumps;
}


#if 0

bool operator<(const general_spiral &s) const
{
  return jumps.size() < s.jumps.size() ||
	(jumps.size() == s.jumps.size() && jumps < s.jumps) ||
        (jumps == s.jumps && spiral < s.spiral);
}


ostream &operator<<(ostream &s, const general_spiral &GS)
{
  return s << make_pair(GS.jumps,GS.spiral); 
}


spiral_nomenclature::spiral_nomenclature(const PlanarGraph &g, const string& atom="",
					 bool compatibility=false) : atom(atom), graph(g),cs(compatibility) {
  if(g.is_triangulation()){
    graph_type     = TRIANGULATION;
    is_a_fullerene = g.dual_graph().is_a_fullerene();
    triangulation  = g;
  } else if(g.is_cubic()){
    graph_type     = CUBIC;
    triangulation  = g.dual_graph();
    is_a_fullerene = g.is_a_fullerene();
  } else {
    graph_type     = GENERAL;
    triangulation  = g.leapfrog_dual();
    is_a_fullerene = false;
  }

  permutation.resize(triangulation.N);
  bool spiral_success = triangulation.get_spiral(GS.spiral,GS.jumps,!cs);
  assert(spiral_success);

  // Find spiral order with respect to triangulation vertices
  permutation.resize(triangulation.N);
  spiral_success = triangulation.get_spiral(GS.spiral,GS.jumps,permutation,!cs);
}

template <typename T> vector<T> read_list
spiral_nomenclature::spiral_nomenclature(const string& name)
{
  // Grammar for nomenclature as described in spiral paper:
  // 
  // MOLECULE_NAME   ::= (POINTGROUP '-')? POLYHEDRON_NAME '-' FORMULA '-' SUFFIX
  // POLYHEDRON_NAME ::= '[' METHOD_SPEC? (JUMP+';')? SPIRAL ']'
  //                 |   '[' METHOD_SPEC? (JUMP+';')? PI     ']'
  //                 |   '[' METHOD_SPEC? (JUMP+';')? FI     "]-(" FACE_SIZES ")" BASE_FACE
  //
  // METHOD_SPEC     ::= ("D"|"LF"|"CS"|"GS")+ ':' 
  // SUFFIX          ::= 'cage' | 'fullerene' | 'fulleroid'
  // SPIRAL          ::= NUMBER_LIST
  // PI              ::= NUMBER_LIST
  // FI              ::= NUMBER_LISTS
  // NUMBER_LIST     ::= number (',' number)*
  // NUMBER_LISTS    ::= NUMBER_LIST (';' NUMBER_LIST)*
  // FORMULA         ::= string
  // POINTGROUP      ::= string

  // Simplified grammar which accepts the above + some more:
  //
  // SEPARATOR  ::= " []-"
  // POLYHEDRON_NAME ::= '[' METHOD_SPEC? NUMBER_LISTS ']'
  //                 |   '[' METHOD_SPEC? NUMBER_LISTS "]-(" NUMBER_LIST ")" number
  //
  // FI              ::= NUMBER_LIST (';' NUMBER_LIST)*
  // FORMULA         ::= string
  // POINTGROUP      ::= string
  
  // The grammar is so simple that we'll just brutally parse by hand.

  for(int i=0;i<name.size();i++){
  }

}

friend ostream& operator<<(ostream &s, const spiral_nomenclature &n){
  string graph_type_string[3] = {"","D,","LF,"};
  s << "["<<graph_type_string[n.graph_type] <<(n.compatibility?"CS":"GS") << ": "
    << (n.GS.jumps.empty()? "": (jumps_to_string(n.GS.jumps)+"; "))
    << (n.is_a_fullerene? spiral_to_rspi_string(n.GS.spiral) : spiral_to_string(n.GS.spiral))
    << "]-" <<n.atom<< n.graph.N <<"-" << (n.is_a_fullerene? "fullerene" : "cage");
  return s;
}

#endif
