#include <utility>
#include <vector>
#include <string>


#include "fullerenes/spiral.hh"
#include "fullerenes/triangulation.hh"

using namespace std;

// TODO: auxiliary.cc
string trim(const string& str, const string& wschars)
{
    size_t first = str.find_first_not_of(wschars);
    if(first == string::npos)
      return "";
    size_t last = str.find_last_not_of(wschars);

    return str.substr(first, (last - first + 1));
}

// Version of split that handles empty strings ("a;;b;c;" -> {"a","","b","c",""} in stead of {"a","b","c"}.
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




string spiral_nomenclature::search_scheme_txt[4]       = {"UNSPECIFIED","CANONICAL_GENERALIZED_SPIRAL","COMPATIBILITY_CANONICAL_SPIRAL"};
string spiral_nomenclature::construction_scheme_txt[4] = {"UNSPECIFIED","CUBIC","TRIANGULATION", "LEAPFROG"};
string spiral_nomenclature::naming_scheme_txt[4]          = {"null","FULLERENE", "FULLEROID", "CAGE"};

spiral_nomenclature::spiral_nomenclature(const string &str) : naming_scheme(CAGE), search_scheme(SS_UNSPECIFIED),
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
    naming_scheme = CAGE;
    cage_constructor(spiral_numbers);
    return;
  }

  if(!(suffix == "fullerene" || suffix == "fulleroid")){
    cerr << "Graph type is \""<<suffix<<"\", must be one of \"cage\", \"fullerene\", or \"fulleroid\".\n";
    abort();
  }

  if (suffix == "fullerene"){
    naming_scheme = FULLERENE;
    base_face_degree = 6;
    face_degrees     = vector<int>(1,5);
  }

  if(suffix == "fulleroid"){
    naming_scheme = FULLEROID;
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

void spiral_nomenclature::cage_constructor(const vector<vector<int>> &spiral_numbers)
{
  assert(spiral_numbers.size() == 1 || spiral_numbers.size() == 2);
  bool has_jumps = (spiral_numbers.size() == 2);

  if(has_jumps){
    for(int i=0; i<spiral_numbers[0].size()/2; i++){
      spiral.jumps.push_back(make_pair(spiral_numbers[0][2*i]-1,spiral_numbers[0][2*i+1])); // -1 because indices start counting at 0
    }
  }
  spiral.spiral_code = spiral_numbers[has_jumps];

  // Set face_degrees to something sensible
  set<int> face_degree_set;
  for(int f: spiral.spiral_code)
    if(f!=base_face_degree) face_degree_set.insert(f);

  face_degrees = vector<int>(face_degree_set.begin(), face_degree_set.end());
}

void spiral_nomenclature::fulleroid_constructor(const vector<vector<int>> &spiral_numbers, vector<int> face_degrees, int base_face_degree)
{
  assert(spiral_numbers.size() ==  face_degrees.size() || spiral_numbers.size() == face_degrees.size()+1);

  int max_index = 0;
  for(const auto &v: spiral_numbers)
    for(const auto &ix: v)
      max_index = max(ix,max_index);

  spiral.spiral_code = vector<int>(2*max_index,6);

  bool has_jumps = (spiral_numbers.size() == face_degrees.size()+1);

  if(has_jumps){
    for(int i=0; i<spiral_numbers[0].size()/2; i++){
      spiral.jumps.push_back(make_pair(spiral_numbers[0][2*i]-1,spiral_numbers[0][2*i+1])); // -1 because indices start counting at 0
    }
  }

  for(int i=0;i<face_degrees.size();i++){
    for(auto ix: spiral_numbers[i+has_jumps]){
      spiral.spiral_code[ix-1] = face_degrees[i];
    }
  }
}

// TODO: Should it be possible to specify base_face_degree?
// TODO: Should it be possible to specify base_face_degree?
spiral_nomenclature::spiral_nomenclature(const PlanarGraph &G, const naming_scheme_t naming_scheme
,					 const construction_scheme_t construction_scheme,
					 bool rarest_special_start) :
  naming_scheme(naming_scheme), search_scheme(rarest_special_start? CANONICAL_GENERALIZED_SPIRAL : COMPATIBILITY_CANONICAL_SPIRAL),
  base_face_degree(6)
{
  Triangulation T;
  if(construction_scheme==CS_NONE)
    T = G.enveloping_triangulation(construction_scheme); // This *writes* to construction_scheme
  else
    T = G.enveloping_triangulation(static_cast<const construction_scheme_t>(construction_scheme)); // This *reads* from construction_scheme
  
  spiral = T.get_general_spiral(rarest_special_start);

  // Which face degrees appear?
  set<int> face_degree_set;
  for(int d: spiral.spiral_code) if(d != base_face_degree) face_degree_set.insert(d);
  face_degrees = vector<int>(face_degree_set.begin(), face_degree_set.end());
}

template <typename T> string riffle(const vector<T>& xs, string delim, string end_if_nonempty="")
{
  string s;
  for(int i=0;i<xs.size();i++) s += to_string(xs[i]) + (i+1<xs.size()? delim : end_if_nonempty);
  return s;
}

string spiral_nomenclature::to_string(bool unpacked) const
{
  if(unpacked){
    ostringstream s;
    s << "<|\n\t"
      << "\"naming_scheme\" -> \""<<spiral_nomenclature::naming_scheme_txt[naming_scheme]<<"\",\n\t"
      << "\"search_scheme\" -> \""<<spiral_nomenclature::search_scheme_txt[search_scheme]<<"\",\n\t"
      << "\"construction_scheme\" -> \""<<spiral_nomenclature::construction_scheme_txt[construction_scheme]<<"\",\n\t"
      << "\"point_group\" -> \""<<(point_group.empty()? "UNSPECIFIED" : point_group) <<"\",\n\t"
      << "\"chemical_formula\" -> \""<<chemical_formula<<"\",\n\t"
      << "\"base_face_degree\" -> "<<base_face_degree<<",\n\t"
      << "\"face_degrees\" -> " << face_degrees << ",\n\t"
      << "\"jumps\" -> " << spiral.jumps << ",\n\t" // indices start counting at 0
      << "\"spiral_code\" -> " << spiral.spiral_code //<< ", (length: " << spiral_code.size() << ") \n\t"
      << "|>";
    return s.str();
  } else {
    // Add point group prefix if present
    string prefix = point_group.empty()? "" : (point_group+"-");

    // Encode construction and search schemes
    string scheme_string;
    vector<string> schemes;
    // Default construction_scheme is CUBIC, default search_scheme is UNSPECIFIED
    if(construction_scheme == LEAPFROG)      schemes.push_back("LF");
    if(construction_scheme == TRIANGULATION) schemes.push_back("T");
    if(search_scheme == CANONICAL_GENERALIZED_SPIRAL)   schemes.push_back("GS");
    if(search_scheme == COMPATIBILITY_CANONICAL_SPIRAL) schemes.push_back("CS");

    scheme_string = riffle(schemes,",",":");
      
    // Encode jumps
    vector<int> jumps_plus_one(spiral.jumps.size()*2);
    for(int i=0;i<spiral.jumps.size();i++){
      jumps_plus_one[2*i]   = spiral.jumps[i].first+1;
      jumps_plus_one[2*i+1] = spiral.jumps[i].second;
    }
    string jump_string = riffle(jumps_plus_one,",","; ");
    
    // Encode spiral and determine suffix
    string spiral_string, suffix;
    switch(naming_scheme){
    case CAGE:
      spiral_string = riffle(spiral.spiral_code,",");
      suffix        = "cage";
      break;
    case FULLERENE:
      suffix = "fullerene";
    case FULLEROID:
      for(int i=0;i<face_degrees.size();i++){
	vector<int> indices;	
	int d = face_degrees[i];

	for(int j=0;j<spiral.spiral_code.size();j++) if(spiral.spiral_code[j] == d) indices.push_back(j+1);
	spiral_string += riffle(indices,",") + (i+1<face_degrees.size()? ";":"");
      }
      if(suffix.empty()) suffix = "("+riffle(face_degrees,",") + ")-fulleroid"; // TODO: non-6 base face
      break;
    default:
      break;			// TODO: Error
    }

    if(!chemical_formula.empty()) suffix = chemical_formula + "-"+suffix;
    
    return prefix + "[" + scheme_string + jump_string + spiral_string + "]-" + suffix;    
  }
  
}

