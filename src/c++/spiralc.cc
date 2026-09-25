#include <utility>
#include <vector>
#include <string>
#include <stdexcept>
#include <sstream>
#include <charconv>
#include <climits>
#include <set>


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


string spiral_nomenclature::search_scheme_txt[4]       = {"UNSPECIFIED","CANONICAL_GENERALIZED_SPIRAL","COMPATIBILITY_CANONICAL_SPIRAL"};
string spiral_nomenclature::construction_scheme_txt[4] = {"UNSPECIFIED","CUBIC","TRIANGULATION", "LEAPFROG"};
string spiral_nomenclature::naming_scheme_txt[4]          = {"null","FULLERENE", "FULLEROID", "CAGE"};

namespace {

// The ints of a ','-separated list (empty for the empty string).  Each entry
// is a decimal numeral without sign that fits an int.
vector<int> parse_int_list(const string& name, const string& list)
{
  vector<int> xs;
  for(const string& tok: split<string>(list, ",")){
    unsigned long long x = 0;
    const char *first = tok.data(), *last = tok.data() + tok.size();
    const auto [end, ec] = std::from_chars(first, last, x);
    if(tok.empty() || ec != std::errc() || end != last || x > (unsigned long long)INT_MAX)
      throw SpiralNameMalformed(name, SpiralNameMalformed::NUMBER,
                                "\"" + tok + "\" is not a decimal int without sign");
    xs.push_back(int(x));
  }
  return xs;
}

// Jumps from the (position, rotations) pairs of a name: 1-based positions
// >= 3 (the windup applies a jump from its third vertex on) and strictly
// increasing (it consumes the jumps in order), rotations >= 1.
jumplist_t parse_jumps(const string& name, const vector<int>& J)
{
  auto refuse = [&](const string& why) -> void {
    throw SpiralNameMalformed(name, SpiralNameMalformed::JUMPS, why);
  };
  if(J.empty() || J.size() % 2 != 0)
    refuse("the jump list holds " + std::to_string(J.size()) + " ints, not a positive even number");
  jumplist_t jumps;
  for(size_t i = 0; i < J.size(); i += 2){
    const int position = J[i], rotations = J[i+1];
    if(position < 3)  refuse("jump position " + std::to_string(position) + " is below 3");
    if(rotations < 1) refuse("jump at " + std::to_string(position) + " rotates by " + std::to_string(rotations));
    if(!jumps.empty() && position - 1 <= jumps.back().first)
      refuse("jump positions are not strictly increasing at " + std::to_string(position));
    jumps.push_back({position - 1, rotations});   // stored 0-based
  }
  return jumps;
}

// Face positions of one degree: every position >= 1, strictly increasing.
void check_positions(const string& name, const vector<int>& P, int degree)
{
  for(size_t i = 0; i < P.size(); i++){
    if(P[i] < 1)
      throw SpiralNameMalformed(name, SpiralNameMalformed::INDEX_RANGE,
                                "face position " + std::to_string(P[i]) + " is below 1");
    if(i > 0 && P[i] <= P[i-1])
      throw SpiralNameMalformed(name, SpiralNameMalformed::INDEX_ORDER,
                                "the positions of the degree-" + std::to_string(degree)
                                + " faces are not strictly increasing at " + std::to_string(P[i]));
  }
}

// The suffix after "]": the graph type, the chemical formula, and for a
// fulleroid the "(d1,...,dk)" face-degree group.
struct NameSuffix {
  string type = "cage", formula, faces;
};

NameSuffix parse_suffix(const string& name, const string& suffix_string)
{
  auto refuse = [&](const string& why) -> void {
    throw SpiralNameMalformed(name, SpiralNameMalformed::SUFFIX, why);
  };
  NameSuffix sfx;
  vector<string> segments = split<string>(suffix_string, "-\u2013");
  if(segments.empty()) return sfx;                              // "[...]": a cage
  if(segments.size() >= 2 && segments[0].empty())               // the leading "-"
    segments.erase(segments.begin());
  sfx.type = segments.back();
  segments.pop_back();
  if(sfx.type != "cage" && sfx.type != "fullerene" && sfx.type != "fulleroid")
    refuse("graph type \"" + sfx.type + "\" is not one of \"cage\", \"fullerene\", \"fulleroid\"");
  for(const string& seg: segments)
    if(seg.empty()) refuse("an empty \"-\"-separated segment");

  // The face-degree group is the segment opening with "(": exactly one for a
  // fulleroid, none otherwise.  It may stand before or after the formula.
  vector<string> formulas;
  for(const string& seg: segments){
    if(seg[0] != '('){ formulas.push_back(seg); continue; }
    if(sfx.type != "fulleroid") refuse("a face-degree group \"" + seg + "\" on a " + sfx.type);
    if(!sfx.faces.empty())      refuse("a second face-degree group \"" + seg + "\"");
    sfx.faces = seg;
  }
  if(sfx.type == "fulleroid" && sfx.faces.empty())
    refuse("a fulleroid names its face degrees \"(d1,...,dk)\"");
  if(formulas.size() > 1)
    refuse(std::to_string(formulas.size()) + " formula segments before the graph type, at most one");
  if(!formulas.empty()) sfx.formula = formulas[0];
  return sfx;
}

// The fulleroid face-degree group "(d1,...,dk)" with an optional base-degree
// tail "6" or "_6": distinct degrees >= 3, none the base degree 6.
vector<int> parse_face_degrees(const string& name, const string& group)
{
  const size_t close = group.find(')');
  if(close == string::npos || group.find('(', 1) != string::npos || group.find(')', close+1) != string::npos)
    throw SpiralNameMalformed(name, SpiralNameMalformed::SUFFIX,
                              "face-degree group \"" + group + "\" is not one \"(\" ... \")\"");
  const string tail = group.substr(close+1);
  if(!(tail.empty() || tail == "6" || tail == "_6"))
    throw SpiralNameMalformed(name, SpiralNameMalformed::SUFFIX,
                              "base degree \"" + tail + "\": only the hexagonal base 6 is supported");

  const vector<int> degrees = parse_int_list(name, group.substr(1, close-1));
  auto refuse = [&](const string& why) -> void {
    throw SpiralNameMalformed(name, SpiralNameMalformed::DEGREES, why);
  };
  if(degrees.empty()) refuse("the face-degree group is empty");
  for(size_t i = 0; i < degrees.size(); i++){
    if(degrees[i] < 3) refuse("face degree " + std::to_string(degrees[i]) + " is below 3");
    if(degrees[i] == 6) refuse("face degree 6 is the base degree");
    for(size_t j = 0; j < i; j++)
      if(degrees[j] == degrees[i]) refuse("face degree " + std::to_string(degrees[i]) + " is listed twice");
  }
  return degrees;
}

} // namespace

spiral_nomenclature::spiral_nomenclature(const string &str) : naming_scheme(CAGE), search_scheme(CANONICAL_GENERALIZED_SPIRAL),
							construction_scheme(CUBIC),
							base_face_degree(6), face_degrees({5})
{
  using M = SpiralNameMalformed;
  const string ws = " \t\r\n";

  //------------------------------ Brackets ------------------------------
  const size_t open = str.find('['), close = str.find(']');
  if(open == string::npos || close == string::npos || close < open
     || str.find('[', open+1) != string::npos || str.find(']', close+1) != string::npos)
    throw M(str, M::BRACKETS, "a name holds exactly one \"[\" followed by one \"]\"");

  //------------------------------ Prefix: the point group ------------------------------
  point_group = trim(str.substr(0, open), "-\u2013" + ws);

  //------------------------------ Scheme tags ------------------------------
  string numbers_spec = trim(str.substr(open+1, close-open-1), ws);
  if(const size_t colon = numbers_spec.find(':'); colon != string::npos){
    if(numbers_spec.find(':', colon+1) != string::npos)
      throw M(str, M::SCHEME_TAG, "more than one \":\"");
    bool construction_given = false, search_given = false;
    for(const string& tag: split<string>(numbers_spec.substr(0, colon), ",")){
      bool is_construction = true;
      if     (tag == "C")  construction_scheme = CUBIC;
      else if(tag == "T")  construction_scheme = TRIANGULATION;
      else if(tag == "LF") construction_scheme = LEAPFROG;
      else if(tag == "GS"){ search_scheme = CANONICAL_GENERALIZED_SPIRAL;   is_construction = false; }
      else if(tag == "CS"){ search_scheme = COMPATIBILITY_CANONICAL_SPIRAL; is_construction = false; }
      else throw M(str, M::SCHEME_TAG, "unknown scheme tag \"" + tag + "\"");
      bool &given = is_construction ? construction_given : search_given;
      if(given) throw M(str, M::SCHEME_TAG, string("a second ") + (is_construction ? "construction" : "search")
                                            + " scheme tag \"" + tag + "\"");
      given = true;
    }
    if(!construction_given && !search_given) throw M(str, M::SCHEME_TAG, "\":\" without a scheme tag");
    numbers_spec = numbers_spec.substr(colon+1);
  }

  // The ';'-separated int lists: optional jumps, then the spiral numbers.
  vector<vector<int>> spiral_numbers;
  for(const string& segment: split<string>(numbers_spec, ";"))
    spiral_numbers.push_back(parse_int_list(str, segment));

  //------------------------------ Suffix ------------------------------
  const NameSuffix sfx = parse_suffix(str, str.substr(close+1));
  chemical_formula = sfx.formula;
  naming_scheme = sfx.type == "cage" ? CAGE : sfx.type == "fullerene" ? FULLERENE : FULLEROID;
  if(naming_scheme == FULLEROID) face_degrees = parse_face_degrees(str, sfx.faces);

  // A cage lists one degree sequence; a fullerene or fulleroid one position
  // list per special degree.  One more list in front is the jump list.
  const size_t n_lists = naming_scheme == CAGE ? 1 : face_degrees.size();
  if(spiral_numbers.size() != n_lists && spiral_numbers.size() != n_lists + 1)
    throw M(str, M::SEGMENT_COUNT, std::to_string(spiral_numbers.size()) + " \";\"-separated lists, expected "
            + std::to_string(n_lists) + " (or " + std::to_string(n_lists+1) + " with jumps)");
  const bool has_jumps = spiral_numbers.size() == n_lists + 1;
  const jumplist_t jumps = has_jumps ? parse_jumps(str, spiral_numbers[0]) : jumplist_t();

  if(naming_scheme == CAGE){
    // A vertex of a triangulated sphere has degree >= 3, and the degrees d
    // satisfy Euler's relation sum(6 - d) = 12.
    const vector<int> &degrees = spiral_numbers[has_jumps];
    long curvature = 0;
    for(int d: degrees){
      if(d < 3) throw M(str, M::DEGREES, "vertex degree " + std::to_string(d) + " is below 3");
      curvature += 6 - d;
    }
    if(curvature != 12)
      throw M(str, M::DEGREES, "the degrees give sum(6 - d) = " + std::to_string(curvature)
              + ", not 12 (Euler's relation for a triangulated sphere)");
    if(!jumps.empty() && jumps.back().first >= int(degrees.size()))
      throw M(str, M::JUMPS, "jump position " + std::to_string(jumps.back().first + 1)
              + " lies past the " + std::to_string(degrees.size()) + "-vertex spiral");
    cage_constructor(spiral_numbers);
    return;
  }

  // Fullerene or fulleroid: positions >= 1, strictly increasing per degree,
  // no position under two degrees, and the faces around hexagons satisfy
  // Euler's relation sum(6 - d) = 12 (for a fullerene: 12 pentagons).
  long curvature = 0;
  set<int> positions;
  int max_number = 0;
  for(size_t i = 0; i < face_degrees.size(); i++){
    const vector<int> &P = spiral_numbers[i + has_jumps];
    check_positions(str, P, face_degrees[i]);
    for(int p: P)
      if(!positions.insert(p).second)
        throw M(str, M::INDEX_ORDER, "face position " + std::to_string(p) + " is listed under two degrees");
    curvature += long(6 - face_degrees[i]) * long(P.size());
  }
  if(naming_scheme == FULLERENE && positions.size() != 12)
    throw M(str, M::INDEX_COUNT, std::to_string(positions.size()) + " pentagon positions, a fullerene has 12");
  if(curvature != 12)
    throw M(str, M::INDEX_COUNT, "the listed faces give sum(6 - d) = " + std::to_string(curvature)
            + ", not 12 (Euler's relation for a sphere tiled around hexagons)");
  // The spiral is padded to twice the largest number written (fulleroid_constructor).
  for(const auto &list: spiral_numbers) for(int x: list) max_number = max(max_number, x);
  if(max_number > INT_MAX / 2)
    throw M(str, M::INDEX_RANGE, "number " + std::to_string(max_number) + " is too large for a spiral length");

  fulleroid_constructor(spiral_numbers, face_degrees, base_face_degree);
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

// Initialises every member except the spiral and its face degrees, which
// name_triangulation sets.
// TODO: Should it be possible to specify base_face_degree?
spiral_nomenclature::spiral_nomenclature(const naming_scheme_t naming_scheme,
					 const construction_scheme_t construction_scheme,
					 bool rarest_special_start) :
  naming_scheme(naming_scheme),
  search_scheme(rarest_special_start? CANONICAL_GENERALIZED_SPIRAL : COMPATIBILITY_CANONICAL_SPIRAL),
  construction_scheme(construction_scheme),
  base_face_degree(6)
{}

void spiral_nomenclature::name_triangulation(const TriangulationView &T, bool rarest_special_start)
{
  spiral = T.get_general_spiral(rarest_special_start);

  // Which face degrees appear?
  set<int> face_degree_set;
  for(int d: spiral.spiral_code) if(d != base_face_degree) face_degree_set.insert(d);
  face_degrees = vector<int>(face_degree_set.begin(), face_degree_set.end());
}

spiral_nomenclature::spiral_nomenclature(const PlanarGraph &G, const naming_scheme_t naming_scheme,
					 const construction_scheme_t construction_scheme,
					 bool rarest_special_start) :
  spiral_nomenclature(naming_scheme, construction_scheme, rarest_special_start)
{
  // CS_NONE: G's shape decides, and the member records the decision.
  if(construction_scheme == CS_NONE) this->construction_scheme = G.enveloping_scheme();
  const Triangulation T = G.enveloping_triangulation(this->construction_scheme);
  name_triangulation(T, rarest_special_start);
}

spiral_nomenclature spiral_nomenclature::fullerene_from_dual(const TriangulationView &dual,
							     bool rarest_special_start)
{
  spiral_nomenclature sn(FULLERENE, CUBIC, rarest_special_start);
  sn.name_triangulation(dual, rarest_special_start);
  return sn;
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
    // Only the non-default schemes are written (spiral.hh, 1. and 2.): the
    // construction scheme unless CUBIC, the search scheme unless
    // CANONICAL_GENERALIZED_SPIRAL.
    if(construction_scheme == LEAPFROG)      schemes.push_back("LF");
    if(construction_scheme == TRIANGULATION) schemes.push_back("T");
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

