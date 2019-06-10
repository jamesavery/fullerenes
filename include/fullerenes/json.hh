#include <map>
#include <vector>
#include <string>
#include <iostream>

using namespace std;

// TODO: Pretty-printing.
// TODO?: Nodes should know their indentation level - distance both from top and bottom.
namespace JSON {


  
  class Scalar {
  public:
    typedef enum { NullType, BoolType, IntegerType, FloatType, StringType } value_t;
    value_t scalar_type;

    bool bool_value;
    long integer_value;
    double float_value;
    string string_value;

    Scalar() : scalar_type(NullType) {}
    Scalar(const bool& x) : scalar_type(BoolType), bool_value(x) {}
    Scalar(const long& x) : scalar_type(IntegerType), integer_value(x) {}
    Scalar(const double& x) : scalar_type(FloatType), float_value(x) {}
    Scalar(const string& x) : scalar_type(StringType), string_value(x) {}
    
    friend ostream& operator << (ostream& out, const Scalar& x){
      switch(x.scalar_type){
      case NullType:
	out << "~"; break;
      case BoolType: 
	out << x.bool_value; break;
      case IntegerType:
	out << x.integer_value; break;
      case FloatType:
	out << x.float_value; break;
      case StringType:
	out << '"' << x.string_value << '"'; break; // TODO: Escape string literal
      }
      return out;
    }
  };


  class Node {
  public:
    typedef enum {  ScalarType, SequenceType, MapType, ErrorType } value_t;
    value_t node_type;

    Scalar           scalar_value;
    vector<Node>     sequence_value;
    map<string,Node> map_value;

    // Add where and what to error type node
    Node() : node_type(ErrorType) {}
    Node(const Scalar& x) : node_type(ScalarType), scalar_value(x) {}
    Node(const vector<Node>& seq) : node_type(SequenceType), sequence_value(seq.begin(),seq.end()) {}
    Node(const map<string,Node>& m) : node_type(MapType), map_value(m.begin(),m.end()) {}

    friend ostream& operator << (ostream& out, const Node& node){
      switch(node.node_type){
      case ScalarType:
	out << node.scalar_value; break;
      case SequenceType: 
	{
	  const vector<Node>& seq(node.sequence_value);
	  out << "[";
	  for(size_t i=0;i<seq.size();i++) out << seq[i] << (i+1<seq.size()?", ":"");
	  out << "]";
	}
	break;
      case MapType: 
	{
	  const map<string,Node>& m(node.map_value);
	  out << "{";
	  for(map<string,Node>::const_iterator e(m.begin()); e!=m.end(); ){
	    out << e->first << ":" << e->second;
	    if(++e != m.end()) out << ", ";
	  }
	  out << "}";
	}
      }
      return out;
	
    }
  };

  class Document : public Node {
  public:

    Node parse_number(const string& input, string::const_iterator& next_char)
    {
      size_t end;
      for(end=0;next_char + end != input.end();end++){
	char c = *(next_char+end);
	if(c < '0' && c > '9' && c != '-' && c != '.' && c != 'e') break;
      }
      string token(next_char, next_char + end);
      const char *start(token.c_str()), *endptr = 0;
      double d = strtod(start,&endptr);
      if(endptr != 0) { next_char += end; return Scalar(double(d)); }
      long i = strtol(start,&endptr);
      if(endptr != 0) { next_char += end; return Scalar(long(i)); }
      
      if(string(next_char,next_char+4) == "null")  { next_char += 4;  return Node(); }
      if(string(next_char,next_char+4) == "true")  { next_char += 4;  return Scalar(bool(true)); }
      if(string(next_char,next_char+5) == "false") { next_char += 5; return Scalar(bool(false)); }

      // REPORT ERROR
      return Node();
    }
    void skip_whitespace(const string& input, string::const_iterator& c)
    {
      while(c!=input.end() && (*c == ' ' || *c == '\t' || *c == '\n' || *c == '\r')) c++;
    }
    bool next_token_is(const string& input, string::const_iterator& next_char, char expected)
    {
	skip_whitespace(input,next_char);
	if(*next_char != ':') return false;
	skip_whitespace(input,++next_char);
	return true;
    }
    Node parse_object(const string& input, string::const_iterator& next_char)
    {
      map<string,Node> elements;
      while(next_char != input.end() && *next_char != ']'){
	skip_whitespace(input,next_char);
	string k(parse_key(input,next_char));
	// Check for error, and report position in array?
	if(!next_token_is(input,next_char,':')) // ERROR
	  ;
	Node v((parse_value(input,next_char)));
	// Check for error, and report position in array?
	skip_whitespace(input,next_char);
	// Check that *(next_char+1) == ']' or *next_char == ','
      }
      return Node(elements);
    }

    Node parse_array(const string& input, string::const_iterator& next_char)
    {
      vector<Node> elements;
      while(next_char != input.end() && *next_char != ']'){
	elements.push_back((parse_value(input,next_char)));
	// Check for error, and report position in array?
	// Check that *(next_char+1) == ']' or *next_char == ','
      }
      return Node(elements);
    }

    Node parse_value(const string& input, string::const_iterator& next_char)
    {
      while(next_char != input.end()){
	switch(*next_char){
	case '{': 
	  return parse_object(input,++next_char);
	case '[':
	  return parse_array(input,++next_char);
	case '"':
	case "'":
	  return parse_string(input,++next_char);

	case ' ': // Skip whitespace
	case '\t':
	case '\n':
	  next_char++;
	  break;
	default:
	  return parse_number(input,next_char);
	}
	next_char++;
      }
    }

    Document(const string& input) : Node(vector<Node>()) {
      list<Node> stack;
      Node current_object;

      string::const_iterator next_char(input.begin());
    }
  };
}
