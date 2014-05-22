%module fullwrap
%include <std_vector.i>
%include <std_string.i>
%include <stdint.i>

%{
#include <vector>
#include <inttypes.h>
#include "FullereneSelect.hh"
%}

%include "FullereneSelect.hh"

class IsomerDB {
 public:
  int N, Nisomers;

  static size_t number_isomers(int N, const string& sym = "Any", bool IPR=true);
  static vector<string> symmetries(int N,bool IPR);
};

/* class Unfolding { */
/*  public: */
/*   Unfolding(const PlanarGraph& dual, bool planar_layout = false); */
/*   Unfolding operator *(const Eisenstein& y) const; */
/*   Unfolding operator +(const Eisenstein& y) const; */
/* } */

%template(entryvector)   std::vector<Entry>;
%template(stringvector)  std::vector<std::string>;
%template(doublevector)  std::vector<double>;
%template(intvector)     std::vector<int>;

// %include polynomial.i
// %include goscinskian.i
