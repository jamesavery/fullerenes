#include "fullerenes/fortran_format.hh"

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <stdexcept>

namespace fortran {

using std::string; using std::string_view; using std::to_string;

string_view field(const string& s, int& pos, int len, const char* descriptor)
{
  if(pos < 0 || len < 0 || size_t(pos) + size_t(len) > s.size())
    throw std::invalid_argument(string("Fortran ") + descriptor + " field at column " +
                                to_string(pos+1) + " runs past the end of the record");
  string_view f(s.data() + pos, len);
  pos += len;
  return f;
}

void read_A(char *result, const string& s, int& pos, int len)
{
  string_view f = field(s, pos, len, "A");
  f.copy(result, len);
}

long read_I(const string& s, int& pos, int len)
{
  string_view f = field(s, pos, len, "I");
  long v = 0; int digits = 0; bool negative = false, sign_seen = false;
  for(char c: f){
    if(c == ' ') continue;
    if((c == '-' || c == '+') && !sign_seen && digits == 0){ sign_seen = true; negative = (c=='-'); continue; }
    if(c >= '0' && c <= '9'){
      if(++digits > 18) throw std::invalid_argument("Fortran I" + to_string(len) + " field '" + string(f) + "' overflows");
      v = 10*v + (c-'0'); continue;
    }
    throw std::invalid_argument("Fortran I" + to_string(len) + " field '" + string(f) + "' is not an integer");
  }
  if(sign_seen && digits == 0)
    throw std::invalid_argument("Fortran I" + to_string(len) + " field '" + string(f) + "' is a sign without digits");
  return negative? -v : v;
}

double read_F(const string& s, int& pos, int len)
{
  string_view f = field(s, pos, len, "F");
  string compact;                          // the field with blanks removed
  for(char c: f) if(c != ' ') compact += c;
  if(compact.empty()) return 0;
  size_t points = 0, digits = 0;
  for(size_t i=0;i<compact.size();i++){
    char c = compact[i];
    if(c == '.') points++;
    else if(c >= '0' && c <= '9') digits++;
    else if(!((c == '-' || c == '+') && i == 0))
      throw std::invalid_argument("Fortran F" + to_string(len) + " field '" + string(f) + "' is not a decimal number");
  }
  if(points != 1 || digits == 0)
    throw std::invalid_argument("Fortran F" + to_string(len) + " field '" + string(f) + "' is not a decimal number");
  return strtod(compact.c_str(), nullptr);
}

string write_I(long v, int len)
{
  string digits = to_string(v);
  if((int)digits.size() > len)
    throw std::invalid_argument("value " + digits + " does not fit a Fortran I" + to_string(len) + " field");
  return string(len - digits.size(), ' ') + digits;
}

string write_F(double v, int len, int decimals)
{
  if(!std::isfinite(v))
    throw std::invalid_argument("value " + to_string(v) + " is not a finite number for a Fortran F field");
  char buf[64];
  int n = snprintf(buf, sizeof buf, "%*.*f", len, decimals, v);
  string s(buf, n < 0? 0 : n);
  if((int)s.size() == len + 1 && s.rfind("-0.", 0) == 0) s.erase(1, 1);   // gfortran: "-.50000"
  if(n < 0 || (int)s.size() > len)
    throw std::invalid_argument("value " + to_string(v) + " does not fit a Fortran F" +
                                to_string(len) + "." + to_string(decimals) + " field");
  return s;
}

}  // namespace fortran
