#include "libgraph/json.hh"

int main()
{
  JSON::Scalar b(bool(true)), i(long(3)), f(3.5), s(string("hello"));
  JSON::Node u(b), v(f);
  JSON::Document d(cin);
  
  cout << b << ", " << i << ", " << f << ", " << s << ", " << u << ", " << v << endl;
  return 0;
}
