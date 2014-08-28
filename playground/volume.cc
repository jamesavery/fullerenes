#include "libgraph/fullerenegraph.hh"
#include "libgraph/polyhedron.hh"

int main(int ac, char **av)
{
  string name(ac<=1? "" : av[1]);
  Polyhedron P;
  cin >> P;
  
  cout << name << "P = " << P << endl;
  cout << name <<"nodes = "<<name<<"P[[1]];\n";
  cout << name <<"xs    = "<<name<<"P[[2]];\n";
  cout << name <<"graph = "<<name<<"P[[3]];\n";
  cout << name << "volume = " << P.volume() << endl;
  cout << name << "volume2 = " << P.volume2() << endl;
  cout << name << "area = " << P.surface_area() << endl;

  return 0;
}
