#include "libgraph/polyhedron.hh"

int main(int ac, char **av)
{
  assert(ac >= 2);

  string filename = av[1];

  Polyhedron P(Polyhedron::from_mol2(filename));

  cout << "P = " << P << endl;

  return 0;
}
