#include <fullerenes/polyhedron.hh>
#include <fullerenes/triangulation.hh>

int main(int ac, char **av)
{
  assert(ac>=2);

  Polyhedron P = Polyhedron::from_file(av[1]);
  
  cout << P.to_povray() << endl;
  

  return 0;
}
