#include <libgraph/polyhedron.hh>
#include <libgraph/triangulation.hh>

int main(int ac, char **av)
{
  assert(ac>=2);

  Polyhedron P(av[1]);
  
  P.move_to_origin();
  P.align_with_axes();

  cout << P.to_povray() << endl;
  

  return 0;
}
