#include <libgraph/polyhedron.hh>
#include <libgraph/triangulation.hh>

int main(int ac, char **av)
{
  assert(ac>=2);

  Polyhedron P(av[1]);
  
  P.move_to_origin();
  matrix3d If(P.inertial_frame());
  P.points = If*P.points;

  cout << P.to_povray() << endl;
  

  return 0;
}
