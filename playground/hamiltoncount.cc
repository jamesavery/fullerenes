#include "libgraph/cubicgraph.hh"

int main()
{
  CubicGraph g(stdin);
  unsigned int hamilton_count = g.hamilton_count();
  printf("%d Hamiltonian cycles out of %ld length-N paths (%ld steps)\n",hamilton_count,g.path_count,g.step_count);
  return 0;
}
