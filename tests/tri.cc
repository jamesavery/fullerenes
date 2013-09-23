#include "libgraph/triangulation.hh"

int main(int ac, char **av)
{
  int default_spiral[12] = {2, 8, 9, 23, 24, 28, 29, 37, 41, 45, 46, 52};

  int N = ac>=13? strtol(av[1],0,0) : 100, M = 2+N/2;
  vector<int> spiral(12), full_spiral(M,6);
  for(int i=0;i<12;i++){
    spiral[i] = (ac>=13? strtol(av[i+2],0,0) : default_spiral[i])-1;
    full_spiral[spiral[i]] = 5;
  }

  Triangulation dual(full_spiral);
  dual.orient_neighbours();
  PlanarGraph      g(dual.dual_graph());

  cout << "dual = " << dual << ";\n";
  cout << "g    = " << g << ";\n";

  printf("dual is %soriented, g is %soriented\n",dual.is_consistently_oriented()? "":"not ",g.is_consistently_oriented()? "":"not ");

  return 0;
}
