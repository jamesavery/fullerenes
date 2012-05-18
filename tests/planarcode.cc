#include "libgraph/cubicgraph.hh"
using namespace std;

int main(int ac, char **av)
{
  int index = ac>1? strtol(av[1],0,0) : 0;

  CubicGraph g(index);

  cout << "g = " << g << endl;

  return 0;
}
