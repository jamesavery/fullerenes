#include "libgraph/cubicgraph.hh"
#include "libgraph/fullerenegraph.hh"

using namespace std;

int main(int ac, char **av)
{

  int n = 192;
//no jump
//  int pentagon_indices[] = {1, 2, 3, 4, 5, 7, 10, 12, 13, 14, 15, 16};
//jump
//  int pentagon_indices[] = {1, 2, 3, 5, 8, 9, 11, 12, 13, 14, 15, 16};
//jump, 380
  int pentagon_indices[] = {1, 2, 3, 137, 138, 147, 148, 157, 158, 166, 173, 180};
  bool ipr=false;
  bool general=true;

  FullereneGraph fg(n, pentagon_indices, ipr, general);

  cout << "fg = " << fg << endl;

  return 0;
}
