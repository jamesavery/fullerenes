#include "libgraph/cubicgraph.hh"
#include "libgraph/fullerenegraph.hh"

using namespace std;

int main(int ac, char **av)
{

  int n = 192;
//no jump
//  int pentagon_indices[] = {1,2,3,4,5,7,10,12,13,14,15,16};
//jump
//  int pentagon_indices[] = {1,2,3,5,8,9,11,12,13,14,15,16};
//jump, 380
  std::vector<int> pentagon_indices;
  pentagon_indices.push_back(1);
  pentagon_indices.push_back(2);
  pentagon_indices.push_back(3);
  pentagon_indices.push_back(137);
  pentagon_indices.push_back(138);
  pentagon_indices.push_back(147);
  pentagon_indices.push_back(148);
  pentagon_indices.push_back(157);
  pentagon_indices.push_back(158);
  pentagon_indices.push_back(166);
  pentagon_indices.push_back(173);
  pentagon_indices.push_back(180);
  bool ipr=false;
  bool general=true;

  FullereneGraph fg(n, pentagon_indices, ipr, general);

  cout << "fg = " << fg << endl;

  return 0;
}
