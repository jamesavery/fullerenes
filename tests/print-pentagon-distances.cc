#include "libgraph/planargraph.hh"
#include "libgraph/cubicgraph.hh"
#include "libgraph/fullerenegraph.hh"

using namespace std;

int main(int ac, char **av)
{
  vector<int> rspi(12);
  FullereneGraph::jumplist_t jumps;

  if(ac != 14) return -1;

  int N = atol(av[1]);
  cout << N << endl;
  for(int i=0;i<12;i++) rspi[i] = strtol(av[2+i],0,0)-1;

  cout << "Attempting to create graph from spiral indices " << rspi << endl;
 
  FullereneGraph g(N, rspi, jumps);
  cout << "done" << endl;

  g.pentagon_distance_mtx();


  return 0;
}
