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
  PlanarGraph dual(g.dual_graph(6));


  vector<int> pentagons;
  for(int i=0;i<dual.N;i++) if(dual.neighbours[i].size() == 5) pentagons.push_back(i);
  cout << "pentagon indices: (normally 0 through 11)" << endl;
  cout << pentagons << endl<<endl;

  // I know the following is very wasteful but it's good enough
  vector<int> all_distances = dual.all_pairs_shortest_paths(INT_MAX);
  vector<int> pentagon_distances(144, 0);

  cout << "pentagon distances: " << endl;
  for(int i=0; i!=12; ++i){
    for(int j=0; j!=12; ++j){
      pentagon_distances[12*i+j] = all_distances[dual.N * pentagons[i] + pentagons[j]];
    }
  }
  
  cout << "{";
  for(int i=0; i!=12; ++i){
    cout << "{";
    for(int j=0; j!=12; ++j){
      cout << pentagon_distances[12*i + j];
      if(j!=11)cout << ", ";
    }
    cout << "}";
    if(i!=11)cout << ",";
    cout << endl;
  }
  cout << "}" << endl;

  return 0;
}
