#include "libgraph/planargraph.hh"
#include "libgraph/cubicgraph.hh"
#include "libgraph/fullerenegraph.hh"

using namespace std;

int main(int ac, char **av)
{

  std::cout << "Attempting to create graph from spiral indices ..." << std::endl;
 
// spiral indices are input starting at 0
//  int n = 20;
//  int pentagon_indices_array[] = {0,1,2,3,4,5,6,7,8,9,10,11};
//  int n = 32;
//  int pentagon_indices_array[] = {0,1,2,3,6,9,10,12,13,15,16,17 };
  int n = 60;
  int pentagon_indices_array[] = {2, 4, 6, 8, 11, 14, 16, 19, 22, 26, 28, 30};
  int jumps_array[][2] = {{171,7}};
  int n_jumps = 0;

  vector<int> pentagon_indices_input(pentagon_indices_array,pentagon_indices_array+12);
  vector<int> pentagon_indices_output;
  vector<int> spiral;
  list<pair<int,int> > jumps;
  
  cout << "N: " << n << endl
       << "Input spiral indices: " << pentagon_indices_input << endl;

  cout << "Input jumps: ";
  for (int i=0; i<n_jumps; ++i){
    jumps.push_back(make_pair(jumps_array[i][0], jumps_array[i][1]));
    cout <<  jumps.rbegin()->first << ", " << jumps.rbegin()->second << "; " ;
  }
  cout << endl;
  cout << endl;

  FullereneGraph fg(n, pentagon_indices_input, jumps);
  PlanarGraph dual(fg.dual_graph(6));

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
