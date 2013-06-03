#include "libgraph/fullerenegraph.hh"

#include <vector>

using namespace std;

int main(int ac, char **av)
{

  cout << "Attempting to create graph from spiral indices ..." << endl;
 
//no jump
//  int n = 28;
//  int pentagon_indices_array[] = {0, 1, 2, 3, 4, 6, 9, 11, 12, 13, 14, 15};
//jump
//  int n = 28;
//  int pentagon_indices_array[] = {0, 1, 2, 4, 7, 8, 10, 11, 12, 13, 14, 15};
//jump, 380
//  int n = 380;
//  int pentagon_indices_array[] = {0, 1, 2, 136, 137, 146, 147, 156, 157, 165, 172, 179};
//  int jumps_array[][2] = {{165,2}, {172,3}, {179,3}};
//  int n_jumps = 3;
//jump, 440
//  int n = 440;
//  int pentagon_indices_array[] = {159, 160, 161, 174, 175, 176, 189, 190, 191, 195, 204, 212};
//  int jumps_array[][2] = {{192,6}, {195,6}, {204,6}, {212,7}};
//  int n_jumps = 4;
//  //0, 1, 9
//jump, 672
  int n = 672;
  int pentagon_indices_array[] = {140, 142, 144, 161, 163, 165, 225, 242, 328, 330, 333, 335};
  int jumps_array[][2] = {{225,5}, {242,6}};
  int n_jumps = 2;
  //0, 115, 151

  vector<int> pentagon_indices_input(12);
  vector<int> pentagon_indices_output;
  list<pair<int,int> > jumps;

  cout << "Input spiral indices: ";
  for (int i=0; i<12; ++i){
    pentagon_indices_input[i] = pentagon_indices_array[i];
    cout <<  pentagon_indices_input[i] << ", ";
  }
  cout << endl;

  cout << "Input jumps: ";
  for (int i=0; i<n_jumps; ++i){
    jumps.push_back(make_pair(jumps_array[i][0], jumps_array[i][1]));
    cout <<  jumps.rbegin()->first << ", " << jumps.rbegin()->second << "; " ;
  }
  cout << endl;

  FullereneGraph fg(n, pentagon_indices_input, jumps);
//  FullereneGraph fg(n, pentagon_indices_input);


//  cout << "fg = " << fg << endl;
//  cout << "number of faces: " << fg.N << endl;
//  cout << "number of edges: " << fg.edge_set.size() << endl;
//  int j=0;
//  for (vector<vector<node_t> >::iterator it(fg.neighbours.begin()); it!=fg.neighbours.end(); ++it, ++j){
//    cout << j+1 << " ";
//    for (vector<node_t>::iterator jt(fg.neighbours[j].begin()); jt!=fg.neighbours[j].end(); ++jt){
//      cout << *jt+1 << " "; 
//    }
//    cout << endl;
//  }

  cout << "Attempting to create spiral from (previously created) graph ..." << endl;

  jumps.clear();
  pentagon_indices_output.clear();
  fg.get_general_spiral_from_fg(0,115,151, pentagon_indices_output, jumps);  

  cout << "Output spiral indices: ";
  for (int i=0; i<12; ++i){
    cout <<  pentagon_indices_output[i] << ", ";
  }
  cout << endl;
  
  cout << "Jumps: ";
  for (list<pair<node_t,int> >::iterator it=jumps.begin(); it!=jumps.end(); ++it){
    cout <<  it->first << ": " << it->second << ", ";
  }
  cout << endl;

  return 0;
}
