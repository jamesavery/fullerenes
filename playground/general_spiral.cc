// compile with lapack and with gsl

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
  int n = 380;
//  int pentagon_indices_array[] = {0, 1, 2, 136, 137, 146, 147, 156, 157, 165, 172, 179};
  int pentagon_indices_array[] = {44, 69, 70, 81, 82, 109, 118, 119, 143, 183, 184, 191};
//  int jumps_array[][2] = {{165,2}, {172,3}, {179,3}};
  int jumps_array[][2] = {{109,2}};
//  int n_jumps = 3;
  int n_jumps = 1;
  //0,1,2
//jump, 440
//  int n = 440;
//  int pentagon_indices_array[] = {159, 160, 161, 174, 175, 176, 189, 190, 191, 195, 204, 212};
//  int jumps_array[][2] = {{192,6}, {195,6}, {204,6}, {212,7}};
//  int n_jumps = 4;
//  //0, 1, 9
//jump, 672
//  int n = 672;
//  int pentagon_indices_array[] = {140, 142, 144, 161, 163, 165, 225, 242, 328, 330, 333, 335};
//  int jumps_array[][2] = {{225,5}, {242,6}};
//  int n_jumps = 2;
//  //0, 115, 151

  vector<int> rspi_in(12);
  vector<int> rspi_out;
  list<pair<int,int> > jumps;

  cout << "Input spiral indices: ";
  for (int i=0; i<12; ++i){
    rspi_in[i] = pentagon_indices_array[i];
    cout <<  rspi_in[i] << ", ";
  }
  cout << endl;

  cout << "Input jumps: ";
  for (int i=0; i<n_jumps; ++i){
    jumps.push_back(make_pair(jumps_array[i][0], jumps_array[i][1]));
    cout <<  jumps.rbegin()->first << ", " << jumps.rbegin()->second << "; " ;
  }
  cout << endl;

  FullereneGraph fg(n, rspi_in, jumps);
  cout << "graph created" << endl;
  fg.layout2d = fg.tutte_layout();
  cout << "layout created" << endl;

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
  rspi_out.clear(); 
  bool pentagon_start = false;
  if(!fg.get_rspi_from_fg(rspi_out, jumps, true, pentagon_start)){
    cerr << "no spiral found, which is not possible. " << endl;
    return 1;
  }

  cout << "Output spiral: ";
  vector<int> out_spiral(n/2+2,6);
  cerr << out_spiral << endl;
  for (int i=0; i<12; ++i){
    out_spiral[rspi_out[i]] = 5;
  }
  cout << out_spiral << endl;
  
  cout << "Output spiral indices: ";
  cout << rspi_out << endl;
  
  cout << "Jumps: ";
  cout << jumps << endl;

  return 0;
}
