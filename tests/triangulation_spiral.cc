#include "libgraph/fullerenegraph.hh"
#include "libgraph/triangulation.hh"

using namespace std;

int main(int ac, char **av)
{

  std::cout << "Attempting to create graph from spiral indices ..." << std::endl;
 
// spiral indices are input starting at 0
  // D3 C1152 -- leapfrog of NS-C384
  //int n = 1152/2+2;
  //int pentagon_indices_array[] = {1,11,13,24,453,496,498,542,563,572,574,577};
  //int n = 4608/2+2;
  //  int pentagon_indices_array[] = {51,108,112,272,1746,1921,2005,2096,2246,2284,2288,2304};
//jump, 380
//  int n = 380;
//  int pentagon_indices_array[] = {0, 1, 2, 136, 137, 146, 147, 156, 157, 165, 172, 179};
//  int jumps_array[][2] = {{165,2}, {172,3}, {179,3}};
//  int n_jumps = 3;
//jump, 384
  int n = 384;
  int pentagon_indices_array[] = {133,134,135, 147,148,149, 161,162,163, 171,179,186};
  int jumps_array[][2] = {{171,7}, {179,7}, {186,8}};
  int n_jumps = 3;
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

  vector<int> pentagon_indices_input(pentagon_indices_array,pentagon_indices_array+12);
  vector<int> pentagon_indices_output,pentagon_indices_output2;
  list<pair<int,int> > jumps_input;
  list<pair<int,int> > jumps_output,jumps_output2;
  
  cout << "N: " << n << endl
       << "Input spiral indices: " << pentagon_indices_input << endl;

  cout << "Input jumps: ";
  for (int i=0; i<n_jumps; ++i){
    jumps_input.push_back(make_pair(jumps_array[i][0], jumps_array[i][1]));
    cout <<  jumps_input.rbegin()->first << ", " << jumps_input.rbegin()->second << "; " ;
  }
  cout << endl;

  FullereneGraph fg(n, pentagon_indices_input, jumps_input);

  fg.get_canonical_general_spiral_from_fg(pentagon_indices_output,jumps_output);

  fg.layout2d = fg.tutte_layout();
  Triangulation dual(fg.dual_graph(6,true));
  
  FullereneDual(dual).get_canonical_fullerene_rspi(pentagon_indices_output2,jumps_output2);

  cout << "FullereneGraph general RSPI: " << jumps_output  << "; " << pentagon_indices_output << endl;
  cout << "FullereneDual  general RSPI: " << jumps_output2 << "; " << pentagon_indices_output2 << endl;

  return 0;
}
