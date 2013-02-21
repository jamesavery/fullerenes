#include "libgraph/geometry.hh"
#include "libgraph/planargraph.hh"
#include "libgraph/cubicgraph.hh"
#include "libgraph/fullerenegraph.hh"

using namespace std;

int main(int ac, char **av)
{

  std::cout << "Attempting to create graph from spiral indices ..." << std::endl;
 
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
//    int pentagon_indices_array[] = {27, 47, 48, 87, 88, 114, 132, 133, 154, 180, 181, 191};
//    int jumps_array[][2] = {{114, 2}};
//    int n_jumps = 1;
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
//  //0, 115, 151

  vector<int> pentagon_indices_input(pentagon_indices_array,pentagon_indices_array+12);
  vector<int> pentagon_indices_output;
  list<pair<int,int> > jumps;
  
  cout << "N: " << n << endl
       << "Input spiral indices: " << pentagon_indices_input << endl;

  cout << "Input jumps: ";
  for (int i=0; i<n_jumps; ++i){
    jumps.push_back(make_pair(jumps_array[i][0], jumps_array[i][1]));
    cout <<  jumps.rbegin()->first << ", " << jumps.rbegin()->second << "; " ;
  }
  cout << endl;

  vector<int> pentagon_indices_output_bak(pentagon_indices_input);
  list<pair<int,int> > jumps_bak(jumps);

  FullereneGraph fg(n, pentagon_indices_input, jumps);
  PlanarGraph dual(fg.dual_graph(6));
  vector<face_t> faces(dual.compute_faces_flat(3));

  for(int i=0;i<faces.size();i++){ 
    int permutations[6][3] = {{0,1,2},{0,2,1},{1,0,2},{1,2,0},{2,0,1},{2,1,0}};
    const face_t& f = faces[i];
    for(int j=0;j<6;j++){

      int f1 = f[permutations[j][0]], f2 = f[permutations[j][1]], f3 = f[permutations[j][2]];

      fg.get_pentagon_indices(f1, f2, f3, pentagon_indices_output, jumps);

      printf("Face %d:%d pentagon_indices(%d,%d,%d)\n",i,j,f1,f2,f3);
      for(list<pair<int,int> >::iterator it(jumps.begin()); it!= jumps.end(); ++it){cout << it->first << ", " << it->second << ", ";}
      for(int i=0; i<12;++i){cout << pentagon_indices_output[i] << ", ";}
      cout << endl;
    }

  if(jumps.size() < jumps_bak.size()){jumps_bak = jumps; pentagon_indices_output_bak = pentagon_indices_output;}

  }

  cout << "finally: " << endl;
  for(list<pair<int,int> >::iterator it(jumps_bak.begin()); it!= jumps_bak.end(); ++it){cout << it->first << ", " << it->second << ", ";}
  for(int i=0; i<12;++i){cout << pentagon_indices_output_bak[i] << ", ";}
  cout << endl;

  return 0;
}
