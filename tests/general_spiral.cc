#include "libgraph/cubicgraph.hh"
#include "libgraph/fullerenegraph.hh"

using namespace std;

int main(int ac, char **av)
{

  std::cout << "Attempting to create graph from spiral indices ..." << std::endl;
 
  int n = 380;
//  int n = 28;
//no jump
//  int pentagon_indices_array[] = {0, 1, 2, 3, 4, 6, 9, 11, 12, 13, 14, 15};
//jump
//  int pentagon_indices_array[] = {0, 1, 2, 4, 7, 8, 10, 11, 12, 13, 14, 15};
//jump, 380
  int pentagon_indices_array[] = {0, 1, 2, 136, 137, 146, 147, 156, 157, 165, 172, 179};
  int jumps_array[][2] = {{165,2}, {172,3}, {179,3}};
  int n_jumps = 3;

  std::vector<int> pentagon_indices_input(12);
  std::vector<int> pentagon_indices_output;
  std::list<pair<int,int> > jumps;
  //bool ipr=false;

  std::cout << "Input spiral indices: ";
  for (int i=0; i<12; ++i){
    pentagon_indices_input[i] = pentagon_indices_array[i];
    std::cout <<  pentagon_indices_input[i] << ", ";
  }
  std::cout << std::endl;

  std::cout << "Input jumps: ";
  for (int i=0; i<n_jumps; ++i){
    jumps.push_back(make_pair(jumps_array[i][0], jumps_array[i][1]));
    std::cout <<  jumps.rbegin()->first << ", " << jumps.rbegin()->second << "; " ;
  }
  std::cout << std::endl;

//  FullereneGraph fg(n, pentagon_indices_input, jumps);
  FullereneGraph fg(n, pentagon_indices_input);

  cout << "fg = " << fg << endl;



  std::cout << "Attempting to create spiral from (previously created) graph ..." << std::endl;

  fg.get_pentagon_indices(0, 1, 2, pentagon_indices_output, jumps);  

  std::cout << "Output spiral indices: ";
  for (int i=0; i<12; ++i){
    std::cout <<  pentagon_indices_output[i] << ", ";
  }
  std::cout << std::endl;
  
  std::cout << "Jumps: ";
  for (std::list<pair<node_t,int> >::iterator it=jumps.begin(); it!=jumps.end(); ++it){
    std::cout <<  it->first << ": " << it->second << ", ";
  }
  std::cout << std::endl;
  


  return 0;
}
