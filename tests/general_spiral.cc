#include "libgraph/cubicgraph.hh"
#include "libgraph/fullerenegraph.hh"

using namespace std;

int main(int ac, char **av)
{

  std::cout << "Attempting to create graph from spiral indices ..." << std::endl;
 
  int n = 192;
//  int n = 16;
//no jump
//  int pentagon_indices_array[] = {1, 2, 3, 4, 5, 7, 10, 12, 13, 14, 15, 16};
//jump
//  int pentagon_indices_array[] = {1, 2, 3, 5, 8, 9, 11, 12, 13, 14, 15, 16};
//jump, 380
  int pentagon_indices_array[] = {1, 2, 3, 137, 138, 147, 148, 157, 158, 166, 173, 180};
  std::vector<int> pentagon_indices_input(12);
  std::vector<int> pentagon_indices_output;
  std::vector<pair<int,int> > jumps;
  bool ipr=false;
  bool general=true;

  std::cout << "Input spiral indices: ";
  for (int i=0; i<12; ++i){
    pentagon_indices_input[i] = pentagon_indices_array[i];
    std::cout <<  pentagon_indices_input[i] << ", ";
  }
  std::cout << std::endl;

  FullereneGraph fg(n, pentagon_indices_input, ipr, general);

  cout << "fg = " << fg << endl;



  std::cout << "Attempting to create spiral from (previously created) graph ..." << std::endl;

  fg.get_pentagon_indices(0, 1, 2, pentagon_indices_output, jumps);  

  std::cout << "Input spiral indices: ";
  for (int i=0; i<12; ++i){
    std::cout <<  pentagon_indices_output[i] << ", ";
  }
  std::cout << std::endl;
  


  return 0;
}
