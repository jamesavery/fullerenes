// usage:  ./app-hog-to-rspi hog-file
// creates one file with all spirals in output/

#include <vector>
#include <fstream>

#include "libgraph/planargraph.hh"
#include "libgraph/triangulation.hh"
#include "libgraph/fullerenegraph.hh"

using namespace std;
typedef list<pair<int,int> > jumplist_t;


int main(int ac, char **av)
{
  string inputfile = av[1];
  FILE *input = fopen(inputfile.c_str(),"r");
  int N = 0;

  vector<int> pentagon_indices;
  jumplist_t jumps;

  // get number of graphs per file
  const int header_size = 15;	
  // Get file size
  fseek(input, 0, SEEK_END);
  size_t file_size = ftell(input);
  //find number of vertices per graph
  //this only works for files with graphs of the equal size
  fseek(input, header_size, SEEK_SET);
  // Read the number N of vertices per graph.
  fread(reinterpret_cast<char*>(&N), 1, 1, input);
  if(N == 0){
    fread(reinterpret_cast<char*>(&N), 2, 1, input);
  }
  //only for files with graphs of the equal size
  unsigned int step;
  if(N<=255)
    {step = N * 4 + 1;}
  else
    {step = N * 8 + 3;}
  unsigned int graphs_per_file = (file_size - header_size ) /step;

  cout << graphs_per_file << " graphs with " << N << " nodes found." << endl;


  for(int i=0; i!=graphs_per_file; i++){
    FullereneGraph fg = FullereneGraph(i, input);

//    cout << 0 << fg.neighbours[0][0] <<  fg.neighbours[0][1] << endl;
    fg.layout2d = fg.tutte_layout(fg.neighbours[0][0],0, fg.neighbours[0][1]);
    //fg.layout2d = fg.tutte_layout();
//    cout << "fg: " << fg << endl;

    fg.get_rspi_from_fg(pentagon_indices, jumps);
//    cout << "pentagon indices: " << pentagon_indices << endl;

    //start indices at 1 for fortran
    for(int i=0;i<12;++i) ++pentagon_indices[i];

    ofstream output;
    output.open(("output/" + to_string(N) + "-rspi").c_str(),ios::out | ios::app);
    output << pentagon_indices << endl;
    output.close();

  }
  
  return 0;
}
