#include <vector>
#include <fstream>

#include "libgraph/planargraph.hh"
#include "libgraph/triangulation.hh"
#include "libgraph/fullerenegraph.hh"

using namespace std;
typedef list<pair<int,int> > jumplist_t;

void get_rspi(const vector<int>& spiral, uint8_t rspi[12])
{
  for(int i=0, j=0;i<spiral.size();i++) if(spiral[i] == 5) rspi[j++] = i;
}


int main(int ac, char **av)
{
  string inputfile = av[1];
  FILE *input = fopen(inputfile.c_str(),"r");
  string outputfile = av[2];

  vector<int> pentagon_indices;
  jumplist_t jumps;

  for(int i=0; i!=1; i++){
    FullereneGraph fg = FullereneGraph(i, input);

    cout << 0 << fg.neighbours[0][0] <<  fg.neighbours[0][1] << endl;
    fg.layout2d = fg.tutte_layout(fg.neighbours[0][0],0, fg.neighbours[0][1]);
    //fg.layout2d = fg.tutte_layout();
    cout << "fg: " << fg << endl;

    fg.get_canonical_general_spiral_from_fg(pentagon_indices, jumps);
    cout << "pentagon indices: " << pentagon_indices << endl;

    ofstream output(("output/" + outputfile).c_str());
    output << pentagon_indices << endl;
    output.close();

  }
  
  return 0;
}
