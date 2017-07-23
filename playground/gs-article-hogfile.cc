#include <limits.h>
#include "libgraph/planargraph.hh"
#include "libgraph/triangulation.hh"

#include <iostream>
#include <chrono>

typedef Triangulation::jumplist_t jumplist_t;

vector<int> rspi_from_fullerene_spiral(const vector<int>& spiral)
{
  vector<int> rspi(12);
  for(int i=0,j=0;i<spiral.size();i++){
    if(spiral[i] == 5) rspi[j++] = i+1;
    assert(j<=12);
  }
  return rspi;
}

int main(int ac, char **av)
{
  const char *filename = av[1];

  FILE *file = fopen(filename,"rb");
  if(!file){
    perror((string("Error opening ")+filename).c_str());
  }
  vector<Graph> graphs(PlanarGraph::read_hog_planarcodes(file));
  fclose(file);

  cerr << "Successfully read " << graphs.size() << " planar graphs.\n";
  
  vector<FullereneDual> duals(graphs.size());

  for(int i=0;i<graphs.size();i++){
    //    cerr << "Dualizing graph g"<<(++i)<<".\n";
    duals[i] = Triangulation(PlanarGraph(graphs[i]).dual_graph(6,false),false);
    //    cerr << "dual"<<i<<" = "<<duals[i-1]<<";\n";
  }
  cerr << "Successfully dualized " << duals.size() << " planar graphs.\n";
  //  cout << "graphs = " << graphs << ";\n";
  //  cout << "duals = " << duals << ";\n";
  //  return 0;
  vector<int> spiral(duals[0].N);
  jumplist_t jumps;

  vector<jumplist_t> all_jumps(duals.size());

  auto start = chrono::steady_clock::now();
  for(int i=0;i<duals.size();i++){
    duals[i].get_spiral(spiral,jumps,true,true,true,true);
    all_jumps[i] = jumps;
  }
  auto end = chrono::steady_clock::now();
  cerr << "Computed all spirals.\n";
  cout << "Nisomers = " << duals.size() << ";\n";
  cout << "timing   = " << chrono::duration<double,std::ratio<1,1>>(end-start).count() << ";\n";
  cout << "alljumps = " << all_jumps << ";\n";
  
  return 0;
}
