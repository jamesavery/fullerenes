#include <limits.h>
#include "libgraph/planargraph.hh"
#include "libgraph/triangulation.hh"

#include <iostream>
#include <chrono>		// For wall-time
#include <ctime>		// For CPU-time

using namespace std;
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
  int N = graphs[0].N, NF = duals[0].N;
  vector<int> spiral(NF);
  jumplist_t jumps;

  vector<jumplist_t> all_jumps(duals.size());

  auto wall_start = chrono::steady_clock::now();
  auto cpu_start  = clock();
  for(int i=0;i<duals.size();i++){
    duals[i].get_spiral(spiral,jumps,true,true,true,true);
    all_jumps[i] = jumps;
  }
  auto cpu_end = clock();
  auto wall_end = chrono::steady_clock::now();
  double walltime = chrono::duration<double,std::ratio<1,1>>(wall_end-wall_start).count();
  double cputime  = (cpu_end-cpu_start) * 1.0 / CLOCKS_PER_SEC;
  cerr << "Computed all spirals.\n";
  cout << "N        = " << N << ";\n";
  cout << "Nisomers = " << duals.size() << ";\n";
  cout << "walltime = " << walltime << ";\n";
  cout << "cputime  = " << cputime << ";\n";
  cout << "rate     = " << cputime/(duals.size()*N) << ";\n";
  cout << "alljumps = " << all_jumps << ";\n";
  
  return 0;
}
