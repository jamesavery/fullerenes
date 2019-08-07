#include <stdlib.h>
#include <fstream>
#include "libgraph/cubicgraph.hh"
#include "libgraph/fullerenegraph.hh"

#include <vector>

using namespace std;

extern "C" void optgraph_(const int *N, const int *iop, const int *iout, 
			  const int *ida,  
			  const int *is, const int *mdist, const int *maxl, 
			  const double *scalePPG, double *xs);

vector<coord2d> OptGraph(const PlanarGraph &g, int method, double scalePPG)
{
  vector<double> xs(2*g.N);
  int iout=3;
  vector<int> A(g.N*g.N), IS(6);	// Check that this is zero filled.

  // Fill in adjacency matrix
  for(set<edge_t>::const_iterator e(g.edge_set.begin()); e!=g.edge_set.end(); e++){
    A[e->first*g.N+e->second] = 1;
    A[e->second*g.N+e->first] = 1;
  }

  // Compute topological distance matrix and maximum distance
  int maxdist = INT_MIN;
  vector<int> dist(g.all_pairs_shortest_paths(g.N));
  for(int i=0;i<g.N*g.N;i++) if(dist[i]>maxdist) maxdist = dist[i];

  cout << "maxdist = " << maxdist << ";\n";

  // Copy outer face to IS
  assert(g.outer_face.size() >= 5);
  for(int i=0;i<g.outer_face.size();i++) IS[i] = g.outer_face[i]+1;

  // Copy coordinates to Fortran order
  for(unsigned int i=0;i<g.N;i++){
    xs[i*2]   = g.layout2d[i].first;
    xs[i*2+1] = g.layout2d[i].second;
  }



  optgraph_(&g.N, &method,&iout,&A[0],&IS[0],&dist[0],&maxdist,&scalePPG,&xs[0]);

  // Copy coordinates back
  vector<coord2d> result(g.N);
  for(int i=0;i<g.N;i++)
    result[i] = coord2d(xs[i*2],xs[i*2+1]);

  cout << "Any changes? " << (result == g.layout2d? "No.":"Yes!") << endl;
  
  return result;
}


int main(int ac, char **av)
{
  vector<int> rspi(12);
  FullereneGraph::jumplist_t jumps;

  if(ac<14) return -1;

  int N = strtol(av[1],0,0);
  for(int i=0;i<12;i++) rspi[i] = strtol(av[2+i],0,0)-1;
  
  int Nleap = ac>=15? strtol(av[14],0,0) : 1;

  cout << "Attempting to create graph from spiral indices " << rspi << endl;
  FullereneGraph g(N, rspi, jumps);
  //  PlanarGraph  leap(g.dual_graph(6));

  g.layout2d = g.tutte_layout();
  //  g.layout2d = OptGraph(g,3,1.9);
  //  g.layout2d = OptGraph(g,1,1.4);

  ofstream mathfile(string("output/leapfrog"+to_string(Nleap)+"-C" + to_string(N)+".m").c_str());
  mathfile << "g = " << g << ";\n";

  cout << "Computing " << Nleap << (Nleap==1?"st":"th") << " leapfrog.\n";
  FullereneGraph leap(g.leapfrog_fullerene(true));
  for(int i=1;i<Nleap;i++) leap = leap.leapfrog_fullerene(true);


  mathfile << "leap = " << leap << ";\n";
  mathfile.close();

  ofstream g_latex("output/lf-g.tex"), leap_latex("output/lf-leap.tex");
  
  g_latex  << g.to_latex(20,20,false,false,false,0,0,0xffffff) << endl;
  leap_latex << leap.to_latex(20,20,false,false,false,0,0,0xffffff) << endl;

  g_latex.close(); 
  leap_latex.close();
  
  return 0;
}
