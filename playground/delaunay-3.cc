#include "libgraph/triangulation.hh"
#include "libgraph/delaunay.hh"
#include "libgraph/polyhedron.hh"
#include "libgraph/fullerenegraph.hh"
#include "libgraph/debug.hh"

#include <fstream>

Triangulation regular_polygon(int n)
{
  vector<tri_t> triangles(n);
  Graph P(n+2,vector<int>((n+2)*(n+2),0));;
  node_t v = n, w = n+1;  
  for(node_t u=0;u<n;u++) P.insert_edge(edge_t(u,(u+1)%n)); // Connect cycle
  for(node_t u=0;u<n;u++){ P.insert_edge(edge_t(u,v)); P.insert_edge(edge_t(u,w)); }   // Connect centers
  
  return Triangulation(P);
}


int testN = 24;
vector<int> testRSPI = {1,2,3,4,5,7,8,10,11,12,13,14};

int main(int ac, char **av) 
{
  int N;
  vector<int> RSPI(12);

  Debug::channel_verbosity["all"] = Debug::INFO1;
  Debug warning("main",Debug::WARNING);
  Debug info("main",Debug::INFO1);

  if(ac==14){
    N = strtol(av[1], 0, 0);
    for (int i = 0; i < 12; i++){
      RSPI[i] = strtol(av[i + 2], 0, 0) - 1;
    }
  } else {
    warning << "Using test fullerene C" << testN << " " << testRSPI << "\n";
    N = testN;
    RSPI = testRSPI;
    for(int i=0;i<12;i++) RSPI[i]--;
  }

  string filename = "output/reduced-graph-C"+to_string<int>(N)+".m";
  ofstream output(filename);

  vector<int> spiral(N/2+2, 6);
  for (int i = 0; i < 12; i++){
    spiral[RSPI[i]] = 5;
  }

  cout << "spiral = " << spiral << endl;

  Triangulation T1(spiral);
  FulleroidDelaunay DY(T1);

  //  Debug::channel_verbosity["Delaunay"] = Debug::INFO3;

  {
    vector<node_t> hole(DY.neighbours[N-1]);
    int n = hole.size();

    vector<double> angles(n);
    for(node_t u=0;u<n;u++) angles[u] = DY.angle_d6y(u,n,(u+1)%n);

    vector<double> distances = DY.new_distances_d6y(n,hole);

    cout << "angles    = " << angles << ";\n";
    cout << "distances = " << distances << ";\n";
    
  }
  

  return 0;
}
