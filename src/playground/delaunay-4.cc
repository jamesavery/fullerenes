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

Polyhedron fullerene_dual_polyhedron(const Triangulation& dg)
{
  PlanarGraph pg(dg.dual_graph());
  cout << "pg = " << pg << endl;

  FullereneGraph g(pg);
  g.layout2d = g.tutte_layout();

  vector<coord3d> points = g.zero_order_geometry();
  points = g.optimized_geometry(points);

  vector<coord3d> dual_points(dg.N);

  vector<face_t> faces(dg.N);
  for(int i=0;i<dg.triangles.size();i++)
    for(int j=0;j<3;j++)
      faces[dg.triangles[i][j]].push_back(i);

  for(int i=0;i<faces.size();i++)
    dual_points[i] = faces[i].centroid(points);

  return Polyhedron(dg, dual_points);
}

int testN = 24;
vector<int> testRSPI = {1,2,3,4,5,7,8,10,11,12,13,14};


int main(int ac, char **av) 
{
  int N, M;
  vector<int> RSPI(12);

  Debug::channel_verbosity["all"] = Debug::INFO1;
  Debug::channel_verbosity["Delaunay"] = Debug::INFO1;
  Debug warning("main",Debug::WARNING);
  Debug info("main",Debug::INFO1);
  
  ofstream steps_file("output/delaunay-steps.m");
  MathematicaDebug::channel_stream["Delaunay"] = &steps_file;
  MathematicaDebug steps("Delaunay",0);

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

  vector<int> spiral(N/2+2, 6);
  for (int i = 0; i < 12; i++){
    spiral[RSPI[i]] = 5;
  }

  cout << "spiral = " << spiral << endl;

  Triangulation T1(spiral);
  FulleroidDelaunay DY(T1);
  cout << "DY: " << DY << endl;
  M = DY.N;

  Polyhedron P0 = fullerene_dual_polyhedron(DY);
  cout << "P0: " << P0 << endl;

  steps << "g[0] = "  << DY << ";\n";
  steps << "points[0] = " << P0.points << ";\n";

  cout << "starting to remove vertices " << endl;
  DY.remove_flat_vertices();

  // Polyhedron P1 = fullerene_dual_polyhedron(DY);

  return 0;
}
