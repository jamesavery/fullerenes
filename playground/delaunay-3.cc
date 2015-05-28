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
  M = DY.N;
  //  Debug::channel_verbosity["Delaunay"] = Debug::INFO3;

  {
    node_t v = M-1;
    vector<node_t> hole(DY.neighbours[v]);
    int n = hole.size();

    cout << "hole = " << hole << ";\n";

    vector<double> angles(n);

    for(node_t u=0;u<n;u++) angles[u] = DY.angle_d6y(hole[u],v,hole[(u+1)%n]);

    vector<double> distances = DY.new_distances_d6y(v,hole);

    cout << "angles    = " << angles << ";\n";
    cout << "distances = " << distances << ";\n";


    output << "g0 = " << DY << ";\n";
    cout   << "degree"<<(v+1)<<" = " << DY.neighbours[v].size();
    DY.remove_flat_vertex(v);
    output << "g1 = " << DY << ";\n";
    //    Debug::channel_verbosity["Delaunay"] = Debug::WARNING;

    while(v>12){
      v--;

      vector<node_t> hole(DY.neighbours[v]);
      int n = hole.size();
      cout << "degree"<<(v+1)<<" = " << DY.neighbours[v].size() <<";\n";
      cout << "hole = " << hole << ";\n";

      vector<double> angles(n);
      for(node_t u=0;u<n;u++) angles[u] = DY.angle_d6y(hole[u],v,hole[(u+1)%n]);

      cout << "angles    = " << angles << ";\n";
      //      cout << "distances = " << distances << ";\n";
      if(v+1==16) Debug::channel_verbosity["Delaunay"] = Debug::INFO3;
      DY.remove_flat_vertex(v);
      if(v+1==16) Debug::channel_verbosity["Delaunay"] = Debug::WARNING;
      output << "g"<<(M-v)<<" = " << DY << ";\n";
      output.flush();
    }



    
  }

  Polyhedron P0 = fullerene_dual_polyhedron(T1);

  output << "points0 = " << P0.points << ";\n";

  output.close();

  return 0;
}
