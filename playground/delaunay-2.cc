#include "libgraph/triangulation.hh"
#include "libgraph/delaunay.hh"
#include "libgraph/polyhedron.hh"
#include "libgraph/fullerenegraph.hh"

#include <fstream>

Polyhedron fullerene_dual_polyhedron(const Triangulation& dg)
{
  PlanarGraph pg(dg.dual_graph());
  cout << "pg = " << pg << endl;

  FullereneGraph g(pg);
  g.layout2d = g.tutte_layout();

  vector<coord3d> points = g.zero_order_geometry();
  points = g.optimised_geometry(points);

  vector<coord3d> dual_points(dg.N);

  vector<face_t> faces(dg.N);
  for(int i=0;i<dg.triangles.size();i++)
    for(int j=0;j<3;j++)
      faces[dg.triangles[i][j]].push_back(i);

  for(int i=0;i<faces.size();i++)
    dual_points[i] = faces[i].centroid(points);

  return Polyhedron(dg, dual_points);
}

int main(int ac, char **av) {
  int N;
  vector<int> RSPI(12);
  if(ac!=14) cout << "questionable number of input parameters: " << ac << endl;
  N = strtol(av[1], 0, 0);
  for (int i = 0; i < 12; i++){
    RSPI[i] = strtol(av[i + 2], 0, 0) - 1;
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

  cout << "DY = " << DY << ";\n"
       << "dDist = " << DY.distances << ";\n";
  output << "DY = " << DY << ";\n"
	 << "dDist = " << DY.distances << ";\n";

//  Polyhedron PT = fullerene_dual_polyhedron(DY);

//  output << "PT = " << PT << "\n";
  output.close();

  DY.remove_flat_vertices();
  cout << "rT = " << DY << ";\n";

  map<edge_t,double> edge_lengths;
  for (int i=0; i<12;i++){
    for (int j=i; j<12;j++){
      if(DY.edge_lengths_d6y(i,j)){
        edge_lengths.insert(make_pair(edge_t(i,j), DY.edge_lengths_d6y(i,j)));
      }
    }
  }
//  cout << edge_lengths << endl;


 
  cout << "creating P from DY" << endl;
  DY.layout2d = DY.tutte_layout();
  Polyhedron P = Polyhedron(DY,DY.zero_order_geometry(),6);
  cout << "created P from DY" << endl;
  cout << "P=" << P << endl;


  cout << "optimizing" << endl;
  P.optimise_other(false, edge_lengths);
  cout << "optimised" << endl;


  string mol2_name = "output/reduced-graph-C"+to_string(N)+".mol2";
  {
    ofstream mol2(mol2_name);
    mol2 << P.to_mol2();
    mol2.close();
  }

  
  return 0;
}
