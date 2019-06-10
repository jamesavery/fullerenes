#include "libgraph/triangulation.hh"
#include "libgraph/polyhedron.hh"
#include "libgraph/isomerdb.hh"


int main(int ac, char **av)
{
//  const int N = 16;
//  int is[N] = {7,5,5,5,5, 5,5,5,5,5, 5,5,5,5,5, 7};
  const int N = 62;
  int is[N] = {5, 6,6,6,6,6, 5,6,5,6,5,6,5,6,5,6, 6,5,6,5,6,5,6,5,6,5, 7,7,7,7,7, 7,7,7,7,7, 5,6,5,6,5,6,5,6,5,6, 6,5,6,5,6,5,6,5,6,5, 6,6,6,6,6, 5};
//  const int N = 86;
//  int is[N] = {6,6,6,6,6,6,6, 5,6,5,6,5,6,5,6,5,6,5,6, 6,5,6,5,6,5,6,5,6,5,6,5, 7,7,7,7,7,7, 5,7,5,7,5,7, 7,5,7,5,7,5, 7,7,7,7,7,7, 5,6,5,6,5,6,5,6,5,6,5,6, 6,5,6,5,6,5,6,5,6,5,6,5, 6,6,6,6,6,6,6};
  vector<int> s(is,is+N);

  string basename("peanut-"+to_string(N));

  Triangulation T(s);
  cout << "T created" << endl;

  CubicGraph g(PlanarGraph(T).dual_graph(3));
  g.layout2d = g.tutte_layout();

  Polyhedron P0(g,g.zero_order_geometry(),6);
  cout << "P0 created" << endl;

   {
     ofstream mol2(("output/"+basename+"-P0.mol2").c_str());
     mol2 << P0.to_mol2();
     mol2.close();
   }

   Polyhedron P(P0);
   P.optimize();

  {
    ofstream mol2(("output/"+basename+".mol2").c_str());
    mol2 << P.to_mol2();
    mol2.close();
  }

//  {
//    P.move_to_origin();
//    matrix3d If(P.inertial_frame());
//    P.points = If*P.points;
//
//    ofstream mol2(("output/"+basename+"-if.mol2").c_str());
//    mol2 << P.to_mol2();
//    mol2.close();
//
//    ofstream pov(("output/"+basename+"-if.pov").c_str());
//    pov << P.to_povray();
//    pov.close();
//  }
//
//   ofstream output(("output/"+basename+".m").c_str());
//
//   facemap_t facemap(g.compute_faces(6,true));
//   output << "g = " << g << ";\n";
//   output << "coordinates0 = " << P0.points << ";\n";
//   output << "coordinates = "  << P.points << ";\n";
//   output << "pentagons = " << facemap[5] << ";\n"
//	  << "hexagons  = " << facemap[6] << ";\n";
//
//   output << "P = " << P << ";\n";
//
//   Polyhedron D(P.dual(6,true));
//   D.layout2d = D.tutte_layout();
//   D.faces    = D.compute_faces_flat(3,true);
//   D.face_max = 3;
//   //   D.optimize();
//   output << "PD = " << D << ";\n";
//  
//  output.close();
//
//
//  {
//    ofstream mol2(("output/"+basename+"-dual.mol2").c_str());
//    mol2 << D.to_mol2();
//    mol2.close() ;
//    ofstream pov(("output/"+basename+"-dual.pov").c_str());
//    pov << D.to_povray();
//    pov.close();
//  }

  return 0;
}
