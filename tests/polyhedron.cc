#include "libgraph/fullerenegraph.hh"
#include "libgraph/polyhedron.hh"
#include "libgraph/isomerdb.hh"

int testN = 80;
int testRSPI[12] = {1,2,3,4,5,6,37,38,39,40,41,42};

int main(int ac, char **av)
{
  int N;
  vector<int> RSPI(12);
  if(ac>=14) {
    N = strtol(av[1],0,0);
    for(int i=0;i<12;i++) RSPI[i] = strtol(av[i+2],0,0)-1;
  } else {
    N = testN;
    for(int i=0;i<12;i++) RSPI[i] = testRSPI[i]-1;
  }

  // int isomer_number = ac>=3? strtol(av[2],0,0) : 1;
  // bool IPR = ac>=4? strtol(av[3],0,0) : false;

  // IsomerDB DB(N,IPR);

  // FullereneGraph g(IsomerDB::makeIsomer(N,DB.getIsomer(N,isomer_number,IPR)));
  string basename("polyhedron-"+to_string(N));


  FullereneGraph g(N,RSPI);
  g.layout2d = g.tutte_layout();

  Polyhedron P0(g,g.zero_order_geometry(),6);

   {
     ofstream mol2(("output/"+basename+"-P0.mol2").c_str());
     mol2 << P0.to_mol2();
     mol2.close();
   }

   Polyhedron P(P0);
   P.optimize();

   ofstream output(("output/"+basename+".m").c_str());

   output << "g = " << g << ";\n";
   output << "coordinates0 = " << P0.points << ";\n";
   output << "coordinates = "  << P.points << ";\n";

   output << "P = " << P << ";\n";

   Polyhedron D(P.dual(6,true));
   
   D.optimize();
   output << "PD = " << D << ";\n";
  
  output.close();


  {
    ofstream mol2(("output/"+basename+".mol2").c_str());
    mol2 << P.to_mol2();
    mol2.close();
  }

  {
    P.move_to_origin();
    coord3d centre_check;
    for(int i=0;i<P.points.size();i++) centre_check += P.points[i];
    centre_check /= double(P.points.size());
    cout << "centre_check = " << centre_check << endl;

    matrix3d If(P.inertial_frame());
    P.points = If*P.points;

    ofstream mol2(("output/"+basename+"-if.mol2").c_str());
    mol2 << P.to_mol2();
    mol2.close();
  }


  {
    ofstream mol2(("output/"+basename+"-dual.mol2").c_str());
    mol2 << D.to_mol2();
    mol2.close();
  }

  return 0;
}
