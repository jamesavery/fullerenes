#include "fullerenes/spiral.hh"
#include "fullerenes/fullerenegraph.hh"
#include "fullerenes/polyhedron.hh"
#include "fullerenes/triangulation.hh"
#include "fullerenes/isomerdb.hh"

int main(int ac, char **av)
{
  int N;
  jumplist_t jumps;
  vector<int> RSPI(12);
  bool from_file = false;

  if(ac==2){			// If only one argument is given, it is a filename to read initial polyhedron from
    from_file = true;
    N = 0;
    
  } else if(ac==3) {            // Pentagon indices in quoted string
    N = strtol(av[1],0,0);
    int n = strlen(av[2]);
    for(int i=0;i<n;i++) if(av[2][i] == ',') av[2][i] = ' ';
    istringstream iss(av[2]);
    for(int i=0;i<12;i++){ iss >> RSPI[i]; RSPI[i]--; }
    cout << "RSPI = " << RSPI << endl;
  } else 
  if(ac>=14) {
    N = strtol(av[1],0,0);
    for(int i=0;i<12;i++) RSPI[i] = strtol(av[i+2],0,0)-1;
  } 


  // int isomer_number = ac>=3? strtol(av[2],0,0) : 1;
  // bool IPR = ac>=4? strtol(av[3],0,0) : false;

  // IsomerDB DB(N,IPR);

  // FullereneGraph g(IsomerDB::makeIsomer(N,DB.getIsomer(N,isomer_number,IPR)));


  vector<int> spiral(N/2+2,6);
  for(int i=0;i<12;i++) spiral[RSPI[i]] = 5;

  Polyhedron P0;
  PlanarGraph g;
  if(from_file){
    P0 = Polyhedron::from_file(av[1]);
    N = P0.N;
    g = P0;
    g.layout2d = g.tutte_layout(-1,-1,-1,8);
  } else {
    Triangulation T(spiral,jumps);
    g = T.dual_graph();
    g.layout2d = g.tutte_layout();
    P0 = Polyhedron(g,g.zero_order_geometry(),6);
  }

  string basename("polyhedron-"+to_string(N));
  Polyhedron::to_file(P0,"output/"+basename+"-P0.mol2");

  printf("P0\n");
  Polyhedron P(P0);
  printf("Optimizing P\n");  
  P.optimize();
  printf("Writing P\n");
  Polyhedron::to_file(P,"output/"+basename+".mol2");

  printf("Aligning P\n");  
  P.move_to_origin();
  P.align_with_axes();

  printf("Aligning P-aligned\n");    
  Polyhedron::to_file(P,"output/"+basename+"-if.mol2");
  Polyhedron::to_file(P,"output/"+basename+"-if.xyz");
  //  Polyhedron::to_file(P,"output/"+basename+"-if.pov");

  ofstream output(("output/"+basename+".m").c_str());

  printf("Writing mathematica version of P\n");  
  vector<face_t> facemap(g.compute_faces());
  output << "g = " << g << ";\n";
  output << "coordinates0 = " << P0.points << ";\n";
  output << "coordinates = "  << P.points << ";\n";
  output << "pentagons = " << facemap[5] << ";\n"
  	  << "hexagons  = " << facemap[6] << ";\n"
  	  << "RSPI = " << RSPI << ";\n";

  output << "P0 = " << P0 << ";\n";
  output << "P = " << P << ";\n";

  Polyhedron D(P.dual());
  D.layout2d = D.tutte_layout();
  D.faces    = D.compute_faces(3,true);
  D.face_max = 3;
  D.optimize();
  output << "PD = " << D << ";\n";
  
  output.close();

  Polyhedron::to_file(P,"output/"+basename+"-dual.mol2");
  //  Polyhedron::to_file(P,"output/"+basename+"-dual.pov");

  return 0;
}
