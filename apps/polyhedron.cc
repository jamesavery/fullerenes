#include "libgraph/fullerenegraph.hh"
#include "libgraph/polyhedron.hh"
#include "libgraph/triangulation.hh"
#include "libgraph/isomerdb.hh"

int testN = 80;
int testRSPI[12] = {1,2,3,4,5,6,37,38,39,40,41,42};

int main(int ac, char **av)
{
  int N;
  jumplist_t jumps;
  vector<int> RSPI(12);

//  bool from_file = false;
//  if(ac==2){
//    from_file = true;
//    N = 0;
//  } else 
  if(ac<14){
    N = testN;
    for(int i=0;i<12;i++) RSPI[i] = testRSPI[i]-1;
  }
  if(ac>=14) {
    N = strtol(av[1],0,0);
    for(int i=0;i<12;i++) RSPI[i] = strtol(av[i+2],0,0)-1;
  } 
  if(ac>14){			// General RSPI: RSPI followed by jumps.
    for(int i=14;i<ac;i+=2)
      jumps.push_back(make_pair(strtol(av[i],0,0)-1,strtol(av[i+1],0,0)));
  }

  // int isomer_number = ac>=3? strtol(av[2],0,0) : 1;
  // bool IPR = ac>=4? strtol(av[3],0,0) : false;

  // IsomerDB DB(N,IPR);

  // FullereneGraph g(IsomerDB::makeIsomer(N,DB.getIsomer(N,isomer_number,IPR)));


  vector<int> spiral(N/2+2,6);
  for(int i=0;i<12;i++) spiral[RSPI[i]] = 5;

  Polyhedron P0;
  PlanarGraph g;
//  if(from_file){
//    P0 = Polyhedron(av[1]);
//    N = P0.N;
//    g = P0;
//    g.layout2d = g.tutte_layout(-1,-1,-1,8);
//  } else
  {
    Triangulation T(spiral,jumps);
    
    g = T.dual_graph();
    g.layout2d = g.tutte_layout();
    P0 = Polyhedron(g,g.zero_order_geometry(),6);
  }

  string basename("polyhedron-"+to_string(N));
  {
    FILE *file = fopen(("output/"+basename+"-P0.mol2").c_str(),"wb");
    Polyhedron::to_mol2(P0, file);
    fclose(file);
  }

  Polyhedron P(P0);
  P.optimize();

  {
    FILE *file = fopen(("output/"+basename+"-P.mol2").c_str(),"wb");
    Polyhedron::to_mol2(P, file);
    fclose(file);
  }

  {
    P.move_to_origin();
    P.align_with_axes();

    FILE *file = fopen(("output/"+basename+"-if.mol2").c_str(),"wb");
    Polyhedron::to_mol2(P, file);
    fclose(file);

    ofstream pov(("output/"+basename+"-if.pov").c_str());
    pov << P.to_povray();
    pov.close();
  }

  ofstream output(("output/"+basename+".m").c_str());

  vector<face_t> faces(g.compute_faces(8,true));
  output << "g = " << g << ";\n";
  output << "coordinates0 = " << P0.points << ";\n";
  output << "coordinates = "  << P.points << ";\n";
  output << "faces = " << faces << ";\n"
  	  << "RSPI = " << RSPI << ";\n";

  output << "P0 = " << P0 << ";\n";
  output << "P = " << P << ";\n";

  // Polyhedron D(P.dual(6,true));
  // D.layout2d = D.tutte_layout();
  // D.faces    = D.compute_faces(3,true);
  // D.face_max = 3;
  //   D.optimize();
  // output << "PD = " << D << ";\n";
  
  output.close();


  // {
  //   FILE *file = fopen(("output/"+basename+"-dual.mol2").c_str(),"wb");
  //   Polyhedron::to_mol2(D, file);
  //   fclose(file);

  //   ofstream pov(("output/"+basename+"-dual.pov").c_str());
  //   pov << D.to_povray();
  //   pov.close();
  // }

  return 0;
}
