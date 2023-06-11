#include "libgraph/fullerenegraph.hh"
#include "libgraph/triangulation.hh"
#include "libgraph/polyhedron.hh"

bool compute_polyhedron = true;
int halma_level = 0;

int testN = 380;

//Canonical GSPI for C380
//int testNjumps = 1;
//int testjump[1][2] = {{110,2}};
//int testRSPI[12] = {45, 70, 71, 82, 83, 110, 119, 120, 144, 184, 185, 192};

// But in fact, we want a spiral with a long-ish starting sequence that starts in one of the corners, i.e.
// in one of the pentagons.
// int RSPIstart[3] = {45, 70, 71}; // We find a 555-spiral by getting the GSPI starting with 45,70,71.
int testNjumps = 3;

// 555 -- too long (almost no visible missing faces)
// int testRSPI[12] = {1, 2, 3, 142, 143, 152, 153, 162, 163, 170, 177, 183};
// int testjump[3][2] = {{170,3},{177,3},{183,4}};

// 556 -- shorter, but still too long.
//int testRSPI[12] = {1, 2, 6, 136, 137, 146, 147, 156, 157, 165, 172, 180};
//int testjump[3][2] = {{165,1},{172,2},{180,1}};

// 565 -- better.
int testRSPI[12] = {1,3,4,127,128,146,147,155,164,165,173,183};
int testjump[3][2] = {{155,1},{173,1},{183,3}};


struct windup_t {
  Triangulation   dual;
  PlanarGraph    graph;
  //  vector<int>     face_spiral; // Is always 0,1,....,N-1
  vector<face_t>  faces;

  windup_t(const vector<int>& spiral_string, const FullereneGraph::jumplist_t& j)
  {
    dual  = Triangulation(spiral_string,j);
    graph = PlanarGraph(dual.dual_graph());
    faces = dual.dual_faces();
  }
};




int main(int ac, char **av)
{
  int N = ac>=2? strtol(av[1],0,0) : testN;
  int E = 3*N/2;
  int F = E-N+2;
  vector<int> spiral_string(F,6);
  FullereneGraph::jumplist_t  jumps;

  if(ac>=14) 
    for(int i=0;i<12;i++) spiral_string[strtol(av[i+2],0,0)-1] = 5;
  else 
    for(int i=0;i<12;i++) spiral_string[testRSPI[i]-1] = 5;

  if(ac>=16) {
    for(int i=0;i<ac-14;i++)
      jumps.push_back(make_pair(strtol(av[i+15],0,0),strtol(av[i+16],0,0)));
  } else
    for(int i=0;i<testNjumps;i++)
      jumps.push_back(make_pair(testjump[i][0]-1,testjump[i][1]));

  cout << "Constructing C"<<N<<" and dual from GSPI "<< make_pair(jumps,spiral_string) << "\n";
  windup_t w(spiral_string,jumps);

  vector<int> rspi;
  FullereneGraph::jumplist_t j;

  // Run once in order to find the spiral we want to use.
  //w.dual.get_spiral(RSPIstart[0]-1,RSPIstart[1]-1,RSPIstart[2]-1,rspi,j);
  //w.dual.get_spiral(0,1,5,rspi,j); // WRT 555...
  w.dual.get_spiral(0,3,2,rspi,j);   // WRT 555
  cout << "RSPI = " << rspi << ";\n"
       << "jumps = " << j << ";\n";

  printf("Computing planar layout\n");
  w.graph.layout2d = w.graph.tutte_layout();
  w.dual.layout2d  = w.dual.tutte_layout();

  vector<face_t> faces(w.graph.compute_faces_flat(6,true));
  vector<int> face_translate(F);

  for(int i=0;i<F;i++) sort(faces[i].begin(),faces[i].end());

  



  FullereneGraph g(w.graph);
  if(halma_level > 0){
    g = g.halma_fullerene(halma_level ,true);
    printf("Computing planar layout of halma transform\n");
    g.layout2d = g.tutte_layout();
  }

  string basename("polyhedron-"+to_string(g.N));
  ofstream output(("output/"+basename+".m").c_str());
  output
    << "graph = " << w.graph << ";\n"
    << "dual  = " << w.dual << ";\n"
    << "faces = " << w.faces << ";\n";

  if(compute_polyhedron){
    printf("Constructing initial geometry for optimization\n");
    Polyhedron P0(g,g.zero_order_geometry(),6);
    
    {
      string filename("output/"+basename+"-P0.mol2");
      ofstream mol2(filename.c_str());
      printf("Writing mol2-file: %s\n",filename.c_str());
      mol2 << P0.to_mol2();
      mol2.close();
    }
    
    Polyhedron P(P0);
    P.optimise();
    
    output << "coordinates0 = " << P0.points << ";\n";
    output << "coordinates = "  << P.points << ";\n";
    output << "P = " << P << ";\n";
  
    Polyhedron D(P.dual(6,true));
    output << "PD = " << D << ";\n";
    output.close();

    {
      string filename("output/"+basename+".mol2");
      ofstream mol2(filename.c_str());
      printf("Writing mol2-file: %s\n",filename.c_str());
      mol2 << P.to_mol2();
      mol2.close();
    }
    {
      string filename("output/"+basename+"-dual.mol2");
      ofstream mol2(filename.c_str());
      printf("Writing mol2-file: %s\n",filename.c_str());
      mol2 << D.to_mol2();
      mol2.close();
    }

    {
      string filename("output/"+basename+".xyz");
      ofstream pov(filename.c_str());
      printf("Writing PoV-file: %s\n",filename.c_str());
      pov << P.to_xyz();
      pov.close();
    }
    {
      string filename("output/"+basename+"-dual.xyz");
      ofstream pov(filename.c_str());
      printf("Writing PoV-file: %s\n",filename.c_str());
      pov << D.to_xyz();
      pov.close();
    }
  }


  return 0;
}
