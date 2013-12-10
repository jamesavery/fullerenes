#include "libgraph/fullerenegraph.hh"
#include "libgraph/triangulation.hh"
#include "libgraph/polyhedron.hh"

bool compute_polyhedron = true;
int halma_level = 0;

int testN = 380;
int testjump[2] = {110,2};
int testRSPI[12] = {45, 70, 71, 82, 83, 110, 119, 120, 144, 184, 185, 192};

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
    jumps.push_back(make_pair(testjump[0]-1,testjump[1]));

  cout << "Constructing C"<<N<<" and dual from GSPI "<< make_pair(jumps,spiral_string) << "\n";
  windup_t w(spiral_string,jumps);

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
    P.optimize();
    
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
