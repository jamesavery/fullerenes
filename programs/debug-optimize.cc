#include "fullerenes/spiral.hh"
#include "fullerenes/fullerenegraph.hh"
#include "fullerenes/polyhedron.hh"
#include "fullerenes/triangulation.hh"
#include "fullerenes/isomerdb.hh"

#include <sys/stat.h>
#include <sys/types.h>

size_t file_size(const char *filename) {
  struct stat st;

  if (stat(filename, &st) == 0)
    return st.st_size;

  return -1;
}

int main(int ac, char **av)
{
  string cubic_neighbours_path = av[1];

  size_t file_length = file_size(cubic_neighbours_path.c_str());
  size_t N = file_length/(sizeof(node_t)*3);	// 4 bytes per node_t, 3 neighbours per node.
  
  vector<node_t> neighbours_flat(3*N);
  FILE *file = fopen(cubic_neighbours_path.c_str(), "rb");
  size_t n_read = fread(&neighbours_flat[0], sizeof(node_t), 3*N, file);
  assert(n_read == 3*N);
  fclose(file);

  PlanarGraph g(N);
  for(node_t u=0;u<N;u++) g.neighbours[u] = {neighbours_flat[3*u],neighbours_flat[3*u+1],neighbours_flat[3*u+2]};
  
  g.layout2d = g.tutte_layout();
  cout << g.to_latex(10, 10) << "\n";
    
  Polyhedron P0(g,g.zero_order_geometry(),6);
  string basename("polyhedron-"+to_string(N));
  Polyhedron::to_file(P0,"output/"+basename+"-P0.mol2");

  return 0;
  
  // printf("P0\n");
  // Polyhedron P(P0);
  // printf("Optimizing P\n");  
  // P.optimize();
  // printf("Writing P\n");
  // Polyhedron::to_file(P,"output/"+basename+".mol2");

  // printf("Aligning P\n");  
  // P.move_to_origin();
  // P.align_with_axes();

  // printf("Aligning P-aligned\n");    
  // Polyhedron::to_file(P,"output/"+basename+"-if.mol2");
  // Polyhedron::to_file(P,"output/"+basename+"-if.xyz");
  // //  Polyhedron::to_file(P,"output/"+basename+"-if.pov");

  // ofstream output(("output/"+basename+".m").c_str());

  // printf("Writing mathematica version of P\n");  
  // vector<face_t> facemap(g.compute_faces());
  // output << "g = " << g << ";\n";
  // output << "coordinates0 = " << P0.points << ";\n";
  // output << "coordinates = "  << P.points << ";\n";
  // output << "pentagons = " << facemap[5] << ";\n"
  // 	  << "hexagons  = " << facemap[6] << ";\n"
  // 	  << "RSPI = " << RSPI << ";\n";

  // output << "P0 = " << P0 << ";\n";
  // output << "P = " << P << ";\n";

  // Polyhedron D(P.dual());
  // D.layout2d = D.tutte_layout();
  // D.faces    = D.compute_faces(3,true);
  // D.face_max = 3;
  // D.optimize();
  // output << "PD = " << D << ";\n";
  
  // output.close();

  // Polyhedron::to_file(P,"output/"+basename+"-dual.mol2");
  // //  Polyhedron::to_file(P,"output/"+basename+"-dual.pov");

  //  return 0;
}
