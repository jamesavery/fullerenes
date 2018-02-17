#include "libgraph/spiral.hh"
#include "libgraph/planargraph.hh"
#include "libgraph/triangulation.hh"

int main(int ac, char **av)
{
  if(ac<2){
    fprintf(stderr,"Syntax: %s <graph-filename> [format=mol2] [index=0] [naming_scheme=cage]\n"
	    "Graph-file formats: %s.\n"
    	    "Index (optional) is index into multi-graph file.\n"
	    "Naming-scheme is one of: {cage, fulleroid, fullerene}.\n\n",
	    av[0], to_string(PlanarGraph::input_formats).c_str());
    return -1;
  }

  string filename = av[1];
  string format   = ac>2? av[2] : filename_extension(filename);
  int index       = ac>3? strtol(av[3],0,0) : 0;
  string naming_scheme = ac>4? av[4] : "cage";

  FILE *file = fopen(filename.c_str(),"rb");
  PlanarGraph G = PlanarGraph::from_file(file,format,index);
  fclose(file);
  
  // TODO: Hack. Make sure the orientation is already present after reading from_file()
  assert(G.is_oriented);

  spiral_nomenclature spiral_name(G);
  
  if(naming_scheme == "fullerene") spiral_name.naming_scheme = spiral_nomenclature::FULLERENE;
  if(naming_scheme == "fulleroid") spiral_name.naming_scheme = spiral_nomenclature::FULLEROID;  

  cout << spiral_name.to_string() << endl;
  
  return 0;
}
