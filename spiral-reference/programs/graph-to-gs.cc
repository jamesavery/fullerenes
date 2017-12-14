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
  string format   = ac>2? av[2] : "mol2";
  int index       = ac>3? strtol(av[3],0,0) : 0;
  string naming_scheme = ac>4? av[4] : "cage";

  PlanarGraph G = PlanarGraph::from_file(filename);

  // TODO: Hack. Make sure the orientation is already present after reading from_file()
  if(!G.is_oriented){
    G.layout2d = G.tutte_layout();
    G.orient_neighbours();
  }

  spiral_nomenclature spiral_name(G);

  if(format == "fullerene") spiral_name.naming_scheme = spiral_nomenclature::FULLERENE;
  if(format == "fulleroid") spiral_name.naming_scheme = spiral_nomenclature::FULLEROID;  

  cout << spiral_name << endl;
  
  return 0;
}
