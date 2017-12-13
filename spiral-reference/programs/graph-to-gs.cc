#include "libgraph/spiral.hh"
#include "libgraph/planargraph.hh"
#include "libgraph/triangulation.hh"

int main(int ac, char **av)
{
  if(ac<2){
    fprintf(stderr,"Syntax: %s <graph-filename> [index=0] [naming_scheme=cage]\n"
	    "Graph-file formats: %s.\n"
    	    "Index (optional) is index into multi-graph file.\n"
	    "Naming-scheme is one of: {cage, fulleroid, fullerene}.\n\n",
	    av[0], to_string(PlanarGraph::input_formats).c_str());
    return -1;
  }
  
  PlanarGraph G = PlanarGraph::from_file(av[1]);

  // TODO: Hack. Make sure the orientation is already present after reading from_file()
  if(!G.is_oriented){
    G.layout2d = G.tutte_layout();
    G.orient_neighbours();
  }

  full_spiral_name spiral_name(G);

  cout << spiral_name << endl;
  
  return 0;
}
