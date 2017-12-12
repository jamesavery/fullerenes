#include "libgraph/spiral.hh"
#include "libgraph/planargraph.hh"
#include "libgraph/triangulation.hh"

int main(int ac, char **av)
{
  if(ac<2){
    fprintf(stderr,"Syntax: %s <graph-filename> [index]\n"
	    "Supported graph formats: %s\n\n", av[0], to_string(PlanarGraph::input_formats).c_str());
    return -1;
  }
  
  PlanarGraph G = PlanarGraph::from_file(av[1]);

  // TODO: Hack. Make sure the orientation is already present after reading from_file()
  if(!G.is_oriented){
    G.layout2d = G.tutte_layout();
    G.orient_neighbours();
  }

  full_spiral_name::construction_scheme_t graph_type;
  
  Triangulation T(G.enveloping_triangulation(graph_type));
  // TODO: Add functionality to compute the full_spiral_name.
  //       Add this method both to Triangulation and to PlanarGraph classes
  vector<int> spiral;
  jumplist_t jumps;

  // TODO: Add option to specify GS or CS for canonical spiral
  T.get_spiral(spiral,jumps);

  // TODO: Output general spiral code
  
  return 0;
}
