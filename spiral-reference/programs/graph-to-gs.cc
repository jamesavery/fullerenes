#include "libgraph/spiral.hh"
#include "libgraph/planargraph.hh"

int main(int ac, char **av)
{
  if(ac<2){
    fprintf(stderr,"Syntax: %s <graph-filename> [index]\n"
	    "Supported graph formats: %s\n\n", av[0], to_string(PlanarGraph::input_formats).c_str());
    return -1;
  }
  
  PlanarGraph G = PlanarGraph::from_file(av[1]);

  return 0;
}
