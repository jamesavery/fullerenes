#include <string>

#include "libgraph/triangulation.hh"
#include "libgraph/spiral.hh"



int main(int ac, char **av)
{
  if(ac<2){
    fprintf(stderr,"Syntax: %s \"spiral_code\" [graph-format]\n",av[0]);
    fprintf(stderr,"Graph output formats: %s\n",to_string(PlanarGraph::output_formats).c_str());
    return -1;
  }
  
  string spiral_code   = av[1];
  string output_format = ac>=3? av[2] : "ascii";

  spiral_nomenclature fn(spiral_code);
  PlanarGraph G(fn);
  
  PlanarGraph::to_file(G,stdout,output_format);
  printf("\n");  

  return 0;
}
