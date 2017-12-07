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

  full_spiral_name fn(spiral_code);
  PlanarGraph G(fn);
  
  PlanarGraph::to_file(G,stdout,output_format);
  
  // for(auto example: spiral_paper_examples){
  //   string name = example.first, gs = example.second;
  //   cout << name << "_name = \"" << gs << "\";\n";
  //   cout << name << " = " << full_spiral_name(gs) << ";\n\n";
  //   cout << name << " = " << PlanarGraph(full_spiral_name(gs)) << ";\n\n";
  // }

  return 0;
}
