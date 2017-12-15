#include "libgraph/spiral.hh"
#include "libgraph/planargraph.hh"
#include "libgraph/triangulation.hh"

int main(int ac, char **av)
{
  if(ac<2){
    fprintf(stderr,"Syntax: %s <graph-filename> [naming_scheme=cage]\n"
	    "Graph-file formats: %s.\n"
    	    "Index (optional) is index into multi-graph file.\n"
	    "Naming-scheme is one of: {cage, fulleroid, fullerene}.\n\n",
	    av[0], to_string(PlanarGraph::input_formats).c_str());
    return -1;
  }

  string filename = av[1];
  string naming_scheme = ac>2? av[2] : "cage";

  FILE *file = fopen(filename.c_str(),"rb");
  size_t graph_count=0, graph_size=0;
  PlanarGraph::read_hog_metadata(file,graph_count,graph_size);

  cerr << "Reading " << graph_count << " planarcode graphs from " << filename << ".\n";
  fseek(file,15,SEEK_SET);  
  for(size_t i=0;i<graph_count;i++){
    PlanarGraph G = PlanarGraph::read_hog_planarcode(file);
    assert(G.is_oriented);
    
    spiral_nomenclature spiral_name(G);
    
    if(naming_scheme == "fullerene") spiral_name.naming_scheme = spiral_nomenclature::FULLERENE;
    if(naming_scheme == "fulleroid") spiral_name.naming_scheme = spiral_nomenclature::FULLEROID;  

    string spiral_code = spiral_name.to_string();
    cout << spiral_code << endl;
    cerr << "Graph " << i << " is of length " << spiral_code.size() << ".\n";
  }
  
  return 0;
}
