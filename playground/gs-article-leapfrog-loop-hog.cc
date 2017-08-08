#include "libgraph/planargraph.hh"
#include "libgraph/polyhedron.hh"
#include "libgraph/triangulation.hh"

typedef Triangulation::jumplist_t jumplist_t;

int main(int ac, char **av)
{
  if(ac <= 2){
    fprintf(stderr,"Syntax: %s <hog-db> <dual?> [index-start] [index-end]",av[0]);
    return -1;
  }

  FILE *file = fopen(av[1],"rb");
  if(!file){
    perror((string("Error opening ")+av[1]).c_str());
    return -2;
  }
  
  vector<Graph> read_graphs(PlanarGraph::read_hog_planarcodes(file));

  bool read_dual = strtol(av[2],0,0);
  int  index_start = ac > 3? strtol(av[3],0,0) : 0;
  int  index_end   = ac > 4? min(strtol(av[4],0,0),long(read_graphs.size())) : read_graphs.size();
  
  vector<PlanarGraph> graphs(read_graphs.size());
  
  for(int i=index_start;i<index_end;i++)
    if(read_dual) {
      graphs[i] = PlanarGraph(read_graphs[i]).dual_graph();
      graphs[i].layout2d = graphs[i].tutte_layout();
      graphs[i].orient_neighbours();
    } else 
      graphs[i] = read_graphs[i];


  vector<Triangulation> LFduals(graphs.size());

  auto LF_start = clock();
  for(int i=0;i<graphs.size();i++)
    LFduals[i] = graphs[i].leapfrog_dual();
  auto LF_end = clock();

  vector<vector<int> > LFspirals(graphs.size());
  vector<jumplist_t>   LFjumps  (graphs.size());

  auto spiral_start = clock();  
  for(int i=0;i<graphs.size();i++){
    LFduals[i].get_spiral(LFspirals[i],LFjumps[i]);
  }
  auto spiral_end = clock();  

  double LFtime = (LF_end-LF_start) * 1.0 / CLOCKS_PER_SEC;
  double LFrate = LFtime/(graphs[0].N*graphs.size());

  double spiraltime = (spiral_end-spiral_start) * 1.0 / CLOCKS_PER_SEC;
  double spiralrate = spiraltime/(graphs[0].N*graphs.size());
  
  cout << "LFspirals  = " << LFspirals << ";\n"
       << "LFjumps    = " << LFjumps  << ";\n"
       << "Ngraphs    = " << graphs.size() << ";\n"
       << "Nv         = " << graphs[0].N << ";\n"
       << "LFtime     = " << LFtime << ";\n"
       << "LFrate     = " << LFrate << "';\n"
       << "spiraltime = " << spiraltime << ";\n"
       << "spiralrate = " << spiralrate << ";\n";


  
  return 0;
}
