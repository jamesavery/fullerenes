#include "libgraph/planargraph.hh"
#include "libgraph/polyhedron.hh"
#include "libgraph/triangulation.hh"

typedef Triangulation::jumplist_t jumplist_t;

bool grspi_from_args(int ac, char **av, int &N, vector<int> &RSPI, jumplist_t &jumps)
{
  if(ac < 14) return false;
  
  N = strtol(av[1],0,0);
  for(int i=0;i<12;i++) RSPI[i] = strtol(av[i+2],0,0)-1;

  if(ac>14){			// General RSPI: RSPI followed by jumps.
    for(int i=14;i<ac;i+=2)
      jumps.push_back(make_pair(strtol(av[i],0,0)-1,strtol(av[i+1],0,0)));
  }

  return true;
}

Graph planargraph_from_args(int ac, char **av)
{
  int N;
  jumplist_t jumps;
  vector<int> RSPI(12);

  if(ac >= 14){
    grspi_from_args(ac,av,N,RSPI,jumps);
    return FullereneDual(N,RSPI,jumps).dual_graph();
  } else if (ac >= 3) { // Read from HoG planarcode file
    FILE *file = fopen(av[1],"rb");
    int index  = strtol(av[2],0,0);
    assert(file != 0);

    vector<Graph> graphs = PlanarGraph::read_hog_planarcodes(file);
    return graphs[index];
  } else {
    fprintf(stderr,"Syntax: %s [<HoG file> <index> | <N> <general spiral>]\n",av[0]);
    abort();
  }
}

int main(int ac, char **av)
{

  PlanarGraph  g(planargraph_from_args(ac,av));
  g.layout2d = g.tutte_layout();
  
  cout << "g  = " << g  << ";\n";

  assert(g.is_consistently_oriented());
  
  Triangulation dG = g.leapfrog_dual();

  assert(dG.is_consistently_oriented());  
  cout << "dG = " << dG << ";\n";

  PlanarGraph G = dG.dual_graph();
  assert(G.is_consistently_oriented());  
  cout << "G = " << G << ";\n";
   
  {
    vector<int> spiral;
    jumplist_t  jumps;
    dG.get_spiral(spiral,jumps);

    cout << "LFspiral = " << spiral << ";\n"
	 << "LFjumps  = " << jumps  << ";\n";
  }

  
  return 0;
}
