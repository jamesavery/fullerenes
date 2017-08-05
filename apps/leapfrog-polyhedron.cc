#include "libgraph/planargraph.hh"
#include "libgraph/polyhedron.hh"
#include "libgraph/triangulation.hh"

typedef Triangulation::jumplist_t jumplist_t;

bool grspi_from_args(int ac, char **av, int &N, jumplist_t &jumps, vector<int> &RSPI)
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


int main(int ac, char **av)
{
  int N;
  jumplist_t jumps;
  vector<int> RSPI(12);

  grspi_from_args(ac,av,N,jumps,RSPI);
  
  FullereneDual dF(N,RSPI,jumps);
  PlanarGraph    F(dF.dual_graph());
  F.layout2d = F.tutte_layout();
  
  cout << "dF = " << dF << ";\n"
       << "F  = " << F  << ";\n";

  assert(dF.is_consistently_oriented());
  assert(F.is_consistently_oriented());
  
  Triangulation G = F.leapfrog_dual();

  assert(G.is_consistently_oriented());  
  cout << "G = " << G << ";\n";

  {
    vector<int> spiral;
    jumplist_t  jumps;
    G.get_spiral(spiral,jumps);

    cout << "LFspiral = " << spiral << ";\n"
	 << "LFjumps  = " << jumps  << ";\n";
  }
  
  return 0;
}
