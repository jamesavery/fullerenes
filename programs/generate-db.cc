#include <stdio.h>
#include <inttypes.h>
#include <libgraph/triangulation.hh>
#include <libgraph/symmetry.hh>

#include "linalg.cc"

typedef enum {PSEUDO,PROPER,META} valence_qualifier_t;

void open_closed(const Graph& g, bool &open, int &qualifier, double &HOMO, double &LUMO)
{
  DenseSqrMatrix A(g.N);
  for(node_t u=0;u<g.N;u++)
    for(int i=0;i<g.neighbours[u].size();i++){
      A(u,g.neighbours[u][i]) = 1;
      A(g.neighbours[u][i],u) = 1;
    }

  vector<double> lambda = A.eigenvalues(g.N/2-1,g.N/2);
  //  vector<double> lambda = A.eigenvalues(0,g.N-1);

  HOMO = lambda[0];
  LUMO = lambda[1];
  open = false;
  
  if(HOMO == LUMO){ open = true; }
  if(!open){
    if(LUMO>0) qualifier = PSEUDO;
    else 
      if (HOMO<0) qualifier = PROPER;
      else qualifier = META;
  } {
    if(LUMO > 0) qualifier = PSEUDO;
    else if (LUMO == 0) qualifier = PROPER;
    else qualifier = META;
  }
}

int main(int ac, char **av)
{
  if(ac<=2) return -1;

  int N = strtol(av[1],0,0);
  string rspipath = string(av[2]);

  FILE *f = fopen(rspipath.c_str(),"rb");
  if(!f){
    perror(rspipath.c_str());
    return -2;
  }
  while(!feof(f)){
    uint8_t rspi[12];
    size_t nread = fread(rspi,1,12,f);
    
    if(nread == 12){
      vector<int> spiral(N/2+2,6);
      for(int i=0;i<12;i++) spiral[rspi[i]] = 5;

      Symmetry S(spiral);
      cout << S.point_group() << ", NMR = " << S.NMR_pattern() << endl;
      //      open_closed(S.dual_graph());
    }
    if(rspi[0] != 0) cout << vector<int>(rspi,rspi+12) << endl;
  }
  return 0;
}
