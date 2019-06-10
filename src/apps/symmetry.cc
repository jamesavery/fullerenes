#include "libgraph/symmetry.hh"

// struct GSpiral {
//   PlanarGraph::jumplist_t   jumps;
//   vector<int>  spiral;
//   GSpiral(const PlanarGraph::jumplist_t& jumps = PlanarGraph::jumplist_t(), 
// 	  const vector<int>& spiral = vector<int>()) : jumps(jumps), spiral(spiral) {}
// };

#include "symmetry-examples.cc"

int main(int ac, char **av)
{
  int N;
  PlanarGraph::jumplist_t jumps;
  vector<int> spiral;

  if(ac>13){
    N = strtol(av[1],0,0);
    spiral = vector<int>(N/2+2,6);
    for(int i=0;i<12;i++) spiral[strtol(av[i+2],0,0)-1] = 5;    
  } else if(ac>1){
    string sym = av[1];
    int number = ac>2? strtol(av[2],0,0) : 0;
    
    cout << "sym = " << sym << "; number = " << number << endl;

    RSPIExample e(SymmetryExample(sym,number));
    N = e.N;
    spiral = vector<int>(N/2+2,6);
    for(int i=0;i<12;i++) spiral[e.RSPI[i]-1] = 5;
  } else {
    N = 20;
    spiral = vector<int>(12,5);
  }

  cout << "n = " << N << ";\n"
       << "spiral = " << spiral << ";\n";

  Symmetry g(spiral);

  printf("Gorder=%d\n",int(g.G.size()));
  vector<int> involutions = g.involutions();
  printf("Ninvolutions=%d\n",int(involutions.size()));
  
  cout << "g = " << g << ";\n";
  // cout << "permutations = " << g.G << ";\n";
  // cout << "permutations = " << g.Gtri << ";\n";
  // cout << "permutations = " << g.Gedge << ";\n";


  vector<int> 
    mF = g.site_symmetry_counts(g.G),
    mV = g.site_symmetry_counts(g.Gtri),
    mE = g.site_symmetry_counts(g.Gedge);
  
  cout << "mF = " << mF << ";\n"
       << "mV = " << mV << ";\n"
       << "mE = " << mE << ";\n";

  vector<bool> reverses(g.G.size());
  for(int i=0;i<g.G.size();i++) reverses[i] = g.reverses_orientation(g.G[i]);
  cout << "reverses = " << reverses << ";\n";
  cout
    << "|Fixed_F| = " << g.group_fixpoints(g.G).size() << "\n"
    << "|Fixed_V| = " << g.group_fixpoints(g.Gtri).size() << "\n"
    << "|Fixed_E| = " << g.group_fixpoints(g.Gedge).size() << "\n";


  cout << "Calculated point group: " << g.point_group() << endl;
  //  cout << "coxmatrix = " << g.coxeter_matrix() << ";\n";
  //  cout << "multtable = " << g.multiplication_table() << ";\n";


  
  return 0;
}

