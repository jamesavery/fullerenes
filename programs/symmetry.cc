#include "fullerenes/spiral.hh"
#include "fullerenes/symmetry.hh"

#include "symmetry-examples.cc"

int main(int argc, char **argv)
{
  int N;
  jumplist_t jumps;
  vector<int> spiral;

  // Parse command line
  if(argc==3){		// Two arguments: Symmetry name and example number
    string sym = argv[1];
    int number = argc>2? strtol(argv[2],0,0) : 0;
    
    cout << "sym = " << sym << "; number = " << number << endl;

    RSPIExample e(SymmetryExample(sym,number));
    N = e.N;
    spiral = vector<int>(N/2+2,6);
    for(int i=0;i<12;i++) spiral[e.RSPI[i]-1] = 5;
  } else if(argc==14){		// 13 arguments: N and 12 pentagon indices (starting from 1)
    N = strtol(argv[1],0,0);
    spiral = vector<int>(N/2+2,6);
    for(int i=0;i<12;i++) spiral[strtol(argv[i+2],0,0)-1] = 5;    
  }  else {
    fprintf(stderr,"Usage:\n"
	    "\t%s <Symmetry name> <example number>                 (See symmetry-examples.cc)\n"
	    "\t%s <N> <12 pentagon indices separated by spaces>\n",
	    argv[0],argv[0]);
    return -1;
  }
  
  cout << "N = " << N << ";\n"
       << "spiral = " << spiral << ";\n";

  // Compute symmetry-information
  Symmetry g(spiral);

  printf("Symmetry-group order=%d\n",int(g.G.size()));

  vector<int> involutions = g.involutions();
  printf("Number of involutions=%d\n",int(involutions.size()));
  
  cout << "g = " << g.neighbours << ";\n";

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

