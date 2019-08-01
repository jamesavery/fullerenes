#include "fullerenes/symmetry.hh"

#include "../apps/symmetry-examples.cc"

int main()
{
  int Nexamples = sizeof(symmetry_examples)/sizeof(RSPIExample);  

  for(int i=0;i<Nexamples;i++){
    const RSPIExample &e(symmetry_examples[i]);
    vector<int> spiral(e.N/2+2,6);
    for(int i=0;i<12;i++) spiral[e.RSPI[i]-1] = 5;

    Symmetry S(spiral);
    // Triangulation T(spiral), U(spiral), V(spiral), W(spiral), Q(spiral);
    // assert(T.N == U.N && U.N == V.N && V.N == W.N);
    if(!(PointGroup(e.sym) == S.point_group())){
      cout << "Example " << i << ", " << e.RSPI << ": " << PointGroup(e.sym) << " != " << S.point_group() << endl;
      cout << "order = " << S.G.size() << ";\n";
      vector<int> 
	mF = S.site_symmetry_counts(S.G),
	mV = S.site_symmetry_counts(S.Gtri),
	mE = S.site_symmetry_counts(S.Gedge);
  
      cout << "mF = " << mF << ";\n"
	   << "mV = " << mV << ";\n"
	   << "mE = " << mE << ";\n";

      abort();
    }
    //    fprintf(stderr,"Example %d OK, detected point group %3s = %3s\n",i,S.point_group().to_string().c_str(),e.sym.c_str());
  }
  return 0;
}
