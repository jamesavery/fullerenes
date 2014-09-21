#include "libgraph/auxiliary.hh"
#include "libgraph/triangulation.hh"

typedef unsigned long int_l;

class minplus_matrix : public vector<int_l> {
public: 
  int m,n;
  static const int_l zero = UINT_MAX;

  minplus_matrix(int m, int n, const vector<int_l>& A) : vector<int_l>(A.begin(),A.end()),m(m), n(n) { assert(A.size() == m*n); }
  minplus_matrix(int m, int n) : vector<int_l>(m*n,zero),m(m),n(n) { }

  minplus_matrix operator+(const minplus_matrix& B)
  {
    const minplus_matrix& A(*this);
    assert(A.m == B.m && A.n == B.n);
    minplus_matrix C(A.m,A.n);

    for(int i=0;i<B.size();i++) C[i] = min(A[i],B[i]);
    return C;
  }

  minplus_matrix operator*(const minplus_matrix& B)
  {
    const minplus_matrix& A(*this);
    assert(A.n == B.m);
    minplus_matrix C(A.m,B.n);

    for(int i=0;i<A.m;i++)
      for(int j=0;j<B.n;j++){
	int_l x = zero;
	for(int k=0;k<A.n;k++) x = min(x, int_l(A[i*A.n+k]+B[k*B.n+j]));
	x = min(x,zero);
	C[i*C.n+j] = x;
      }
    return C;
  }

  friend ostream& operator<<(ostream& S, const minplus_matrix& A)
  {
    vector< vector<int_l> > VV(A.m, vector<int_l>(A.n,zero));
    for(int i=0;i<A.m;i++) 
      for(int j=0;j<A.n;j++)
	if(A[i*A.n+j] != zero) VV[i][j] = A[i*A.n+j];

    S << VV;
    return S;
  }
};

#include <math.h>

minplus_matrix APSP_unweighted(const Graph& g)
{
  minplus_matrix A(g.N,g.N);

  for(node_t u=0;u<g.N;u++){
    A[u*(A.n+1)] = 0;
    for(int i=0;i<g.neighbours[u].size();i++){
      node_t v = g.neighbours[u][i];
      A[u*A.n+v] = 1;
    }
  }
  
  int count = ceil(log2(g.N));
  

  for(int i=0;i<count;i++){
    //    cout << "step " << i << ": " << A << endl;
    A = A*A;
  }

  return A;
}

int main(int ac, char **av)
{
  assert(ac>13);
  int N = strtol(av[1],0,0);
  vector<int> spiral(N/2+2,6);
  for(int i=0;i<12;i++) spiral[strtol(av[2+i],0,0)-1] = 5;

  cout << "spiral = " << spiral << ";\n";

  Triangulation   dg(spiral);
  PlanarGraph      g(dg.dual_graph());

  if(0){
    cout << "\n\nMatrix-APSP:\n";
    minplus_matrix Da(dg.N,dg.N);
    for(int i=0;i<500;i++)
      Da = APSP_unweighted(dg);
    cout << "Da = " << Da << endl;
} else {
    cout << "Dijkstra-APSP:\n";
    vector<int> Dt(dg.N*dg.N);
    for(int i=0;i<500;i++)
      Dt = dg.all_pairs_shortest_paths(dg.N);
  
    minplus_matrix Db(dg.N,dg.N,vector<int_l>(Dt.begin(),Dt.end()));
    cout << "Db = " << Db << endl;
}
    //cout << ((Da ==Db)? "OK" : "Error") << endl;
  

  return 0;
}
