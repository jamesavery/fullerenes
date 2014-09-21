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

minplus_matrix APSP_weighted(const minplus_matrix &A0)
{
  int count = ceil(log2(A0.m));

  minplus_matrix A(A0);
  for(int i=0;i<count;i++) A = A*A;

  return A;
}

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
  return APSP_weighted(A);
}


node_t findEndOfPath(const Triangulation& G, int u, int v, int a, int b)
{
  int t[3],t_next[3];

  // First triangle
  t[0] = u;
  t[1] = v;
  t[2] = G.nextCCW(dedge_t(u,v));
  
  if(a == 1 && b == 0) return t[1];
  if(a == 0 && b == 1) return t[2];
  
  t_next[0] = t[2];
  t_next[1] = t[1];
  t_next[2] = G.nextCCW(dedge_t(t_next[0],t_next[1]));

  int p[3][2] = {{0,0},{1,0},{0,1}};
}

minplus_matrix findSurfaceDistances(const Triangulation& G)
{
  minplus_matrix H = APSP_unweighted(G);  // Shortest path u--v is upper bound to d(u,v)
  int M = *max_element(H.begin(),H.end()); // M is diameter of G
  
  for(int i=0;i<H.size();i++) H[i] *= H[i]; // Work with square distances, so that all distances are integers.
  

  for(int u=0;u<G.N;u++) 
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

  assert(dg.is_consistently_oriented());

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
