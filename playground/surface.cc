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

  int_l& operator()(int i, int j)      { return (*this)[i*n+j]; }
  int_l operator()(int i, int j) const { return (*this)[i*n+j]; }

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


vector<int> draw_path(int major, int minor)
{
  int slope = major/minor, slope_remainder = major%minor, slope_accumulator = 0;

  vector<int> runs(minor);

  for(int i=0; i<minor; i++){
    runs[i] = slope;
    slope_accumulator += slope_remainder;

    runs[i] += (slope_accumulator >= minor);
    slope_accumulator %= minor;
  }

  return runs;
}

// Given start node u0 and adjacent face F_i, lay down triangles along the the straight
// line to Eisenstein number (a,b), and report what the final node is.
//
// Assumes a,b >= 1.
//
node_t end_of_the_line(const Triangulation& G, node_t u0, int i, int a, int b)
{
  node_t q,r,s,t;		// Current square

  auto go_north = [&](){ 
    const node_t S(s), T(t); // From old square
    q = S; r = T; s = G.nextCCW(dedge_t(S,T)); t = G.nextCCW(dedge_t(s,r)); 
  };
  
  auto go_east = [&](){
    const node_t R(r), T(t); // From old square
    q = R; s = T; r = G.nextCCW(dedge_t(s,q)); t = G.nextCCW(dedge_t(s,r)); 
  };

  // Square one
  q = u0; 			// (0,0)
  r = G.neighbours[u0][i];	// (1,0)
  s = G.nextCCW(dedge_t(q,r));	// (0,1)
  t = G.nextCCW(dedge_t(s,r));	// (1,1)

  vector<int> runlengths = draw_path(max(a,b), min(a,b));

  for(int i=0;i<runlengths.size();i++){
    int L = runlengths[i];

    if(a>=b){			// a is major axis
      for(int j=0;j<L-1;j++)    go_north();
      if(i+1<runlengths.size()) go_east();
    } else {			// b is major axis
      for(int j=0;j<L-1;j++)    go_east();
      if(i+1<runlengths.size()) go_north();
    }
  }
  
  return t;			// End node is upper right corner.
}


minplus_matrix semisimple_distances(const minplus_matrix& Hinit, const Triangulation& G)
{
  minplus_matrix H(Hinit);  
  int M = *max_element(H.begin(),H.end());      // M is upper bound to path length
  
  for(int i=0;i<H.size();i++) H[i] *= H[i]; // Work with square distances, so that all distances are integers.

  for(node_t u=0;u<G.N;u++)
    for(int i=0;i<G.neighbours[u].size();i++){

      // Note: All Eisenstein numbers of the form (a,0) or (0,b) yield same lengths
      //       as graph distance, and are hence covered by initial step. So start from 1.
      //       M is upper bound for distance, so only need to do a^2+ab+b^2 strictly less than M.
      for(int a=1; a<M;a++)
	for(int b=1; a*a + a*b + b*b < M*M; b++){
	  // Check: if(gcd(a,b) != 1) continue.
	  const node_t v = end_of_the_line(G,u,i,a,b);
	  H(u,v) = min(H(u,v), int_l(a*a + a*b + b*b));
	}
    }
  return H;
}

minplus_matrix surface_distances(const Triangulation& G)
{
  minplus_matrix Hinit(APSP_unweighted(G)); // Shortest path u--v is upper bound to d(u,v)
  
  minplus_matrix Hsimple = semisimple_distances(Hinit,G);
  minplus_matrix H = APSP_weighted(Hsimple);

  return H;
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
