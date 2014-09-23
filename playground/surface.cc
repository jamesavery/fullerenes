#include "libgraph/auxiliary.hh"
#include "libgraph/triangulation.hh"

typedef unsigned long int_l;

template <typename T> class matrix : public vector<T> {
public:
  int m,n; 

  matrix(int m, int n, const vector<T>& A) : vector<T>(A.begin(), A.end()), m(m), n(n) { assert(A.size()==m*n); }
  matrix(int m, int n, const T& zero = 0) : vector<T>(m*n,zero), m(m), n(n) {}

  // Convert from compatible type
  template <typename S> matrix(const matrix<S>& A) : vector<T>(A.begin(),A.end()), m(A.m), n(A.n) {}

  T& operator()       (int i, int j){ return (*this)[i*n+j]; }
  T  operator()(int i, int j) const { return (*this)[i*n+j]; }

  static matrix minplus_multiply(const matrix& A, const matrix& B, const T& infty_value = USHRT_MAX)
  {
    assert(A.n == B.m);
    matrix C(A.m,B.n);

    for(int i=0;i<A.m;i++)
      for(int j=0;j<B.n;j++){
	T x = infty_value;
	for(int k=0;k<A.n;k++) x = min(x, T(A[i*A.n+k]+B[k*B.n+j]));
	x = min(x,infty_value);
	C[i*C.n+j] = x;
      }
    return C;    
  }

  matrix APSP() const {
    int count = ceil(log2(m));
    matrix A(*this);
    for(int i=0;i<count;i++) A = minplus_multiply(A,A);
    
    return A;
  }

  matrix operator*(const matrix& B)
  {
    assert(n == B.m);
    matrix C(m,B.n);

    for(int i=0;i<m;i++)
      for(int j=0;j<B.n;j++){
	T x = 0;
	for(int k=0;k<n;k++) x += (*this)[i*n+k]*B[k*B.n+j];
	C[i*C.n+j] = x;
      }
    return C;    
  }

  friend ostream& operator<<(ostream& S, const matrix& A)
  {
    vector< vector<T> > VV(A.m, vector<T>(A.n));
    for(int i=0;i<A.m;i++) 
      for(int j=0;j<A.n;j++)
	VV[i][j] = A[i*A.n+j];

    S << VV;
    return S;
  }
};



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
      for(int j=0;j<L-1;j++)    go_east();
      if(i+1<runlengths.size()) go_north();
    } else {			// b is major axis
      for(int j=0;j<L-1;j++)    go_north();

      if(i+1<runlengths.size()) go_east();
    }
  }
  
  return t;			// End node is upper right corner.
}


matrix<int> semisimple_distances(const matrix<int>& Hinit, const Triangulation& G)
{
  matrix<int> H(Hinit);  
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
	  H(u,v) = min(H(u,v), a*a + a*b + b*b);
	}
    }
  return H;
}

matrix<double> surface_distances(const Triangulation& G)
{
  matrix<int> Hinit(G.N,G.N,G.all_pairs_shortest_paths()); // Shortest path u--v is upper bound to d(u,v)
  
  matrix<int> Hsimple_sqr = semisimple_distances(Hinit,G);
  matrix<float> Hsimple(Hsimple_sqr);
  for(int i=0;i<Hsimple.size();i++) Hsimple[i] = sqrt(Hsimple[i]);

  matrix<double> H = Hsimple.APSP();

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

  matrix<int> Hinit(dg.N,dg.N,dg.all_pairs_shortest_paths());
  cout << "Hinit = " << Hinit << "];\n\n";

  matrix<int> Hsimple(semisimple_distances(Hinit,dg));
  
  cout << "Hsimple = Sqrt[" << Hsimple << "];\n\n";

  matrix<double> H = surface_distances(dg);
  
  cout << "H = " << H << ";\n\n";

  return 0;
}
