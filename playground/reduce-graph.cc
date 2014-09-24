#include <math.h>
#include <fstream>

#include "libgraph/fullerenegraph.hh"
#include "libgraph/triangulation.hh"
#include "libgraph/polyhedron.hh"

void insert_before(vector<int> &v, const int before_what, const int value){
  vector<int>::iterator pos = std::find(v.begin(), v.end(), before_what);
  assert(pos != v.end());
  v.insert(pos, value);
}

Triangulation Delaunayify(Triangulation T, double distances[12][12]){

//         C                     C       
//     c / | \ b               /   \      .
//     /   |   \             /       \    .
//   D     |     B   --->  D-----------B 
//     \   |   /             \       /   
//     d \ | / a               \   /     
//         A                     A       

  assert(T.is_consistently_oriented());

  vector< vector<double> > Pdist;
  for(int i=0;i<12;i++) Pdist.push_back(vector<double>(distances[i],distances[i]+12));
  cout << "Pdist = " << Pdist << ";\n\n";

  int A, B, C, D;
  unsigned int delaunay_flips = 1;
  vector< vector<node_t> > moves;
  vector<Triangulation>    steps(1,T);
  

  auto flip = [&](){
    delaunay_flips++;

    T.remove_edge(edge_t(A,C));

    insert_before(T.neighbours[B],A,D);
    insert_before(T.neighbours[D],C,B);

    moves.push_back({A+1,B+1,C+1,D+1});    
    steps.push_back(T);
  };


  ofstream debug("output/delaunayify.m");
  int step = 0;

  while(delaunay_flips != 0){
    delaunay_flips=0;

    double total_angles = 0;
    for(node_t u=0; u<T.N; ++u){
      for(int j=0; j<T.neighbours[u].size(); ++j){

        A = u;
        B = T.neighbours[u][j];
        C = T.neighbours[u][(j+1)%T.neighbours[u].size()];
        D = T.neighbours[u][(j+2)%T.neighbours[u].size()];
        const double a = distances[A][B],
                     b = distances[B][C],
                     c = distances[C][D],
                     d = distances[D][A],
                     AC = distances[A][C];
        const double 
	  beta  = acos((a*a + b*b - AC*AC) / (2*a*b)),
	  delta = acos((c*c + d*d - AC*AC) / (2*c*d));
	//	printf("beta + delta = %g+%g = %g: ",beta,delta,beta+delta);

	total_angles += beta+delta;
        if(beta + delta > M_PI){
	  //	  printf("flip!\n");
	  flip();
	} else { 
	  //	  printf("don't flip!\n");
	  //	  printf("a^2 + b^2 - AC^2 = %g^2 + %g^2 - %g^2 = %g\n",a,b,AC,a*a+b*b-AC*AC);
	}

      }
    }
    printf("Number of flips: %d; Total_angles: %g\n",delaunay_flips, total_angles);
    if(++step > 20) break;
  }

  debug << "moves = " << moves << ";\n"
	<< "steps = " << steps << ";\n";
  debug.close();

  return T;
}


Triangulation reduce_triangulation(Triangulation T){

  for( int k=T.N; k>12; --k){
  //for( ; T.N>12; --T.N){

    // cout << "n1" << T.neighbours << endl;
    // cout << "n1" << T.neighbours << endl;
    // cout << "N: " << T.N << endl;    
    // cout << "T1=" << T << endl;
    // cout << "deg N:" << T.neighbours[T.N-1].size() << endl;

    // find a deg-6-node (if the vertices are sorted, then just take the last)
  
    // save its neighbours
    vector<int> hole(T.neighbours.back());
 
    // remove node and neighbor entries
    for(vector<int>::iterator it=T.neighbours.back().begin(), to=T.neighbours.back().end(); it!=to; ++it){
      vector<int>::iterator to_erase = std::find(T.neighbours[*it].begin(), T.neighbours[*it].end(), T.N-1);

      assert(to_erase != T.neighbours[*it].end());
      T.neighbours[*it].erase(to_erase);
    }
    T.neighbours.pop_back();
    --T.N;
  
    // patch the hole:
    // phase one: find nodes with degree-2 because they require a new connection
    for(int i=0; i<hole.size(); ++i){
      if(T.neighbours[hole[i]].size() == 2){ // if hole[k] is deg-2,
        // cout << hole << "" << hole.size()<< endl;
        // cout << "deg2 found, connecting " << hole[i] << " and " << hole[(i+2) % hole.size()] << endl;

        // then connect hole[k] with hole[k+2i%hole.size] and remove hole[k+1] from hole
        insert_before(T.neighbours[hole[i]],                 hole[(i-1+hole.size())%hole.size()], hole[(i+2)%hole.size()]);
        insert_before(T.neighbours[hole[(i+2)%hole.size()]], hole[(i+1)%hole.size()],             hole[i]);
        hole.erase(hole.begin() + ((i+1)%hole.size()));
      }
    }
  
    // phase two: triangulate the remaining hole in a fan-like manner:
    while(hole.size() > 3){

      int shift = 0;
      // check if 0 and 2 are connected already
      if(T.neighbours[hole[0]].end() != std::find(T.neighbours[hole[0]].begin(), T.neighbours[hole[0]].end(), hole[2])){
        shift = 1;
        // printf("%i and %i were connected already, connecting %i and %i instead.\n", hole[0], hole[2], hole[1], hole[3%hole.size()]);
      }

      // connect hole[k] with hole[k+2i%hole.size] and remove hole[k+1] from hole
      // (this is not a very good retriangulation but valid and delaunay on 12
      // points is fast anyway)
      insert_before(T.neighbours[hole[0+shift]],                 hole[(0+shift-1+hole.size())%hole.size()], hole[(0+shift+2)%hole.size()]);
      insert_before(T.neighbours[hole[(0+shift+2)%hole.size()]], hole[(0+shift+1)%hole.size()],             hole[0+shift]);
      hole.erase(hole.begin() + (0+shift+1)%hole.size());
    }
//    cout << "n3" << T.neighbours << endl;
//    cout << "T3=" << T << endl;
//  if(! T.is_consistently_oriented()) cout << "not consistently oriented" << endl;
//    cout << "------------------------" << endl;
  }

  return T;
}



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

  vector<int> paths(minor+1,0), runs(minor);

  for(int i=0; i<minor; i++){
    slope_accumulator += slope_remainder;

    paths[i+1] = paths[i] + slope + (slope_accumulator != 0);

    if((i+1<minor) && (slope_accumulator >= minor || slope_remainder == 0)){
      paths[i+1]++;
      slope_accumulator %= minor;
    }

    runs[i]    = paths[i+1]-paths[i];
  }

  //  cout << make_pair(major,minor) << " runlengths is " << runs << endl;

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

  if(a==1 && b==1) return t;

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

	  // printf("min(H(%d,%d),|(%d,%d)|^2)  = min(%d,%d)\n",
	  // 	 u,v,a,b,H(u,v), a*a+a*b+b*b);
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

matrix<int> fullerene_surface_distances_sqr(const Triangulation& G)
{
  matrix<int> Hinit(G.N,G.N,G.all_pairs_shortest_paths()); // Shortest path u--v is upper bound to d(u,v)
  return semisimple_distances(Hinit,G);
}


Polyhedron fullerene_dual_polyhedron(const Triangulation& dg)
{
  FullereneGraph g(dg.dual_graph());
  g.layout2d = g.tutte_layout();

  vector<coord3d> points = g.zero_order_geometry();
  points = g.optimized_geometry(points);

  vector<coord3d> dual_points(dg.N);

  vector<face_t> faces(dg.N);
  for(int i=0;i<dg.triangles.size();i++)
    for(int j=0;j<3;j++)
      faces[dg.triangles[i][j]].push_back(i);

  for(int i=0;i<faces.size();i++)
    dual_points[i] = faces[i].centroid(points);

  return Polyhedron(dg, dual_points);
}



int main(int ac, char **av) {
  int N;
  vector<int> RSPI(12);
  N = strtol(av[1], 0, 0);
  for (int i = 0; i < 12; i++)
    RSPI[i] = strtol(av[i + 2], 0, 0) - 1;

  // int isomer_number = ac>=3? strtol(av[2],0,0) : 1;
  // bool IPR = ac>=4? strtol(av[3],0,0) : false;

  // IsomerDB DB(N,IPR);

  // FullereneGraph
  // g(IsomerDB::makeIsomer(N,DB.getIsomer(N,isomer_number,IPR)));

  string filename = "output/reduce-graph-C"+to_string<int>(N)+".m";
  ofstream output(filename);

  vector<int> spiral(N / 2 + 2, 6);
  for (int i = 0; i < 12; i++)
    spiral[RSPI[i]] = 5;

  Triangulation T1(spiral);
  Triangulation T(T1.sort_nodes());

  Polyhedron PT = fullerene_dual_polyhedron(T);

  cout << "number vertices in T: " << T.N << endl;
  cout << "neighbours in T: " << T.neighbours << endl;

  output << "T = " << T << ";\n"
	 << "PT =" << PT<< ";\n";
  Triangulation rT = reduce_triangulation(T);

  cout << "number vertices in rT: " << rT.N << endl;
  cout << "neighbours in rT: " << rT.neighbours << endl;
  cout << "rT=" << rT << endl;

  output << "rT = " << rT << ";\n";

  double Pdist[12][12];
  
  matrix<double> D = surface_distances(T);
  matrix<int>    Dsqr = fullerene_surface_distances_sqr(T);
 
  output << "Dist = " << D << ";\n";
  output << "iDist = Sqrt[" << Dsqr << "];\n";
  
  for(int i=0;i<12;i++)
    for(int j=0;j<12;j++)
      Pdist[i][j] = D(i,j);

  rT.layout2d = rT.tutte_layout();

  // return 0; /* NB: Comment out this line to inspect previous results in output/reduce-graph-CN.m if Delaunayify gets stuck. */
  Triangulation dT = Delaunayify(rT,Pdist);



  output << "dT = " << dT << ";\n";

  output.close();
  return 0;
}

