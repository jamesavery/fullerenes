#include <limits.h>
#include "fullerenes/buckygen-wrapper.hh"
#include "fullerenes/triangulation.hh"
#include "fullerenes/polyhedron.hh"

#include <iostream>
#include <array>


// TODO: Sensible class structure for clustering

namespace clustering {
  constexpr int n_nodes_max=14;
  constexpr int distance_max=65535;  
  typedef uint16_t  distance_t;
  typedef uint16_t  csr_offset_t; 	// Up to 2^16-1 edges  
  typedef uint8_t   node_t;	
  typedef uint16_t  bitset_t;   	

  int popcount(bitset_t v){
    // From Stanford's bit-twiddling hacks, counts up to 14 bits in 3 machine instructions. Change for n_nodes>14.
    // https://graphics.stanford.edu/~seander/bithacks.html#CountBitsSet64
    return (uint64_t(v) * 0x200040008001ULL & 0x111111111111111ULL) % 0xf;
  }
 
  template <typename t> struct static_stack: public vector<t> {
    int i;

    static_stack(int n) : vector<t>(n), i(0) {}

    t pop(){
      assert(i>0);
      i--;
      return (*this)[i];
    }

    int push(const t &x) {
      assert(i<this->size());
      (*this)[i++] = x;
      return i;
    }

    bool empty(){ return (i==0); }
  };

  struct csr_adjacency {

    vector<csr_offset_t> row_starts;
    vector<node_t>       neighbours;

    int  n_neighbours(node_t u)              const { return row_starts[u+1]-row_starts[u]; }
    const node_t operator()(node_t u, int i) const { return neighbours[row_starts[u]+i];   }
    node_t& operator()(node_t u, int i)            { return neighbours[row_starts[u]+i];   }
  
    csr_adjacency(int n_nodes, int n_edges) : row_starts(n_nodes+1,0), neighbours(n_edges) {}
    csr_adjacency(int n_nodes, const vector<pair<node_t,node_t>> &edges) : row_starts(n_nodes+1), neighbours(2*edges.size()) {
      node_t counts[n_nodes];

      // Count how many neighbours each node has
      for(int i=0;i<n_nodes;i++) counts[i] = 0;    
      for(auto &e: edges){
	counts[e.first]++;
	counts[e.second]++;
      }

      // Cumsum into row_starts
      for(int i=0;i<n_nodes;i++)
	row_starts[i+1] = row_starts[i]+counts[i];

      // Fill in neighbour lists
      for(auto &e: edges){
	node_t u = e.first, v = e.second;
	neighbours[row_starts[u]+(--counts[u])] = v;
	neighbours[row_starts[v]+(--counts[v])] = u;      
      }
    }
  
  };

  
  struct dendrogram_node {
    uint16_t distance;
    uint8_t left, right;

    friend ostream &operator<<(ostream &s, const dendrogram_node &n) {
      s << vector<int>{n.distance,n.left,n.right};
      return s;
    }
  };

  struct dendrogram: public vector<dendrogram_node> {
    dendrogram(int capacity=12)          { reserve(capacity); }
    void merge(const dendrogram_node& n) { push_back({n.distance, min(n.left,n.right), max(n.left,n.right)}); }

    // Collect edge list for k-cluster dendrogram into CSR sparse matrix representation of graph
    csr_adjacency sparse_graph(int k) const {
      const vector<dendrogram_node> &dendro_edges(*this);
      const int N = dendro_edges.size()+1;

      vector<pair<node_t,node_t>> edges(N-k);
      for(int i=0;i<N-k; i++) edges[i] = make_pair(dendro_edges[i].left,dendro_edges[i].right);

      return csr_adjacency(N,edges);
    }
  
    // Separate into the two classes of the most distant clusters.
    // 1-bits represent the class containing first left element, 0-bits the one that doesn't
    vector<bitset_t> cluster_classes(int k) const {
      const vector<dendrogram_node> &edges(*this);
      const int N = edges.size()+1;
    
      bitset_t             visited = 0;
      vector<bitset_t>     clusters(k);
      static_stack<node_t> work_stack(N);

      csr_adjacency graph = sparse_graph(k);
      // cout << "\ndendrogram = " << *this << ";\n";
      // cout << "\ngraph = {" << graph.neighbours << "," << graph.row_starts << "};\n";
      // For each cluster
      for(int c=0;c<k;c++){

	// 1. Let u be the smallest node that has not yet been visited 
	node_t u=0;
	for(;u<N;u++) if(~visited & (1<<u)) break;
	work_stack.push(u);
      
	// 2. Depth first traversal of component c
	while(!work_stack.empty()){
	  u = work_stack.pop();
	  visited     |= (1<<u);	// Mark u as visited
	  clusters[c] |= (1<<u);	// Output u to current cluster

	  int n_v = graph.n_neighbours(u);
	  for(int i=0;i<n_v;i++){
	    node_t v = graph(u,i);
	    if(~visited & (1<<v)) work_stack.push(v);
	  }
	}
      }
      return clusters;    
    }
  };

  // Symmetric k x k-matrix of distances between clusters (minimum distance and maximum distance).
  // Diagonal elements in min-matrix measures minimal distance to neighbour within cluster (i.e., is not 0),
  // and diagonal elements in max-matrix is maximum distance within cluster.
  pair<matrix<distance_t>,matrix<distance_t>>
  cluster_distances(const vector<bitset_t> &clusters,
		    const matrix<distance_t> &P)
  {
    int k = clusters.size();
    int n = P.n;
    matrix<distance_t> min_dist(k,k), max_dist(k,k);

    for(int c1=0;c1<k;c1++)
      for(int c2=c1;c2<k;c2++){
	bitset_t   C1 = clusters[c1], C2 = clusters[c2];      
	distance_t mn = distance_max, mx = 0;

	// For each pair
	for(int i=0;i<n;i++)
	  if(C1 & (1<<i))
	    for(int j=i+1;j<n;j++)
	      if(C2 & (1<<j)) {
		mn = min(mn, P(i,j));
		mx = max(mx, P(i,j));
	      }

	min_dist(c1,c2) = mn;
	min_dist(c2,c1) = mn;
	max_dist(c1,c2) = mx;
	max_dist(c2,c1) = mx;      
      }
    return make_pair(min_dist,max_dist);
  }
  


  // TODO:
  //  1. Færdiggør debugging
  //  2. Halver memory-footprint med pakket symmetrisk matrix
  dendrogram hierarchical_clustering(const matrix<distance_t>& P)
  {
    size_t N = P.n;    
    matrix<uint8_t> dist = P;
    dendrogram class_tree(N-1);

    uint8_t order[N], row[N];
    for(int i=0;i<N;i++) order[i] = i;

    for(int h=0;h<=N-2;h++){
      int min_length = 0xffff;

      //uint8_t min_length = 0xff;
      int A=-1,B=-1;

      // Find smallest distance between clusters
      for(uint8_t i=0;i<N-h;i++)
	for(uint8_t j=i+1;j<N-h;j++)
	  if(dist(i,j) != 0 && dist(i,j) <  min_length)
	    min_length = dist(i,j), A = i, B = j;
    
      /* 
	 assert(A < B);
	 assert(A != 0xff);
	 assert(B != 0xff);

	 for(int i=0;i<N-h;i++)
	 for(int j=0;j<N-h;j++)
	 if(dist(i,j) != dist(j,i)) abort();

	 // A = min(A,B), B = max(A,B) per konstruktion
	 cout << "dist"<<h<<" = " << dist <<";\n";
      */
      //    printf("# merge (%d,%d) at %d\n",A,B,dist(A,B));

    
      // Merge equivalence classes
      class_tree.merge({dist(A,B),order[A],order[B]});
    
      // Update distance matrix.
    
      // 1. Set dist[A,:] = maximum( dist[A,:], dist[B,:] )
      //        dist[:,A] = maximum( dist[:,A], dist[:,B] )
      //
      // Copy
      for(uint8_t i=0;i<N;i++) row[i] = (i==A || i==B)? 0 : max(dist(A,i),dist(B,i)); 
      // Update
      for(uint8_t i=0;i<N;i++){
	dist(A,i) = row[i];
	dist(i,A) = row[i];
	dist(B,i) = row[i];
	dist(i,B) = row[i];
      }
      //    cout << "dist"<<h<<"b = " << dist <<";\n";

    
      // 2. Reduce dimension: Swap last row/col into position B.
      for(uint8_t i=0;i<N;i++) row[i] = (i!=B)? dist(N-h-1,i) : 0;     
      swap(order[B], order[N-h-1]);

      for(uint8_t i=0;i<N-h-1;i++){
	dist(B,i) = row[i];
	dist(i,B) = row[i];
      }    
      //    cout << "dist"<<h<<"c = " << dist <<";\n";    

    }
    return class_tree;
  }
}

matrix<double> self_distances(const Triangulation &G, const vector<node_t> &nodes)
{
  vector<int> nodes_inverse(G.N,-1);
  for(node_t U=0;U<nodes.size();U++){
    node_t     u = nodes[U];
    nodes_inverse[u] = U;
  }

  // Initialize H to graph distances, which are upper bound to surface distances,
  matrix<int>          Dg(nodes.size(),nodes.size(),G.all_pairs_shortest_paths(nodes));
  vector<double>       D(nodes.size());
  vector<int>          M(nodes.size(),0);// M[u] = max_v(d_g(u,v)) is upper bound to surface distance from u

  for(node_t U=0; U<nodes.size();U++){
    D[U] = FLOAT_MAX;
    for(node_t V=0;V<nodes.size();V++){
      M[U]   = max(M[U], 2*Dg(U,V));
    }

    for(node_t u:nodes){
      node_T U = nodes_inverse[u];

      // We want to explore all triangular slices with (1,0) along each edge of u
      for(int axis=0;axis<neighbours[u].size();axis++){
	int a = 0, b = 1, c = 1, d = M[U];

	// The simple geodesics correspond to the coprime pairs (a,b), a>=b.
	// Since the graph diameter M is an upper bound to shortest path length,
	// shortest geodesics must also have a^2 + ab + b^2 <= M^2.
	// We generate the coprime pairs obeying this inequality as the M-Farey sequence,
	// and treat (a,0) separately
	while(c <= M[U]){
	  // For each coprime pair (a,b) we get geodesics from u to v=end_of_theline(u,axis,(a,b)*n) for n = 1 until ||(a,b)*n||<=M
	  // BUT: if deg(v) != 6, there is a split, as the line can be continued in multiple ways. For deg(v)=5, it splits in two
	  auto simple_path = draw_path(a,b);

	  stack<pair<dedge_t,Eisenstein>> workset;
	  workset.push({{u,axis},{a,b}});
	  while(!workset.empty()){                                //   y---z     s---v
	    auto uxab = workset.pop();                            //  / \ / ... / \ /
	    dedge_t sv = end_of_the_line(uxab.first,uxab.second); // u---x     t---r  
	    Eisenstein dut = uwab.second
	      
	      switch(G.degree(v)){
	      case 5:
		workset.push({v,CW(wv,3)});
		workset.push({v,CCW(wv,3)});	
	      case 6:
		workset.push({v,CW(wv,3)/*How do we get the right axis?*/})
	      case 7:
		workset.push({v,CW(wv,3)});
		workset.push({v,CCW(wv,3)});
	      default:
		fprintf(stderr,"Geodesics not implemented for degree-%d curvature",G.degree(v));
		abort();
	      }
	    }

	  // Generate next (a,b).
	  // TODO: Since we generate (a,b) in order of ascending fraction a/b, we can do some
	  // dynamic programming to reduce the run-time by a full order (floor(a/b)=1, floor(a/b)=2,...)
	  int k = (Nmax+b) / d;
	  int a_next = c, b_next = d, c = k*c - a, d = k*d - b;
	  a = a_next, b = b_next;
	}
      }
    }
    
  for(node_t u: nodes){
    for(int i=0;i<neighbours[u].size();i++){
      node_t U  = nodes_inverse[u];

      for(int a=1; a<M[U]; a++){	
	for(int b=0; a*a + a*b + b*b < M[U]*M[U]; b++){
	  int d = gcd(a,b);
	  int ad = a/d, bd = b/d;

	  node_t w = u;
	  for(int i=0;i<d;i++){
	    const node_t v = end_of_the_line(w,i,ad,bd);

	  if(nodes_inverse[v] == u){ // Endpoint v is in nodes

	    if(d_sqr < H(U,V)){
	      //	      cout << u << "->" << vector<int>{{a,b,d_sqr}} << "->" << v <<endl;
	      H(U,V) = d_sqr;
	      G(U,V) = simple_geodesic(a,b,i);
	    }
	  }
	}
      }
    }
  }
  //  cout << "Hend = " << H << endl;    
  return G;
  
}

// Graph distance:
//  Square of all-pairs shortest paths matrix between pentagon nodes.
//  (Squared in order to be comparable to the surface distance matrix, which tracks the square distances as integers)
matrix<clustering::distance_t> pentagon_graph_distance(const Triangulation &G)
{
  vector<int> pentagon_indices(12);
  for(int u=0, i=0;u<G.N;u++) if(G.neighbours[u].size() == 5) pentagon_indices[i++] = u;

  auto D = G.all_pairs_shortest_paths(pentagon_indices);
  for(int i=0;i<D.m*D.n;i++) D[i] *= D[i];

  return D;
}

matrix<int> pentagon_surface_distance(const Triangulation& G)
{
  vector<int> pentagon_indices(12);
  for(int u=0, i=0;u<G.N;u++) if(G.neighbours[u].size() == 5) pentagon_indices[i++] = u;
								
  auto Hsqr = G.simple_square_surface_distances(pentagon_indices,true);
  cerr << "Hsqr = " << Hsqr << ";\n";  
  auto H    = G.surface_distances(pentagon_indices,true);
  cerr << "H    = " << H << ";\n";

  return H;
}

auto pentagon_geodesics(const Triangulation& G)
{
  vector<int> pentagon_indices(12);
  for(int u=0, i=0;u<G.N;u++) if(G.neighbours[u].size() == 5) pentagon_indices[i++] = u;
								
  return G.simple_geodesics(pentagon_indices,true);
}

vector<vector<tri_t> > tri_runs(const Triangulation& dG, const node_t u,const Triangulation::simple_geodesic geo)
{
  vector<vector<node_t>> quad_runs = dG.quads_of_the_line(u,geo.axis,geo.g.first, geo.g.second);
  vector<vector<tri_t>>  tri_runs;
  
  for(auto run: quad_runs){
    int n = run.size()/2;	// Number of columns 
    vector<tri_t> tri_run((n-1)*2);

    // Iterate over quads      
    //   s---t
    //  / \ /
    // q---r      
    for(int i=0;i+1<n;i++){
      node_t q = run[i*2+0], s = run[i*2+1], r = run[i*2+2], t = run[i*2+3];
      tri_run[i*2  ] = {q,r,s};
      tri_run[i*2+1] = {s,r,t};
    }
    tri_runs.push_back(tri_run);
  }
  return tri_runs;
}

vector< pair<dedge_t,float> > crossings(const Triangulation& dG, const node_t u, const Triangulation::simple_geodesic geo)
{
  vector<pair<dedge_t,float>> result;  
  vector<vector<node_t>> quad_runs = dG.quads_of_the_line(u,geo.axis,geo.g.first, geo.g.second);

  int a = geo.g.first, b = geo.g.second, axis = geo.axis;
  int minor = min(a,b), major = max(a,b);
  int Nruns = quad_runs.size();

  // No matter what, we start in node u. 
  node_t q = u, r = dG.neighbours[u][geo.axis];
  result.push_back({{q,r},0});
  
  // Deal with degenerate case first
  // 1. Edge-aligned path
  if(minor == 0){
    assert(quad_runs.size()==1);
    auto run = quad_runs[0];
    int    n = run.size()/2;
    
    for(int i=0;i+1<n;i++){
      node_t q = run[i*2+0], r = run[i*2+2];

      result.push_back( {{q,r},1} );
    }
    return result;
  }
  // 2. Diagonal path
  if(minor == major){
    for(int I=0,J=0;J<Nruns;I++,J++){
      auto run = quad_runs[I];
      //      cerr << "run size:" << (run.size()/2-1) << endl;
      node_t q = run[0], s = run[1], r = run[2], t = run[3];
   
      result.push_back( {{s,r},0.5} );       // s-r crossed half-way
      result.push_back( {{r,t},  1} );       //   t hit spot-on
    }
    return result;
  }

  // Everything else is a line that 
  //  - In the final run:  ends in the upper right corner of the final quad
  //  - In previous runs:  crossing the upper edge (s--t) in the final quad
  for(int I=0, J=0;J<Nruns;J++){
    auto     run = quad_runs[J];
    int        n = run.size()/2-1; // Number of quads in Ith run 
    float lambda = 0;
    
    // Iterate over quads      
    //   s---t
    //  / \ /
    // q---r

    // Point on line-segment (x0,y0) + \lambda (x1,y1) determined by
    //     \lambda = -(a*x0 - b*y0) / (a*x1 - b*y1)    
    auto intersect = [&](int x0,int y0,int x1,int y1) -> float
      { return (minor*x0 - major*y0)/float(minor*(x0-x1)-major*(y0-y1)); };

    // For each quad in the run, the positions in the Eisenstein-strip of the vertices are:
    // q: (I,J) r: (I+1,J), s: (I,J+1), t: (I+1,J+1)
    for(int i=0;i<n;i++,I++){
      node_t q = run[i*2+0], s = run[i*2+1], r = run[i*2+2], t = run[i*2+3];

      // Which edges are crossed, and where?
      
      // 1. Edge s,r is always crossed
      lambda = intersect(/*s*/I,J+1, /*r*/I+1,J); 
      result.push_back( {{s,r}, lambda} );
      
      if(i+1<n){
	// If this is not the final quad of the run, we exit to the East:
	// 2. Edge r,t
	lambda = intersect(/*r*/I+1,J, /*t*/I+1,J+1);
	result.push_back( {{r,t}, lambda} );
      } else // If this is the final quad of the run, we	
	if(I+1<Nruns){ // Exit to the north, if this is not the final run
	  // 3a. Edge s,t
	  lambda = intersect(/*s*/I,J+1, /*t*/I+1,J+1);
	  result.push_back( {{s,t}, lambda} );	  
	} else        // Hit the NE corner, if this is the final run
	  // 3b. Vertecx t
	  result.push_back({{r,t},1});	  
    }
  }

  return result;
}


vector<coord3d> line_points(const vector<pair<dedge_t,float>> &crossings, const vector<coord3d> &points)
{
  vector<coord3d> result(crossings.size());
  
  for(int i=0;i<crossings.size();i++){
    auto c = crossings[i];
    node_t s = c.first.first, t = c.first.second;
    float lambda = c.second;

    coord3d x0 = points[s], x1 = points[t];
    result[i] = x0+(x1-x0)*lambda;
  }
  return result;
}


int main(int ac, char **argv)
{
  int N                = ac>=2? strtol(argv[1],0,0) : 20;     // Argument 1: Number of vertices N
  
  if(N<20 || N==22 || N&1){
    fprintf(stderr,"Syntax: %s <N:int> [output_dir] [IPR:0|1] [only_nontrivial:0|1]\n",argv[0]);
    return -1;
  } 
 
  int pick_number     = ac>=3? strtol(argv[2],0,0) : 1;
  
  string output_dir   = ac>=4? argv[3] : "output";    // Argument 2: directory to output files to
  int IPR             = ac>=5? strtol(argv[4],0,0):0; // Argument 3: Only generate IPR fullerenes?
  int only_nontrivial = ac>=6? strtol(argv[5],0,0):0; // Argument 4: Only generate fullerenes with nontrivial symmetry group?

  ofstream failures((output_dir+"/failures.txt").c_str()); // output/failures.txt contains list of any fullerenes that failed optimization

  int i=0;  
  FullereneDual dualG;
  BuckyGen::buckygen_queue Q = BuckyGen::start(N,IPR,only_nontrivial);  
  while(BuckyGen::next_fullerene(Q,dualG)) if((++i) == pick_number) {      
      cerr << "isomer "<< i << endl;
      dualG.update();
      
      Polyhedron P  = Polyhedron::fullerene_polyhedron(dualG.dual_graph());
      Polyhedron dP = P.dual();

      Polyhedron::to_file(P, "output/P.mol2");
      Polyhedron::to_file(dP,"output/dP.mol2");      
      
      ofstream output_file("output/geodesics.m");
      LIST_OPEN = '{'; LIST_CLOSE = '}';
      output_file << "neighbours = " << dP.neighbours << "+1;\n"
      		  << "faces      = " << dP.faces      << "+1;\n"
		  << "points     = " << dP.points     << ";\n";

      dualG = FullereneDual(dP);

      auto Ds = pentagon_surface_distance(dualG); 
      //      auto Gs = pentagon_geodesics(dualG);

      output_file << "Ds = " << Ds << ";\n";
      
      // for(int i=0;i<12;i++){
      // 	Eisenstein g(Gs(i,i).g);
      // 	int axis = Gs(i,i).axis;
      // 	cout << make_pair(i,i) << " shortest self-geodesic";
      // 	if(g != Eisenstein{0,0}){
      // 	  cout << ": " << make_pair(g,axis) << "\n";
      //     cout << "\t Path: " << dualG.quads_of_the_line(i,axis,g.first,g.second) << endl;										 
      // 	} else {
      // 	  cout << " passes through cone-point.\n"; 
      // 	}
      // }

      // //      Gs(0,0).g = {3,2};
      // auto crossings00 = crossings(dualG,0,Gs(0,0));
      // auto run00       = tri_runs(dualG,0,Gs(0,0));
      // auto linepts3D00 = line_points(crossings00,dP.points);
      // output_file << "geodesic00 = " <<Gs(0,0).g << ";\n";
      // output_file << "run00 = " << run00 << "+1;\n";
      // output_file << "crossings00 = " << crossings00 << ";\n";
      // output_file << "linepts3D00 = " << linepts3D00 << ";\n";      
      
  }
  failures.close();
  
  return 0;
}
