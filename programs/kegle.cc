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
								
  return G.simple_square_surface_distances(pentagon_indices);
}

auto pentagon_geodesics(const Triangulation& G)
{
  vector<int> pentagon_indices(12);
  for(int u=0, i=0;u<G.N;u++) if(G.neighbours[u].size() == 5) pentagon_indices[i++] = u;
								
  return G.simple_geodesics(pentagon_indices,true);
}


int main(int ac, char **argv)
{
  int N                = ac>1? strtol(argv[1],0,0) : 20;     // Argument 1: Number of vertices N
  
  if(N<20 || N==22 || N&1){
    fprintf(stderr,"Syntax: %s <N:int> [output_dir] [IPR:0|1] [only_nontrivial:0|1]\n",argv[0]);
    return -1;
  } 

  string output_dir   = ac>=3? argv[2] : "output";    // Argument 2: directory to output files to
  int IPR             = ac>=4? strtol(argv[3],0,0):0; // Argument 3: Only generate IPR fullerenes?
  int only_nontrivial = ac>=5? strtol(argv[4],0,0):0; // Argument 4: Only generate fullerenes with nontrivial symmetry group?
  
  ofstream failures((output_dir+"/failures.txt").c_str()); // output/failures.txt contains list of any fullerenes that failed optimization

  int i=0;  
  FullereneDual dualG;
  BuckyGen::buckygen_queue Q = BuckyGen::start(N,IPR,only_nontrivial);  

  while(BuckyGen::next_fullerene(Q,dualG)){ // Generate all appropriate C_N isomer duals 
    i++;

    //    if(i%100000 == 0)
      cerr << "isomer "<< i << endl;
    //    dualG.update();		            // Update metadata

    auto Ds = pentagon_surface_distance(dualG); //
    auto Gs = pentagon_geodesics(dualG);

    cout << "Ds = " << Ds << endl;

    // cout << "Gs = [";
    // for(int i=0;i<12;i++){
    //   cout << "[";
    //   for(int j=0;j<12;j++)
    // 	cout << make_pair(Gs(i,j).g,Gs(i,j).axis) << (j+1<Gs.n? ", ":"]\n");
    // }
    // cout << "]\n";

    for(int j=0;j<12;j++){
      Eisenstein g(Gs(j,j).g);
      int axis = Gs(j,j).axis;
      cout << j<<"-"<<j<<" shortest self-geodesic: " << make_pair(g,axis) << "\n";
      cout << "\t Path: " << dualG.quads_of_the_line(j,axis,g.first,g.second) << endl;
    }
    
    //auto Dg = pentagon_graph_distance(dualG);

    auto D = Ds;
    auto hierarchy = clustering::hierarchical_clustering(D); 
    auto Cs        = hierarchy.cluster_classes(2);
    auto Dc        = clustering::cluster_distances(Cs,D);
    auto Dc_min = Dc.first, Dc_max = Dc.second;
    cout << "Ds = " << Ds << endl;
    cout << "Dc_min = " << Dc_min << "\nDc_max = " << Dc_max << "\n\n";

    
    // // bool is_nanotube(const Triangulation &G)
    // // Rules:
    // //        1. Two pentagon clusters, each with 6 pentagons (caps, Gauss curvature 2pi each)
    // //        2. The two clusters are further away from each other than their own diameter (e.g.: at least 1*, 3/2*, 2* - which?)
    // if(clustering::popcount(Cs[0]) == 6 && clustering::popcount(Cs[1]) == 6 &&
    //    Dc_min(0,1) >= 1.5*Dc_max(0,0) &&
    //    Dc_min(0,1) >= 1.5*Dc_max(1,1)) {
      
    //   cout << "Isomer " << i << " is a nanotube of sorts: ";
    //   general_spiral spiral = dualG.get_rspi();
    //   cout << (spiral.spiral+1) << endl;
    // }

    // bool is_nanocone(int p, int q, const Triangulation& G)
    // Rules:
    //         1. Two pentagon clusters, one with p pentagons and one with q, p<q
    //         2. Pointy:    p-cluster has smaller diameter than q-cluster
    // Either
    //         3a. Elongated: inter-cluster distance (min) is greater than diameter of q-cluster
    // or
    //         3b. Diameters of clusters constrained instead
    int p = clustering::popcount(Cs[0]), q = clustering::popcount(Cs[1]);
    int p_cluster = 0, q_cluster = 1;
    if(p>q){
      swap(p,q);
      swap(p_cluster, q_cluster);
    }
    
    if(clustering::popcount(Cs[0]) == 6 && clustering::popcount(Cs[1]) == 6 &&
       Dc_min(0,1) >= 1.5*Dc_max(0,0) &&
       Dc_min(0,1) >= 1.5*Dc_max(1,1)) {
      
      cout << "Isomer " << i << " is a nanotube of sorts: ";
      general_spiral spiral = dualG.get_rspi();
      cout << (spiral.spiral_code+1) << endl;
    }
    
    
#if 0

    //    auto C = cluster_sizes(pentagon_distance(dualG), 2, 4);
    if(C.first){
      cout << "C"<<i<<" = " << C.second<< ";\n";
      
      // if(0)
      //   if(max_length(pentagon_distance(dualG)) >= 10){
      //	cout << i << endl;
	
      dualG.update();		            // Update metadata
      FullereneGraph G = dualG.dual_graph();  // Construct fullerene graph
      vector<int> RSPI;
      jumplist_t jumps;
      G.get_rspi_from_fg(RSPI,jumps);
      cout << i << ", RSPI="<<(RSPI+1) << ", jumps="<<jumps<<endl;
    }
#endif
  }
  cout << i << " graphs.\n";
  failures.close();
  
  return 0;
}
