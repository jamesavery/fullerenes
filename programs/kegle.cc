#include <limits.h>
#include <getopt.h>
#include <set>
#include "fullerenes/buckygen-wrapper.hh"
#include "fullerenes/triangulation.hh"

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

	// For each pair (i,j) with i<j belonging to C1 x C2 in either order
	for(int i=0;i<n;i++)
	  for(int j=i+1;j<n;j++)
	    if(((C1 & (1<<i)) && (C2 & (1<<j))) ||
	       ((C2 & (1<<i)) && (C1 & (1<<j)))) {
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


  dendrogram hierarchical_clustering(const matrix<distance_t>& P)
  {
    size_t N = P.n;
    matrix<distance_t> dist(P);
    dendrogram class_tree(N-1);

    node_t order[N];
    distance_t row[N];
    for(int i=0;i<N;i++) order[i] = i;

    for(int h=0;h<=N-2;h++){
      int min_length = 0xffff;

      int A=-1,B=-1;

      // Find smallest distance between clusters
      for(uint8_t i=0;i<N-h;i++)
	for(uint8_t j=i+1;j<N-h;j++)
	  if(dist(i,j) != 0 && dist(i,j) <  min_length)
	    min_length = dist(i,j), A = i, B = j;

      // Merge equivalence classes
      class_tree.merge({dist(A,B),order[A],order[B]});

      // Update distance matrix (complete linkage).

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

      // 2. Reduce dimension: Swap last row/col into position B.
      for(uint8_t i=0;i<N;i++) row[i] = (i!=B)? dist(N-h-1,i) : 0;
      swap(order[B], order[N-h-1]);

      for(uint8_t i=0;i<N-h-1;i++){
	dist(B,i) = row[i];
	dist(i,B) = row[i];
      }
    }
    return class_tree;
  }
}

// Collect the vertex indices of the 12 pentagons (degree-5 vertices) in the dual
vector<int> pentagon_indices(const Triangulation& G)
{
  vector<int> indices;
  indices.reserve(12);
  for(int u=0;u<G.N;u++) if(G.degree(u) == 5) indices.push_back(u);
  return indices;
}

// Graph distance:
//  Square of all-pairs shortest paths matrix between pentagon nodes.
//  (Squared in order to be comparable to the surface distance matrix, which tracks the square distances as integers)
matrix<clustering::distance_t> pentagon_graph_distance(const Triangulation &G, const vector<int>& pent_idx)
{
  auto D = G.all_pairs_shortest_paths(pent_idx);
  for(int i=0;i<D.m*D.n;i++) D[i] *= D[i];
  return D;
}

matrix<int> pentagon_surface_distance(const Triangulation& G, const vector<int>& pent_idx)
{
  return G.simple_square_surface_distances(pent_idx);
}

auto pentagon_geodesics(const Triangulation& G, const vector<int>& pent_idx)
{
  return G.simple_geodesics(pent_idx, true);
}


// Squared surface-distance diameter of the tip cluster including its degree-6 boundary.
// The boundary consists of all degree-6 vertices adjacent to any tip pentagon.
// This gives a non-degenerate diameter even for single-pentagon tips (p=1).
int tip_boundary_diameter(const Triangulation& G,
			  const vector<int>& pent_idx,
			  clustering::bitset_t tip_cluster)
{
  // Collect tip pentagon vertex indices
  vector<int> tip_vertices;
  for(int i=0;i<12;i++)
    if(tip_cluster & (1<<i))
      tip_vertices.push_back(pent_idx[i]);

  // Collect boundary: degree-6 vertices adjacent to any tip pentagon
  set<int> seen(tip_vertices.begin(), tip_vertices.end());
  vector<int> extended = tip_vertices;

  for(int u : tip_vertices)
    for(auto v : G[u])
      if(G.degree(v) == 6 && seen.find(v) == seen.end()) {
	extended.push_back(v);
	seen.insert(v);
      }

  // Surface squared distances for extended tip set
  auto D = G.simple_square_surface_distances(extended);

  // Return max distance (diameter)
  int diameter = 0;
  int m = D.n;
  for(int i=0;i<m;i++)
    for(int j=i+1;j<m;j++)
      diameter = max(diameter, D(i,j));
  return diameter;
}


// Classify a fullerene dual as a (p,q)-nanocone based on pentagon clustering.
//
// Returns (p, q) with p <= q if the isomer is a nanocone, or (0, 0) if not.
// Also outputs the cluster assignments and distances.
//
// A (p,q)-nanocone (p+q=12) has p pentagons at the tip and q=12-p at the base.
// Detection criteria:
//   1. Hierarchical clustering (k=2) gives clusters of size p and q
//   2. Tip boundary diameter (including degree-6 ring) <= max_tip_diam
//   3. Inter-cluster min graph distance >= alpha * tip boundary diameter (separation)
//
// Special case: (6,6) is a nanotube (both caps equivalent).
pair<int,int> classify_nanocone(const Triangulation& G,
				const vector<int>& pent_idx,
				const matrix<clustering::distance_t>& Ds,
				double alpha,
				int max_tip_diam,
				// outputs:
				vector<clustering::bitset_t>& clusters,
				matrix<clustering::distance_t>& Dc_min_out,
				matrix<clustering::distance_t>& Dc_max_out,
				int& tip_bd_diam,
				int& separation_out)
{
  auto hierarchy = clustering::hierarchical_clustering(Ds);
  clusters = hierarchy.cluster_classes(2);
  auto Dc = clustering::cluster_distances(clusters, Ds);
  Dc_min_out = Dc.first;
  Dc_max_out = Dc.second;

  int p = clustering::popcount(clusters[0]), q = clustering::popcount(clusters[1]);
  int tip = 0, base = 1;
  if(p > q) { swap(p, q); tip = 1; base = 0; }

  // Compute tip boundary diameter in surface squared distance (includes degree-6 ring)
  tip_bd_diam = tip_boundary_diameter(G, pent_idx, clusters[tip]);

  // Criterion 1: tip is compact
  if(tip_bd_diam > max_tip_diam) { separation_out = 0; return {0, 0}; }

  // Criterion 2: tip is well-separated from base
  // Dc_min_out is already computed from the surface distance matrix Ds
  separation_out = Dc_min_out(tip, base);

  if(separation_out < alpha * tip_bd_diam) return {0, 0};

  return {p, q};
}


const char *cone_type_name(int p, int q)
{
  if(p == 6 && q == 6) return "nanotube";
  if(p == 5 && q == 7) return "nanocone(5,7)";
  if(p == 4 && q == 8) return "nanocone(4,8)";
  if(p == 3 && q == 9) return "nanocone(3,9)";
  if(p == 2 && q == 10) return "nanocone(2,10)";
  if(p == 1 && q == 11) return "nanocone(1,11)";
  return "nanocone";
}


void usage(const char *progname)
{
  fprintf(stderr,
	  "Usage: %s [options] <N>\n"
	  "\n"
	  "Scan all C_N fullerene isomers for nanocones and nanotubes.\n"
	  "\n"
	  "Options:\n"
	  "  -q          Quiet: only output detected structures\n"
	  "  -v          Verbose: repeat for more detail (-vv, -vvv)\n"
	  "  -I          IPR isomers only\n"
	  "  -S          Only nontrivial symmetry groups\n"
	  "  -a <float>  Separation threshold alpha (default: 1.5)\n"
	  "  -d <int>    Max tip boundary diameter in surface squared distance (default: N/2)\n"
	  "\n"
	  "Verbosity levels:\n"
	  "  default     One line per detection + summary counts\n"
	  "  -v          Add cluster distances for detected structures\n"
	  "  -vv         Add distance matrix + cluster info for all isomers\n"
	  "  -vvv        Add self-geodesics and paths for all isomers\n",
	  progname);
}


int main(int ac, char **argv)
{
  int verbose = 0;
  int IPR = 0;
  int only_nontrivial = 0;
  double alpha = 1.5;
  int max_tip_diam = 0; // 0 = auto

  int opt;
  while((opt = getopt(ac, argv, "qvISa:d:h")) != -1) {
    switch(opt) {
    case 'q': verbose = -1; break;
    case 'v': verbose++; break;
    case 'I': IPR = 1; break;
    case 'S': only_nontrivial = 1; break;
    case 'a': alpha = strtod(optarg, 0); break;
    case 'd': max_tip_diam = strtol(optarg, 0, 0); break;
    case 'h': usage(argv[0]); return 0;
    default:  usage(argv[0]); return -1;
    }
  }

  if(optind >= ac) { usage(argv[0]); return -1; }
  int N = strtol(argv[optind], 0, 0);

  if(N<20 || N==22 || N&1){
    fprintf(stderr,"Error: N must be an even integer >= 20, and N != 22.\n");
    return -1;
  }

  if(max_tip_diam == 0) max_tip_diam = N/2; // auto: scales with fullerene size

  int n_isomers = 0;
  int n_cones[7] = {}; // n_cones[p] counts (p, 12-p)-nanocones; p=6 is nanotube

  Graph dualG_buf;
  BuckyGen::buckygen_queue Q = BuckyGen::start(N, IPR, only_nontrivial);

  while(BuckyGen::next_fullerene(Q, dualG_buf)){
    FullereneDual dualG(dualG_buf);
    n_isomers++;

    if(verbose >= 0 && n_isomers % 100000 == 0)
      cerr << "isomer " << n_isomers << endl;

    auto pent_idx = pentagon_indices(dualG);
    auto Ds = pentagon_surface_distance(dualG, pent_idx);

    // Geodesic computation (level 3 only)
    if(verbose >= 3) {
      auto Gs = pentagon_geodesics(dualG, pent_idx);
      cout << "Ds = " << Ds << endl;
      for(int j=0;j<12;j++){
	Eisenstein g(Gs(j,j).g);
	int axis = Gs(j,j).axis;
	cout << j<<"-"<<j<<" shortest self-geodesic: " << make_pair(g,axis) << "\n";
	cout << "\t Path: " << dualG.quads_of_the_line(j,axis,g.first,g.second) << endl;
      }
    }

    // Classify
    vector<clustering::bitset_t> Cs;
    matrix<clustering::distance_t> Dc_min(2,2), Dc_max(2,2);
    int tip_bd_diam, separation;
    auto [p, q] = classify_nanocone(dualG, pent_idx, Ds, alpha, max_tip_diam,
				    Cs, Dc_min, Dc_max, tip_bd_diam, separation);

    // Level 2: all isomers get cluster info
    if(verbose >= 2) {
      int p2 = clustering::popcount(Cs[0]), q2 = clustering::popcount(Cs[1]);
      cout << "Isomer " << n_isomers << ": clusters (" << min(p2,q2) << "," << max(p2,q2) << ")"
	   << " tip_bd=" << tip_bd_diam
	   << " Dc_min=" << Dc_min << " Dc_max=" << Dc_max << "\n";
    }

    if(p == 0) continue; // not a nanocone

    n_cones[p]++;
    general_spiral spiral = dualG.get_rspi();

    // Level 0 (default): one line per detection
    cout << "Isomer " << n_isomers << ": " << cone_type_name(p, q)
	 << " spiral=" << (spiral.spiral_code+1) << "\n";

    // Level 1: add cluster distances and separation
    if(verbose >= 1) {
      cout << "  tip_bd=" << tip_bd_diam
	   << " separation=" << separation
	   << " ratio=" << (tip_bd_diam > 0 ? (double)separation/tip_bd_diam : 0)
	   << " Dc_min=" << Dc_min << " Dc_max=" << Dc_max << "\n";
    }
  }

  // Summary
  cout << "\n" << n_isomers << " isomers scanned.\n";
  for(int p=1;p<=6;p++)
    if(n_cones[p] > 0)
      cout << "  " << cone_type_name(p, 12-p) << ": " << n_cones[p] << "\n";

  return 0;
}
