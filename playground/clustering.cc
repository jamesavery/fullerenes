#include <limits.h>
#include "fullerenes/buckygen-wrapper.hh"
#include "fullerenes/triangulation.hh"
#include "fullerenes/polyhedron.hh"

#include <iostream>
#include <array>

namespace clustering {
  constexpr int n_nodes_max=16;
  constexpr int distance_max=65535;  
  typedef uint16_t  distance_t;
  typedef uint16_t  csr_offset_t; 	// Up to 2^16-1 edges  
  typedef uint8_t   node_t;	
  typedef uint16_t  bitset_t;   	
  
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

  int  n_neighbours(node_t u)            const { return row_starts[u+1]-row_starts[u]; }
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
  // Max 16 elements - increase size to accommodate more elements

  
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
	  for(int j=0;j<n;j++)
	    if(C2 & (1<<j)) {
	      mn = ::min(mn, P(i,j));
	      mx = ::max(mx, P(i,j));
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

matrix<uint8_t> dist1d(vector<int> data)
{
  int n = data.size();
  matrix<uint8_t> dist(n,n);
  for(int i=0;i<n;i++)
    for(int j=0;j<n;j++)
      dist(i,j) = abs(data[i]-data[j]);
  return dist;
}
}

int main()
{
  using namespace clustering;
  //vector<int> names{{7, 10, 20, 28, 35}};
  //vector<int> names{{35, 20, 7, 28, 10}};  
  //vector<int> names{{1,2,5,10}};
  vector<int> names{{12, 20, 9, 13, 17, 14, 1, 5, 8, 3, 16, 18}};
  matrix<uint8_t> P(dist1d(names));
  
  dendrogram clusters;
  vector<bitset_t> Cs;  
  //  for(int i=0;i<1000000;i++){
  clusters = hierarchical_clustering(P);
  //    Cs = clusters.cluster_classes(2);
  //  if(i%100000==0) cout << Cs << endl;
    
  //  }

  Cs = clusters.cluster_classes(3);
  
  cout << "clusters_raw = " << clusters << ";\n";
  cout << "names = " << names << ";\n";

  printf("Two clusters at distance %d:\n",clusters[clusters.size()-2].distance);


  for(auto c: Cs){
    printf("\n%x: ",c);
    for(int i=0;i<P.n;i++) if(c & (1<<i)) printf("%d ",names[i]);
  }
  printf("\n");

  for(int i=0;i<clusters.size();i++){
    clusters[i].left  = names[clusters[i].left];
    clusters[i].right = names[clusters[i].right];
  }

  printf("Distances:\n");
  auto Dc = cluster_distances(Cs,P);

  cout << "Dmin = " << Dc.first  << endl;
  cout << "Dmxa = " << Dc.second << endl;  
  
  return 0;
}


