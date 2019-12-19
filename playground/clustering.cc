#include <limits.h>
#include "fullerenes/buckygen-wrapper.hh"
#include "fullerenes/triangulation.hh"
#include "fullerenes/polyhedron.hh"

#include <iostream>

struct dendrogram_node {
  uint16_t distance;
  uint8_t left, right;

  friend ostream &operator<<(ostream &s, const dendrogram_node &n) {
    s << vector<int>{n.distance,n.left,n.right};
    return s;
  }
};

struct dendrogram: public vector<dendrogram_node> {
  dendrogram(int capacity=12) { reserve(capacity); }
  void merge(const dendrogram_node& n) { push_back(n); }
};



// TODO:
//  1. Færdiggør debugging
//  2. Halver memory-footprint med pakket symmetrisk matrix
dendrogram hierarchical_clustering(const matrix<uint8_t>& P)
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
    printf("# merge (%d,%d) at %d\n",A,B,dist(A,B));
    */
    
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

int main()
{
  //  vector<int> names{{7, 10, 20, 28, 35}};
  //vector<int> names{{35, 20, 7, 28, 10}};  
  //vector<int> names{{1,2,5,10}};
  vector<int> names{{12, 20, 9, 13, 17, 14, 1, 5, 8, 3, 16, 18}};
  matrix<uint8_t> P(dist1d(names));
  
  dendrogram clusters;
  for(int i=0;i<1000000;i++)
    clusters = hierarchical_clustering(P);
    
  cout << "clusters_raw = " << clusters << ";\n";
  cout << "names = " << names << ";\n";
  for(int i=0;i<clusters.size();i++){
    clusters[i].left  = names[clusters[i].left];
    clusters[i].right = names[clusters[i].right];
  }
  cout << "clusters = " << clusters << ";\n";

  return 0;
}
