#include "fullerenes/triangulation.hh"

using namespace std;

typedef vector< pair<int,int> > path_t;
typedef vector<vector<path_t>> pattern_t;
typedef vector<vector<int>> distance_pattern_t;

bool debug = false;

#include <queue>

struct candidate {
  double error;
  vector<int> rspi;

  bool operator<(const candidate &b) const { return error < b.error; }
};

template <typename T> struct k_smallest : public priority_queue<T> {
public:
  size_t k;
  
  k_smallest(size_t k) : k(k) {}

  bool insert(const T& x){
    if(this->size() < k || x < this->top()){
      if(this->size() == k) this->pop();
      this->push(x);
      return true;
    }
    return false;
  }

  vector<T>& as_vector() {
    return (*this).*&k_smallest::c;
  }
  
};


int incident_edge(const Graph &G, node_t u, node_t v)
{
  for(int i=0;i<G.degree(v);i++) if(G.neighbours[v][i] == u) return i;
  return -1;			// Not a neighbour
}

Triangulation reverse(const Triangulation &G)
{
  Triangulation Gr = G;
  for(node_t u=0;u<G.N;u++) std::reverse(Gr.neighbours[u].begin(),Gr.neighbours[u].end());
  return Gr;
}

template <typename T, typename t> bool contains(const T& c, const t& x)
{
  return find(c.begin(),c.end(),x) != c.end();
}

template <typename T> vector<int> sort_permutation(const vector<T>& xs)
{
  vector<int> pi(xs.size());
  for(size_t i=0; i<xs.size(); i++) pi[i] = i;

  sort(pi.begin(), pi.end(), [&](int i, int j) { return (xs[i] < xs[j]); });

  return pi;
}

// void sort_by_row(const matrix<int>& A, int i, matrix<int>& R)
// {
//   vector<int> pi = sort_permutation(A(i));
//   for(int i=0;i<A.m;i++)
//     for(int j=0;j<A.m;j++)
//       R(pi[i],pi[j]) = A(i,j);
// }

struct nanocone {
  int nU, dminU, dmaxU;
  int nV, dminV, dmaxV;
  int dUV;

  // U = pentagons[0], V = 
  int pentagons[12];
};

bool match_nanocone(const Triangulation &G, nanocone &C, int dminUV, int dmax_all)
{
  // Fill out remaining pentagons
  vector<int> pentagons;
  for(node_t u=0;u<G.N;u++) if(G.degree(u)==5) pentagons.push_back(u);

  //  matrix<int> D = G.all_pairs_shortest_paths(pentagons,dmax_all); // Alternativt: G.convex_square_surface_distances(pentagons)
  matrix<int> D = G.simple_square_surface_distances(pentagons);

  // Find vertex 1 and 12: the furthest from each other.
  int &U = C.pentagons[0], &V = C.pentagons[11];
  C.dUV=0;
  for(int i=0;i<12;i++)
    for(int j=i+1;j<12;j++)
      if(D(i,j) > C.dUV){
	U = i;
	V = j;	
	C.dUV = D(U,V);
      }
  
  if(C.dUV < dminUV){
    // cerr << "Too compact: " << C.dUV << " < " << dminUV << endl;
    // cerr << "D = " << D << ";\n";
    return false; // Too compact to be called a nanocone
  } else {
       // printf("Long enough: %d >= %d\n",C.dUV, dminUV);
       // printf("U = %d, V =%d\n",U,V);
  }
  
  // In a nanocone, the caps should *at the very least* each be smaller than half the cone length.
  int dcap_max = C.dUV/2;
  C.nU = 0;          C.nV = 0;
  C.dminU = INT_MAX; C.dminV = INT_MAX;
  C.dmaxU = 0;       C.dmaxV = 0;
  for(int k=0;k<12;k++){
    if(D(U,k) <= dcap_max && k!=U){
      //      printf("Pentagon %d (node %d) belongs to U-cap\n",k,pentagons[k]);
      C.dminU  = min(C.dminU,D(U,k));
      C.dmaxU  = max(C.dmaxU,D(U,k));
      
      C.nU++;
      C.pentagons[C.nU] = k;
    }
    if(D(V,k) <= dcap_max && k!=V){
      C.dminV  = min(C.dminV,D(V,k));
      C.dmaxV  = max(C.dmaxV,D(V,k));
      C.nV++;
      //      printf("Pentagon %d (node %d) at distance %d is %dth node of V-cap (dmin=%d)\n",k,pentagons[k],D(V,k),C.nV,C.dminV);            
      C.pentagons[12-C.nV] = k;
    }
  }

  // We define U to be the cap with the smallest number of pentagons, 
  if(C.nU > C.nV){
    // Swap caps
    reverse(C.pentagons, C.pentagons+12);
    swap(C.nU,   C.nV);
    swap(C.dminU,C.dminV);
    swap(C.dmaxU,C.dmaxV);
  }

  return true;
}


struct filter_nanocones {
  int dminUV, dmax_all, dminU, dmaxU, dminV, dmaxV, nU, nV;
  filter_nanocones(int dminUV, int dmax_all, int dminU, int dmaxU, int dminV, int dmaxV, int nU=-1, int nV=-1) :
    dminUV(dminUV), dmax_all(dmax_all), dminU(dminU), dmaxU(dmaxU), dminV(dminV), dmaxV(dmaxV), nU(nU), nV(nV) {}

  bool operator()(const Triangulation &G) const
  {
    nanocone C;
    bool matched = match_nanocone(G,C,dminUV,dmax_all);
    
    if(matched){
      //      printf("U: %d vertices between length %d and %d.\n",C.nU,C.dminU,C.dmaxU);
      //      printf("V: %d vertices between length %d and %d.\n",C.nV,C.dminV,C.dmaxV);      
      matched &= C.dminU >= dminU && C.dmaxU <= dmaxU && C.dminV >= dminV && C.dmaxV <= dmaxV;
      if(nU > 0) matched &= C.nU==nU;
      if(nV > 0) matched &= C.nV==nV;

      return matched;
    } else return false;
  }
    
};



int main(int ac, char **av)
{
  size_t N           = ac>=2? strtol(av[1],0,0) : 100;
  size_t N_chunks    = ac>=3? strtol(av[2],0,0) : 1;
  size_t chunk_index = ac>=4? strtol(av[3],0,0) : 0;
 
  vector<general_spiral> isomers = FullereneDual::isomer_search(filter_nanocones(14*14,30,1,3*3,2*2,10*10,5,7), N,100000,false,false,N_chunks,chunk_index);

  
  
  return 0;
}
