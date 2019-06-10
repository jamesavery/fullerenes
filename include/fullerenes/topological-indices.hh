#ifndef TOPOLOGICAL_INDICES_HH
# define TOPOLOGICAL_INDICES_HH

#include <math.h>

class TopologicalIndices {
public:
  const FullereneGraph& g;
  const vector<int> D;
  const int N;
  vector<int> dmax;
  int diameter, radius;

  TopologicalIndices(const FullereneGraph& g) : g(g), D(g.all_pairs_shortest_paths()), N(g.N),
						dmax(int(N))
  { 
    for(int u=0;u<N;u++)
      for(int v=0;v<N;v++) if(dmax[u] < D[u*N+v]) dmax[u] = D[u*N+v];
    
    diameter = 0, radius = INT_MAX;
    for(int u=0;u<N;u++){
      if(diameter < dmax[u]) diameter = dmax[u];
      if(radius   > dmax[u]) radius   = dmax[u];
    }
  }
  
  int Wiener() const {
    int w = 0;
    vector<int> wi(N);

    for(int u=0;u<N;u++)
      for(int v=u+1;v<N;v++) w += D[u*N+v];

    return w;
  }

  int hyperWiener() const {
    int ww = 0;
    for(int u=0;u<N;u++)
      for(int v=u+1;v<N;v++) ww += D[u*N+v] + D[u*N+v]*D[u*N+v];

    return ww;
  }

  int reverseWiener() const {
    return diameter*N*(N-1)/2 - Wiener();
  }

  int Szeged() const {
    int S = 0;

    for(set<edge_t>::const_iterator e(g.edge_set.begin()); e!=g.edge_set.end();e++){
      int i=e->first, j = e->second, ni = 0, nj = 0;

      for(int k = 0;k<N;k++){
	if(D[i*N+k] < D[j*N+k]) ni++;
	if(D[i*N+k] > D[j*N+k]) nj++;
      }
      S += ni*nj;
    }
    return S;
  }

  double Balaban() const {
    vector<int> Wi(N);
    // Wi[u] = \sum_v dist(u,v)
    for(int u=0;u<N;u++)
      for(int v=0;v<N;v++) Wi[u] += D[u*N+v];
    
    double B = 0;
    for(set<edge_t>::const_iterator e(g.edge_set.begin()); e!=g.edge_set.end();e++){
      int u = e->first, v = e->second;
      B += 1/sqrt(Wi[u]*Wi[v]);
    }

    return B*N/(N/2+2.0);
  }

  double Estrada() const {
    double E = 0;

    // Needs eigenvalues -- should diagonalization be done at
    // initialization, or only if Estrada or bipartivity is needed?
    //    for(int u=0;u<N;u++) E += exp(lambda[u]);
    return E;
  }
  
  double bipartivity() const {
    double E = Estrada(), B = 0;
    // Needs eigenvalues
    //    for(int u=0;u<N;u++) B += cosh(lambda[u]);
    return B/E;
  }

  bool check_all() const {
    // TODO: check the various relations and invariants hold
    return true;
  }
};

#endif
