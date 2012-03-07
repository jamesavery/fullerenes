#include "fullerenegraph.hh"
 // Exported interface
extern "C" {
  typedef FullereneGraph* fullerene_graph_ptr;
  typedef Graph* graph_ptr;

  // Graph construction and destruction
  fullerene_graph_ptr new_fullerene_graph_(const int *nmax, const int *N, const int *adjacency);
  void delete_fullerene_graph_(fullerene_graph_ptr*);

  graph_ptr new_graph_(const int *nmax, const int *N, const int *adjacency);
  void delete_graph_(graph_ptr*);

  // General graph operations -- look in fullerenegraph/graph.hh for others to potentially add
  graph_ptr dual_graph_(const graph_ptr *);
  int hamiltonian_count_(const graph_ptr *);
  
  // Fullerene graph generation 
  fullerene_graph_ptr halma_fullerene_(const fullerene_graph_ptr *g, const unsigned int n);
  fullerene_graph_ptr leapfrog_fullerene_(const fullerene_graph_ptr *g, const bool);

  // Planar and spherical graph layout methods 
  void tutte_layout_(fullerene_graph_ptr* g, double *LAYOUT);
  void spherical_layout_(const fullerene_graph_ptr* g, double *LAYOUT3D);


};


//------------------ Wrapper functions --------------------------------
#include <math.h>
#include <string.h>



graph_ptr new_graph_(const int *nmax, const int *n, const int *adjacency){
  const int Nmax(*nmax), N(*n);
  vector<int> A(N*N);

  for(unsigned int i=0;i<N;i++)
    for(int j=0;j<N;j++)
      A[i*N+j] = adjacency[i*Nmax+j];

  return new Graph(N,A);
}

void delete_graph_(graph_ptr *g){  delete (*g); }

fullerene_graph_ptr new_fullerene_graph_(const int *nmax, const int *n, const int *adjacency){
  const int Nmax(*nmax), N(*n);
  vector<int> A(N*N);

  for(unsigned int i=0;i<N;i++)
    for(int j=0;j<N;j++)
      A[i*N+j] = adjacency[i*Nmax+j];

  return new FullereneGraph(Graph(N,A));
}

void delete_fullerene_graph_(fullerene_graph_ptr *g){  delete (*g); }

graph_ptr dual_graph_(const graph_ptr *g){  return new Graph((*g)->dual_graph()); }
int hamiltonian_count_(const graph_ptr *g){ return (*g)->hamiltonian_count(); }

fullerene_graph_ptr halma_fullerene_(const fullerene_graph_ptr *g, const unsigned int n)
{
  return new FullereneGraph((*g)->halma_fullerene(n,false));
}

fullerene_graph_ptr leapfrog_fullerene_(const fullerene_graph_ptr *g)
{
  return new FullereneGraph((*g)->leapfrog_fullerene(false));
}

void tutte_layout_(fullerene_graph_ptr *g, double *LAYOUT)
{
  FullereneGraph &G(**g);
  G.layout2d = G.tutte_layout();
  for(unsigned int i=0;i<G.N;i++){
    LAYOUT[i*2]   = G.layout2d[i].first;
    LAYOUT[i*2+1] = G.layout2d[i].second;
  }
}

void spherical_layout_(const fullerene_graph_ptr *g, double *LAYOUT3D)
{
  const FullereneGraph &G(**g);
  vector<coord2d> layout2d(G.tutte_layout());
  vector<coord2d> angles(G.spherical_projection(layout2d));
  
  for(int i=0;i<G.N;i++){
    double theta = angles[i].first, phi = angles[i].second;
    fprintf(stderr,"%f,%f\n",theta,phi);
    double x = cos(theta)*sin(phi), y = sin(theta)*sin(phi), z = cos(phi);
    LAYOUT3D[i*3] = x;
    LAYOUT3D[i*3+1] = y;
    LAYOUT3D[i*3+2] = z;
  }
}
