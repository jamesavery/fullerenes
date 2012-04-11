#include "fullerenegraph.hh"
#include "polyhedron.hh"
 // Exported interface
extern "C" {
  typedef FullereneGraph* fullerene_graph_ptr;
  typedef PlanarGraph* graph_ptr;
  typedef Polyhedron* polyhedron_ptr;

  // Graph construction and destruction
  fullerene_graph_ptr new_fullerene_graph_(const int *nmax, const int *N, const int *adjacency);
  fullerene_graph_ptr read_fullerene_graph_(const char *path);
  void delete_fullerene_graph_(fullerene_graph_ptr*);

  polyhedron_ptr new_polyhedron_(const graph_ptr *g, const double *points);
  polyhedron_ptr read_polyhedron_(const char *path);
  void delete_polyhedron_(polyhedron_ptr*);

  graph_ptr new_graph_(const int *nmax, const int *N, const int *adjacency);
  void delete_graph_(graph_ptr*);

  // General graph operations -- look in fullerenegraph/graph.hh for others to potentially add
  graph_ptr dual_graph_(const graph_ptr *);
  int hamiltonian_count_(const graph_ptr *);
  void all_pairs_shortest_path_(const graph_ptr *g, const int *max_depth, const int *outer_dim, int *D);
  void adjacency_matrix_(const graph_ptr *g, const int *outer_dim, int *adjacency);
  int graph_is_a_fullerene_(const graph_ptr *);

  void print_graph_(const char *name, const graph_ptr *);

  // Fullerene graph generation 
  fullerene_graph_ptr halma_fullerene_(const fullerene_graph_ptr *g, const unsigned int *n);
  fullerene_graph_ptr leapfrog_fullerene_(const fullerene_graph_ptr *g, const int *n_leaps);


  // Planar and spherical graph layout methods 
  void tutte_layout_(graph_ptr* g, double *LAYOUT);
  void tutte_layout_b_(graph_ptr* g, int *s, int *t, int *r, double *LAYOUT);
  void spherical_layout_(const graph_ptr* g, double *LAYOUT3D);

  //  double *get_layout2d(const graph_ptr *g);
  //  void    set_layout2d(const graph_ptr *g, double *layout2d);

  double shortest_planar_distance_(const graph_ptr *g);

  // Operations on polyhedra
  double surface_area_(polyhedron_ptr *P);
  double volume2_(polyhedron_ptr *P);
  // graph_ptr triangulation_(polyhedron_ptr *P)
  // polyhedron_ptr convex_hull_(polyhedron_ptr *P)
};


//------------------ Wrapper functions --------------------------------
#include <math.h>
#include <string.h>
#include <stdio.h>

void adjacency_matrix_(const graph_ptr *g, const int *outer_dim, int *adjacency)
{
  Graph G(*(*g));
  unsigned int N = *outer_dim;

  memset(adjacency,0,N*N*sizeof(int));

  for(node_t u=0;u<G.N;u++){
    const vector<node_t> &ns(G.neighbours[u]);
    for(unsigned int i=0;i<ns.size();i++)
      adjacency[u*N+ns[i]] = 1;
  }
}

graph_ptr new_graph_(const int *nmax, const int *n, const int *adjacency){
  const int Nmax(*nmax), N(*n);
  set<edge_t> edge_set;

  for(unsigned int i=0;i<N;i++)
    for(int j=i+1;j<N;j++)
      if(adjacency[i*Nmax+j])  edge_set.insert(edge_t(i,j));

  return new PlanarGraph(edge_set);
}

void delete_graph_(graph_ptr *g){  delete (*g); }

fullerene_graph_ptr new_fullerene_graph_(const int *nmax, const int *n, const int *adjacency){
  const int Nmax(*nmax), N(*n);
  set<edge_t> edge_set;

  for(node_t i=0;i<N;i++)
    for(node_t j=i+1;j<N;j++)
      if(adjacency[i*Nmax+j]) edge_set.insert(edge_t(i,j));

  return new FullereneGraph(Graph(edge_set));
}

fullerene_graph_ptr read_fullerene_graph_(const char *path){
  fullerene_graph_ptr g;
  FILE *f = fopen(path,"r");
  if(!f){
    fprintf(stderr,"Cannot open file %s for reading: ",path);
    perror(path);
    return NULL;
  }
  g = new FullereneGraph(f);
  fclose(f);
  return g;
}

void delete_fullerene_graph_(fullerene_graph_ptr *g){  delete (*g); }

polyhedron_ptr new_polyhedron_(const graph_ptr *g, const double *points)
{
  PlanarGraph G(*(*g));
  vector<coord3d> vertex_points(G.N);
  for(size_t i=0;i<G.N;i++) vertex_points[i] = coord3d(&points[3*i]);
  return new Polyhedron(G,vertex_points,6); 
}
polyhedron_ptr read_polyhedron_(const char *path)
{
  polyhedron_ptr P = new Polyhedron(6);
  ifstream file(path);
  file >> *P;
  file.close();
  return P;
}
void delete_polyhedron_(polyhedron_ptr *P){ delete *P; }

graph_ptr dual_graph_(const graph_ptr *g){  return new PlanarGraph((*g)->dual_graph()); }
int hamiltonian_count_(const graph_ptr *g){ return (*g)->hamiltonian_count(); }

fullerene_graph_ptr halma_fullerene_(const fullerene_graph_ptr *g, const unsigned int *n)
{
  return new FullereneGraph((*g)->halma_fullerene(*n,false));
}

fullerene_graph_ptr leapfrog_fullerene_(const fullerene_graph_ptr *g, const int *n_leaps)
{
  FullereneGraph frog((*g)->leapfrog_fullerene(false));
  for(int i=1;i<*n_leaps;i++) 
    frog = frog.leapfrog_fullerene(false);

  return new FullereneGraph(frog);
}

void tutte_layout_(graph_ptr *g, double *LAYOUT)
{
  int null = 0, one = 1;
  //  fprintf(stderr,"tutte_layout_a()\n");
  tutte_layout_b_(g,&one,&null,&null,LAYOUT);
}
#include <fstream>
void tutte_layout_b_(graph_ptr *g, int *s, int *t, int *r, double *LAYOUT)
{
  //  fprintf(stderr,"tutte_layout_b(%d,%d,%d)\n",*s,*t,*r);
  PlanarGraph &G(**g);
  G.layout2d = G.tutte_layout(*s-1,*t-1,*r-1);
  for(unsigned int i=0;i<G.N;i++){
    LAYOUT[i*2]   = G.layout2d[i].first;
    LAYOUT[i*2+1] = G.layout2d[i].second;
  }
  //  ofstream f("output/tutte.m");
  //  f << "g = " << G << ";\n";
  //  f.close();
}

void spherical_layout_(const graph_ptr *g, double *LAYOUT3D)
{
  PlanarGraph G(**g,(*g)->tutte_layout()); // TODO: Reducer antallet af gange, tl kaldes.
  vector<coord2d> angles(G.spherical_projection());
  
  for(int i=0;i<G.N;i++){
    double theta = angles[i].first, phi = angles[i].second;
    //    fprintf(stderr,"%f,%f\n",theta,phi);
    double x = cos(theta)*sin(phi), y = sin(theta)*sin(phi), z = cos(phi);
    LAYOUT3D[i*3] = x;
    LAYOUT3D[i*3+1] = y;
    LAYOUT3D[i*3+2] = z;
  }
}

void all_pairs_shortest_path_(const graph_ptr *g, const int *max_depth, const int *outer_dim, int *D)
{
  unsigned int M=(*g)->N, N = *outer_dim;
  
  vector<unsigned int> distances((*g)->all_pairs_shortest_paths(*max_depth));
  for(unsigned int i=0;i<M;i++)
    for(unsigned int j=0;j<M;j++)
      D[i*N+j] = distances[i*M+j];
}

double surface_area_(polyhedron_ptr *P){ return (*P)->surface_area(); }
double volume2_(polyhedron_ptr *P){ return (*P)->volume(); }

int graph_is_a_fullerene_(const graph_ptr *g){ return (*g)->this_is_a_fullerene(); }

double shortest_planar_distance_(const graph_ptr *g)
{
  const PlanarGraph& G(*(*g));

  const vector<double> lengths(G.edge_lengths());
  
  return *min_element(lengths.begin(),lengths.end());
}

void print_graph_(const char *name, const graph_ptr *g)
{
  const PlanarGraph& G(*(*g));
  cout << name << " = " << G << ";\n";
}
