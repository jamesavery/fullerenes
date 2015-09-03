#include "fullerenegraph.hh"
#include "polyhedron.hh"
 // Exported interface
extern "C" {
  typedef FullereneGraph* fullerene_graph_ptr;
  typedef PlanarGraph* graph_ptr;
  typedef Polyhedron* polyhedron_ptr;

  // Graph construction and destruction
  fullerene_graph_ptr new_fullerene_graph_(const int *nmax, const int *N, const int *adjacency);
  fullerene_graph_ptr read_fullerene_graph_(const char *f_path);
  fullerene_graph_ptr read_fullerene_graph_hog_(const unsigned int *index, const char *f_path);
  fullerene_graph_ptr windup_general_(const int *n, const int indices[12], const int jumps_array[10]);
  void delete_fullerene_graph_(fullerene_graph_ptr*);

  polyhedron_ptr new_polyhedron_(const graph_ptr *g, const double *points);
  polyhedron_ptr read_polyhedron_(const char *path);
  polyhedron_ptr new_c20_();	
  void delete_polyhedron_(polyhedron_ptr*);

  graph_ptr new_graph_(const int *nmax, const int *N, const int *adjacency);
  void delete_graph_(graph_ptr*);

  // General graph operations -- look in fullerenegraph/graph.hh for others to potentially add
  int nvertices_(const graph_ptr *);
  int nedges_(const graph_ptr *);
  void edge_list_(const graph_ptr *, int *edges /* 2x|E| */, int *length);
  void adjacency_list_(const fullerene_graph_ptr *, const int *k, int *neighbours /* kxN */);
  void adjacency_matrix_(const graph_ptr *g, const int *outer_dim, int *adjacency);
  void get_arc_face_(const graph_ptr *g, const int *u, const int *v, int *face, int *l);

  int hamiltonian_count_(const graph_ptr *);
  int perfect_match_count_(const graph_ptr *);
  void all_pairs_shortest_path_(const graph_ptr *g, const int *max_depth, const int *outer_dim, int *D);
  void vertex_depth_(const graph_ptr *g, const int *outer_face, const int *of_length, 
		     int *depths, int *max_depth);
  int graph_is_a_fullerene_(const graph_ptr *);
  void print_graph_(const graph_ptr *);
  void draw_graph_(const graph_ptr *g, const char *filename, const char *format, const int *show_dual, const int *show_labels, const double *dimensions,
		   const int *edge_colour, const int *vertex_colour, const double *edge_width, const double *vertex_diameter);
  void draw_graph_with_path_(const graph_ptr *g, const char *filename, const char *format, const double *dimensions,
			     const int *edge_colour,const int *path_colour, const int *vertex_colour, const double *edge_width,
			     const double *path_width, const double *vertex_diameter,const int *Npath, int *path);
  graph_ptr dual_graph_(const graph_ptr *);

  void get_general_spiral_(const fullerene_graph_ptr*, int rspi_a[12], int jumps_a[10]);


  // Fullerene graph generation 
  fullerene_graph_ptr halma_fullerene_(const fullerene_graph_ptr *g, const int *n);
  fullerene_graph_ptr leapfrog_fullerene_(const fullerene_graph_ptr *g, const int *n_leaps);
  fullerene_graph_ptr goldberg_coxeter_(const fullerene_graph_ptr *g, const int *k, const int *l);

  // Fullerene graph operations
  void compute_fullerene_faces_(const fullerene_graph_ptr *, int *pentagons, int *hexagons);
  void get_pentagon_distance_mtx_(const fullerene_graph_ptr *, int *pentagon_distances);
  void get_face_distance_mtx_(const fullerene_graph_ptr *, int *face_distances);


  // Planar and spherical graph layout methods 
  void tutte_layout_(graph_ptr* g, double *LAYOUT);
  void tutte_layout_b_(graph_ptr* g, int *s, int *t, int *r, double *LAYOUT);
  void spherical_layout_(const graph_ptr* g, double *LAYOUT3D);
  void get_layout2d_(const graph_ptr *p, double *layout2d);
  void get_face_(const graph_ptr *g, const int *s, const int *t, const int *r, const int *fmax, int *face, int *l);

  void set_layout2d_(graph_ptr *g, const double *layout2d, const int *layout_is_spherical);

  double shortest_planar_distance_(const graph_ptr *g);

  // Operations on polyhedra
  double get_surface_area_(const polyhedron_ptr *P);
  double get_volume_(const polyhedron_ptr *P);
  // graph_ptr triangulation_(polyhedron_ptr *P)
  polyhedron_ptr convex_hull_(const polyhedron_ptr *P);

  void draw_polyhedron_(const polyhedron_ptr *P, const char *filename_, const char *format, 
			const int *edge_colour, const int *node_colour, const int *face_colour,
			const double *edge_width, const double *node_diameter, const double *face_opacity);

  void get_layout3d_(const polyhedron_ptr *p, double *layout3d);
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

string fortran_string(const char *s, int max)
{
  int i;
  for(i=0;i<max && s[i] != ' '; i++) ;
  return string(s,i);
}

fullerene_graph_ptr read_fullerene_graph_(const char *f_path){
  fullerene_graph_ptr g;

  string path(fortran_string(f_path,50));

  FILE *f = fopen(path.c_str(),"r");
  if(!f){
    fprintf(stderr,"Cannot open file %s for reading: ",path.c_str());
    perror(path.c_str());
    return NULL;
  }

  g = new FullereneGraph(f);
  fclose(f);
  return g;
}

fullerene_graph_ptr read_fullerene_graph_hog_(const unsigned int *index, const char *f_path){
  fullerene_graph_ptr g;

// strip the whitespace (fortran passes an array of 50 characters)
  string path(fortran_string(f_path,50));

  FILE *f = fopen(path.c_str(),"r");
  if(!f){
    fprintf(stderr,"Cannot open file %s for reading: ",path.c_str());
    perror(path.c_str());
    return NULL;
  }

  g = new FullereneGraph(*index,f);
  fclose(f);
  return g;
}

fullerene_graph_ptr windup_general_(const int *n, const int spiral_indices_array[12], const int jumps_array[10]){
  
  vector<int> spiral_indices(12);
  for(int i=0; i<12; ++i){
    spiral_indices[i] = spiral_indices_array[i]-1;
  //  cout << spiral_indices[i] << endl;
  }

  list<pair<int,int> > jumps;
  for(int i=0; i<10; ++i,++i){
    if(jumps_array[i]==0) break;
    jumps.push_back(make_pair(jumps_array[i]-1,jumps_array[i+1]));
  //  cout << jumps.back().first << ", " << jumps.back().second << endl;
  }

  return new FullereneGraph(*n, spiral_indices, jumps);
}


void delete_fullerene_graph_(fullerene_graph_ptr *g){  delete (*g); }

polyhedron_ptr new_polyhedron_(const graph_ptr *g, const double *points)
{
  PlanarGraph G(*(*g));
  vector<coord3d> vertex_points(G.N);
  for(size_t i=0;i<G.N;i++) vertex_points[i] = coord3d(&points[3*i]);
  return new Polyhedron(G,vertex_points,6); 
}

polyhedron_ptr read_polyhedron_(const char *path) { return new Polyhedron(path); }

void delete_polyhedron_(polyhedron_ptr *P){ delete *P; }

// The Fortran call to dual_graph() assumes that g is either a fullerene
// or a fullerene dual.
graph_ptr dual_graph_(const graph_ptr *g){  
  (*g)->layout2d = (*g)->tutte_layout();
  bool is_fullerene = (*g)->neighbours[0].size() == 3;
  graph_ptr pg = new PlanarGraph((*g)->dual_graph(is_fullerene? 6 : 3));
  pg->layout2d = pg->tutte_layout();
  return pg;
}
//int hamiltonian_count_(const graph_ptr *g){ return (*g)->hamiltonian_count(); }
int perfect_match_count_(const graph_ptr *g){ return (*g)->count_perfect_matchings(); }

fullerene_graph_ptr halma_fullerene_(const fullerene_graph_ptr *g, const int *n)
{
  bool has_layout = (*g)->layout2d.size() == (*g)->N;
  return new FullereneGraph((*g)->halma_fullerene(*n,has_layout));
}

fullerene_graph_ptr leapfrog_fullerene_(const fullerene_graph_ptr *g, const int *n_leaps)
{
  FullereneGraph frog((*g)->leapfrog_fullerene(false));
  for(int i=1;i<*n_leaps;i++) 
    frog = frog.leapfrog_fullerene(false);

  return new FullereneGraph(frog);
}

fullerene_graph_ptr goldberg_coxeter_(const fullerene_graph_ptr *g, const int *k, const int *l)
{
  FullereneGraph fg(**g);
  fg.layout2d = fg.tutte_layout(); // FIXME remove, and pass argument to GCtransform?
  fg.layout_is_spherical = false;
  return new FullereneGraph(fg.GCtransform(*k,*l));
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
  
  vector<int> distances((*g)->all_pairs_shortest_paths(*max_depth));
  for(unsigned int i=0;i<M;i++)
    for(unsigned int j=0;j<M;j++)
      D[i*N+j] = distances[i*M+j];
}

double get_surface_area_(const polyhedron_ptr *P){ return (*P)->surface_area(); }
double get_volume_(const polyhedron_ptr *P){ return (*P)->volume(); }

int graph_is_a_fullerene_(const graph_ptr *g){ return (*g)->is_a_fullerene(); }

double shortest_planar_distance_(const graph_ptr *g)
{
  const PlanarGraph& G(*(*g));

  const vector<double> lengths(G.edge_lengths());
  
  return *min_element(lengths.begin(),lengths.end());
}

polyhedron_ptr convex_hull_(const polyhedron_ptr *p)
{
  return new Polyhedron((*p)->convex_hull());
}


void print_graph_(const graph_ptr *g)
{
  const PlanarGraph& G(*(*g));
  cout << G << ";\n";
}

void get_layout3d_(const polyhedron_ptr *p, double *points)
{
  const Polyhedron& P(*(*p));
  for(node_t u=0;u<P.N;u++){
    const coord3d &x(P.points[u]);
    for(int i=0;i<3;i++)
      points[u*3+i] = x[i];
  }
}

void get_layout2d_(const graph_ptr *g, double *points)
{
  const PlanarGraph& G(*(*g));
  for(node_t u=0;u<G.N;u++){
    const coord2d &x(G.layout2d[u]);
      points[u*2+0] = x.first;
      points[u*2+1] = x.second;
  }
}

void set_layout2d_(graph_ptr *g, const double *layout2d, const int *layout_is_spherical)
{
  PlanarGraph& G(*(*g));

  G.layout2d.resize(G.N);

  for(node_t u=0;u<G.N;u++){
    G.layout2d[u] = coord2d(layout2d[u*2],layout2d[u*2+1]);
  }
  G.orient_neighbours();
  G.layout_is_spherical = *layout_is_spherical;
}





polyhedron_ptr new_c20_(){  return new Polyhedron(Polyhedron::C20()); }

int nvertices_(const graph_ptr *g){ return (*g)->N; }
int nedges_(const graph_ptr *g){ 
  set<edge_t> edge_set = (*g)->undirected_edges(); // TODO: Better solution to this.
  return edge_set.size(); 
}

void draw_graph_(const graph_ptr *g, const char *filename_, const char *format, const int *show_dual, const int *show_labels, const double *dimensions,
		 const int *line_colour, const int *vertex_colour, const double *line_width, const double *vertex_diameter)
{
  string fmt(format,3), filename;
  filename = fortran_string(filename_,50)+"-2D."+fmt;

//   printf("draw_graph({\n"
//   	 "\tformat:     '%3s',\n"
//   	 "\tfilename:   '%s',\n"
//   	 "\tshow_dual:  %d,\n"
//   	 "\tdimensions: %gcm x %gcm,\n"
//   	 "\tline_colour: #%.6x,\n"
//   	 "\tnode_colour: #%.6x,\n"
//   	 "\tline_width: %.2gmm,\n"
//   	 "\tnode_diameter: %.2gmm,\n"
//   	 "})\n\n",
//   	 format, filename.c_str(), *show_dual, dimensions[0], dimensions[1],
//   	 *line_colour, *vertex_colour,
//   	 *line_width, *vertex_diameter);
//
//   cout << endl;
//   cout << **g << endl;

  ofstream graph_file(filename.c_str(),ios::out | ios::binary);
  if        (fmt == "tex"){
    graph_file <<  (*g)->to_latex(dimensions[0],dimensions[1],*show_dual,*show_labels,true,*line_colour,0,*vertex_colour,*line_width,0,*vertex_diameter,0,0);
  } else if (fmt == "pov"){
    graph_file << (*g)->to_povray(dimensions[0],dimensions[1],*line_colour,*vertex_colour,*line_width,*vertex_diameter);
  }
//  cout << "write done" << endl;
  graph_file.close();
}

void draw_graph_with_path_(const graph_ptr *g, const char *filename_, const char *format, const double *dimensions,
			   const int *edge_colour,const int *path_colour, const int *vertex_colour, const double *edge_width,
			   const double *path_width, const double *vertex_diameter, const int *Npath, int *path_)
{
  string fmt(format,3), filename;
  filename = fortran_string(filename_,50)+"-2D."+fmt;

  int path[*Npath];// Change from Fortran 1-indexing to C/C++ 0-indexing
  for(int i=0;i<*Npath;i++) path[i] = path_[i]-1;

  ofstream graph_file(filename.c_str(),ios::out | ios::binary);
  if        (fmt == "tex"){
    graph_file <<  (*g)->to_latex(dimensions[0],dimensions[1],false,true,true,*edge_colour,*path_colour,
				  *vertex_colour,*edge_width,*path_width,*vertex_diameter,*Npath,path);
  } 
  graph_file.close();  
}
void draw_polyhedron_(const polyhedron_ptr *P, const char *filename_, const char *format, 
		      const int *edge_colour, const int *node_colour, const int *face_colour,
		      const double *edge_width, const double *node_diameter, const double *face_opacity)
{
  string fmt(format,3), filename;
  filename = fortran_string(filename_,50)+"."+fmt;

  ofstream polyhedron_file(filename.c_str(),ios::out|ios::binary);

  if     (fmt == "pov"){
    polyhedron_file << (*P)->to_povray(10.0,10.0,*edge_colour,*node_colour,*face_colour,*edge_width,*node_diameter,*face_opacity);
  }
}

void edge_list_(const graph_ptr *g, int *edges, int *length)
{
  const set<edge_t>& edge_set((*g)->undirected_edges());
  *length = edge_set.size();
  int i=0;
  for(set<edge_t>::const_iterator e(edge_set.begin());e!=edge_set.end();e++,i++){
    edges[2*i+0] = e->first;
    edges[2*i+1] = e->second;
  }
}


// Assumes k-regular graph, since fortran handles non-flat data structures poorly.
void adjacency_list_(const fullerene_graph_ptr *g, const int *k, int *neighbours)
{
  const neighbours_t& ns((*g)->neighbours);

  for(node_t u=0;u<ns.size();u++)
    for(int i=0;i<*k;i++)
      neighbours[(*k)*u+i] = ns[u][i]+1;
}

void compute_fullerene_faces_(const fullerene_graph_ptr *g, int *pentagons, int *hexagons)
{
  assert((*g)->layout2d.size() == (*g)->N);

  facemap_t faces((*g)->compute_faces_oriented());
  int NH = (*g)->N/2-10;
  
  assert(faces[5].size() == 12);
  assert(faces[6].size() == NH);

  const vector<face_t> p(faces[5].begin(),faces[5].end()), 
                       h(faces[6].begin(),faces[6].end());

  for(int i=0;i<12;i++)
    for(int j=0;j<5;j++)
      pentagons[i*5+j] = p[i][j];

  for(int i=0;i<NH;i++)
    for(int j=0;j<6;j++)
      hexagons[i*6+j] = h[i][j];
}

void get_pentagon_distance_mtx_(const fullerene_graph_ptr *fg, int *pentagon_distances){
  vector<int> mtx_v = (*fg)->pentagon_distance_mtx();
//  assert(mtx_v.size()==144);
  std::copy(mtx_v.begin(), mtx_v.end(), pentagon_distances);
}

void get_face_distance_mtx_(const fullerene_graph_ptr *fg, int *face_distances){

  PlanarGraph dual((*fg)->dual_graph(6));
  
  vector<int> all_distances = dual.all_pairs_shortest_paths(INT_MAX);

  int i=0;
  for (vector<int>::iterator it=all_distances.begin();it!= all_distances.end(); ++it, ++i){
    face_distances[i] = all_distances[i];
  }
}

// rspi_a and jumps_a start counting at 1
void get_general_spiral_(const fullerene_graph_ptr* fg, int rspi_a[12], int jumps_a[10]){
//  12 will always be 12, 10 is just an arbitrary magic number
  assert((*fg)->layout2d.size() == (*fg)->N);
  vector<int> rspi_v;
  FullereneGraph::jumplist_t jumps_v;
  const bool canonical=true, general=true;
  (*fg)->get_rspi_from_fg(rspi_v, jumps_v, canonical, general);

  for(int i=0; i!=12; i++){
    rspi_a[i] = rspi_v[i] +1;//start counting at 1
  }
  int j=0;
  std::fill(jumps_a,jumps_a+10, 0);
  for(list<pair<int,int> >::iterator it(jumps_v.begin()); it!=jumps_v.end(); it++){
  	jumps_a[j++]=it->first +1;//start counting at 1
  	jumps_a[j++]=it->second;
  }
}


// Only works for trivalent graphs. 
void get_face_(const graph_ptr *g, const int *s, const int *t, const int *r, const int *fmax, int *face, int *l)
{
  face_t cycle((*g)->shortest_cycle(*s-1,*t-1,*r-1,*fmax));
  for(int i=0;i<cycle.size();i++) face[i] = cycle[i]+1;
  //  cerr << "get_face(" << *s << "," << *t << "," << *r << ") = " << cycle << endl; 
  *l = cycle.size();
}

// This version works for all planar graphs.
void get_arc_face_(const graph_ptr *g, const int *u, const int *v, int *face, int *l)
{
  assert((*g)->layout2d.size() == (*g)->N);
  face_t f((*g)->get_face_oriented(*u,*v));
  
  for(int i=0;i<f.size();i++) face[i] = f[i];
  *l = f.size();
}

void vertex_depth_(const graph_ptr *g, const int *outer_face, const int *of_length, int *depths, int *max_depth)
{
  const int N = (*g)->N;

  
  vector<node_t> outer(*of_length);
  for(int i=0;i<outer.size();i++) outer[i] = outer_face[i]-1;

  vector<int> D((*g)->multiple_source_shortest_paths(outer,vector<bool>(N*(N-1)/2),vector<bool>(N)));
  
  *max_depth = -1;
  for(node_t u=0;u<N;u++){
    depths[u] = D[u];
    if(depths[u] > *max_depth) *max_depth = depths[u];
  }
}
