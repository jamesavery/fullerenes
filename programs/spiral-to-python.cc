#include <string>

#include "fullerenes/triangulation.hh"
#include "fullerenes/spiral.hh"
#include "fullerenes/polyhedron.hh"

using namespace std;

int face_size(const Graph &g, node_t u, node_t v)
{
  int d = 1;
  node_t u0 = u;
  while(v != u0){
    node_t w = v;
    v = g.next_on_face(u,v);
    u = w;
    d++;
  }
  return d;
}


int main(int ac, char **av)
{
  string spiral_name = av[1];
  spiral_nomenclature fsn(spiral_name);
  Triangulation t(fsn);
  //  t = t.sort_nodes();
  
  PlanarGraph   g = t.dual_graph();
  g.layout2d = g.tutte_layout();
  Polyhedron P0(g,g.zero_order_geometry());
  Polyhedron P = P0;
  P.optimize();

  int N = g.N, Nf = t.N;
  
  vector<face_t> pentagons(12), hexagons(Nf-12);
  neighbours_t next_on_face(N,vector<node_t>(3));
  neighbours_t prev_on_face(N,vector<node_t>(3));
  neighbours_t face_right(N,vector<node_t>(3));

  
  for(node_t u=0;u<N;u++){
    for(int j=0;j<3;j++){
      node_t v = g.neighbours[u][j];
      next_on_face[u][j] = g.next_on_face(u,v);
      prev_on_face[u][j] = g.prev_on_face(u,v);
      face_right[u][j]  = face_size(g,u,v);
    }
  }

  vector<face_t> faces = g.compute_faces_oriented(6);
  vector<int> degrees(Nf);
  
  for(int f=0,p=0,h=0;f<Nf;f++){
    int df = faces[f].size(), dd = t.neighbours[f].size();
    // if(df != dd){
    //   printf("Face %d has size %ld, but dual node %d has degree %d\n",
    // 	     f,faces[f].size(),
    // 	     f,t.degree(f));
    //   abort();
    // }

    if(df == 5) pentagons[p++] = faces[f];
    if(df == 6) pentagons[h++] = faces[f];    
  }
  
  cerr << "from numpy import array, nan\n\n";
  cerr << "name = \"" << spiral_name <<"\";\n";

  cerr << "# Cubic graph, its faces, 3D embedding, and 2D Tutte-embedding\n";
  cerr << "cubic_neighbours  = array(" << g.neighbours << ");\n\n";
  cerr << "pentagons    = array(" << pentagons  << ");\n\n"; // TODO
  cerr << "hexagons     = array(" << hexagons   << ");\n\n"; // TODO
  cerr << "next_on_face = array(" << next_on_face << ");\n\n";
  cerr << "prev_on_face = array(" << prev_on_face << ");\n\n";
  cerr << "face_right   = array(" << face_right << ");\n\n"; 
  cerr << "points_start = array(" << P0.points << ");\n\n";
  cerr << "points_opt   = array(" << P.points << ");\n\n";   
  cerr << "tutte_layout = array(" << g.layout2d << ");\n\n";
  cerr << "# Dual graph and its faces\n";
  cerr << "dual_neighbours   = array(" << t.neighbours << ");\n\n";
  //  cerr << "next_on_tri       = array(" << next_on_tri << ");\n\n";  
  cerr << "triangles         = array(" << t.compute_faces() << ");\n\n";
  cerr << "# prev_on_tri is the same as next_on_tri\n";
  
    
  return 0;
}
