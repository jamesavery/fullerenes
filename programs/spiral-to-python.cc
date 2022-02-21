#include <string>

#include "fullerenes/triangulation.hh"
#include "fullerenes/spiral.hh"
#include "fullerenes/polyhedron.hh"
#include "fullerenes/symmetry.hh"

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

  if(ac>3){
    int k = strtol(av[2],0,0),
        l = strtol(av[3],0,0);

    t = t.GCtransform(k,l);
    general_spiral spiral = t.get_general_spiral();
    fsn.jumps  = spiral.jumps;	// TODO: general_spiral in spiral_nomenclature
    fsn.spiral_code = spiral.spiral;
  }
  
  
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

  vector<face_t> faces = t.cubic_faces(); 
  vector<int> degrees(Nf);

  
  for(int f=0,p=0,h=0;f<Nf;f++){
    int df = faces[f].size(), dd = t.neighbours[f].size();
    //    if(f<12) assert(dd==5); // NEJ! MUST BE IN SPIRAL ORDER, OR SYMMETRY CALC BREAKS

    if(df == 5) pentagons[p++] = faces[f];
    if(df == 6) hexagons[h++] = faces[f];    
  }

  // Dual node coordinates
  vector<coord3d> dual_coords(Nf);
  
  assert(t.triangles.size() == N);

  for(node_t u=0;u<N;u++){
    tri_t triangle = t.triangles[u];    
    coord3d triangle_center = P.points[u];
    
    for(int i=0;i<3;i++)
      dual_coords[triangle[i]] += triangle_center;
  }
  for(node_t f=0;f<Nf;f++)
    dual_coords[f] /= t.degree(f);

  Polyhedron T(t,dual_coords);
  Polyhedron::to_file(T,"dual.mol2");  
  
  //  Symmetry S(fsn.spiral_code);
  Symmetry S(t);
  
  cerr << "from numpy import array, nan\n\n";
  cerr << "name = \"" << fsn <<"\";\n";
  cerr << "# Symmetry information\n";
  cerr << "point_group = \"" << S.point_group() << "\";\n";
  cerr << "equivalent_nodes = " << S.equivalence_classes(S.Gtri) << ";\n";
  cerr << "equivalent_faces = " << S.equivalence_classes(S.G) << ";\n";
  //  cerr << "equivalent_edges = " << S.equivalence_classes(S.Gedge) << ";\n";    // TODO: Needs to be arcs in correct order
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
  cerr << "dual_neighbours   = " << t.neighbours << ";\n\n";
  //  cerr << "next_on_tri       = array(" << next_on_tri << ");\n\n";  
  cerr << "triangles         = array(" << t.triangles << ");\n\n";
  cerr << "dual_points     = array(" << dual_coords << ");\n\n";  
  cerr << "# prev_on_tri is the same as next_on_tri\n";


    
  return 0;
}
