#include <fullerenes/polyhedron.hh>
#include <fullerenes/triangulation.hh>

vector<int> spiral_to_rspi(const vector<int>& spiral)
{
  vector<int> RSPI;
  for(int i=0;i<spiral.size();i++)
    if(spiral[i] == 5) RSPI.push_back(i+1);
  return RSPI;
}

int main(int ac, char **av)
{
  assert(ac>=2);

  Polyhedron P = Polyhedron::from_file(av[1]);
  PlanarGraph dual =  static_cast<PlanarGraph>(P).dual_graph();
  assert(dual.is_triangulation());
  Triangulation Tdual(dual);

  jumplist_t jumps;
  vector<int> spiral;
  Tdual.get_spiral(spiral, jumps);

  vector<face_t>
    pentagons(vector<face_t>(12)),
    hexagons (vector<face_t>(P.faces.size()-12));

  for(int j=0,npent=0,nhex=0;j<P.faces.size();j++)
    if      (P.faces[j].size() == 5) pentagons[npent++] = P.faces[j];
    else if (P.faces[j].size() == 6) hexagons [nhex++]  = P.faces[j];  

  cerr << "from numpy import array\n\n";
  cerr << "N = " << P.N << ";\n";
  cerr << "Nf = " << P.faces.size() << ";\n\n";  
  cerr << "spiral       = array(" << spiral << ")\n\n";
  cerr << "spiral_jumps = array(" << jumps << ");\n\n";
  cerr << "neighbours   = array(" << P.neighbours << ").reshape(1,N,3);\n\n";
  cerr << "pentagons    = array(" << pentagons << ").reshape(1,12,5);\n\n";      
  cerr << "hexagons     = array(" << hexagons  << ").reshape(1,Nf-12,6);\n\n";
  cerr << "points_opt   = array(" << P.points     << ").reshape(1,N,3);\n\n";
  
  
  return 0;
}
