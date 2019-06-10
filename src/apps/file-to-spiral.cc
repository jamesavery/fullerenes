#include <libgraph/polyhedron.hh>
#include <libgraph/triangulation.hh>

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

  Polyhedron P(av[1]);
  PlanarGraph dual =  static_cast<PlanarGraph>(P).dual_graph();
  dual.layout2d = dual.tutte_layout();
  
  cout << "dual = " << dual << endl;

  assert(dual.is_triangulation());

  Triangulation Tdual(dual);

  Triangulation::jumplist_t jumps;
  vector<int> spiral;
  Tdual.get_spiral(spiral, jumps, false);
  
  cout << "jumps = " << jumps << ";\n"
       << "spiral = " << spiral << ";\n";

  if(P.is_a_fullerene())
    cout << "RSPI = " << spiral_to_rspi(spiral) << ";\n";
  
  return 0;
}
