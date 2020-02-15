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
  Polyhedron::to_file(P,av[2]);
  
  return 0;
}
