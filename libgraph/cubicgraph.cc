#include "cubicgraph.hh"



ostream& operator<<(ostream& s, const CubicGraph& g)
{
  s.precision(30);
  s << fixed << endl;
  s << "{" << static_cast<PlanarGraph>(g); 
  
  if(g.spherical_layout.size() == g.N){
    s << g.name << ", {";
    for(unsigned int i=0;i<g.N;i++){
      coord2d angles(g.spherical_layout[i]);
      s << angles << (i+1<g.N? ", " : "}");
    }
  }
  s << "}";
  return s;
}
