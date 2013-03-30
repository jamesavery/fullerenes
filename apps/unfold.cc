#include <stdlib.h>
#include <fstream>
#include "libgraph/cubicgraph.hh"
#include "libgraph/fullerenegraph.hh"

#include <vector>

using namespace std;

class Eisenstein: public pair<int,int> {
public:
  Eisenstein(int a=0, int b=0) : pair<int,int>(a,b) {}
  Eisenstein operator*(const Eisenstein& y) const { return Eisenstein(first*y.first,second*y.second); }
  Eisenstein operator+(const Eisenstein& y) const { return Eisenstein(first+y.first,second+y.second); }
  Eisenstein operator-(const Eisenstein& y) const { return Eisenstein(first-y.first,second-y.second); } 
  Eisenstein& operator+=(const Eisenstein& y) { first += y.first; second += y.second; return *this; }
  Eisenstein& operator-=(const Eisenstein& y) { first -= y.first; second -= y.second; return *this; }

  //  
  // (-1,1)   \ /  (1,1)
  // (-1,0)  --x-- (1,0)
  // (-1,-1)  / \  (1,-1)
  // 
  static Eisenstein nextCW(const Eisenstein& x, const Eisenstein& y){
    Eisenstein d(y-x);
    switch(d.second + 10*d.first){
    case  10+0: /*( 1 ,0)*/  return Eisenstein(1,-1);
    case  10-1: /*( 1,-1)*/  return Eisenstein(0,-1); 
    case   0-1: /*( 0,-1)*/  return Eisenstein(-1,0);
    case -10+0: /*(-1, 0)*/  return Eisenstein(-1,1);
    case -10+1: /*(-1, 1)*/  return Eisenstein(0,1);
    case   0+1: /*( 0, 1)*/  return Eisenstein(1,0);
    default:
      cerr << "nextCW(): " << x << " does not neighbour " << y << "(difference: "<<(y-x)<<")\n";
      abort();
    }
  }

  coord2d coord() const { return coord2d(1,0)*first + coord2d(0.5,0.8660254037844386)*second; }
};


// Preconditions: Triangles are oriented consistently
map<Eisenstein,node_t> unfold(const vector<tri_t> &triangulation)
{
#define set_dedge(u,v,ux,vx) {	          \
  dedge_t uv(u,v), vu(v,u);               \
  printf("set_dedge(%d,%d)\n",u,v);       \
  dedge_done[uv] = true;                  \
  workset.erase(uv);                      \
  if(!dedge_done[vu]){                    \
    workset.insert(vu);			  \
    dedge_position[vu] = make_pair(vx,ux);\
  }                                       \
}

  map<dedge_t,node_t> nextNode;
  for(int i=0;i<triangulation.size();i++){
    const tri_t &t(triangulation[i]);
    for(int j=0;j<3;j++)
      nextNode[dedge_t(t[j],t[(j+1)%3])] = t[(j+2)%3];
  }

  map<dedge_t,bool> dedge_done;
  set<dedge_t>   workset;
  map<dedge_t, pair<Eisenstein,Eisenstein> > dedge_position;

  map<Eisenstein,node_t> grid;
  Eisenstein zero(0,0), veci(1,0), vecj(0,1);

  // 1. Place first triangle. 
  tri_t t(triangulation[0]);
  cout << "Place " << t[0] << " at " << zero << endl;
  cout << "Place " << t[1] << " at " << veci << endl;
  cout << "Place " << t[2] << " at " << (veci-vecj) << endl;
  grid[zero]      = t[0];
  grid[veci]      = t[1];
  grid[veci-vecj] = t[2];

  set_dedge(t[0],t[1],zero,veci);
  set_dedge(t[1],t[2],veci,veci-vecj);
  set_dedge(t[2],t[0],veci-vecj,zero);

  while(!workset.empty()){
    dedge_t uv(*workset.begin());
    cout << "Next unused dedge is " << uv << endl;
    // set_triangle(uv)
    node_t u(uv.first), v(uv.second), w(nextNode[uv]);

    pair<Eisenstein,Eisenstein> uvpos(dedge_position[uv]);
    Eisenstein ux(uvpos.first), vx(uvpos.second), wx(ux+Eisenstein::nextCW(ux,vx));

    cout << "Place " << w << " at " << wx << endl;
    grid[wx] = w;

    set_dedge(u,v,ux,vx);
    set_dedge(v,w,vx,wx);
    set_dedge(w,u,wx,ux);

  }
  return grid;
}


int main(int ac, char **av)
{
  vector<int> rspi(12);
  FullereneGraph::jumplist_t jumps;

  if(ac<13) return -1;

  int N = strtol(av[1],0,0);
  for(int i=0;i<12;i++) rspi[i] = strtol(av[2+i],0,0)-1;

  cout << "Attempting to create graph from spiral indices " << rspi << endl;
 

  FullereneGraph g(N, rspi, jumps);
  PlanarGraph dual(g.dual_graph(6));

  g.layout2d = g.tutte_layout();
  dual.layout2d = dual.tutte_layout();

  cout << "g = "  << g << ";\n";
  cout << "dg = " << dual << ";\n";

  vector<face_t> faces(dual.compute_faces_flat(3,true));
  vector<tri_t>  triangles(faces.begin(),faces.end());

  map<Eisenstein,node_t> grid(unfold(triangles));

  for(map<Eisenstein,node_t>::const_iterator i(grid.begin()); i!=grid.end();i++){
    Eisenstein x(i->first);
    node_t u(i->second);
    cout << u << " at " << x << endl;
  }
  
  return 0;
}
