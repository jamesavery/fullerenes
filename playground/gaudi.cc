#include <stdlib.h>
#include <fstream>
#include <vector>

#include "libgraph/cubicgraph.hh"
#include "libgraph/fullerenegraph.hh"
#include "libgraph/eisenstein.hh"
#include "libgraph/unfold.hh"
#include "libgraph/triangulation.hh"
#include "libgraph/polyhedron.hh"

using namespace std;

PlanarGraph fold(vector< pair<Eisenstein, node_t> > &outline);

typedef pair<Eisenstein,Eisenstein> dedgecoord_t;




int turn_direction(const Eisenstein& xi,const Eisenstein& xj,const Eisenstein& xk) 
{
  Eisenstein dx1(xj-xi), dx2(xk-xj);
  return sgn(dx2.first * dx1.second - dx2.second * dx1.first);
}


Eisenstein tfm(const Eisenstein& x, const Eisenstein& x0, const Eisenstein& w, const Eisenstein& x0p)
{
  return (x-x0)*w + x0p;
}


PlanarGraph GCTransform(const PlanarGraph& dual, int K=1, int L=0)
{
  Unfolding U(dual);
  Folding F(U*Eisenstein(K,L));
  return F.fold();
}


vector< pair<Eisenstein, node_t> > GCDreduce(const vector< pair<Eisenstein, node_t> > &outline)
{
  vector<Eisenstein> segments(outline.size());

  for(int i=0;i<outline.size();i++) segments[i] = outline[(i+1)%outline.size()].first - outline[i].first;

  cout << "segments  = " << segments << ";\n";

  Eisenstein d(Eisenstein::gcd(segments)); // TODO: Only do GCD between pentagon nodes.

  cout << "GCD = " << d << endl;
  for(int i=0;i<segments.size();i++) segments[i] = segments[i].div(d);

  vector< pair<Eisenstein,node_t> > new_outline(outline);
  for(int i=0;i+1<outline.size();i++) new_outline[i+1].first = new_outline[i].first+segments[i];

  return new_outline;
}


// the cube
Graph cube()
{
  const int N = 8;
  neighbours_t neighbours(N,vector<node_t>(3));

  for(int i=0; i<4; i++){
    neighbours[i][0] = (i+1)%4;
    neighbours[i][1] = (i-1+4)%4;
    neighbours[i][2] = (i+4)%4 + 4;

    neighbours[i+4][0] = (i+1)%4 + 4;
    neighbours[i+4][1] = (i-1+4)%4 + 4;
    neighbours[i+4][2] = (i+4)%4;
  }
  cout << neighbours << endl;
  return Graph(neighbours);
}


// smallest polyhedron with only pentagons and heptagons
Graph example1(){
  const int M=7, N=28;
  neighbours_t neighbours(N,vector<node_t>(3));

  for(int i=0; i!=7; ++i){
    neighbours[i][0] = (i+1)%M;
    neighbours[i][1] = (i-1+M)%M;
    neighbours[i][2] = M+i;
    
    neighbours[1*M+i][0] =     i;
    neighbours[1*M+i][1] = 2*M+i;
    neighbours[1*M+i][2] = 2*M+(i-1+M)%M;
    
    neighbours[2*M+i][0] = 1*M+i%M;
    neighbours[2*M+i][1] = 1*M+(i+1)%M;
    neighbours[2*M+i][2] = 3*M+i;
    
    neighbours[3*M+i][0] = 2*M+i;
    neighbours[3*M+i][1] = 3*M+(i+1)%M;
    neighbours[3*M+i][2] = 3*M+(i-1+M)%M;
  }

  return Graph(neighbours);
}


Graph examples[2] = {cube(), example1()};



int main(int ac, char **av)
{
  vector<int> rspi(12);
  FullereneGraph::jumplist_t jumps;

  if(ac!=4) {cout << "three arguments required" << endl;}
  const int index = strtol(av[1],0,0) - 1;
  const int K = strtol(av[2],0,0);
  const int L = strtol(av[3],0,0);
  cout << "index, K, L: " << index << ", " << K << ", " <<  L << endl;


  PlanarGraph g(examples[index]);
  cout << "planar graph created" << endl;
  cout << g << endl;
  g.layout2d = g.tutte_layout(0,-1,-1,4);
  cout << "layout created" << endl;

  const int N = g.N;
  ofstream output(("output/C"+to_string(N)+"-unfold.m").c_str());

  PlanarGraph dual(g.dual_graph(6));
  cout << "dual graph created" << endl;
  dual.layout2d = dual.tutte_layout();
  cout << "layout created" << endl;

  output << "g = "  << g << ";\n";
  output << "dg = " << dual << ";\n";
//  cout << "Need to place 2x"<<dual.edge_set.size()<<" edges.\n";

  Unfolding unfld(dual,true);

  
  cout << "Placed " << unfld.edgecoords.size() << " edges.\n";

  output << "dedges   = " << get_keys(unfld.edgecoords) << ";\n";
  output << "dedgepos = " << get_values(unfld.edgecoords) << ";\n";
  output << "outline = "  << unfld.outline << ";\n";
  output << "outlinecoords = " << get_keys(unfld.outline) << ";\n";

  Unfolding gct_unfld = unfld * Eisenstein(K,L);

  output << "gctoutline = " << get_keys(gct_unfld.outline) << ";\n";

  Folding fld(gct_unfld);
  PlanarGraph gctdual = fld.fold(), gct = gctdual.dual_graph(3,false);
  cout << "gctdual = " << gctdual << ";\n"
	 << "gct     = " << gct << ";\n";


  gct.layout2d = gct.tutte_layout();
  Polyhedron P0 = Polyhedron(gct,gct.zero_order_geometry(),6);

  string basename("polyhedron-"+to_string(N));
  {
    ofstream mol2(("output/"+basename+"-P0.mol2").c_str());
    mol2 << P0.to_mol2();
    mol2.close();
  }

   Polyhedron P(P0);
   P.optimize_other();

  cout << P << endl;
  
  {
    ofstream mol2(("output/"+basename+".mol2").c_str());
    mol2 << P.to_mol2();
    mol2.close();
  }


  output.close();
  return 0;
}
