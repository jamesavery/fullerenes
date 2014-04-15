#include <stdlib.h>
#include <fstream>
#include "libgraph/cubicgraph.hh"
#include "libgraph/fullerenegraph.hh"
#include "libgraph/eisenstein.hh"
#include "libgraph/unfold.hh"
#include "libgraph/triangulation.hh"
#include <vector>

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


int N_default = 28;
int rspi_default[12] = {1,2,3,5,7,9,10,11,12,13,14,15};

int main(int ac, char **av)
{
  vector<int> rspi(12);
  FullereneGraph::jumplist_t jumps;

  int N = N_default, K=1, L=0;
  if(ac<13){
    rspi = vector<int>(rspi_default,rspi_default+12);
    for(int i=0;i<12;i++) rspi[i]--;
  }

  else {
    N = strtol(av[1],0,0);
    for(int i=0;i<12;i++) rspi[i] = strtol(av[2+i],0,0)-1;
  }
  if(ac>=15){
    K = strtol(av[14],0,0);
    L = strtol(av[15],0,0);
  }

  cout << "Attempting to create graph from spiral indices " << rspi << endl;
 
  ofstream output(("output/C"+to_string(N)+"-unfold.m").c_str());


  FullereneGraph g(N, rspi, jumps);
  PlanarGraph dual(g.dual_graph(6));

  g.layout2d = g.tutte_layout();
  dual.layout2d = dual.tutte_layout();

  output << "g = "  << g << ";\n";
  output << "dg = " << dual << ";\n";
  cout << "Need to place 2x"<<dual.edge_set.size()<<" edges.\n";

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
  output << "gctdual = " << gctdual << ";\n"
	 << "gct     = " << gct << ";\n";

  vector<int> gct_rspi(12);
  FullereneGraph::jumplist_t gct_jumps;
  FullereneDual Tgctdual(gctdual);
  
  Tgctdual.get_rspi(gct_rspi,gct_jumps,false,false);

  output << "RSPI = "    << (rspi+1) << ";\n"
         << "gctRSPI = " << (gct_rspi+1) << ";\n";
  
  cout << "P = " << fld.P << ";\n";
  
  polygon::scanline scans(fld.P.scanConvert());
  
  output << "scans[\"minY\"]      = " << scans.minY << ";\n";
  output << "scans[\"xs\"]        = " << scans.xs << ";\n";

  output.close();
  return 0;
}
