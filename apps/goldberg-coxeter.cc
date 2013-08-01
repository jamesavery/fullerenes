#include <stdlib.h>
#include <fstream>
#include "libgraph/cubicgraph.hh"
#include "libgraph/fullerenegraph.hh"
#include "libgraph/eisenstein.hh"
#include <vector>

using namespace std;

PlanarGraph fold(vector< pair<Eisenstein, node_t> > &outline);

typedef pair<Eisenstein,Eisenstein> dedgecoord_t;


struct dedge_sort : public std::binary_function<dedge_t, dedge_t, bool>
{
    bool operator()(const dedge_t &x, const dedge_t &y) const
    {   
      int maxx = max(x.first,x.second), maxy = max(y.first,y.second);
      return maxx < maxy || (maxx == maxy && min(x.first,x.second) < min(y.first,y.second));
    }
};



// This function unfolds a triangulation and lays it out on an equilateran
// triangular grid, such that if one were to cut along the outline
// and glue together the nodes with the same labels, one would obtain
// again the original fullerene dual.
//
// Preconditions: Triangles are oriented consistently, i.e. CW or CCW.
// TODO: A more intricate directed edge selection scheme could lead to
// more compact unfoldings (i.e. with more interior points).
map<dedge_t,dedgecoord_t> unfold(const vector<tri_t> &triangulation)
{
#define set_dedge(u,v,ux,vx) {	          \
  dedge_t uv(u,v), vu(v,u);               \
  dedge_done[uv] = true;                  \
  workset.erase(uv);                      \
  dedge_position[vu] = make_pair(vx,ux);  \
  if(!dedge_done[vu])                     \
    workset.insert(vu);			  \
}

  // A single directed edge uniquely defines the third node in the oriented triangle
  map<dedge_t,node_t> nextNode;
  for(int i=0;i<triangulation.size();i++){
    const tri_t &t(triangulation[i]);
    for(int j=0;j<3;j++)
      nextNode[dedge_t(t[j],t[(j+1)%3])] = t[(j+2)%3];
  }

  map<dedge_t,bool> dedge_done;
  set<dedge_t, dedge_sort>   workset;
  map<dedge_t, dedgecoord_t > dedge_position
;
  map<Eisenstein,node_t> grid;
  Eisenstein zero(0,0), veci(1,0), vecj(0,1);

  // 1. Place first triangle. 
  tri_t t(triangulation[0]);

  set_dedge(t[0],t[1],zero,veci);
  set_dedge(t[1],t[2],veci,veci-vecj);
  set_dedge(t[2],t[0],veci-vecj,zero);

  while(!workset.empty()){
    dedge_t uv(*workset.rbegin()); // Next placeable directed edge 
    // set_triangle(uv)
    node_t u(uv.first), v(uv.second), w(nextNode[uv]);

    dedgecoord_t uvpos(dedge_position[uv]);
    Eisenstein ux(uvpos.first), vx(uvpos.second), wx(ux+(vx-ux).nextCW());

    set_dedge(u,v,ux,vx);
    set_dedge(v,w,vx,wx);
    set_dedge(w,u,wx,ux);

  }
  return dedge_position;
}


// Given the output of unfold(), this function efficiently computes the polygon outlining
// the unfolded triangulation and returns it in clockwise order. 
vector< pair<Eisenstein,node_t> > get_outline(const map<dedge_t,dedgecoord_t>& edgecoords)
{
  map<Eisenstein,node_t>    label;
  map<Eisenstein,Eisenstein> next;

  // Collect the directed edges u->v whose positions do not coincide with the reverse edge v->u. 
  // These form the outline of the polygon. 
  for(map<dedge_t,dedgecoord_t>::const_iterator i(edgecoords.begin()); i!= edgecoords.end(); i++){
    const dedge_t &uv(i->first), vu(uv.second,uv.first);
    const dedgecoord_t &uvpos(i->second), vupos(edgecoords.find(vu)->second);

    if(uvpos != make_pair(vupos.second,vupos.first)){
      next[uvpos.first]   = uvpos.second; 
      label[uvpos.first]  = uv.first;
    }
  }
  
   // Now we're ready to find a CW walk through the polygon outline coordinates 
  // and to assign the corresponding nodes in the triangulation graph as labels.
  vector< pair<Eisenstein,node_t> > outline(next.size());
  Eisenstein nextpos = next.begin()->first;

  for(int i=0;i<outline.size();i++, nextpos = next[nextpos]){
    //    cout << i << ": " << label[nextpos] << " at " << nextpos << endl;
    outline[i] = make_pair(nextpos,label[nextpos]);
  }

  // If the outline doesn't form a closed loop, something is wrong.
  assert(nextpos == next.begin()->first);

  reverse(outline.begin(),outline.end());
  return outline;
}



PlanarGraph GCTransform(const PlanarGraph& dual, int K=1, int L=0)
{
  vector<face_t> faces(dual.compute_faces_flat(3,true));
  vector<tri_t>  triangles(faces.begin(),faces.end());

  map<dedge_t,dedgecoord_t>  dgrid(unfold(triangles));
  vector< pair<Eisenstein,node_t> > outline(get_outline(dgrid));

  for(int i=0;i<outline.size();i++) 
    outline[i].first = outline[i].first * Eisenstein(K,L);

  return fold(outline);
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


void transform_line(const dedgecoord_t& l1, const dedgecoord_t& l2, Eisenstein& x0, Eisenstein& x0p, Eisenstein& w)
{
  Eisenstein Duv(l1.second-l1.first), Dvu(l2.first-l2.second), Tuvvu((Duv.invertn()*Dvu)/Dvu.norm2());

  x0  = l1.first;
  x0p = l2.second;
  w   = Tuvvu;
}

Eisenstein tfm(const Eisenstein& x, const Eisenstein& x0, const Eisenstein& w, const Eisenstein& x0p)
{
  return (x-x0)*w + x0p;
}

int turn_direction(const Eisenstein& xi,const Eisenstein& xj,const Eisenstein& xk) 
{
  Eisenstein dx1(xj-xi), dx2(xk-xj);
  return sgn(dx2.first * dx1.second - dx2.second * dx1.first);
}


vector<edge_t> connect_cross(vector< pair<Eisenstein, node_t> > &outline, const IDCounter<Eisenstein>& grid, const Eisenstein& w)
{
  vector<edge_t> edges;

  // Register reverse arcs
  map<dedge_t,dedgecoord_t> reverse_arc;
  vector< pair<Eisenstein,Eisenstein> > EC, ECp;

  for(int i=0;i<outline.size();i++){
    const Eisenstein xu(outline[i].first), xv(outline[(i+1)%outline.size()].first);
    const node_t u(outline[i].second), v(outline[(i+1)%outline.size()].second);
    reverse_arc[dedge_t(v,u)] = dedgecoord_t(xu,xv);
  }

  Eisenstein iw(w.invertn());
  for(int i=0;i<outline.size();i++){
    Eisenstein X0(outline[i].first*w), X1(outline[(i+1)%outline.size()].first*w);
    node_t     U(outline[i].second), V(outline[(i+1)%outline.size()].second);
    dedgecoord_t Xuv(X0,X1), Xvu(reverse_arc[dedge_t(U,V)]);
    Xvu.first  = Xvu.first *w;
    Xvu.second = Xvu.second*w;

    Eisenstein x0,x0p,T;
    transform_line(Xuv,Xvu, x0,x0p, T);

    vector<Eisenstein> segment(polygon::draw_line(X0,X1)), 
      revsegment(polygon::draw_line(X1,X0));
    reverse(revsegment.begin(),revsegment.end());
    assert(segment.size() == revsegment.size());

    for(int j=0;j<segment.size();j++){
      const Eisenstein& x(segment[j]), y(revsegment[j]);
      if(x != y){
	Eisenstein xp((x-x0)*T+x0p);
	Eisenstein yp((y-x0)*T+x0p);
	node_t u = grid(x*iw), v = grid(yp*iw);
	assert(u>=0 && v>=0);
	edges.push_back(edge_t(u,v));
	
	EC.push_back(make_pair(x,y));
	ECp.push_back(make_pair(xp,yp));
      }
      
    }
  }

  polygon P(w*get_keys(outline));
  ofstream output("output/connect.m");
  output
    << "points = " << P.allpoints() << ";\n"
    << "outlinecoords = " << P.outline << ";\n" 
    << "EC = " << EC << ";\n"
    << "ECp = " << ECp << ";\n"
    ;
  output.close();

  return edges;
}

vector<edge_t> connect_polygon(vector< pair<Eisenstein, node_t> > &outline, const IDCounter<Eisenstein>& grid, const Eisenstein& w)
{
  Eisenstein iw(w.invertn());
  vector<Eisenstein> outline_coords(w*get_keys(outline));
  polygon P(outline_coords);
  polygon::scanline S(P.scanConvert());
  
  vector<edge_t> edges;
  for(int i=0;i<S.xs.size();i++){
    for(int j=0;j<S.xs[i].size()/2;j++){
      int start = S.xs[i][2*j], end = S.xs[i][2*j+1];
      for(int x=start;x<end;x++){
	node_t u = grid(Eisenstein(x,i+S.minY)*iw), v = grid(Eisenstein(x+1,i+S.minY)*iw);
	if(u>=0 && v>=0) edges.push_back(edge_t(u,v));
      }
    }
  }
  return edges;
}

vector<edge_t> connect(vector< pair<Eisenstein, node_t> > &outline, const IDCounter<Eisenstein>& grid, const Eisenstein& w)
{
  vector<edge_t> 
    polygon_edges(connect_polygon(outline,grid,w)), 
    cross_edges(connect_cross(outline,grid,w)),
    edges;
  
  copy(polygon_edges.begin(),polygon_edges.end(),back_inserter(edges));
  copy(cross_edges.begin(),cross_edges.end(),back_inserter(edges));

  return edges;
}

vector<edge_t> connect(vector< pair<Eisenstein, node_t> > &outline, const IDCounter<Eisenstein>& grid)
{
  vector<edge_t> 
    e(connect(outline,grid,Eisenstein(1,0))),
    eCW(connect(outline,grid,Eisenstein(1,-1))),
    eCCW(connect(outline,grid,Eisenstein(0,1))), edges;
  
  copy(e.begin(),   e.end(),back_inserter(edges));
  copy(eCW.begin(), eCW.end(),back_inserter(edges));
  copy(eCCW.begin(),eCCW.end(),back_inserter(edges));

  return edges;
}


vector<int> identify_nodes(vector< pair<Eisenstein, node_t> > &outline, const IDCounter<Eisenstein>& grid)
{
  set<edge_t> same_as;

  // Register reverse arcs
  map<dedge_t,dedgecoord_t> reverse_arc;

  for(int i=0;i<outline.size();i++){
    const Eisenstein xu(outline[i].first), xv(outline[(i+1)%outline.size()].first);
    const node_t u(outline[i].second), v(outline[(i+1)%outline.size()].second);
    reverse_arc[dedge_t(v,u)] = dedgecoord_t(xu,xv);
  }

  for(int i=0;i<outline.size();i++){
    Eisenstein X0(outline[i].first), X1(outline[(i+1)%outline.size()].first);
    node_t     U(outline[i].second), V(outline[(i+1)%outline.size()].second);
    dedgecoord_t Xuv(X0,X1), Xvu(reverse_arc[dedge_t(U,V)]);
    Xvu.first  = Xvu.first ;
    Xvu.second = Xvu.second;

    Eisenstein x0,x0p,T;
    transform_line(Xuv,Xvu, x0,x0p, T);

    //TODO: Handle horizontal lines.
    vector<Eisenstein> segment(polygon::draw_line(X0,X1)), 
      revsegment(polygon::draw_line(X1,X0));
    reverse(revsegment.begin(),revsegment.end());
    assert(segment.size() == revsegment.size());

    for(int j=0;j<segment.size();j++){
      const Eisenstein& x(segment[j]), y(revsegment[j]);
      if(x == y){
	Eisenstein xp((x-x0)*T+x0p);
	node_t u = grid(x), v = grid(xp);
	assert(u>=0 && v>=0);
	same_as.insert(edge_t(u,v));
      }
    }
  }


  // Find connected components
  vector<int> same(grid.size());
  for(int i=0;i<grid.size();i++) same[i] = i;

  Graph S(same_as);
  list< list<node_t> > components(S.connected_components());

  //  cout << "samecomponents = " << components << ";\n";

  for(list< list<node_t> >::const_iterator s(components.begin()); s!=components.end(); s++){
    node_t canonical = *s->begin();

    for(list<node_t>::const_iterator t(s->begin()); t!=s->end(); t++){
      same[*t] = canonical;
    }
  }

  return same;
}


PlanarGraph fold(vector< pair<Eisenstein, node_t> > &outline)

{
  IDCounter<Eisenstein> grid;
  vector<Eisenstein> outline_coords(get_keys(outline));
  polygon P(outline_coords);
  set<Eisenstein> allpoints(P.allpoints());

  for(set<Eisenstein>::const_iterator x(allpoints.begin()); x!=allpoints.end();x++) grid.insert(*x);


  vector<edge_t> edge_list(connect(outline,grid));

  vector<int> same_as = identify_nodes(outline,grid);



  if(1) for(int i=0;i<edge_list.size();i++){
      node_t u = edge_list[i].first, v = edge_list[i].second;
      edge_list[i].first = same_as[u];
      edge_list[i].second = same_as[v];
  }
  

  // Compactify node names
  set<edge_t> edges;
  IDCounter<int> new_names;

  if(1) {
    for(int i=0;i<edge_list.size();i++){
      edge_t e = edge_list[i];
      node_t u = e.first, v = e.second;

      int newu = new_names.insert(u);
      int newv = new_names.insert(v);

      edges.insert(edge_t(newu,newv));
    }
  } else {
    edges = set<edge_t>(edge_list.begin(),edge_list.end());
    for(int i=0;i<grid.size();i++) new_names.insert(i);
  }

  int N = new_names.size();
  vector<coord2d> layout2d(N);
  for(IDCounter<Eisenstein>::const_iterator xi(grid.begin()); xi!=grid.end(); xi++){
    Eisenstein xy(xi->first);
    node_t      u(new_names(xi->second));

    if(u>=0) layout2d[u] = coord2d(xy.first,xy.second);
  }

  PlanarGraph G(edges,layout2d);
  return G;
}


int main(int ac, char **av)
{
  vector<int> rspi(12);
  FullereneGraph::jumplist_t jumps;

  if(ac<13) return -1;

  int N = strtol(av[1],0,0), K=1, L=0;
  for(int i=0;i<12;i++) rspi[i] = strtol(av[2+i],0,0)-1;
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

  vector<face_t> faces(dual.compute_faces_flat(3,true));
  vector<tri_t>  triangles(faces.begin(),faces.end());

  map<dedge_t,dedgecoord_t>         dgrid(unfold(triangles));

  cout << "Placed " << dgrid.size() << " edges.\n";

  output << "dedges   = " << get_keys(dgrid) << ";\n";
  output << "dedgepos = " << get_values(dgrid) << ";\n";

  vector< pair<Eisenstein,node_t> > outline(get_outline(dgrid));
  output << "outline = " << outline << ";\n";
  output << "outlinecoords = " << get_keys(outline) << ";\n";

  vector<pair<Eisenstein,node_t> > gct_outline(outline);
  for(int i=0;i<outline.size();i++) gct_outline[i].first = outline[i].first * Eisenstein(K,L);

  output << "gctoutline = " << get_keys(gct_outline) << ";\n";
  cout << "gctoutline = " << get_keys(gct_outline) << ";\n";

  output << "gct = " << GCTransform(dual,K,L) << ";\n";

  //  vector<pair<Eisenstein,node_t> > reduced_outline(GCDreduce(outline));

  //  cout << "outline  = " << outline << ";\n"
  //       << "routline = " << reduced_outline << ";\n";
  //  output << "reduced = " << fold(reduced_outline) << ";\n";

  //  output << "gctb = " << fold(gct_outline) << ";\n";



  //  Graph gct = GCTransform(outline, dgrid, K,L);
  //  cout << "gct = " << gct << ";\n";
  
  vector<Eisenstein> outline_coords(get_keys(outline));
  typedef pair<int,int> coord;
  polygon outline_polygon(convert_vector<Eisenstein, coord>(Eisenstein(K,L)*outline_coords));
  
  cout << "P = " << outline_polygon << ";\n";
  
  polygon::scanline scans(outline_polygon.scanConvert());
  
  output << "scans[\"minY\"]      = " << scans.minY << ";\n";
  output << "scans[\"xs\"]        = " << scans.xs << ";\n";

  output.close();
  return 0;
}
