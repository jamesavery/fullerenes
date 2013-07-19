#include <stdlib.h>
#include <fstream>
#include "libgraph/cubicgraph.hh"
#include "libgraph/fullerenegraph.hh"
#include "libgraph/eisenstein.hh"
#include <vector>

using namespace std;


typedef pair<Eisenstein,Eisenstein> dedgecoord_t;
#if 0
Graph GCTransform(const vector< pair<Eisenstein,node_t> > &outline, const map<dedge_t,dedgecoord_t> &edgecoords, int K=1, int L=0)
{
#define insert_node(ijx) \
  if(grid.find(ijx) == grid.end()) /* First time looking at this? */	\
    grid[ijx] = new_node++ \

  map<Eisenstein, node_t> grid;

  const EOp GC(EOp::GC(K,L)), iGC(EOp::GCInverse(K,L)), Delta10(EOp::Delta(Eisenstein(1,0)));
  node_t new_node = 0; 	// Set to number of nodes in original dual

  // Rasterize the standard GC triangle (0,0) -- (k,l) -- (-l,k+l)
  vector<Eisenstein> GCtriangle;
  Eisenstein a(0,0), b(K,L),c(-L,K+L);
  double
    sab = K/double(L),
    sac = -L/double(K+L),
    sbc = (-L-K)/double(K);

  cout << " a = " << a << "; b = " << b << "; c = " << c << endl;

  double xleft = 0, xright = 0;
  for(int j=0;j<L;j++,xleft+=sac,xright+=sab){
    int ileft = ceil(xleft), iright = floor(xright);
    for(int i=ileft;i<=iright;i++) GCtriangle.push_back(Eisenstein(i,j));
    printf("%d: %g -- %g :: %d -- %d\n",j,xleft,xright,ileft,iright);

  }

  xright = b.first;
  for(int j=L;j<=K+L;j++,xleft+=sac,xright+=sbc){
    int ileft = ceil(xleft), iright = floor(xright);
    for(int i=ileft;i<=iright;i++) GCtriangle.push_back(Eisenstein(i,j));
    printf("%d: %g -- %g :: %d -- %d\n",j,xleft,xright,ileft,iright);
  }

  // Begin by computing all vertices inside GC-transformed graph
  for(map<dedge_t,dedgecoord_t>::const_iterator i(edgecoords.begin()); i!= edgecoords.end(); i++){
    const dedgecoord_t &vupos(i->second);

    const Eisenstein ux(vupos.second), vx(vupos.first); // directed edge coordinates are stored in reverse order.

    EOp T(GC*EOp::Delta(vx-ux)*iGC);
    // Send triangle points into their final destinations
    for(int k=0;k<GCtriangle.size();k++)
      insert_node(T*GCtriangle[k] + GC*ux);
  }

  
  // Now transform all exterior triangles to their matching position
  // and identify the nodes with each other
  map<Eisenstein, node_t> grid2;
  map<Eisenstein, node_t>::const_iterator g,g1,g2;
  int No = outline.size();
  for(int ii=0;ii<No;ii++){
    const node_t u = outline[ii].second, v = outline[(ii+1)%No].second;
    const dedgecoord_t &uvx(edgecoords.find(dedge_t(u,v))->second);
    const dedgecoord_t &vux(edgecoords.find(dedge_t(v,u))->second);
    const Eisenstein ux0(uvx.second), vx0(uvx.first);
    const Eisenstein ux1(vux.second), vx1(vux.first);
    
    EOp  T0(GC*EOp::Delta(vx0-ux0)*iGC),
         T1(GC*EOp::Delta(ux1-vx1)*iGC);
    
    for(int k=0;k<GCtriangle.size();k++){
      g = grid.find(T0*GCtriangle[k]+GC*ux0);
      if(g != grid.end())
	grid2[T1*GCtriangle[k]+GC*vx1] = g->second;
    }
  }
  
  // Now connect every node to all of its neighbours or transformed neighbours
  set<edge_t> edge_set, new_edge_set;
  for(g = grid.begin(); g!=grid.end(); g++){
    const Eisenstein x(g->first);
    const node_t u(g->second);

    cerr << "Connecting node " << u << " at " << x << endl;
    
    Eisenstein d(1,0);
    for(int i=0;i<6;i++,d=d.nextCW()){
      g1 = grid.find(x+d);
      g2 = grid2.find(x+d);
      
      if     (g2 != grid2.end()) edge_set.insert(edge_t(u,g2->second));
      else if(g1 != grid.end())  edge_set.insert(edge_t(u,g1->second));
    }
  }
  // write out lines corresponding to edges.
  Graph G(edge_set);

  set<edge_t> same_as;
  Graph same;

  // Finally merge nodes that are supposed to be the same
  for(g2 = grid2.begin(); g2 != grid2.end(); g2++){
    g1 = grid.find(g2->first);

    if(g1 != grid.end() && g1->second != g2->second)
      same_as.insert(edge_t(g1->second,g2->second));
  }


  // New nodes are the connected components of the graph same_as.
  list< list<node_t> > components(Graph(same_as).connected_components());
  cout << "components = " << components << endl;
  map<node_t,node_t> id;
  int i=0;
  for(list< list<node_t> >::const_iterator c(components.begin()); c!=components.end();c++,i++)
    for(list<node_t>::const_iterator u(c->begin()); u!=c->end();u++){
      id[*u] = i;
      printf("new id %d -> %d\n",*u,i);
    }

  // Build graph with new names
  for(set<edge_t>::const_iterator e(edge_set.begin()); e!=edge_set.end();e++){
    const node_t u(id[e->first]), v(id[e->second]);
    if(u!=v) new_edge_set.insert(edge_t(u,v));
  }
    
  ofstream f("misc/grid.m");
  f << "outline = " << outline << ";\n";
  f << "grid  = " << get_keys(grid) << ";\n";
  f << "grid2 = " << get_keys(grid2) << ";\n";
  f << "gct = " << Graph(edge_set) << ";\n";
  f << "gctnew = " << Graph(new_edge_set) << ";\n";
  f.close();

  return Graph(new_edge_set);
}
#endif

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
    cout << i << ": " << label[nextpos] << " at " << nextpos << endl;
    outline[i] = make_pair(nextpos,label[nextpos]);
  }

  // If the outline doesn't form a closed loop, something is wrong.
  assert(nextpos == next.begin()->first);

  return outline;
}



// (Not finished) Output a LaTeX/TiKZ figure of the unfolded triangulation, optionally GC(K,L)-transformed.
void latex_GCunfold(ostream& latexfile, const vector< pair<Eisenstein,node_t> > &outline, const map<dedge_t,dedgecoord_t> &dedge_positions, int K=1, int L=0,
		    bool equilateralp=false, int label_vertices=1, bool include_headers=false)
{
  if(include_headers)
    latexfile << 
"\\documentclass{standalone}\n\
\\usepackage{tikz}\n\
\\begin{document}\n\
\\definecolor{darkgreen}{rgb}{.4, .7, .2}\n\
\\tikzstyle{outline}=[draw=black, ultra thick, fill opacity=.2, fill=darkgreen]\n\
\\tikzstyle{vertex}=[circle, draw, inner sep=0, fill=white, minimum width=4.00000mm]\n\
\n\
\\pgfmathsetmacro{\\xcoord}{cos(60)}\n\
\\pgfmathsetmacro{\\ycoord}{sin(60)}\n\
\n\
\\newcommand{\\drawEGrid}{\n\
    \\path[clip,preaction = {draw=black}] (0,0) -- (\\last,0) -- (\\cols,\\rows) --(\\first,\\rows) -- cycle;\n\
    \\draw (\\first,0) grid (\\last,\\rows);\n\
    \\foreach \\x in {1,2,...,\\total}\n\
        \\draw (-\\x,\\x*2) -- (\\x,0);\n\
}\n\
\\newcommand{\\drawRGrid}{\n\
    \\path[clip,preaction = {draw=black}] (0,0) -- (\\last,0) -- (\\cols,\\rows) --(0,\\rows) -- cycle;\n\
    \\draw (0,0) grid (\\cols,\\rows);\n\
    \\foreach \\x in {1,2,...,\\total}\n\
        \\draw (-\\x,\\x*2) -- (\\x,0);\n\
}\n\
";

  vector<Eisenstein> outline_gc(outline.size());
  for(int i=0;i<outline.size();i++) outline_gc[i] = outline[i].first.GCtransform(K,L);

  // Extract (I,J)-bounds
  int imin = INT_MAX, imax = INT_MIN, jmin = INT_MAX, jmax = INT_MIN;
  for(int i=0;i<outline_gc.size();i++){
    Eisenstein x(outline_gc[i]);

    if(equilateralp){
      coord2d X(x.coord());
      x.first  = floor(X.first);
      x.second = floor(X.second);
    }

    if(x.first < imin) imin = x.first;
    if(x.first > imax) imax = x.first;
    if(x.second < jmin) jmin = x.second;
    if(x.second > jmax) jmax = x.second;

  }
  Eisenstein gcmin(imin-1,jmin-1), gcmax(imax+1,jmax+1);
  if(equilateralp){
    gcmin = Eisenstein(coord2d(gcmin.first-1,gcmin.second));
    gcmax = Eisenstein(coord2d(gcmax.first,gcmax.second+1));
  }

  latexfile << "\\begin{tikzpicture}"<<(K==1 && L==0?"[scale=2.5]":"")<<"\n";
  // Define bounds
  Eisenstein D(gcmax-gcmin);
  latexfile << "\\newcommand*{\\cols}{"<<D.first<<"}\n"
	    << "\\newcommand*{\\rows}{"<<D.second<<"}\n"
	    << "\\newcommand*{\\total}{"<<(D.first+D.second)<<"}\n"
	    << "\\newcommand*{\\first}{-\\rows/2}\n"
	    << "\\newcommand*{\\last}{\\cols+\\rows/2}\n";

  // Draw E-grid
  latexfile << "\\bgroup\n";
  if(equilateralp){
    // -- as equilateral triangles
    latexfile << "\\pgftransformcm{1}{0}{\\xcoord}{\\ycoord}{\\pgfpointorigin}\n"
	      << "\\drawEGrid{}\n";
  } else // -- or as a regular grid
    latexfile << "\\drawRGrid{}\n";

  
  // Draw outline polygon
    latexfile << "\\draw[outline] ";
  for(int i=0;i<outline.size();i++){
    const Eisenstein &x((outline_gc[i]-gcmin));
    latexfile << "(" << x.first << "," << x.second << ") -- ";
  }
  latexfile << "cycle;\n"
	    << "\\egroup\n\n";

  // Place vertex labels according to scheme chosen in parameter 'label_vertices':
  latexfile << "\\foreach \\place/\\name/\\lbl in {";
  switch(label_vertices){
  case 0: break; // Don't label vertices at all.
  case 1:        // Only label vertices on polygon outline.
    if(label_vertices == 1)	
      for(int i=0;i<outline.size();i++){
	const Eisenstein &ij(outline_gc[i]-gcmin);

	coord2d x;
	if(equilateralp) x = ij.coord();
	else             x = coord2d(ij.first,ij.second);

	const node_t &u(outline[i].second);
	latexfile << "{(" << x.first << "," << x.second << ")/"<<i<<"/"<<u<<(i+1<outline.size()?"},":"}");
      }
    break;
  case 2: // Label all original vertices, including internal ones
    {
      int i=0;
      for(map<dedge_t,dedgecoord_t>::const_iterator it(dedge_positions.begin()); it!=dedge_positions.end(); i++){
	node_t u(it->first.first);
	Eisenstein ij(it->second.first.GCtransform(K,L)-gcmin);

	coord2d x;
	if(equilateralp) x = ij.coord();
	else             x = coord2d(ij.first,ij.second);

	latexfile << "{(" << x.first << "," << x.second << ")/"<<i<<"/"<<u<<(++it != dedge_positions.end()? "},":"}");
      }
    }
    break;
  case 3: // Label ALL vertices, both old and new.
    break;
  }
  latexfile << "}\n"
	    << "\t \\node[vertex] (\\name) at \\place {\\lbl};\n\n"
	    << "\\end{tikzpicture}\n";
  if(include_headers) 
    latexfile << "\\end{document}\n";
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
 
  ofstream output("output/C"+to_string(N)+"-unfold.m");


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

  output.close();
  ofstream latex_output("output/C"+to_string(N)+"-GC"+to_string(K)
			+"x"+to_string(L)+"-unfold.tex");
  latex_GCunfold(latex_output,outline,dgrid,K,L,false,2,true);
  latex_output.close();

  //  Graph gct = GCTransform(outline, dgrid, K,L);
  //  cout << "gct = " << gct << ";\n";
  
  vector<Eisenstein> outline_coords(get_keys(outline));
  //  const vector< set<int> > raster(Eisenstein::rasterize_polygon(EOp::GC(K,L)*outline_coords));
  //  vector<int> raster(Eisenstein::rasterize_line(Eisenstein(0,0),Eisenstein(2,3)));

  //  cout << "raster = " << raster << endl;  

  return 0;
}
