#include <stdlib.h>
#include <fstream>
#include "libgraph/cubicgraph.hh"
#include "libgraph/fullerenegraph.hh"
#include "libgraph/eisenstein.hh"
#include <vector>

using namespace std;


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

vector<edge_t> scanlineedges(const IDCounter<Eisenstein>& inner_nodes, const polygon::scanline& S, int N = 0, const Eisenstein& transform = Eisenstein(1,0))
{
  vector<edge_t> edges;
  
  for(int i=0;i<S.xs.size();i++){
    const vector<int> &row(S.xs[i]);
    for(int j=0;j<row.size()/2;j++)
      for(int x=row[2*j];x+1<=row[2*j+1];x++){
	Eisenstein x0(x,S.minY+i), x1(x+1,S.minY+i);
	node_t u = inner_nodes(x0*transform)+N, v = inner_nodes(x1*transform)+N;
	// cout << "0: inner_nodes("<<(x0*transform)<<") = " << u << endl;
	// cout << "1: inner_nodes("<<(x1*transform)<<") = " << v << endl;
	edges.push_back(edge_t(u,v));
      }
  }
  return edges;
}

// TODO: Get rid of GC-transform, rename to "fold", and make this work for any polygon from unfolded triangulation?
PlanarGraph GCTransform(const PlanarGraph& dual, int K=1, int L=0)
{
  vector<face_t> faces(dual.compute_faces_flat(3,true));
  vector<tri_t>  triangles(faces.begin(),faces.end());
  map<dedge_t,dedgecoord_t>  dgrid(unfold(triangles));
  vector< pair<Eisenstein,node_t> > outline(get_outline(dgrid));


  // TODO: Split into ConnectInner() and FoldOutline()
  IDCounter<Eisenstein> inner_nodes, outer_nodes;
  vector<Eisenstein> outline_coords(get_keys(outline));
  polygon outline_polygon = convert_vector<Eisenstein,pair<int,int> >(Eisenstein(K,L)*outline_coords);
  polygon::scanline scans = outline_polygon.scanConvert();

  // Register all nodes on interior of polygon.
  for(int i=0;i<scans.xs.size();i++){
    const vector<int> &row(scans.xs[i]);
    for(int j=0;j<row.size()/2;j++)
      for(int x=row[2*j];x<=row[2*j+1];x++){
	Eisenstein xy(x,i+scans.minY);
	inner_nodes.insert(xy);
      }
  }
  
  // Register all nodes on boundary of polygon
  for(int i=0;i<scans.edge_xs.size();i++){
    const vector<int> &pts(scans.edge_xs[i]);
    for(int j=0;j<pts.size();j++) outer_nodes.insert(Eisenstein(pts[j],i+scans.minY));
  }
  
  cout << "inner_nodes.keys   = " << get_keys(inner_nodes) << endl;
  cout << "inner_nodes.values = " << get_values(inner_nodes) << endl;

  // Now connect internal node grid
  polygon outline_polygonCW = convert_vector<Eisenstein,pair<int,int> >(Eisenstein(K,L).nextCW()*outline_coords),
    outline_polygonCCW = convert_vector<Eisenstein,pair<int,int> >(Eisenstein(K,L).nextCCW()*outline_coords);

  polygon::scanline 
    scansA = outline_polygonCW.scanConvert(),
    scansB = outline_polygonCCW.scanConvert();

  vector<edge_t> 
    e (scanlineedges(inner_nodes,scans,0)), 
    eA(scanlineedges(inner_nodes,scansA,0,Eisenstein(0,1))), 
    eB(scanlineedges(inner_nodes,scansB,0,Eisenstein(1,-1)));
 
  // Connect outline:
  // 1. Sort edge_xs CW
  // 2. Traverse outline while connecting edge nodes

  // Connect outline to interior:
  // For each outline node u, if x_u + (1,0) or x_u + (1,-1) is in interior, make edge.

  // Identify identical nodes on outline:
  // Problem: How to connect nodes on interior across cuts?
  // Idea: 
  //  1. Make bounding box grid Rab of a--b and Rba of b--a, using nodes from inner_nodes and outer_nodes.
  //  2. Transform Rab with T(ab->ba). 
  //  3. For each node u in Rab
  //  3.1 If x(ab,u) in outer_nodes, identify u with v = Rba(x(ba,u)) (replace v by u everywhere)
  //  4. For each edge u--v in Rab
  //  4.1 If Rab(T(u)) is inner, Rab(T(v)) is undefined, and Rba(T(v)) is inner, add edge (u,v)


  for(int i=0;i<scans.edge_xs.size();i++){
    const vector<int> &pts(scans.edge_xs[i]);
    for(int j=0;j<pts.size();j++) outer_nodes.insert(Eisenstein(pts[j],i+scans.minY));
  }

  // Join up all the edges
  vector<edge_t> inner_edges(e.begin(),e.end());
  inner_edges.reserve(e.size()+eA.size()+eB.size());
  copy(eA.begin(),eA.end(),inserter(inner_edges,inner_edges.end()));
  copy(eB.begin(),eB.end(),inserter(inner_edges,inner_edges.end()));

  printf("%ld + %ld + %ld = %ld\n",e.size(),eA.size(),eB.size(),inner_edges.size());

  Graph dualGC(set<edge_t>(inner_edges.begin(), inner_edges.end()));
  
  int N = inner_nodes.nextid;
  vector<coord2d> layout2d(N);
  for(IDCounter<Eisenstein>::const_iterator ei(inner_nodes.begin()); ei!=inner_nodes.end(); ei++){
    Eisenstein xy(ei->first);
    node_t      u(ei->second);
    layout2d[u] = coord2d(xy.first,xy.second);
  }

  return PlanarGraph(dualGC,layout2d);
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

  output << "gct = " << GCTransform(dual,K,L) << ";\n";

  ofstream latex_output(("output/C"+to_string(N)+"-GC"+to_string(K)
			 +"x"+to_string(L)+"-unfold.tex").c_str());
  latex_GCunfold(latex_output,outline,dgrid,K,L,true,2,true);
  latex_output.close(); 

  //  Graph gct = GCTransform(outline, dgrid, K,L);
  //  cout << "gct = " << gct << ";\n";
  
  vector<Eisenstein> outline_coords(get_keys(outline));
  typedef pair<int,int> coord;
  polygon outline_polygon(convert_vector<Eisenstein, coord>(Eisenstein(K,L)*outline_coords));
  
  cout << "P = " << outline_polygon << ";\n";
  
  polygon::scanline scans(outline_polygon.scanConvert());
  
  output << "scans[\"minY\"]      = " << scans.minY << ";\n";
  output << "scans[\"xs\"]        = " << scans.xs << ";\n";
  output << "scans[\"edge_xs\"]   = " << scans.edge_xs << ";\n";

  //  const vector< set<int> > raster(Eisenstein::rasterize_polygon(EOp::GC(K,L)*outline_coords));
  //  vector<int> raster(Eisenstein::rasterize_line(Eisenstein(0,0),Eisenstein(2,3)));

  //  cout << "raster = " << raster << endl;  

  output.close();
  return 0;
}
