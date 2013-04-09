#include <stdlib.h>
#include <fstream>
#include "libgraph/cubicgraph.hh"
#include "libgraph/fullerenegraph.hh"

#include <vector>

using namespace std;

class Eisenstein: public pair<int,int> {
public:
  Eisenstein(int a=0, int b=0) : pair<int,int>(a,b) {}
  Eisenstein(const coord2d& x) : pair<int,int>(round(x.first-x.second/sqrt(3)), round(2*x.second/sqrt(3)))
  { }
  Eisenstein operator*(const Eisenstein& y) const { return Eisenstein(first*y.first,second*y.second); }
  Eisenstein operator+(const Eisenstein& y) const { return Eisenstein(first+y.first,second+y.second); }
  Eisenstein operator-(const Eisenstein& y) const { return Eisenstein(first-y.first,second-y.second); } 
  Eisenstein& operator+=(const Eisenstein& y) { first += y.first; second += y.second; return *this; }
  Eisenstein& operator-=(const Eisenstein& y) { first -= y.first; second -= y.second; return *this; }

  Eisenstein GCtransform(int k, int l) const {
    return Eisenstein(k*first - l*second, l*first + (k+l)*second);
  }
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

  // Given a line segment, return all the grid points either directly on the line
  // or of distance < 1 to it on the right hand side.
  static vector<Eisenstein> rightofline(const Eisenstein &x0, const Eisenstein &x1)
  {
    vector<Eisenstein> result;

    // Deal with degenerate cases first
    if(x0.first == x1.first){
      result.reserve(x1.second - x0.second + 1);
      for(int y=x0.second; y<=x1.second; y++) result.push_back(Eisenstein(x0.first,y));
    } 

    if(x0.second == x1.second){
      result.reserve(x1.first - x0.first + 1);
      for(int x=x0.first; x<=x1.first; x++) result.push_back(Eisenstein(x,x0.second));
    }

    // Line segment is not vertical or horizontal
    int Dx = x1.first - x0.first, Dy = x1.second - x0.second;
    double slope = Dx/Dy, sign = slope/abs(slope);

    for(int i=0,y=x0.second; y<=x1.second; y++,i++){
      
    }

    return result;
  }

  coord2d coord() const { return coord2d(1,0)*first + coord2d(0.5,0.8660254037844386)*second; }
};

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
  printf("set_dedge(%d,%d)\n",u,v);       \
  dedge_done[uv] = true;                  \
  workset.erase(uv);                      \
  dedge_position[vu] = make_pair(vx,ux);  \
  if(!dedge_done[vu]){                    \
    workset.insert(vu);			  \
  } \
}

  map<dedge_t,node_t> nextNode;
  for(int i=0;i<triangulation.size();i++){
    const tri_t &t(triangulation[i]);
    for(int j=0;j<3;j++)
      nextNode[dedge_t(t[j],t[(j+1)%3])] = t[(j+2)%3];
  }

  map<dedge_t,bool> dedge_done;
  set<dedge_t,dedge_sort>   workset;
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
    dedge_t uv(*workset.rbegin());
    //    cout << "Next unused dedge is " << uv << endl;
    // set_triangle(uv)
    node_t u(uv.first), v(uv.second), w(nextNode[uv]);

    dedgecoord_t uvpos(dedge_position[uv]);
    Eisenstein ux(uvpos.first), vx(uvpos.second), wx(ux+Eisenstein::nextCW(ux,vx));

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

  map<dedge_t,dedgecoord_t>         grid(unfold(triangles));

  cout << "Placed " << grid.size() << " edges.\n";

  output << "dedges   = " << get_keys(grid) << ";\n";
  output << "dedgepos = " << get_values(grid) << ";\n";

  vector< pair<Eisenstein,node_t> > outline(get_outline(grid));
  output << "outline = " << outline << ";\n";

  output.close();
  ofstream latex_output("output/C"+to_string(N)+"-GC"+to_string(K)
			+"x"+to_string(L)+"-unfold.tex");
  latex_GCunfold(latex_output,outline,grid,K,L,true,2,true);
  latex_output.close();
  
  return 0;
}
