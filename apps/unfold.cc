#include <stdlib.h>
#include <fstream>
#include "libgraph/cubicgraph.hh"
#include "libgraph/fullerenegraph.hh"
#include "libgraph/eisenstein.hh"
#include <vector>

using namespace std;


typedef pair<Eisenstein,Eisenstein> dedgecoord_t;

int right_of(const pair<Eisenstein,Eisenstein>& line, const Eisenstein& x)
{
  Eisenstein Dx(line.second - line.first), xr(x-line.first);
  return sgn(xr.first * Dx.second - xr.second * Dx.first);
}


vector<Eisenstein> bresenham_line(const Eisenstein& x0, const Eisenstein& x1)
{
  vector<Eisenstein> result;

  Eisenstein D(x1-x0);
  int Dx = D.first, Dy = D.second;

  if(abs(Dx) > abs(Dy)) { // Draw transposed line
    printf("X is major axis; transposing.\n");
    result = bresenham_line(x0.transpose(),x1.transpose());
    for(int i=0;i<result.size();i++) result[i] = result[i].transpose();
    return result;
  }
  // Now Y is major axis

  int sx = sgn(Dx), sy = sgn(Dy);
  Eisenstein xy(x0);

  if(sx == 0) 
    while(xy.second != x1.second){ result.push_back(xy); xy.second += sy; }

  int epsilon = 0;
  do{
    result.push_back(xy);

    if(2*(epsilon+Dx) < Dy) {
      epsilon += sx*Dx;
    } else {
      xy.first+=sx;
      epsilon += sx*Dx-sy*Dy;
    } 
    xy.second+=sy;
  } while(xy.second-sy != x1.second);

  return result;
}

int peak(const Eisenstein& a,const Eisenstein& b,const Eisenstein& c)
{
  return sgn((b.second - a.second)*(c.second - a.second));
}

struct dedge_sort : public std::binary_function<dedge_t, dedge_t, bool>
{
    bool operator()(const dedge_t &x, const dedge_t &y) const
    {   
      int maxx = max(x.first,x.second), maxy = max(y.first,y.second);
      return maxx < maxy || (maxx == maxy && min(x.first,x.second) < min(y.first,y.second));
    }
};


struct sort_ccw_eisenstein {
  const coord2d centre;
  sort_ccw_eisenstein(const coord2d& centre) : centre(centre) { }
  
  bool operator()(const Eisenstein& x, const Eisenstein& y) const {
    return angle(x) >= angle(y);
  }

  double angle(const Eisenstein& x) const {
    coord2d dx(x-centre);
    return atan2(dx.first,dx.second);
  }
  double operator()(const Eisenstein& x) const { return angle(x); }
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

vector<edge_t> scanlineedges(const IDCounter<Eisenstein>& inner_nodes, const IDCounter<Eisenstein>& grid, const polygon::scanline& S, int N = 0, const Eisenstein& transform = Eisenstein(1,0))
{
  vector<edge_t> edges;
  
  for(int i=0;i<S.xs.size();i++){
    const vector<int> &row(S.xs[i]);
    for(int j=0;j<row.size()/2;j++)
      for(int x=row[2*j];x+1<=row[2*j+1];x++){
	Eisenstein x0(x,S.minY+i), x1(x+1,S.minY+i);
	node_t u = grid(x0*transform), v = grid(x1*transform);

	edges.push_back(edge_t(u,v));
      }
  }

  return edges;
}

set<edge_t> connect_corners(const vector<Eisenstein>& edge_coords, const IDCounter<Eisenstein>& outer_points)
{
  set<edge_t> edges;
  for(int i=0;i<edge_coords.size();i++){
    int j = (i+1)%edge_coords.size();
    int k = (i+2)%edge_coords.size();
    int trn = Eisenstein::turn(edge_coords[i],edge_coords[j],edge_coords[k]);
    //    cout << "Node " << (outer_points(edge_coords[j])+1) << " at " << edge_coords[j] << " has turn " << trn << endl;
    if(trn == 1 && (edge_coords[i]-edge_coords[k]).isUnit()) // j is right turning corner
      edges.insert(edge_t(outer_points(edge_coords[i]), outer_points(edge_coords[k])));
  }
  return edges;
}


PlanarGraph fold(vector< pair<Eisenstein, node_t> > &outline);

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




vector< vector<int> >  arcgrid(const Eisenstein& Duv)
{
  Eisenstein D(Duv.abs());
  vector<vector<int> > G(D.first+1, vector<int>(D.second+1));
  
  for(int i=0;i<=D.first;i++)
    for(int j=0;j<=D.second;j++){
      int side = D.second*i-D.first*j;
      G[i][j] = side < 0? -1 : (side == 0? 0 : 1);
    }
  return G;
}

vector< pair<Eisenstein, node_t> > GCDreduce(const vector< pair<Eisenstein, node_t> > &outline)
{
  vector<Eisenstein> segments(outline.size());

  for(int i=0;i<outline.size();i++) segments[i] = outline[(i+1)%outline.size()].first - outline[i].first;

  cout << "segments = " << segments << ";\n";

  Eisenstein d(Eisenstein::gcd(segments)); // TODO: Only do GCD between pentagon nodes.

  cout << "GCD = " << d << endl;
  for(int i=0;i<segments.size();i++) segments[i] = segments[i].div(d);

  vector< pair<Eisenstein,node_t> > new_outline(outline.size());
  new_outline[0] = outline[0];
  for(int i=0;i+1<outline.size();i++) new_outline[i+1] = make_pair(new_outline[i].first+segments[i],outline[i].second);

  return new_outline;
}


PlanarGraph GCTransformTCG(const PlanarGraph& dual, int K=1, int L=0)
{
  vector<face_t> faces(dual.compute_faces_flat(3,true));
  vector<tri_t>  triangles(faces.begin(),faces.end());

  map<dedge_t,dedgecoord_t>  dgrid(unfold(triangles));
  vector< pair<Eisenstein,node_t> > outline(get_outline(dgrid));

  cout << "outline = " << outline << endl;
  for(int i=0;i<outline.size();i++) 
    outline[i].first = outline[i].first * Eisenstein(K,L);
  cout << "GC(outline) = " << outline << endl;

  vector< pair<Eisenstein,node_t> > reduced_outline(GCDreduce(outline));

  cout << "reduce(GC(outline)) = " << reduced_outline << endl;

  return fold(reduced_outline);
}

int gridnode(const dedgecoord_t& xuv,  const Eisenstein& x, const IDCounter<Eisenstein>& grid, const Eisenstein& xu, const Eisenstein& Tuvvu)
{
  switch(right_of(xuv,x)){
  case 1:
    return grid((x-xuv.first)*Tuvvu+xu);
  case -1:
  case  0:
    return grid(x);
  }
  return -2;
}


set<edge_t> connect_edge(const vector< pair<Eisenstein,node_t> >& outline, const Eisenstein& w, 
			  map<dedge_t,dedgecoord_t>& reverse, 
			  const IDCounter<Eisenstein> &grid,const IDCounter<Eisenstein> &inner_nodes, const IDCounter<Eisenstein> &outer_nodes)
{
  set<edge_t> edges;

  cout << "outline = " << outline << ";\n";

  //  const Eisenstein neighbours[3] = { Eisenstein(1,0), Eisenstein(1,-1),Eisenstein(0,-1)};
  const Eisenstein neighbours[6] = { Eisenstein(1,0),Eisenstein(0,1),Eisenstein(1,-1),
				     Eisenstein(-1,0),Eisenstein(0,-1),Eisenstein(-1,1) };

  for(IDCounter<Eisenstein>::const_iterator p(inner_nodes.begin()); p!=inner_nodes.end(); p++){
    Eisenstein x(p->first);
    int u = p->second;

    for(int i=0;i<6;i++){
      Eisenstein xn(x+neighbours[i]);
      int v = outer_nodes(xn);
      if(v>=0) edges.insert(edge_t(u+outer_nodes.nextid,v));
    }
  }			 
  return edges;
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

// set<edge_t> connect_outline(const vector< pair<Eisenstein,node_t> >& outline, 
// 			    const IDCounter<Eisenstein> &inner_nodes, const IDCounter<Eisenstein>& outer_nodes)
// {
//   map<dedge_t,dedgecoord_t>& reverse;
//   for(int i=0;i<outline.size();i++){
//     const Eisenstein xu(outline[i].first), xv(outline[(i+1)%outline.size()].first);
//     const node_t u(outline[i].second), v(outline[(i+1)%outline.size()].second);
//     reverse[dedge_t(v,u)] = dedgecoord_t(xu,xv);
//   }

  
  
// }

vector<int> identify_nodes(const vector<pair<Eisenstein,node_t> >& outline, const IDCounter<Eisenstein>& outer_nodes, map<dedge_t,dedgecoord_t>& reverse_arc)
{
  vector<int> same(outer_nodes.size(),-1);
  set<edge_t> same_as;

  for(int i=0;i<outline.size();i++){
    node_t U = outline[i].second,  V = outline[(i+1)%outline.size()].second;
    Eisenstein X0 = outline[i].first, X1 = outline[(i+1)%outline.size()].first;

    dedgecoord_t Xuv(X0,X1), Xvu(reverse_arc[dedge_t(U,V)]);
    Eisenstein x0,x0p,T;

    transform_line(Xuv,Xvu, x0,x0p, T);

    Eisenstein Duv(X1-X0), delta(Duv/abs(gcd(Duv.first,Duv.second)));

    for(Eisenstein x(X0); x!=X1; x += delta){
      Eisenstein xp((x-x0)*T + x0p);
 
      same_as.insert(edge_t(outer_nodes(x),outer_nodes(xp)));;
    }
  }
  
  Graph S(same_as);
  list< list<node_t> > components(S.connected_components());

  cout << "samecomponents = " << components << ";\n";

  for(list< list<node_t> >::const_iterator s(components.begin()); s!=components.end(); s++){
    node_t canonical = *s->begin();

    for(list<node_t>::const_iterator t(s->begin()); t!=s->end(); t++){
      same[*t] = canonical;
    }
  }

  return same;
}

/*
set<edge_t> connect_mathias(const vector<pair<Eisenstein,node_t> >& outline, const IDCounter<Eisenstein> &inner_nodes, const IDCounter<Eisenstein> &outer_nodes, map<dedge_t,dedgecoord_t>& reverse_arc, const Eisenstein& w)
{
  set<edge_t> edges;
  vector<Eisenstein> outline_coords(w*get_keys(outline));
  Eisenstein iw = w.invertn();

  vector< pair<Eisenstein,Eisenstein> > smap;

  polygon p = convert_vector<Eisenstein,pair<int,int> >(outline_coords);
  polygon::inneredge e(p.getInnerEdgePoints());

  reverse( outline_coords.begin(), outline_coords.end() );
  map<Eisenstein,node_t> lavpaenere(outline.begin(),outline.end());
  vector< pair<Eisenstein,Eisenstein> > EC,ECp,O,Ot;



  for(int i=0;i<e.edgePoints.size();i++){
    vector<Eisenstein> xs(e.innerEdgePoints[i].begin(),e.innerEdgePoints[i].end());
    Eisenstein X1 = outline_coords[i], X0 = outline_coords[(i+1)%outline_coords.size()];
    node_t      U = lavpaenere[X0*iw], V = lavpaenere[X1*iw];
    dedgecoord_t Xuv(X0,X1), Xvu(reverse_arc[dedge_t(U,V)]);
    Xvu.first  = Xvu.first *w;
    Xvu.second = Xvu.second*w;

    Eisenstein x0,x0p,T;

    transform_line(Xuv,Xvu, x0,x0p, T);
    //    cout << "Transform " << Xuv << " to " << Xvu << ": (x-" << x0 << " )*"<<T<<" + " << x0p << endl;

    for(int j=0;j<xs.size();j++){
      int d = X1.second < X0.second? -1 : 1;

      Eisenstein x  = xs[j];
      Eisenstein y  = x+Eisenstein(d,0);
      Eisenstein xp = tfm(x,x0,T,x0p);
      Eisenstein yp = tfm(y,x0,T,x0p);

      int u = inner_nodes(x*iw), v = inner_nodes(yp*iw);

      EC.push_back(make_pair(x,y));
      ECp.push_back(make_pair(xp,yp));
      O.push_back(Xuv);
      Ot.push_back(Xvu);
      edges.insert(edge_t(u+outer_nodes.nextid,v+outer_nodes.nextid));
    }
  }

  vector<int> same = identify_nodes(outline,outer_nodes,reverse_arc);
  map< int,set<Eisenstein> > samemap; 
  map< int,set<int> > samenodesmap; 
  for(int i=0;i<same.size();i++){ samemap[same[i]].insert(w*outer_nodes.invert(i));  samenodesmap[same[i]].insert(i); }

  ofstream output("output/connect.m");
  output << "innernodes = " << w*get_keys(inner_nodes) << ";\n"
	 << "outernodes = " << w*get_keys(outer_nodes) << ";\n"
	 << "innerEdgePoints = " << e.innerEdgePoints << ";\n"
	 << "edgePoints = " << e.edgePoints << ";\n"
	 << "outlinecoords = " << outline_coords << ";\n"
	 << "EC = " << EC << ";\n"
	 << "ECp = " << ECp << ";\n"
	 << "OO  = " << O << ";\n"
	 << "samenodes = " << get_values(samenodesmap) << ";\n"
	 << "samecoords = " << get_values(samemap) << ";\n"
	 << "OOt = " << Ot << ";\n";
  output.close();
  return edges;
}
*/



// TODO: Get rid of GC-transform, rename to "fold", and make this work for any polygon from unfolded triangulation?
PlanarGraph fold(vector< pair<Eisenstein, node_t> > &outline)

{
  //  for(int i=0;i<outline.size();i++) outline[i].first = outline[i].first * Eisenstein(0,1);

  // TODO: Split into ConnectInner() and FoldOutline()
  IDCounter<Eisenstein> inner_nodes, outer_nodes;
  vector<Eisenstein> outline_coords(get_keys(outline));
  polygon outline_polygon = convert_vector<Eisenstein,pair<int,int> >(outline_coords);
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
  vector<Eisenstein> edge_coords;
  for(int i=0;i<scans.edge_xs.size();i++){
    const vector<int> &pts(scans.edge_xs[i]);
    for(int j=0;j<pts.size();j++){
      Eisenstein xy(pts[j],i+scans.minY);
      outer_nodes.insert(xy);
      edge_coords.push_back(xy);
    }
  }

  // Register reverse arcs
  map<dedge_t,dedgecoord_t> reverse;
  for(int i=0;i<outline.size();i++){
    const Eisenstein xu(outline[i].first), xv(outline[(i+1)%outline.size()].first);
    const node_t u(outline[i].second), v(outline[(i+1)%outline.size()].second);
    reverse[dedge_t(v,u)] = dedgecoord_t(xu,xv);
  }
  // Register node union
  IDCounter<Eisenstein> grid;
  for(IDCounter<Eisenstein>::const_iterator i(outer_nodes.begin()); i!=outer_nodes.end(); i++)
    grid[i->first] = i->second;

  for(IDCounter<Eisenstein>::const_iterator i(inner_nodes.begin()); i!=inner_nodes.end(); i++)
    grid[i->first] = i->second+outer_nodes.size();
  grid.nextid = grid.size();


  // Now connect internal node grid
  polygon outline_polygonCW = convert_vector<Eisenstein,pair<int,int> >(Eisenstein(1,-1) * outline_coords),
          outline_polygonCCW = convert_vector<Eisenstein,pair<int,int> >(Eisenstein(0,1) * outline_coords);

  polygon::scanline 
    scansCW  = outline_polygonCW.scanConvert(),
    scansCCW = outline_polygonCCW.scanConvert();

  int Nouter = outer_nodes.size();
  vector<edge_t> 
    e   (scanlineedges(inner_nodes,grid,scans,Nouter)), 
    eCW (scanlineedges(inner_nodes,grid,scansCW,Nouter,Eisenstein(0,1))), 
    eCCW(scanlineedges(inner_nodes,grid,scansCCW,Nouter,Eisenstein(1,-1))),
    eoutline(Nouter);
 
  // Connect outline:
  // 1. Sort edge_xs CW
  sort_ccw_eisenstein CCW(Eisenstein::average(edge_coords));
  sort(edge_coords.begin(),edge_coords.end(),CCW);     
  // cout << "ctr3 = " << CCW.centre << ";\n";
  // cout << "outer3 = " << edge_coords << ";\n";
  // vector<double> angles(edge_coords.size());
  // transform(edge_coords.begin(),edge_coords.end(),angles.begin(),CCW);
  // cout << "angles3 = " << angles << ";\n";

  // 2. Traverse outline while connecting edge nodes
  for(int i=0;i<Nouter;i++){
    node_t u(outer_nodes(edge_coords[i])), v(outer_nodes(edge_coords[(i+1)%Nouter]));
    eoutline[i] = edge_t(u,v);
  }

  // Connect boundary nodes to inner nodes and identify identical nodes
  map<edge_t,bool> done;
  set<edge_t> eedge,ecross,ecrossCW,ecrossCCW,ecorner;

  eedge   = connect_edge(outline, Eisenstein(1,0), reverse, grid,inner_nodes,outer_nodes);
  ecorner = connect_corners(edge_coords,outer_nodes);
  ecross    = connect_mathias(outline, inner_nodes,outer_nodes,reverse,Eisenstein(1,0));
  ecrossCW  = connect_mathias(outline, inner_nodes,outer_nodes,reverse,Eisenstein(0,1));
  ecrossCCW = connect_mathias(outline, inner_nodes,outer_nodes,reverse,Eisenstein(1,-1));

  // Join up all the edges
  vector<edge_t> edge_list;
  edge_list.reserve(eoutline.size()+e.size()*3+ecross.size()*3+eedge.size()+ecorner.size());

  copy(e.begin(),e.end(),back_inserter(edge_list));
  copy(eCW.begin(),eCW.end(),back_inserter(edge_list));
  copy(eCCW.begin(),eCCW.end(),back_inserter(edge_list));

  copy(ecross.begin(),   ecross.end(),back_inserter(edge_list));
  copy(ecrossCW.begin(), ecrossCW.end(),back_inserter(edge_list));
  copy(ecrossCCW.begin(),ecrossCCW.end(),back_inserter(edge_list));

  //  copy(eoutline.begin(),eoutline.end(),back_inserter(edge_list));
  copy(eedge.begin(),    eedge.end(),  back_inserter(edge_list));
  copy(ecorner.begin(),ecorner.end(),  back_inserter(edge_list)); // Should be replaced by scanlines of outer nodes, -,CW,CCW.


    
  cout << "edgelist = " << edge_list << ";\n";

  // Identify nodes
  vector<int> same_as = identify_nodes(outline,outer_nodes,reverse);
  
  if(1) for(int i=0;i<edge_list.size();i++){
      node_t u = edge_list[i].first, v = edge_list[i].second;
      if(u<outer_nodes.size()) edge_list[i].first = same_as[u];
      if(v<outer_nodes.size()) edge_list[i].second = same_as[v];
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

      printf("(%d,%d) -> (%d,%d)\n",u,v,newu,newv);

      edges.insert(edge_t(newu,newv));
    }
  } else {
    edges = set<edge_t>(edge_list.begin(),edge_list.end());
    for(int i=0;i<outer_nodes.size()+inner_nodes.size();i++) new_names.insert(i);
  }
  cout << "edgelist = " << edge_list << ";\n";
  cout << "edges = " << edges << ";\n";

  int N = new_names.nextid;
  printf("N = %d (%ld)\n",N,inner_nodes.size()+outer_nodes.size());
  Graph dualGC(edges);

  vector<coord2d> layout2d(N);
  for(IDCounter<Eisenstein>::const_iterator ei(inner_nodes.begin()); ei!=inner_nodes.end(); ei++){
    Eisenstein xy(ei->first);
    node_t      u(new_names(ei->second+outer_nodes.nextid));
    if(u>=0) layout2d[u] = coord2d(xy.first,xy.second);
  }
  for(IDCounter<Eisenstein>::const_iterator ei(outer_nodes.begin()); ei!=outer_nodes.end(); ei++){
    Eisenstein xy(ei->first);
    node_t      u(new_names(ei->second));
    if(u>=0) layout2d[u] = coord2d(xy.first,xy.second);
  }

  for(int i=0;i<N;i++){
    printf("d(%d -> %d) = %ld: ",new_names.invert(i), i,dualGC.neighbours[i].size());
    cout << dualGC.neighbours[i] << endl;
  }



  PlanarGraph PG(dualGC,layout2d);

  return PG;
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

  vector<pair<Eisenstein,node_t> > reduced_outline(GCDreduce(outline));

  output << "reduced = " << fold(reduced_outline) << ";\n";
  output << "gct = " << GCTransform(dual,K,L) << ";\n";

  vector<pair<Eisenstein,node_t> > gct_outline(outline);
  for(int i=0;i<outline.size();i++) gct_outline[i].first = outline[i].first * Eisenstein(K,L);
  output << "gctb = " << fold(gct_outline) << ";\n";



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
