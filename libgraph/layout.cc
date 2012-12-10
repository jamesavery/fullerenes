#include "graph.hh"
#include "cubicgraph.hh"
#include <sstream>

struct ToleranceLess {
  const double tolerance;
  ToleranceLess(const double tolerance) : tolerance(tolerance) {}
  bool operator()(const coord2d& x,const coord2d& y) const { return x<y && (y-x).norm() > tolerance; }
};

vector<coord2d> PlanarGraph::tutte_layout(node_t s, node_t t, node_t r, unsigned int face_max) const
{
  if(!is_cubic())
    printf("tutte_layout called for non-cubic graph. Tutte embedding is only guaranteed planar for cubic graphs.\n");

  if(s<0) s = 0;
  if(t<0){
    //    fprintf(stderr,"t = %d\n",t);
    face_t c(shortest_cycle(s,face_max));
    t = c[1];
    r = c[2];
    //    fprintf(stderr,"s,t,r ~> %d,%d,%d\n",s,t,r);
  } else if(r < 0) {
    //    fprintf(stderr,"r = %d\n",r);
    face_t c(shortest_cycle(s,t,face_max));
    r = c[2];
    //    fprintf(stderr,"s,t,r ~> %d,%d,%d\n",s,t,r);
  }
  //  fprintf(stderr,"tutte_layout(%d,%d,%d)\n",s,t,r);
  outer_face = shortest_cycle(s,t,r,face_max);

  unsigned int Nface = outer_face.size();
  vector<coord2d> outer_coords(Nface);
  for(unsigned int i=0;i<Nface;i++){
    outer_coords[i] = coord2d(sin(i*2*M_PI/double(Nface)),cos(i*2*M_PI/double(Nface)));
  }

  return tutte_layout_direct(outer_face,outer_coords);
}

#include "../contrib/mgmres.hpp"
vector<coord2d> PlanarGraph::tutte_layout_direct(const face_t& outer_face, const vector<coord2d>& outer_coords) const
{
  vector<coord2d> result(N);

  // Construct right hand side: rhs_i = {0,0} everywhere except for outer-face nodes with fixed solutions 
  double rhs[N*2], x[N*2];

  memset(rhs,0,N*2*sizeof(double));
  for(int i=0;i<outer_face.size();i++){
    rhs[    outer_face[i] ] = outer_coords[i].first;
    rhs[N + outer_face[i] ] = outer_coords[i].second;
  }
  memcpy(x,rhs,N*2*sizeof(double));

  // Construct matrix I-1/3*A in triplet form, where A is adjacency, except for rows A_i=\delta_{ij} for i an outer-face node
  double A[4*N];
  int IA[N+1], JA[4*N];
  int nz = 0;
  {
    double Afull[N*N];
    memset(Afull,0,N*N*sizeof(double));
    for(node_t u=0;u<N;u++){
      Afull[u*(N+1)] = 1.0;
      for(int i=0;i<neighbours[u].size();i++){
	const node_t& v(neighbours[u][i]);
	Afull[u*N+v] = -1.0L/3.0;
      }
    }
    for(int i=0;i<outer_face.size();i++){
      const node_t& u(outer_face[i]);
      for(node_t v=0;v<N;v++) Afull[u*N+v] = (u==v)? 1.0 : 0.0;
    }

    for(node_t u=0;u<N;u++){
      IA[u] = nz;
      for(node_t v=0;v<N;v++)
	if(Afull[u*N+v]!=0){
	  A[nz] = Afull[u*N+v];
	  JA[nz] = v;
	  nz++;
	}
    }
  }
  IA[N] = nz;
  // Solve sparse linear system for x-coordinates and y-coordinates 
  pmgmres_ilu_cr(N,nz,IA,JA,A,x,  rhs,  50000,N-1,1e-15,1e-15);
  pmgmres_ilu_cr(N,nz,IA,JA,A,x+N,rhs+N,50000,N-1,1e-15,1e-15);
  //  mgmres_st(N,nz,IA,JA,A,x,  rhs,  50000,N-1,1e-10,1e-10);
  //  mgmres_st(N,nz,IA,JA,A,x+N,rhs+N,50000,N-1,1e-10,1e-10);

  for(node_t u=0;u<N;u++) result[u] = coord2d(x[u],x[u+N]);

  return result;
}

vector<coord2d> PlanarGraph::tutte_layout_iterative(const face_t& outer_face, const vector<coord2d>& outer_coords) const
{
  vector<coord2d> xys(N), newxys(N);
  vector<bool> fixed(N);

  cerr << "tutte_layout: Outer face: " << outer_face << endl;

  unsigned int Nface = outer_face.size();
  for(unsigned int i=0;i<Nface;i++){
    fixed[outer_face[i]] = true;
    xys[outer_face[i]] = outer_coords[i];
  }
    
  bool converged = false;
  const unsigned int TUTTE_MAX_ITERATION = 500000;
  const double TUTTE_CONVERGENCE = 5e-7;
  unsigned int i;
  double max_change;
  for(i=0;!converged && i<TUTTE_MAX_ITERATION; i++){

    max_change = 0;
#pragma omp parallel for
    for(node_t u=0;u<N;u++)
      if(fixed[u]){
	newxys[u] = xys[u];
      } else {
	const vector<node_t>& ns(neighbours[u]);
	coord2d neighbour_sum(0.0);

	for(int i=0;i<ns.size();i++) neighbour_sum += xys[ns[i]];
	newxys[u] = xys[u]*0.15 + (neighbour_sum/ns.size())*0.85;

	// Did the solution converge yet?
	double neighbour_dist = 0;
	for(size_t i=0;i<ns.size();i++) neighbour_dist += (xys[u]-xys[ns[0]]).norm()/ns.size();
	double relative_change = (xys[u]-newxys[u]).norm()/neighbour_dist;
	if(relative_change > max_change) max_change = relative_change;
      }
      
    if(max_change <= TUTTE_CONVERGENCE) converged = true;
    xys = newxys;
  }
  if(i>=TUTTE_MAX_ITERATION){
    printf("Planar Tutte embedding failed to converge. Increase TUTTE_MAX_ITERATION. ");
    abort();
  }
  //  cerr << "Tutte layout of "<<N<<" vertices converged after " << i << " iterations, with maximal relative change " << max_change << endl;
  // Test that points are distinct
  ToleranceLess lt(0.0);
  set<coord2d,ToleranceLess> point_set(xys.begin(),xys.end(),lt);
  if(point_set.size() != N){
    fprintf(stderr,"Tutte layout failed: only %d unique coordinates out of %d vertices (up to tolerance %g).\n",
	    int(point_set.size()),N,0.0);
  }
  return xys;
}

vector<coord2d> PlanarGraph::spherical_projection() const
{
  vector<node_t> outer_face(find_outer_face());
  vector<unsigned int> vertex_depth(multiple_source_shortest_paths(outer_face,vector<bool>(N*(N-1)/2),vector<bool>(N)));

  // Step 1. Sort nodes wrt. vertex depth; partition into dmax sets V[d]. 
  unsigned int dmax = 0;
  map<int,list<node_t> > V;
  for(node_t u=0;u<N;u++){
    V[vertex_depth[u]].push_back(u);
    dmax = max(dmax,vertex_depth[u]);
  }

  // Step 2. Calculate the centroid for vertices grouped by vertex depth.
  vector<coord2d> centroids(dmax+1);
  for(unsigned int d=0;d<=dmax;d++){
    const list<node_t>& Vd(V[d]);
    coord2d c(0);
    for(list<node_t>::const_iterator u(Vd.begin());u!=Vd.end();u++)
      c += layout2d[*u];
    centroids[d] = c/V[d].size();
  }

  // Step 3. Lay out the vertices in order of the distance from the outer face.
  // The angle, when seen from above, is the angle in the flat layout. The
  // center used for calculating the angle is the centroid of the vertices
  // at the same depth.
  double dtheta = M_PI/(dmax+1.0);

  coord2d centroid(centre2d(layout2d));
  vector< coord2d > spherical_layout(N);
  for(unsigned int d=0;d<=dmax;d++){
    double phi = dtheta*(d+0.5);
    const list<node_t>& Vd(V[d]);
    const coord2d& centroid(centroids[d]);

    for(list<node_t>::const_iterator ui(Vd.begin());ui!=Vd.end();ui++){
      const node_t u(*ui);
      coord2d xy(layout2d[u]-centroid);
      double theta = atan2(xy.first,xy.second);

      spherical_layout[u] = coord2d(theta,phi);
    }
  }
  return spherical_layout;
}



// Orient all neighbours clockwise according to a given 2D layout.
void Graph::orient_neighbours(const vector<coord2d>& layout)
{
  for(node_t u=0;u<N;u++){
    vector<node_t>& ns(neighbours[u]);

    sort_ccw_point CW(layout,layout[u]);
    sort(ns.begin(), ns.end(), CW);
  }

}

coord2d Graph::centre2d(const vector<coord2d>& layout) const {
  coord2d centre(0,0);
  for(node_t u=0;u<N;u++) centre += layout[u];
  return centre/N;
}

coord3d Graph::centre3d(const vector<coord3d>& layout) const {
  coord3d centre(0,0,0);
  for(node_t u=0;u<N;u++) centre += layout[u];
  return centre/N;
}
// TODO: 
// * Move layout to member variable
// * Check if layout is planar before allowing it (this function crashes if it is not).


string PlanarGraph::to_latex(double w_cm, double h_cm, bool show_dual, bool number_vertices, bool include_latex_header,
		  int line_colour, int vertex_colour, double line_width, double vertex_diameter) const
{
  ostringstream s;
  s << fixed;
  // If we're outputting a stand-alone LaTeX file, spit out a reasonable header.
  if(include_latex_header)
    s << "\\documentclass{standalone}\n"
      "\\usepackage{tikz}\n"
      "\\begin{document}\n"
      "\\definecolor{vertexcolour}{RGB}{"<<(vertex_colour>>16)<<","<<((vertex_colour>>8)&0xff)<<","<<(vertex_colour&0xff)<<"}\n"
      "\\definecolor{edgecolour}{RGB}{"<<(line_colour>>16)<<","<<((line_colour>>8)&0xff)<<","<<(line_colour&0xff)<<"}\n"
      "\\definecolor{dualvertexcolour}{RGB}{205,79,57}\n"
      "\\definecolor{dualedgecolour}{RGB}{0,0,0}\n"
      "\\tikzstyle{vertex}=[circle, draw, inner sep="<<(number_vertices?"1pt":"0")<<", fill=vertexcolour, minimum width="<<vertex_diameter<<"mm]\n"
      "\\tikzstyle{dualvertex}=[circle, draw, inner sep="<<(number_vertices?"1pt":"0")<<", fill=dualvertexcolour, minimum width="<<vertex_diameter<<"mm]\n"
      "\\tikzstyle{edge}=[draw,color=edgecolour,line width="<<line_width<<"mm]\n"
      "\\tikzstyle{dualedge}=[dotted,draw,color=dualedgecolour,line width="<<line_width<<"mm]\n"
      "\\tikzstyle{invisible}=[draw=none,inner sep=0,fill=none,minimum width=0pt]\n"
      ;

  // Find "internal" width and height of layout and scale to w_cm x h_cm
  coord2d wh(width_height());
  double xscale = w_cm/wh.first;
  double yscale = h_cm/wh.second;


  s << "\\begin{tikzpicture}\n";
  for(node_t u=0;u<N;){
    s << "\\foreach \\place/\\name/\\lbl in {";
    for(node_t u_=0;u_<100 && u<N;u++,u_++){
      const coord2d xs(layout2d[u]*coord2d(xscale,yscale));
    s << "{(" << xs.first << "," << xs.second << ")/v" << u << "/$" << u << "$}" << ((u+1<N && u_+1<100)
? ", ":"}\n\t");
    }
    s << "\\node[vertex] (\\name) at \\place {"<<(number_vertices?"\\lbl":"")<<"};\n";
  }

  
  
  for(set<edge_t>::const_iterator e(edge_set.begin()); e!=edge_set.end();){
    s << "\\foreach \\u/\\v in {";
    for(int i=0;i<100 && e!=edge_set.end();){
      s << "{v"<<e->first<<"/v"<<e->second<<"}";
      if(++e != edge_set.end() && i+1<100) s << ", ";
    }
    s << "}\n\t\\draw[edge] (\\u) -- (\\v);\n";
  }


  if(show_dual){
    PlanarGraph dual(dual_graph());	
    s << "\\foreach \\place/\\name/\\lbl in {";
    for(node_t u=0;u<dual.N;u++){
      const coord2d xs(dual.layout2d[u]*coord2d(xscale,yscale));
      s << "{(" << xs.first << "," << xs.second << ")/v" << u << "/$" << u << "$}" << (u+1<dual.N? ", ":"}\n\t");
    }    
    s << "\\node[dualvertex] (\\name) at \\place {"<<(number_vertices?"\\lbl":"")<<"};\n";
    s << "\\foreach \\u/\\v in {";
    for(set<edge_t>::const_iterator e(dual.edge_set.begin()); e!=dual.edge_set.end();){
      s << "{v"<<e->first<<"/v"<<e->second<<"}";
      if(++e != dual.edge_set.end()) s << ", ";
    }
    s << "}\n\t\\draw[dualedge] (\\u) -- (\\v);\n";
  }

  s<<"\\end{tikzpicture}\n";
  if(include_latex_header)
    s << "\\end{document}\n";

  return s.str();
}

#define byte0(x) ((x)&0xff)
#define byte1(x) (((x)>>8)&0xff)
#define byte2(x) (((x)>>16)&0xff)
#define byte3(x) (((x)>>24)&0xff)

string PlanarGraph::to_povray(double w_cm, double h_cm, 
			      int line_colour, int vertex_colour, double line_width, double vertex_diameter) const
{
  ostringstream s;
  s << fixed;

  s << "#declare Nvertices="<<N<<";\n";
  s << "#declare Nedges="<<edge_set.size()<<";\n";
  s << "#declare edgecolour=color rgb <" << byte2(line_colour)/256. << "," << byte1(line_colour)/256. << "," << byte0(line_colour)/256. << ">;\n";
  s << "#declare nodecolour=color rgb <" << byte2(vertex_colour)/256. << "," << byte1(vertex_colour)/256. << "," << byte0(vertex_colour)/256. << ">;\n";
  s << "#declare edgewidth="<<line_width/10.<<";\n";
  s << "#declare nodediameter="<<vertex_diameter/10.<<";\n";

  if(layout2d.size() == N){
    coord2d wh(width_height());
    double xscale = w_cm/wh.first;
    double yscale = h_cm/wh.second;
    s << "#declare layout2D=array[Nvertices][2]{"; for(int i=0;i<N;i++) s<<layout2d[i]*coord2d(xscale,yscale)<<(i+1<N?",":"}\n\n");
  }
  s << "#declare edges   =array[Nedges][2]{"; 
  for(set<edge_t>::const_iterator e(edge_set.begin());e!=edge_set.end();){ s << *e; s <<(++e != edge_set.end()? ",":"}\n\n"); }

  return s.str();
}

