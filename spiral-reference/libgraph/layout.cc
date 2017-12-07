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
  if(count_edges() == 0)       // empty graph
    return vector<coord2d>(); // -> empty layout
  
  if(s<0) s = 0;
  if(t<0) t = neighbours[s][0];

  if(is_oriented){
    outer_face = get_face_oriented({s,t},face_max);
  } else {
    if(r<0) r = next(t,s);     
    outer_face = shortest_cycle(s,t,r,face_max);
  }

  return tutte_layout(outer_face);
}

vector<coord2d> PlanarGraph::tutte_layout(const face_t& outer_face) const
{
  unsigned int Nface = outer_face.size();
  vector<coord2d> initial_coords(N);
  for(unsigned int i=0;i<Nface;i++){
    // the outer face is layed out CW:
    initial_coords[outer_face[i]] = coord2d(sin(i*2*M_PI/double(Nface)),cos(i*2*M_PI/double(Nface)));
  }

// TODO: There seems to currently be a bug in the direct solver! When did this appear? Need to investigate and fix.
//  initial_coords = tutte_layout_direct(outer_face,initial_coords);
  return tutte_layout_iterative(outer_face,initial_coords);
}


vector<coord2d> PlanarGraph::tutte_layout_iterative(const face_t& outer_face, const vector<coord2d>& initial_coords) const
{
  vector<coord2d> xys(initial_coords.begin(), initial_coords.end()), newxys(N);
  vector<bool> fixed(N);

  //  cerr << "tutte_layout: Outer face: " << outer_face << endl;

  unsigned int Nface = outer_face.size();
  for(unsigned int i=0;i<Nface;i++)
    fixed[outer_face[i]] = true;

  
  bool converged = false;
  const unsigned int TUTTE_MAX_ITERATION = 10000000;
  const double TUTTE_CONVERGENCE = 5e-4;
  unsigned int i;
  double max_change;
  for(i=0;!converged && i<TUTTE_MAX_ITERATION; i++){

    max_change = 0;

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
        for(size_t i=0;i<ns.size();i++) neighbour_dist += (xys[u]-xys[ns[i]]).norm()/ns.size();
        if(neighbour_dist > 0.0){ // let's not divide by zero
          double relative_change = (xys[u]-newxys[u]).norm()/neighbour_dist;
          if(relative_change > max_change) max_change = relative_change;
        }
      }
      
    if(max_change <= TUTTE_CONVERGENCE) converged = true;
    xys = newxys;
  }
  if(i>=TUTTE_MAX_ITERATION){
    printf("Planar Tutte embedding failed to converge. Increase TUTTE_MAX_ITERATION. ");
    cout << "layout = " << xys << ";\n";
    abort();
  }
  //  cerr << "Tutte layout of "<<N<<" vertices converged after " << i << " iterations, with maximal relative change " << max_change << endl;
  // Test that points are distinct
  ToleranceLess lt(0.0);
  set<coord2d,ToleranceLess> point_set(xys.begin(),xys.end(),lt);
  if(point_set.size() != N){
    fprintf(stderr,"Tutte layout failed: only %d unique coordinates out of %d vertices (up to tolerance %g).\n",
	    int(point_set.size()),N,0.0);
    abort();
  }
  return xys;
}

vector<coord2d> PlanarGraph::spherical_projection() const
{
  vector<node_t> outer_face(find_outer_face());
  vector<int> vertex_depth(multiple_source_shortest_paths(outer_face,vector<bool>(N*(N-1)/2),vector<bool>(N)));

  // Step 1. Sort nodes wrt. vertex depth; partition into dmax sets V[d]. 
  int dmax = 0;
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

    sort_ccw_point CCW(layout,layout[u]);
    sort(ns.begin(), ns.end(), CCW);
    reverse(ns.begin(),ns.end()); // reverse to get CW
  }
}

coord2d Graph::centre2d(const vector<coord2d>& layout) const {
  coord2d centre(0,0);
  for(node_t u=0;u<layout.size();u++) centre += layout[u];
  return centre/double(layout.size());
}

coord3d Graph::centre3d(const vector<coord3d>& layout) const {
  coord3d centre(0,0,0);
  for(node_t u=0;u<layout.size();u++) centre += layout[u];
  return centre/double(layout.size());
}


// TODO: 
// * Move layout to member variable
// * Check if layout is planar before allowing it (this function crashes if it is not).
string PlanarGraph::to_latex(double w_cm, double h_cm, bool show_dual, bool print_numbers, bool include_latex_header,
			     int edge_colour, int path_colour, int vertex_colour, double edge_width, double path_width,
			     double vertex_diameter, int Npath, int *path) const
{
  string str;
  ostringstream s(str);
  s << std::fixed;

  if(show_dual && !layout_is_crossingfree()) {
    s << "Get a crossing free layout first.  For example by optimising the layout or using a different algorithm to create it." << endl;
    cerr << "Get a crossing free layout first.  For example by optimising the layout or using a different algorithm to create it." << endl;
    return s.str();
  }

// If we're outputting a stand-alone LaTeX file, spit out a reasonable header.
  if(include_latex_header)
    s << "\\documentclass{standalone}\n"
      "\\usepackage{tikz}\n"
      "\\begin{document}\n"
      "\\definecolor{vertexcolour}{RGB}{"<<(vertex_colour>>16)<<","<<((vertex_colour>>8)&0xff)<<","<<(vertex_colour&0xff)<<"}\n"
      "\\definecolor{edgecolour}{RGB}{"<<(edge_colour>>16)<<","<<((edge_colour>>8)&0xff)<<","<<(edge_colour&0xff)<<"}\n"
      "\\definecolor{pathcolour}{RGB}{"<<(path_colour>>16)<<","<<((path_colour>>8)&0xff)<<","<<(path_colour&0xff)<<"}\n"
      "\\definecolor{dualvertexcolour}{RGB}{205,79,57}\n"
      "\\definecolor{dualedgecolour}{RGB}{0,0,0}\n"
      "\\tikzstyle{vertex}=[circle, draw, inner sep="<<(print_numbers?"1pt":"0pt")<<", fill=vertexcolour, minimum width="<<vertex_diameter<<"mm]\n"
      "\\tikzstyle{dualvertex}=[circle, draw, inner sep="<<(print_numbers?"1pt":"0pt")<<", fill=dualvertexcolour, minimum width="<<vertex_diameter<<"mm]\n"
      "\\tikzstyle{edge}=[draw,color=edgecolour,line width="<<edge_width<<"mm]\n"
      "\\tikzstyle{pth}=[draw,color=pathcolour,line width="<<path_width<<"mm]\n"
      "\\tikzstyle{dualedge}=[dotted,draw,color=dualedgecolour,line width="<<edge_width<<"mm]\n"
      "\\tikzstyle{invisible}=[draw=none,inner sep=0pt,fill=none,minimum width=0pt]\n"
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
      s << "{(" << xs.first << "," << xs.second << ")/v" << u << "/$" << (u+1) << "$}" << ((u+1<N && u_+1<100)
? ", ":"}\n\t");
    }
    s << "\\node[vertex] (\\name) at \\place {"<<(print_numbers?"\\lbl":"")<<"};\n";
  }

  set<edge_t> edge_set = undirected_edges();

  for(set<edge_t>::const_iterator e(edge_set.begin()); e!=edge_set.end();){
    s << "\\foreach \\u/\\v in {";
    for(int i=0;i<100 && e!=edge_set.end();){
      s << "{v"<<e->first<<"/v"<<e->second<<"}";
      if(++e != edge_set.end() && i+1<100) s << ", ";
    }
    s << "}\n\t\\draw[edge] (\\u) -- (\\v);\n";
  }


  // Draw path if non-empty
  if(Npath){
    s << "\\foreach \\u/\\v in {";  
    for(int i=0;i+1<Npath;i++)
      s << "{v"<<path[i]<<"/v"<<path[i+1]<<"}, ";
    s << "{v"<<path[Npath-1]<<"/v"<<path[0]<<"}";
    s << "}\n\t\\draw[pth] (\\u) -- (\\v);\n";
  }    

  if(show_dual){
    PlanarGraph dual(dual_graph(6,true));
    s << "\\foreach \\place/\\name/\\lbl in {";
    for(node_t u=0;u<dual.N;u++){
      const coord2d xs(dual.layout2d[u]*coord2d(xscale,yscale));
      s << "{(" << xs.first << "," << xs.second << ")/v" << u << "/$" << (u+1) << "$}" << (u+1<dual.N? ", ":"}\n\t");
    }    
    s << "\\node[dualvertex] (\\name) at \\place {"<<(print_numbers?"\\lbl":"")<<"};\n";
    s << "\\foreach \\u/\\v in {";

    set<edge_t> dual_edges = dual.undirected_edges();
    for(set<edge_t>::const_iterator e(dual_edges.begin()); e!=dual_edges.end();){
      s << "{v"<<e->first<<"/v"<<e->second<<"}";
      if(++e != dual_edges.end()) s << ", ";
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
			      int edge_colour, int vertex_colour, double edge_width, double vertex_diameter) const
{
  ostringstream s;
  s << fixed;

  set<edge_t> edge_set = undirected_edges();

  s << "#declare Nvertices="<<N<<";\n";
  s << "#declare Nedges="<<edge_set.size()<<";\n";
  s << "#declare edgecolour=color rgb <" << byte2(edge_colour)/256. << "," << byte1(edge_colour)/256. << "," << byte0(edge_colour)/256. << ">;\n";
  s << "#declare nodecolour=color rgb <" << byte2(vertex_colour)/256. << "," << byte1(vertex_colour)/256. << "," << byte0(vertex_colour)/256. << ">;\n";
  s << "#declare edgewidth="<<edge_width/10.<<";\n";
  s << "#declare nodediameter="<<vertex_diameter/10.<<";\n\n";

  vector<int> degrees(N);
  for(node_t u=0;u<N;u++) degrees[u] = neighbours[u].size();
  s << "#declare vertexdegree=array["<<N<<"]" << degrees << ";\n";

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

