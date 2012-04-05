#include "graph.hh"
#include "cubicgraph.hh"
#include <sstream>

struct ToleranceLess {
  const double tolerance;
  ToleranceLess(const double tolerance) : tolerance(tolerance) {}
  bool operator()(const coord2d& x,const coord2d& y) const { return x<y && (y-x).norm() > tolerance; }
};

vector<coord2d> PlanarGraph::tutte_layout(const node_t s, node_t t, node_t r) const
{
  if(!t) t = neighbours[s][0];
  if(!r) {
	r = neighbours[t][0];
	for(int i=1;i<neighbours[t].size();i++) if(r==s) r = neighbours[t][i];
  }
  face_t outer_face(shortest_cycle(s,t,r,6));
  vector<coord2d> xys(N), newxys(N);
  vector<bool> fixed(N);

  cout << "tutte_layout: Outer face: " << outer_face << endl;

  unsigned int Nface = outer_face.size();
  for(unsigned int i=0;i<Nface;i++){
    fixed[outer_face[i]] = true;
    xys[outer_face[i]] = coord2d(sin(i*2*M_PI/double(Nface)),cos(i*2*M_PI/double(Nface)));
  }
    
  bool converged = false;
  const unsigned int TUTTE_MAX_ITERATION = 50000;
  const double TUTTE_CONVERGENCE = 5e-6;
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
	newxys[u] = xys[u]*0.2 + (neighbour_sum/ns.size())*0.8;

	// Did the solution converge yet?
	double neighbour_dist = 0;
	for(size_t i=0;i<ns.size();i++) neighbour_dist += (xys[u]-xys[ns[0]]).norm()/ns.size();
	double relative_change = (xys[u]-newxys[u]).norm()/neighbour_dist;
	if(relative_change > max_change) max_change = relative_change;
      }
      
    if(max_change <= TUTTE_CONVERGENCE) converged = true;
    xys = newxys;
  }
  cout << "Tutte layout converged after " << i << " iterations, with maximal relative change " << max_change << endl;
  // Test that points are distinct
  ToleranceLess lt(0.0);
  set<coord2d,ToleranceLess> point_set(xys.begin(),xys.end(),lt);
  if(point_set.size() != N){
    fprintf(stderr,"Tutte layout failed: only %d unique coordinates out of %d vertices (up to tolerance %g).\n",
	    int(point_set.size()),N,0.0);
  }
  return xys;
}

vector<coord2d> PlanarGraph::spherical_projection(const vector< coord2d >& layout2d) const
{
  vector<node_t> outer_face(shortest_cycle(0,neighbours[0][0],6));
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
facemap_t PlanarGraph::compute_faces_oriented() const 
{
  facemap_t facemap;

  set<dedge_t> workset;
  for(set<edge_t>::const_iterator e(edge_set.begin()); e!= edge_set.end(); e++){
    const node_t s = e->first, t = e->second;
    workset.insert(dedge_t(s,t));
    workset.insert(dedge_t(t,s));
  }

  // Outer face must exist and be ordered CW, rest of faces CCW
  assert(outer_face.size() > 0);

  cerr << "Outer face: " << outer_face << endl;

  coord2d centre(centre2d(layout2d));

  // Add outer face to output, remove directed edges from work set
  facemap[outer_face.size()].insert(outer_face);
  for(unsigned int i=0;i<outer_face.size();i++){
    const node_t u = outer_face[i], v = outer_face[(i+1)%outer_face.size()];
    //    printf("Removing directed edge (%d,%d)\n",u,v);
    workset.erase(dedge_t(u,v));
  }

  // Now visit every other edge once in each direction.
  while(!workset.empty()){
    dedge_t e = *workset.begin(); workset.erase(workset.begin());

    face_t face;
    face.push_back(e.first);
    face.push_back(e.second);

    //    printf("%d->%d\n",e.first,e.second);
    while(e.second != face[0]){
      const node_t u = e.first, v = e.second;
      const vector<node_t>& ns(neighbours[v]);

      coord2d vu(layout2d[u]-layout2d[v]);
      double angle_min = -M_PI;

      node_t w=0;
      for(unsigned int i=0;i<ns.size();i++) {
	//	printf("%d : %d (%d->%d) angle %g\n",i,ns[i],u,v,vu.line_angle(layout[ns[i]]-layout[v]));
	if(ns[i] != u) { // Find and use first unvisited edge in order of angle to u->v
	  set<dedge_t>::iterator ei(workset.find(dedge_t(v,ns[i])));

	  if(ei != workset.end()){ // directed edge is not yet visited
	    coord2d vw(layout2d[ns[i]]-layout2d[v]);
	    double angle = vu.line_angle(vw);

	    if(angle>= angle_min){
	      angle_min = angle;
	      w = ns[i];
	    } // else {
	    //   fprintf(stderr,"\t[%d->%d already used.]\n",v,ns[i]);
	    // }
	  }
	}
      }
      e = dedge_t(v,w);
      workset.erase(e);

      if(e.second != face[0]) face.push_back(e.second);
    }
    //    cout << "face = " << face << endl;
    facemap[face.size()].insert(face);
  }
  return facemap;
}


string PlanarGraph::to_latex(double w_cm, double h_cm, bool show_dual, bool number_vertices, bool include_latex_header) const 
{
  ostringstream s;
  s << fixed;
  // If we're outputting a stand-alone LaTeX file, spit out a reasonable header.
  if(include_latex_header)
    s << "\\documentclass{article}\n"
         "\\usepackage{fullpage,fourier,tikz}\n"
         "\\begin{document}\n"
      "\\tikzstyle{vertex}=[circle, draw, inner sep="<<(number_vertices?"1pt":"0")<<", fill=blue!20, minimum width=4pt]\n"
      "\\tikzstyle{dualvertex}=[circle, draw, inner sep="<<(number_vertices?"1pt":"0")<<", fill=red!40, minimum width=2pt]\n"
      "\\tikzstyle{invisible}=[draw=none,inner sep=0,fill=none,minimum width=0pt]\n"
      "\\tikzstyle{dualedge}=[dotted,draw]\n"
      "\\tikzstyle{edge}=[draw]\n"
      ;

  // Find "internal" width and height of layout and scale to w_cm x h_cm
  double xmin=INFINITY,xmax=-INFINITY,ymin=INFINITY,ymax=-INFINITY, xscale, yscale;
  for(node_t u=0;u<N;u++){
    double x = layout2d[u].first, y = layout2d[u].second;
    if(x<xmin) xmin = x;
    if(x>xmax) xmax = x;
    if(y<ymin) ymin = y;
    if(y>ymax) ymax = y;
  }
  xscale = w_cm/(xmax-xmin);
  yscale = h_cm/(ymax-ymin);


  s << "\\begin{tikzpicture}[xscale="<<xscale<<",yscale="<<yscale<<"]\n";
  s << "\\foreach \\place/\\name/\\lbl in {";
  for(node_t u=0;u<N;u++){
    const coord2d& xs(layout2d[u]);
    s << "{(" << xs.first << "," << xs.second << ")/v" << u << "/$" << u << "$}" << (u+1<N? ", ":"}\n\t");
  }
  s << "\\node[vertex] (\\name) at \\place {"<<(number_vertices?"\\lbl":"")<<"};\n";
  s << "\\foreach \\u/\\v in {";
  for(set<edge_t>::const_iterator e(edge_set.begin()); e!=edge_set.end();){
    s << "{v"<<e->first<<"/v"<<e->second<<"}";
    if(++e != edge_set.end()) s << ", ";
  }
  s << "}\n\t\\draw[edge] (\\u) -- (\\v);\n";

  if(show_dual){
    PlanarGraph dual(dual_graph(6));	// TODO: This breaks for everything else than fullerenes
    s << "\\foreach \\place/\\name/\\lbl in {";
    for(node_t u=0;u<dual.N;u++){
      const coord2d& xs(dual.layout2d[u]);
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


