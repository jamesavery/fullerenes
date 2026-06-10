#include <sstream>
#include <math.h>
#include "fullerenes/layout2d.hh"
#include "fullerenes/cubicgraph.hh"

#ifdef HAS_GSL
#include <gsl/gsl_vector.h>
#include <gsl/gsl_multimin.h>
#endif

using namespace std;

namespace layout2d {

// Orient a planar graph by computing a Tutte embedding and sorting
// neighbours by angle. For 3-connected planar graphs (fullerenes),
// this produces a consistent planar embedding.
//
// Returns true if the graph is planar, false otherwise.
bool planar_orient(GraphView& G)
{
  const int N = G.N;
  if(N == 0) return true;

  // Find a face of the graph using BFS shortest cycle through an edge.
  // For 3-connected planar graphs, the shortest cycle through any edge
  // is a face boundary.
  node_t s = 0, t = G.nbrs(s)[0];
  face_t face = G.shortest_cycle({s, t}, N);
  if(face.empty()) return false;

  // Use this face as the outer boundary for Tutte layout
  PlanarGraph PG(G);
  vector<coord2d> layout = PG.tutte_layout(face);
  if(layout.empty()) return false;

  // Orient neighbours by CCW angle sort
  orient_neighbours(G, layout);

  return G.is_consistently_oriented();
}

void orient_neighbours(GraphView& G, const vector<coord2d>& layout)
{
  for(node_t u=0;u<G.N;u++){
    auto ns = G[u];
    sort_ccw_point CCW(layout,layout[u]);
    sort(ns.begin(), ns.end(), CCW);
    // CCW sort matches the codebase convention:
    // neighbours are ordered counterclockwise in the 2D embedding.
  }
}

bool layout_is_crossingfree(const PlanarGraphView& G, const vector<coord2d>& layout)
{
  assert(layout.size() == G.N);
  vector<edge_t> es = G.undirected_edges();
  for (edge_t e1: es){
    for (edge_t e2: es){
      if (e1.first == e2.first || e1.second == e2.first || e1.first == e2.second || e1.second == e2.second) continue;
      const double e1ax = layout[e1.first].first,
                   e1ay = layout[e1.first].second,
                   e1bx = layout[e1.second].first,
                   e1by = layout[e1.second].second,
                   e2ax = layout[e2.first].first,
                   e2ay = layout[e2.first].second,
                   e2bx = layout[e2.second].first,
                   e2by = layout[e2.second].second;
      const double a1 = (e1ay - e1by)/(e1ax - e1bx);
      const double b1 = e1ay - a1 * e1ax;
      if ((e2ay > a1*e2ax+b1 && e2by > a1*e2bx+b1) || (e2ay < a1*e2ax+b1 && e2by < a1*e2bx+b1)) continue;
      const double a2 = (e2ay - e2by)/(e2ax - e2bx);
      const double b2 = e2ay - a2 * e2ax;
      if ((e1ay > a2*e1ax+b2 && e1by > a2*e1bx+b2) || (e1ay < a2*e1ax+b2 && e1by < a2*e1bx+b2)) continue;
      cerr << "edges " << e1 << " and " << e2 << " intersect." << endl;
      return false;
    }
  }
  return true;
}


face_t find_outer_face(const PlanarGraphView& G, const vector<coord2d>& layout)
{
  assert(layout.size() == G.N);

  vector<double> radii(G.N);

  node_t u_farthest = 0;
  double rmax = 0;
  for(node_t u=0;u<G.N;u++){
    radii[u] = layout[u].norm();
    if(radii[u] > rmax){ rmax = radii[u]; u_farthest = u; }
  }

  face_t outer_face;
  int i = 0;
  for(node_t t = u_farthest, u = u_farthest, v = -1; v != u_farthest && i <= G.N; i++){
    auto ns = G.nbrs(u);
    double r = 0;
    for(int i=0;i<ns.size();i++)
      if(ns[i] != t && ns[i] != u && radii[ns[i]] > r){ r = radii[ns[i]]; v = ns[i]; }
    outer_face.push_back(u);
    t = u;
    u = v;
  }

  assert(i<G.N);

  sort_ccw_point CCW(layout,outer_face.centroid(layout));
  sort(outer_face.begin(),outer_face.end(),CCW); // sort CCW
  reverse(outer_face.begin(),outer_face.end()); // reverse to get CW

  return outer_face;
}

vector<double> edge_lengths(const PlanarGraphView& G, const vector<coord2d>& layout)
{
  assert(layout.size() == G.N);
  vector<edge_t> edges = G.undirected_edges();

  vector<double> lengths(edges.size());

  for(int i=0;i<edges.size();i++)
    lengths[i] = (layout[edges[i].first]-layout[edges[i].second]).norm();

  return lengths;
}

coord2d width_height(const vector<coord2d>& layout)
{
  double xmin=INFINITY,xmax=-INFINITY,ymin=INFINITY,ymax=-INFINITY;
  for(size_t u=0;u<layout.size();u++){
    double x = layout[u].first, y = layout[u].second;
    if(x<xmin) xmin = x;
    if(x>xmax) xmax = x;
    if(y<ymin) ymin = y;
    if(y>ymax) ymax = y;
  }
  return coord2d(xmax-xmin,ymax-ymin);
}

void scale(vector<coord2d>& layout, const coord2d& x) {
  for(size_t u=0;u<layout.size();u++) layout[u] *= x;
}

void move(vector<coord2d>& layout, const coord2d& x) {
  for(size_t u=0;u<layout.size();u++) layout[u] += x;
}


vector<coord2d> spherical_projection(const PlanarGraphView& G, const vector<coord2d>& layout)
{
  face_t outer_face_nodes = find_outer_face(G, layout);
  vector<int> vertex_depth(G.multiple_source_shortest_paths(outer_face_nodes));

  // Step 1. Sort nodes wrt. vertex depth; partition into dmax sets V[d].
  int dmax = 0;
  map<int,list<node_t> > V;
  for(node_t u=0;u<G.N;u++){
    V[vertex_depth[u]].push_back(u);
    dmax = max(dmax,vertex_depth[u]);
  }

  // Step 2. Calculate the centroid for vertices grouped by vertex depth.
  vector<coord2d> centroids(dmax+1);
  for(unsigned int d=0;d<=dmax;d++){
    coord2d c = {0,0};

    for(node_t u: V[d])  c += layout[u];
    centroids[d] = c/V[d].size();
  }

  // Step 3. Lay out the vertices in order of the distance from the outer face.
  double dtheta = M_PI/(dmax+1.0);

  vector< coord2d > spherical_layout(G.N);
  for(unsigned int d=0;d<=dmax;d++){
    double phi = dtheta*(d+0.5);

    for(node_t u: V[d]){
      coord2d xy(layout[u]-centroids[d]);
      double theta = atan2(xy.first,xy.second);

      spherical_layout[u] = coord2d(theta,phi);
    }
  }
  return spherical_layout;
}


string to_latex(const PlanarGraphView& G, const vector<coord2d>& layout,
                double w_cm, double h_cm, bool show_dual, bool print_numbers, bool include_latex_header,
                int edge_colour, int path_colour, int vertex_colour,
                double edge_width, double path_width, double vertex_diameter,
                int Npath, int *path)
{
  string str;
  ostringstream s(str);
  s << std::fixed;

  if(show_dual && !layout_is_crossingfree(G, layout)) {
    s << "Get a crossing free layout first.  For example by optimising the layout or using a different algorithm to create it." << endl;
    cerr << "Get a crossing free layout first.  For example by optimising the layout or using a different algorithm to create it." << endl;
    return s.str();
  }

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

  coord2d wh(width_height(layout));
  double xscale = w_cm/wh.first;
  double yscale = h_cm/wh.second;

  s << "\\begin{tikzpicture}\n";
  for(node_t u=0;u<G.N;){
    s << "\\foreach \\place/\\name/\\lbl in {";
    for(node_t u_=0;u_<100 && u<G.N;u++,u_++){
      const coord2d xs(layout[u]*coord2d(xscale,yscale));
      s << "{(" << xs.first << "," << xs.second << ")/v" << u << "/$" << (u+1) << "$}" << ((u+1<G.N && u_+1<100)
? ", ":"}\n\t");
    }
    s << "\\node[vertex] (\\name) at \\place {"<<(print_numbers?"\\lbl":"")<<"};\n";
  }

  vector<edge_t> edges = G.undirected_edges();

  if(!edges.empty()){
    s << "\\foreach \\u/\\v in {";
    for(int i=0;i<edges.size();i++){
      s << "{v"<<edges[i].first<<"/v"<<edges[i].second<<"}"<<(i+1<edges.size()?", ":"");
    }
    s << "}\n\t\\draw[edge] (\\u) -- (\\v);\n";
  }

  if(Npath){
    s << "\\foreach \\u/\\v in {";
    for(int i=0;i+1<Npath;i++)
      s << "{v"<<path[i]<<"/v"<<path[i+1]<<"}, ";
    s << "{v"<<path[Npath-1]<<"/v"<<path[0]<<"}";
    s << "}\n\t\\draw[pth] (\\u) -- (\\v);\n";
  }

  if(show_dual){
    PlanarGraph dual(G.dual_graph(6));
    vector<coord2d> dual_layout = dual.tutte_layout();
    s << "\\foreach \\place/\\name/\\lbl in {";
    for(node_t u=0;u<dual.N;u++){
      const coord2d xs(dual_layout[u]*coord2d(xscale,yscale));
      s << "{(" << xs.first << "," << xs.second << ")/v" << u << "/$" << (u+1) << "$}" << (u+1<dual.N? ", ":"}\n\t");
    }
    s << "\\node[dualvertex] (\\name) at \\place {"<<(print_numbers?"\\lbl":"")<<"};\n";
    s << "\\foreach \\u/\\v in {";

    vector<edge_t> dual_edges = dual.undirected_edges();
    for(int i=0;i<dual_edges.size();i++){
      edge_t e = dual_edges[i];
      s << "{v"<<e.first<<"/v"<<e.second<<"}" << (i+1<dual_edges.size()?", ":"");
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

string to_povray(const PlanarGraphView& G, const vector<coord2d>& layout,
                 double w_cm, double h_cm,
                 int edge_colour, int vertex_colour, double edge_width, double vertex_diameter)
{
  ostringstream s;
  s << fixed;

  vector<edge_t> edges = G.undirected_edges();

  s << "#declare Nvertices="<<G.N<<";\n";
  s << "#declare Nedges="<<edges.size()<<";\n";
  s << "#declare edgecolour=color rgb <" << byte2(edge_colour)/256. << "," << byte1(edge_colour)/256. << "," << byte0(edge_colour)/256. << ">;\n";
  s << "#declare nodecolour=color rgb <" << byte2(vertex_colour)/256. << "," << byte1(vertex_colour)/256. << "," << byte0(vertex_colour)/256. << ">;\n";
  s << "#declare edgewidth="<<edge_width/10.<<";\n";
  s << "#declare nodediameter="<<vertex_diameter/10.<<";\n\n";

  vector<int> degrees(G.N);
  for(node_t u=0;u<G.N;u++) degrees[u] = G.degree(u);
  s << "#declare vertexdegree=array["<<G.N<<"]" << degrees << ";\n";

  if(layout.size() == G.N){
    coord2d wh(width_height(layout));
    double xscale = w_cm/wh.first;
    double yscale = h_cm/wh.second;
    s << "#declare layout2D=array[Nvertices][2]{"; for(int i=0;i<G.N;i++) s<<layout[i]*coord2d(xscale,yscale)<<(i+1<G.N?",":"}\n\n");
  }
  s << "#declare edges   =array[Nedges][2]{";
  for(int i=0;i<edges.size();i++){ s << edges[i]; s <<(i+1<edges.size()? ",":"}\n\n"); }

  return s.str();
}

#undef byte0
#undef byte1
#undef byte2


// ==================== optimize_layout (from layout-optimize.cc) ====================

#ifdef HAS_GSL

struct optlayout_params_t
{
  PlanarGraphView *graph;
  face_t outer_face;
  vector<double> *zero_values_dist;
  vector<double> *k_dist;
  vector<double> *k_angle;
  vector<double> *k_area;
};


static double optlayout_pot(const gsl_vector* coordinates, void* parameters)
{
  optlayout_params_t &params = *static_cast<optlayout_params_t*>(parameters);
  PlanarGraphView &graph = *params.graph;
  const face_t &outer_face = params.outer_face;
  vector<double> &zero_values_dist = *params.zero_values_dist;
  vector<double> &k_dist = *params.k_dist;
  vector<double> &k_angle = *params.k_angle;
  vector<double> &k_area = *params.k_area;

  vector<edge_t> edge_set = graph.undirected_edges();
  const int n_faces = 2 + edge_set.size() - graph.N;

  assert(zero_values_dist.size() == edge_set.size());
  assert(k_dist.size() == edge_set.size());
  assert(k_angle.size() == 2*edge_set.size());
  assert(k_area.size() == n_faces);

  // DISTANCE TERM
  double log_sum_edge_length=0;

  for(const edge_t &e: edge_set){
    const double ax = gsl_vector_get(coordinates, 2 * e.first);
    const double ay = gsl_vector_get(coordinates, 2 * e.first +1);
    const double bx = gsl_vector_get(coordinates, 2 * e.second);
    const double by = gsl_vector_get(coordinates, 2 * e.second +1);
    log_sum_edge_length += log(coord2d(ax-bx,ay-by).norm());
  }
  log_sum_edge_length /= edge_set.size();

  double potential_energy = 0.0;
  int i=0;
  for(const edge_t &e: edge_set){
    vector<node_t>::const_iterator it1, it2;
    it1 = find (outer_face.begin(), outer_face.end(), e.first);
    it2 = find (outer_face.begin(), outer_face.end(), e.second);
    if (it1 != outer_face.end() && it2 != outer_face.end() && ( it1 == it2+1 || it1 == it2-1 || (it1 == outer_face.begin() && it2 == outer_face.end()-1) || (it1 == outer_face.end()-1 && it2 == outer_face.begin()))){
      continue; // edge is part of outer face
    }

    const double ax = gsl_vector_get(coordinates, 2 * e.first);
    const double ay = gsl_vector_get(coordinates, 2 * e.first +1);
    const double bx = gsl_vector_get(coordinates, 2 * e.second);
    const double by = gsl_vector_get(coordinates, 2 * e.second +1);
    potential_energy += 0.5 * k_dist[i] * pow(coord2d(ax-bx,ay-by).norm() - zero_values_dist[i], 2);
    i++;
  }

  // ANGLE TERM
  vector<face_t> faces(graph.compute_faces());
  for(const face_t &f: faces)
    for(int i=0;i<f.size();i++){
      int d = f.size();

      const double ax = gsl_vector_get(coordinates, 2*f[(i+d-1) % d]  );
      const double ay = gsl_vector_get(coordinates, 2*f[(i+d-1) % d]+1);
      const double bx = gsl_vector_get(coordinates, 2*f[i]  );
      const double by = gsl_vector_get(coordinates, 2*f[i]+1);
      const double cx = gsl_vector_get(coordinates, 2*f[(i+1) % d]  );
      const double cy = gsl_vector_get(coordinates, 2*f[(i+1) % d]+1);

      const double angle_beta = coord3d::angle(coord3d(ax,ay,0) - coord3d(bx,by,0),
                                               coord3d(cx,cy,0) - coord3d(bx,by,0));

      potential_energy += 0.5 * k_angle[f[i]] * pow(angle_beta - M_PI*(1.0-2.0/d),2);
    }

  // AREA TERM
  double A_tot=0;
  const double bx = gsl_vector_get(coordinates, 2* outer_face[0]);
  const double by = gsl_vector_get(coordinates, 2* outer_face[0] +1);
  for (int i=1; i<outer_face.size()-1; ++i){
    const double ax = gsl_vector_get(coordinates, 2* outer_face[i]  );
    const double ay = gsl_vector_get(coordinates, 2* outer_face[i]+1);
    const double cx = gsl_vector_get(coordinates, 2* outer_face[i+1]  );
    const double cy = gsl_vector_get(coordinates, 2* outer_face[i+1]+1);

    A_tot += ((ax-bx)*(cy-by) - (ay-by)*(cx-bx))/2;
  }

  const double A_av = abs(A_tot)/(n_faces-1);

  for(int i=0;i<faces.size();i++){
    const face_t &f(faces[i]);
    double A=0;
    const double bx = gsl_vector_get(coordinates, 2*f[0]);
    const double by = gsl_vector_get(coordinates, 2*f[0]+1);
    for (int j=1; j+1<f.size();j++){
      const double ax = gsl_vector_get(coordinates, 2*f[j]  );
      const double ay = gsl_vector_get(coordinates, 2*f[j]+1);
      const double cx = gsl_vector_get(coordinates, 2*f[j+1]  );
      const double cy = gsl_vector_get(coordinates, 2*f[j+1]+1);

      A += ((ax-bx)*(cy-by) - (ay-by)*(cx-bx))/2;
    }
    potential_energy += 0.5 * k_area[i] * (abs(A) - A_av)*(abs(A)-A_av);
  }

  return potential_energy;
}


static void optlayout_grad(const gsl_vector* coordinates, void* parameters, gsl_vector* gradient)
{
  optlayout_params_t &params = *static_cast<optlayout_params_t*>(parameters);
  PlanarGraphView &graph = *params.graph;
  const face_t &outer_face = params.outer_face;
  vector<double> &zero_values_dist = *params.zero_values_dist;
  vector<double> &k_dist = *params.k_dist;
  vector<double> &k_angle = *params.k_angle;
  vector<double> &k_area = *params.k_area;

  vector<edge_t> edge_set = graph.undirected_edges();
  const int n_faces = 2+edge_set.size() - graph.N;

  assert(zero_values_dist.size() == edge_set.size());
  assert(k_dist.size() == edge_set.size());
  assert(k_angle.size() == 2*edge_set.size());
  assert(k_area.size() == n_faces);

  vector<coord2d> derivatives(graph.N);

  double log_sum_edge_length=0;
  for(const edge_t &e: edge_set){
    const double ax = gsl_vector_get(coordinates, 2 * e.first);
    const double ay = gsl_vector_get(coordinates, 2 * e.first +1);
    const double bx = gsl_vector_get(coordinates, 2 * e.second);
    const double by = gsl_vector_get(coordinates, 2 * e.second +1);
    log_sum_edge_length += log(coord2d(ax-bx,ay-by).norm());
  }
  log_sum_edge_length /= edge_set.size();

  int i=0;
  for(const edge_t &e: edge_set){
    const double ax = gsl_vector_get(coordinates, 2 * e.first);
    const double ay = gsl_vector_get(coordinates, 2 * e.first +1);
    const double bx = gsl_vector_get(coordinates, 2 * e.second);
    const double by = gsl_vector_get(coordinates, 2 * e.second +1);
    derivatives[e.first]  += coord2d::dnorm(coord2d(ax-bx,ay-by))
                           * k_dist[i] * (coord2d(ax-bx,ay-by).norm() - zero_values_dist[i]);
    derivatives[e.second] -= coord2d::dnorm(coord2d(ax-bx,ay-by)) * k_dist[i]
                           * (coord2d(ax-bx,ay-by).norm() - zero_values_dist[i]);
    i++;
  }

  // ANGLE TERM
  vector<face_t> faces(graph.compute_faces());
  for(const face_t &f: faces)
    for (int i=0; i<f.size(); i++){
      const int d = f.size();

      const double ax = gsl_vector_get(coordinates, 2*f[(i+d-1) % d]  );
      const double ay = gsl_vector_get(coordinates, 2*f[(i+d-1) % d]+1);
      const double bx = gsl_vector_get(coordinates, 2*f[i]  );
      const double by = gsl_vector_get(coordinates, 2*f[i]+1);
      const double cx = gsl_vector_get(coordinates, 2*f[(i+1) % d]  );
      const double cy = gsl_vector_get(coordinates, 2*f[(i+1) % d]+1);

      coord3d a(coord3d(ax,ay,0) - coord3d(bx,by,0)), c(coord3d(cx,cy,0) - coord3d(bx,by,0)), da, dc;
      coord3d::dangle(a, c, da, dc);

      const double angle_beta = coord3d::angle(coord3d(ax,ay,0) - coord3d(bx,by,0),
                                               coord3d(cx,cy,0) - coord3d(bx,by,0));

      derivatives[f[(i+d-1) % d]] +=  coord2d(da[0],da[1])*(angle_beta-M_PI*(1.0-2.0/d)) * k_angle[f[i]];
      derivatives[f[i]] += -coord2d((da+dc)[0],(da+dc)[1])*(angle_beta-M_PI*(1.0-2.0/d)) * k_angle[f[i]];
      derivatives[f[(i+1) % d]]   += coord2d(dc[0],dc[1]) *(angle_beta-M_PI*(1.0-2.0/d)) * k_angle[f[i]];
    }

  // AREA TERM
  double A_tot=0;
  const double bx = gsl_vector_get(coordinates, 2* outer_face[0]);
  const double by = gsl_vector_get(coordinates, 2* outer_face[0] +1);
  for (int i=1; i+1<outer_face.size(); ++i){
    const double ax = gsl_vector_get(coordinates, 2* outer_face[i  ]  );
    const double ay = gsl_vector_get(coordinates, 2* outer_face[i  ]+1);
    const double cx = gsl_vector_get(coordinates, 2* outer_face[i+1]  );
    const double cy = gsl_vector_get(coordinates, 2* outer_face[i+1]+1);

    A_tot += ((ax-bx)*(cy-by) - (ay-by)*(cx-bx))/2;
  }
  const double A_av = abs(A_tot)/(n_faces-1);

  i=0;
  for(const face_t &f: faces){
      double A=0;
      int d = f.size();

      const double bx = gsl_vector_get(coordinates, 2*f[0]  );
      const double by = gsl_vector_get(coordinates, 2*f[0]+1);
      for (int j=1; j+1<f.size(); j++){
        const double ax = gsl_vector_get(coordinates, 2*f[j  ]  );
        const double ay = gsl_vector_get(coordinates, 2*f[j  ]+1);
        const double cx = gsl_vector_get(coordinates, 2*f[j+1]  );
        const double cy = gsl_vector_get(coordinates, 2*f[j+1]+1);

        A += ((ax-bx)*(cy-by) - (ay-by)*(cx-bx))/2;
      }
      const double sign = A/abs(A);

      for (int j=0; j<f.size(); j++){
        const double ax = gsl_vector_get(coordinates, 2*f[(j+d-1) % d]);
        const double ay = gsl_vector_get(coordinates, 2*f[(j+d-1) % d] +1);
        const double cx = gsl_vector_get(coordinates, 2*f[(j+1)   % d]);
        const double cy = gsl_vector_get(coordinates, 2*f[(j+1)   % d] +1);

        derivatives[f[j]] += coord2d(cy-ay, ax-cx)/2 * k_area[i] * (abs(A) - A_av) * sign;
      }
      i++;
  }

  // fix outer face
  for(vector<node_t>::iterator it = outer_face.begin(); it != outer_face.end(); ++it){
    derivatives[*it].first = 0;
    derivatives[*it].second = 0;
  }
  for(int i=0; i < graph.N; ++i) {
    gsl_vector_set(gradient, 2*i, derivatives[i].first);
    gsl_vector_set(gradient, 2*i+1, derivatives[i].second);
  }
}


static void optlayout_pot_grad(const gsl_vector* coordinates, void* parameters, double* potential, gsl_vector* gradient)
{
  *potential = optlayout_pot(coordinates, parameters);
  optlayout_grad(coordinates, parameters, gradient);
}


bool optimize_layout(PlanarGraphView& G, vector<coord2d>& layout,
                     double zv_dist_inp, double k_dist_inp, double k_angle_inp, double k_area_inp)
{
  if(layout.size()!=G.N){
    layout = G.tutte_layout();
  }

  const double stepsize = 1e-2;
  const double terminate_gradient = 1e-2;
  const double tol = 1e-1;
  const int max_iterations = 5000;

  vector<edge_t> edge_set = G.undirected_edges();
  const int n_faces = 2 + edge_set.size() - G.N;

  vector<double> zero_values_dist(edge_set.size(), zv_dist_inp);
  vector<double> k_dist(edge_set.size(), k_dist_inp);
  vector<double> k_angle(3*G.N, k_angle_inp);
  vector<double> k_area(n_faces, k_area_inp);

  optlayout_params_t params;
  params.graph = &G;
  params.outer_face = find_outer_face(G, layout);
  params.zero_values_dist = &zero_values_dist;
  params.k_dist = &k_dist;
  params.k_angle = &k_angle;
  params.k_area = &k_area;

  gsl_multimin_function_fdf potential_function;
  potential_function.n = 2*G.N;
  potential_function.f = &optlayout_pot;
  potential_function.df = &optlayout_grad;
  potential_function.fdf = &optlayout_pot_grad;
  potential_function.params = static_cast<void*>(&params);

  gsl_vector* coordinates = gsl_vector_alloc(potential_function.n);
  for(int i=0; i<G.N; ++i){
    gsl_vector_set(coordinates, 2*i, layout[i].first);
    gsl_vector_set(coordinates, 2*i+1, layout[i].second);
  }

  const gsl_multimin_fdfminimizer_type *fT = gsl_multimin_fdfminimizer_conjugate_fr;
  gsl_multimin_fdfminimizer *s = gsl_multimin_fdfminimizer_alloc(fT, potential_function.n);
  gsl_multimin_fdfminimizer_set(s, &potential_function, coordinates, stepsize, tol);
  size_t iter = 0;
  int status=0;
  do {
    ++iter;
    status = gsl_multimin_fdfminimizer_iterate(s);
    if(status) break;
    status = gsl_multimin_test_gradient(s->gradient, terminate_gradient);
  } while (status == GSL_CONTINUE && iter < max_iterations);

  for(int i=0; i<G.N; ++i){
    layout[i].first = gsl_vector_get(s->x, 2*i);
    layout[i].second = gsl_vector_get(s->x, 2*i+1);
  }

  gsl_multimin_fdfminimizer_free(s);
  gsl_vector_free(coordinates);

  return status==0 ? true : false;
}

#else
bool optimize_layout(PlanarGraphView& G, vector<coord2d>& layout,
                     double, double, double, double){
  cerr << "Optimizing layouts is only available through GSL." << endl;
  return 0;
}
#endif


} // namespace layout2d
