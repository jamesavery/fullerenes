#include <iomanip>
#include <limits>

#include "fullerenes/triangulation.hh"
#include "fullerenes/polyhedron.hh"


double Polyhedron::diameter() const {
  double dmax = -INFINITY;
  for(int i=0;i<N;i++)
    for(int j=i+1;j<N;j++){
      double d = (points[i]-points[j]).norm();
      if(d > dmax) dmax = d;
    }
  return dmax;
};

double Polyhedron::surface_area() const {
  double A = 0;

  // vector<tri_t> tris(triangulation(faces));
  
  // for(size_t i=0;i<tris.size();i++){
  //   const tri_t& tri(tris[i]);
  //   Tri3D T(points[tri[0]],points[tri[1]],points[tri[2]]);
  //   A += T.area();
  // }
  // return A;

  vector<tri_t> tris(centroid_triangulation(faces));
  vector<coord3d> centroid_points(points.begin(),points.end());
  for(int i=0;i<faces.size();i++) centroid_points.push_back(faces[i].centroid(points));
  
  for(size_t i=0;i<tris.size();i++){
    const tri_t& tri(tris[i]);
    Tri3D T(centroid_points[tri[0]],centroid_points[tri[1]],centroid_points[tri[2]]);
    A += T.area();
  }
  return A;
}

double Polyhedron::volume_tetra() const {
  vector<tri_t>   tris(centroid_triangulation(faces));
  vector<coord3d> centroid_points(points.begin(),points.end());
  for(int i=0;i<faces.size();i++) centroid_points.push_back(faces[i].centroid(points));

  double V = 0,Vm=0,Vp=0;

  // Now generate tetrahedra and either add or subtract volume according to which direction the face is pointing
  coord3d zero(0,0,0);
  for(size_t i=0;i<tris.size();i++){
    const face_t& t(tris[i]);
    Tri3D T(centroid_points[t[0]],centroid_points[t[1]],centroid_points[t[2]]);
    double dV = Tetra3D(T.a,T.b,T.c,zero).volume();
    V += (T.back_face(zero)? 1 : -1)*dV;
    if(T.back_face(zero)) Vp += dV;
    else Vm += dV;
    //    if(!T.back_face(zero))
      //      cerr << "Tri " << t << " / " << T << " is a front face.\n";
  }
  fprintf(stderr,"V = %f - %f = %f\n",Vp,Vm,V);
  return fabs(V);
}

double Polyhedron::volume_divergence() const {
  vector<tri_t>   tris(centroid_triangulation(faces));
  vector<coord3d> centroid_points(points.begin(),points.end());
  for(int i=0;i<faces.size();i++) centroid_points.push_back(faces[i].centroid(points));

  //  cerr << "points = {"; for(int i=0;i<centroid_points.size();i++) cerr << centroid_points[i] << (i+1<centroid_points.size()? ", ":"};\n");

  double V = 0;

  // Now generate tetrahedra and either add or subtract volume according to which direction the face is pointing
  for(size_t i=0;i<tris.size();i++){
    const face_t& t(tris[i]);
    Tri3D T(centroid_points[t[0]],centroid_points[t[1]],centroid_points[t[2]]);

    V += ((T.a).dot(T.n))*T.area()/T.n.norm();
  }
  return fabs(V/3.0);
}

Polyhedron Polyhedron::incremental_convex_hull() const {
  list<tri_t> output;
  typedef list<tri_t>::iterator triit;
  list<node_t> work_queue;
  srandom(42); // Seed random numbers with constant for reproducible behaviour

  // 1. Create initial tetrahedron.
  // 1.1 Find 4 non-coplanar points
  Tri3D T(points,tri_t(0,1,2));
  
  for(node_t u=3;u<N;u++) work_queue.push_front(u);

  // Get the point furthest from the (0,1,2)-plane
  double distmax = 0;
  list<node_t>::iterator v(work_queue.begin());
  for(list<node_t>::iterator u(work_queue.begin());u!=work_queue.end();u++){
    double dist = T.distance(points[*u]);
    if(dist > distmax){
      distmax = dist;
      v = u;
    }
  }

  // 1.2 Add faces to output
  // cerr << "// 1. Create initial tetrahedron: [0, 1, 2, " << *v << "], volume = "<<Tetra3D(points[0],points[1],points[2],points[*v]).volume() << ". \n";
  output.push_back(tri_t(0,1,2));
  output.push_back(tri_t(0,1,*v));
  output.push_back(tri_t(0,2,*v));
  output.push_back(tri_t(1,2,*v));

  coord3d c((points[0]+points[1]+points[2]+points[*v])/4.0);
  work_queue.erase(v);

  for(triit t(output.begin());t!=output.end();t++){
    // Make sure all faces point away from the centroid. 
    if(!Tri3D(points,*t).back_face(c)) t->flip(); 
  }
    


  // 2. For each remaining vertex u
  // cerr << "// 2. For each remaining vertex u\n";
  for(list<node_t>::const_iterator u(work_queue.begin());u!=work_queue.end();u++){
    long r = random();
    const coord3d perturbation(r&0xff,(r>>8)&0xff,(r>>16)&0xff);

    // Perturb p randomly
    coord3d p(points[*u]);
    p *= (coord3d(1,1,1)+perturbation*1e-13);

    // 2.1 Find all faces visible from p ( (f.centroid() - p).dot(f.n) > 0 ) 
    list<triit> visible;
    map<dedge_t,bool> is_visible;
    coord3d centre; // Centre of visible faces
    for(triit t(output.begin());t!=output.end();t++){
      if(!Tri3D(points,*t).back_face(p)) { 
        visible.push_back(t);
        for(int i=0;i<3;i++) 
          is_visible[dedge_t(t->u(i),t->u((i+1)%3))] = true; 
        centre += t->centroid(points);
      }
    }
    if(visible.size() != 0) centre /= visible.size();

    // 2.2 Build set of horizon edges: each edge e in visible faces that has f_a visible, f_b invisible
    list<edge_t> horizon;
    for(list<triit>::const_iterator tvi(visible.begin()); tvi!=visible.end(); tvi++){
      const tri_t& tv(**tvi);

      for(int j=0;j<3;j++){
        const dedge_t e(tv[j],tv[(j+1)%3]);

        if( (is_visible[e] && !is_visible[dedge_t(e.second,e.first)]) || (!is_visible[e] && is_visible[dedge_t(e.second,e.first)]) )
          horizon.push_back(edge_t(e));
      }
      // 2.3 Delete visible faces from output set. 
      output.erase(*tvi);
    }

    // 2.4 For each e in horizon, add tri_t(u,e[0],e[1]) to output set. 
    for(list<edge_t>::const_iterator e(horizon.begin()); e!=horizon.end(); e++){
      tri_t t(*u,e->first,e->second);

      //        Make sure new faces point outwards. 
      if(!Tri3D(points,t).back_face(centre)) t.flip();


      triit ti = output.insert(output.end(),t);
      //      for(int j=0;j<3;j++)
        //        edgetri[dedge_t(t[j],t[(j+1)%3])] = ti;
    }
    if(output.size() > N*N*10){
      fprintf(stderr,"Something went horribly wrong in computation of convex hull:\n");
      fprintf(stderr,"Data sizes: output(%ld), visible(%ld), is_visible(%ld), horizon(%ld), horizon-visible: %ld\n",
              output.size(),visible.size(),is_visible.size(),horizon.size(),horizon.size()-visible.size());
    }
  }
    
  // 3. Finally, construct the graph and the output polyhedron object
  set<node_t> used_nodes;
  for(triit t(output.begin()); t!=output.end(); t++)
    for(int i=0;i<3;i++)
      used_nodes.insert(t->u(i));
  map<node_t,node_t> nodemap;
  vector<coord3d> remaining_points(used_nodes.size());
  node_t i=0;
  for(set<node_t>::const_iterator u(used_nodes.begin()); u!=used_nodes.end(); u++,i++){
    nodemap[*u] = i;
    remaining_points[i] = points[*u];
  }

  set<edge_t> edges;
  for(triit t(output.begin()); t!=output.end(); t++)
    for(int i=0;i<3;i++)
      edges.insert(edge_t(nodemap[t->u(i)],nodemap[t->u((i+1)%3)]));
    
  PlanarGraph g(edges);
  cout << "Polyhedron is "<< (g.N != N?"not ":"") << "equal to convex hull.\n"; 
  vector<face_t> faces;
  for(list<tri_t>::const_iterator o(output.begin());o!=output.end();o++){
    const tri_t& f(*o);
    faces.push_back(tri_t(nodemap[f[0]],nodemap[f[1]],nodemap[f[2]]));
  }
  g.outer_face = faces[0];
  return Polyhedron(g,remaining_points,3,faces);
}

struct sort_ccw_coord3d {
  const vector<coord3d> &points;
  coord3d X, Y, n;

  // Points are only neighbour displacements from origin-node
  sort_ccw_coord3d(const vector<coord3d>& points)
    : points(points) {

    // TODO: More numerically robust method    
    coord3d xc(0,0,0);
    for(const coord3d &p: points) xc += p;
    xc /= points.size();

    n = /*x0*/-xc; n /= n.norm();	//

    coord3d x1 = points[0] /* - x0*/;
    X = x1 - n*x1.dot(n); X /= X.norm();
    Y = n.cross(x1);
  }

  bool operator()(const node_t& s, const node_t& t) const {
    coord3d xs(points[s]/*-x0*/), xt(points[t]/*-x0*/);
    coord2d Xs(X.dot(xs), Y.dot(xs)), Xt(X.dot(xt), Y.dot(xt));
    
    double angs = atan2(Xs.first,Xs.second), angt = atan2(Xt.first,Xt.second);
    return angs <= angt;
  } 
};



void Polyhedron::orient_neighbours()
{
  
  if(layout2d.size() != N)
    layout2d = tutte_layout();
  PlanarGraph::orient_neighbours();

  // TODO: Sort this out
  //  else if(points.size() == N){
  // Orient neighbours locally (CW or CCW depending on luck)
  // for(node_t u=0;u<N;u++){
  //   vector<node_t> &ns(neighbours[u]);
  //   size_t degree = ns.size();
      
  //   int ns_index[degree];
  //   for(int i=0;i<degree;i++) ns_index[i] = i;
				     
  //   vector<coord3d> neighbour_points(degree);
  //   const coord3d &x0 = points[u];

  //   for(int i=0;i<degree;i++) neighbour_points[i] = points[ns[i]] - x0;
  //   sort_ccw_coord3d CCW(neighbour_points);

  //   sort(ns_index,ns_index+degree,CCW);
  //   vector<node_t> ns_sorted(degree);
  //   for(int i=0;i<degree;i++) ns_sorted[i] = ns[ns_index[i]];
  //   ns = ns_sorted;
  // }


  // // Choose that first node is correct, then flip to consistency
  // // TODO: How? For now, cheat slowly.
  // if(is_consistently_oriented())
  //   is_oriented = true;
  // else {
  //   layout2d = tutte_layout();
  //   PlanarGraph::orient_neighbours();
  // }
  //}
  
  // Calculate volume
  double V=0;
  for(node_t u=0;u<N;u++){
    const face_t nu(neighbours[u]);
    const coord3d ux(points[u]);

    for(int i=0;i<nu.size();i++){
      Tri3D T(ux, points[nu[i]],points[nu[(i+1)%nu.size()]]);
      V += ((T.a).dot(T.n))*T.area()/T.n.norm();
    }
  }

  if(V<0){ // Calculated normals are pointing inwards - reverse order.
    //    printf("Inverted normals - reversing neighbours lists.\n");
    for(node_t u=0;u<N;u++) reverse(neighbours[u].begin(), neighbours[u].end());
  }
}

Polyhedron::Polyhedron(const PlanarGraph& G, const vector<coord3d>& points_, const int face_max_, const vector<face_t> faces_) : 
  PlanarGraph(G), face_max(face_max_), points(points_), faces(faces_)
{
  if(!is_oriented) orient_neighbours();

  //  for(node_t u=0;u<N;u++) points[u] = points_[u];
  
  if(faces.size() == 0){
    faces = compute_faces(face_max);
    assert(outer_face.size() <= face_max);
    face_max = 0;
    for(int i=0;i<faces.size();i++) if(faces[i].size() > face_max) face_max = faces[i].size();
  } 

  if(!is_oriented) orient_neighbours();
}

Polyhedron::Polyhedron(const vector<coord3d>& xs, double tolerance) 
{
  double bondlength = INFINITY;

  for(int i=0;i<xs.size();i++){
    for(int j=i+1;j<xs.size();j++){
      double d = (xs[i]-xs[j]).norm();
      if(d < bondlength) bondlength = d;
    }
  }
     
  set<edge_t> edges;
  for(int i=0;i<xs.size();i++){
    for(int j=i+1;j<xs.size();j++){
      double d = (xs[i]-xs[j]).norm();
      if(d <= bondlength*tolerance) {
        edges.insert(edge_t(i,j));
      }
    }
  }
  
  (*this) = Polyhedron(PlanarGraph(edges), xs);
}


matrix3d Polyhedron::inertia_matrix() const
{
  matrix3d I;

  for(int k=0;k<points.size();k++){
    const coord3d& x(points[k]);
    const long double xx(x.dot(x));
    for(int i=0;i<3;i++){
      I(i,i) += xx;

      for(int j=0;j<3;j++)
        I(i,j) -= x[i]*x[j];
    }
  }
  return I;
}

matrix3d Polyhedron::principal_axes() const
{
  const matrix3d I(inertia_matrix());
  pair<coord3d,matrix3d> ES(I.eigensystem());

  matrix3d Id;
  Id(0,0) = 1; 
  Id(1,1) = 1; 
  Id(2,2) = 1; 
/*
  cerr << "Inertial frame:\n " 
       << " inertia_matrix = " << I << ";\n"
       << " lambda  = " << ES.first << ";\n"
       << " vectors = " << ES.second << ";\n";
*/

  for(int i=0;i<3;i++) 
    if(std::isnan(ES.first[i])){
      cerr << "Warning: Inertial frame returned NaN. Setting inertial frame transformation to identity.\n";
      return Id;
    }
  
  if((ES.second*ES.second.transpose() - Id).norm() > 1e-2){
    cerr << "Warning: Inertial frame transform is not unitary. Setting inertial frame transformation to identity.\n";
    return Id;
  }

  return ES.second;
}

coord3d Polyhedron::width_height_depth() const {
  double xmin=INFINITY,xmax=-INFINITY,ymin=INFINITY,ymax=-INFINITY,zmin=INFINITY,zmax=-INFINITY;
  for(node_t u=0;u<N;u++){
    const coord3d& x(points[u]);
    if(x[0]<xmin) xmin = x[0];
    if(x[0]>xmax) xmax = x[0];
    if(x[1]<ymin) ymin = x[1];
    if(x[1]>ymax) ymax = x[1];
    if(x[2]<zmin) zmin = x[2];
    if(x[2]>zmax) zmax = x[2];
  }
  return coord3d(xmax-xmin,ymax-ymin,zmax-zmin);
}



Polyhedron Polyhedron::dual() const 
{
  PlanarGraph d(dual_graph());

  vector<coord3d> coordinates(d.N);
  for(node_t u=0;u<d.N;u++){
    const face_t& f = faces[u];
    coord3d avg;
    for(auto v: f) avg += points[v];
    coordinates[u] = avg/double(f.size());
  }
 
  return Polyhedron(d,coordinates);
}


Polyhedron Polyhedron::leapfrog_dual() const 
{
  assert(is_oriented);
  size_t Nf = faces.size();
  
  Polyhedron Plf(Graph(N+Nf,true));
  Plf.points.reserve(N+Nf);
   
  // Start with all the existing nodes
  for(node_t u=0;u<N;u++){
    Plf.neighbours[u] = neighbours[u];
    Plf.points[u]     = points[u];
  }

  // The result is a deltahedron: a polyhedron consisting of
  // only triangles. The first Nv points are the original vertices,
  // the last Nf are the midpoints of the original faces.
  int n_tris = 0;
  for(auto f: faces) n_tris += f.size();
  Plf.faces = vector<face_t>(n_tris);

  // Now connect new face-center nodes in oriented order
  for(int i=0;i<faces.size();i++){
    const face_t &f  = faces[i];
    node_t c = N+i;                // Face-center node
    
    // cerr << "new node " << c << " at face " << f << "\n";
    coord3d xc = {0,0,0};
    size_t   d = f.size();
    for(int j=0;j<d;j++){
      node_t u = f[j], v = f[(j+1)%f.size()];

      // Center node position is middle of face
      xc += points[u]/d;

      // Add edge mumble mumble
      Plf.insert_edge(dedge_t{v,c},u,-1);

      // Add triangle
      Plf.faces[c] = tri_t{u,v,c};
    }
    Plf.points[c] = xc;
  }

  return Plf;
}

Polyhedron Polyhedron::fullerene_polyhedron(FullereneGraph G)
{
  if(G.layout2d.empty())
    G.layout2d       = G.tutte_layout();
  
  Polyhedron P(G,G.zero_order_geometry(),6);
  P.points = G.optimized_geometry(P.points);

  P.move_to_origin();		// Center of mass at (0,0,0)
  P.align_with_axes();		// Align with principal axes
  
  return P;
}

bool Polyhedron::optimize(int opt_method, double ftol)
{
  if(is_a_fullerene()){
    FullereneGraph g(*this);
    points = g.optimized_geometry(points,opt_method,ftol);
    return true;
  } if(is_cubic()) {
    bool optimize_angles = true;
    return optimize_other(optimize_angles);
  } else if(is_triangulation()) {
    bool optimize_angles = false;
    return optimize_other(optimize_angles);
  }else{
     Triangulation LFD = leapfrog_dual();

    // inverse_tranges for faster lookup:
    // indices in 'triangles' at which there are triangles containing this vertex
    vector<vector<int>> inverse_triangle_list(LFD.N);
    for (int i=0; i<LFD.triangles.size(); i++){
      const tri_t& tri = LFD.triangles[i];
      inverse_triangle_list[tri[0]].push_back(i);
      inverse_triangle_list[tri[1]].push_back(i);
      inverse_triangle_list[tri[2]].push_back(i);
    }

    // generate and optimise LF of initial polyhedron
    PlanarGraph LF = LFD.dual_graph();
    LF.layout2d = LF.tutte_layout();
    Polyhedron P(LF,LF.zero_order_geometry());
    bool optimize_angles = true;
    bool opt_success = P.optimize_other(optimize_angles);

    // for each face in LF which corresponds to a vertex in the initial graph,
    // find the average coordinates of all vertices (ie the face centre)
    for (int i=0; i<N; i++){
      const int face_size(inverse_triangle_list[i].size());
      coord3d face_centre;
      for (int j=0; j<face_size; j++){
        face_centre += P.points[inverse_triangle_list[i][j]];
      }
      face_centre /= face_size;
      points[i] = face_centre/sqrt(3);
    }
    
    return opt_success;
  }
}

bool Polyhedron::is_triangulation() const {
  for(int i=0;i<faces.size();i++) if(faces[i].size()!=3) return false;
  return true;
}

// TODO: Add function for checking if forcefield convergence is achieved

bool Polyhedron::is_invalid() const {
  bool has_nans = false;
  for(auto p: points){
    if(std::isnan(p[0])||std::isnan(p[1])||std::isnan(p[1])) has_nans = true;
  }
  return has_nans;
}

pair<coord3d,coord3d> Polyhedron::bounding_box() const {
  coord3d xmin{INFINITY,INFINITY,INFINITY},xmax{-INFINITY,-INFINITY,-INFINITY};
  for(const auto& p: points)
    for(int i=0;i<3;i++){
      if(p[i] < xmin[i]) xmin[i] = p[i];
      if(p[i] > xmax[i]) xmax[i] = p[i];
    }
  return make_pair(xmin,xmax);
}

double Polyhedron::C20_points[20][3] = {{-1.376381920471174,0,0.2628655560595668},{1.376381920471174,0,-0.2628655560595668},{-0.4253254041760200,-1.309016994374947,0.2628655560595668},{-0.4253254041760200,1.309016994374947,0.2628655560595668},{1.113516364411607,-0.8090169943749474,0.2628655560595668},{1.113516364411607,0.8090169943749474,0.2628655560595668},{-0.2628655560595668,-0.8090169943749474,1.113516364411607},{-0.2628655560595668,0.8090169943749474,1.113516364411607},{-0.6881909602355868,-0.5000000000000000,-1.113516364411607},{-0.6881909602355868,0.5000000000000000,-1.113516364411607},{0.6881909602355868,-0.5000000000000000,1.113516364411607},{0.6881909602355868,0.5000000000000000,1.113516364411607},{0.8506508083520399,0,-1.113516364411607},{-1.113516364411607,-0.8090169943749474,-0.2628655560595668},{-1.113516364411607,0.8090169943749474,-0.2628655560595668},{-0.8506508083520399,0,1.113516364411607},{0.2628655560595668,-0.8090169943749474,-1.113516364411607},{0.2628655560595668,0.8090169943749474,-1.113516364411607},{0.4253254041760200,-1.309016994374947,-0.2628655560595668},{0.4253254041760200,1.309016994374947,-0.2628655560595668}};
