#include "polyhedron.hh"

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

  work_queue.erase(v);

  coord3d c((points[0]+points[1]+points[2]+points[*v])/4.0);
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
    coord3d centre;		// Centre of visible faces
    for(triit t(output.begin());t!=output.end();t++){
      if(!Tri3D(points,*t).back_face(p)) { 
	visible.push_back(t);
	for(int i=0;i<3;i++) 
	  is_visible[dedge_t(t->u(i),t->u((i+1)%3))] = true; 
	centre += t->centroid(points);
      }
    }
    centre /= visible.size();

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

      //	Make sure new faces point outwards. 
      if(!Tri3D(points,t).back_face(centre)) t.flip();


      triit ti = output.insert(output.end(),t);
      //      for(int j=0;j<3;j++)
	//	edgetri[dedge_t(t[j],t[(j+1)%3])] = ti;
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

string Polyhedron::to_latex(bool show_dual, bool number_vertices, bool include_latex_header) const 
{
  ostringstream s;
  s.precision(2);
  s << fixed;
  if(include_latex_header)
    s << "\\documentclass{article}\n"
         "\\usepackage{fullpage,fourier,tikz}\n"
         "\\usetikzlibrary{calc,3d}"
         "\\begin{document}\n"
      "\\tikzstyle{vertex}=[circle, draw, inner sep="<<(number_vertices?"1pt":"0")<<", fill=blue!20, minimum width=3mm]\n"
      "\\tikzstyle{dualvertex}=[circle, draw, inner sep="<<(number_vertices?"1pt":"0")<<", fill=red!40, minimum width=2mm]\n"
      "\\tikzstyle{invisible}=[draw=none,inner sep=0,fill=none,minimum width=0pt]\n"
      "\\tikzstyle{edge}=[line width=1mm,brown]\n"
      "\\tikzstyle{dualedge}=[dotted,draw]\n"
      ;

  s << "\\begin{tikzpicture}\n";
  s << "\\foreach \\place/\\name/\\lbl in {";
  for(node_t u=0;u<N;u++){
    const coord3d& xs(points[u]);
    s << "{(" << xs[0] << "," << xs[1] << "," << xs[2] << ")/v" << u << "/$" << u << "$}" << (u+1<N? ", ":"}\n\t");
  }
  s << "\\node[vertex] (\\name) at \\place {"<<(number_vertices?"\\lbl":"")<<"};\n";
  s << "\\foreach \\u/\\v in {";
  for(set<edge_t>::const_iterator e(edge_set.begin()); e!=edge_set.end();){
    s << "{v"<<e->first<<"/v"<<e->second<<"}";
    if(++e != edge_set.end()) s << ", ";
  }
  s << "}\n\t\\draw[edge] (\\u) -- (\\v);\n";
#if 0
  vector<face_t> faces(compute_faces_flat(face_max));
  for(vector<face_t>::const_iterator f(faces.begin());f!=faces.end();f++){
    s << "\\fill[red!"<<50*(-points[(*f)[0]][0]+1)<<"]" ;
    for(size_t i=0;i<f->size();i++){
      coord3d xs(points[(*f)[i]]);
      s << "(" << xs[0] << "," << xs[1] << "," << xs[2] << ") -- " << (i+1<f->size()?"":"cycle;\n");
    }
  }
#endif


  if(show_dual){
    PlanarGraph dual(dual_graph(face_max));	// TODO: This breaks for everything else than fullerenes
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

Polyhedron::Polyhedron(const PlanarGraph& G, const vector<coord3d>& points_, const int face_max, const vector<face_t> faces_) : 
  PlanarGraph(G), face_max(face_max), points(points_), centre(centre3d(points)), faces(faces_)
{
  //cerr << "New polyhedron has " << N << " points. Largest face is "<<face_max<<"-gon.\n";
  if(faces.size() == 0){
    if(layout2d.size() != N){
      if(is_cubic()){
	layout2d = tutte_layout(-1,-1,-1,face_max);
	faces = compute_faces_flat(face_max,true);
	assert(outer_face.size() <= face_max);
      } else {
	layout2d = polar_angles(); 
	layout_is_spherical = true;
	faces = compute_faces_flat(face_max,true);
      }
    } else faces = compute_faces_flat(face_max,true);
  }
  
  //  assert(faces[0] == outer_face);
  //   }

  // cerr << "points = {"; for(int i=0;i<points.size();i++) cerr << points[i] << (i+1<points.size()? ", ":"};\n");
  // cerr << "faces  = {"; for(int i=0;i<faces.size();i++) cerr << faces[i] << (i+1<faces.size()? ", ":"};\n");
  // cerr << "G = " << static_cast<PlanarGraph>(*this) << endl;
  // cerr << "Layout has " << layout2d.size() << " points.\n";

  //  if(points.size() != N) 
  //    points = polar_mapping(spherical_projection());




  // cerr << "Found " << faces.size() << " faces.\n";
  // cerr << "Volume divergence: " << volume_divergence() << "\n";
  // cerr << "Volume tetra:      " << volume_tetra() << "\n";
  // cerr << "P = " << *this << endl;
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


string Polyhedron::to_povray(double w_cm, double h_cm, 
		   int line_colour, int vertex_colour, int face_colour,
		   double line_width, double vertex_diameter, double face_opacity) const 
{
  coord3d whd(width_height_depth());
  double xscale = w_cm/whd[0];

  ostringstream s;
  s << "#declare facecolour=color rgb <"<<((face_colour>>16)&0xff)/256.<<","<<((face_colour>>8)&0xff)/256.<<","<<(face_colour&0xff)/256.<<">;\n";
  s << "#declare faceopacity="<<face_opacity<<";\n";

  s << PlanarGraph(*this).to_povray(w_cm,h_cm,line_colour,vertex_colour,line_width,vertex_diameter);
  s << "#declare layout3D=array["<<N<<"][3]{"; for(int i=0;i<N;i++) s<<(points[i]*xscale)<<(i+1<N?",":"}\n\n"); 
  s << "#declare faces   =array["<<faces.size()<<"]["<<(face_max+1)<<"]{"; 
  for(int i=0;i<faces.size();i++) {
    const face_t& f(faces[i]);
    s << "{";
    for(int j=0;j<f.size();j++) s << f[j] << ",";
    for(int j=f.size();j<face_max;j++) s << "-1,";
    s << "-1}" << (i+1<faces.size()? ",":"}\n\n");
  }
  s << "#declare facelength=array["<<faces.size()<<"]{";for(int i=0;i<faces.size();i++) s<< faces[i].size() << (i+1<faces.size()?",":"}\n\n");

  vector<tri_t>   tris(centroid_triangulation(faces));
  vector<coord3d> centroid_points(points.begin(),points.end());
  for(int i=0;i<faces.size();i++) centroid_points.push_back(faces[i].centroid(points));

  s << "#declare tris = array["<<tris.size()<<"][3]{";
  for(int i=0;i<tris.size();i++) s <<  "{" << tris[i][0] << "," << tris[i][1] << "," << tris[i][2] << "}" << (i+1<tris.size()?",":"}\n\n");

  s << "#declare cpoints=array["<<centroid_points.size()<<"][3]{"; 
  for(int i=0;i<centroid_points.size();i++) s<<(centroid_points[i]*xscale)<<(i+1<centroid_points.size()?",\n":"}\n\n"); 

  s << "#include \"drawpolyhedron.pov\"\n\n";
  return s.str();
}

Polyhedron Polyhedron::dual(int Fmax, bool planar_layout) const 
{
  PlanarGraph d(dual_graph(Fmax,planar_layout));

  printf("|faces| = %ld\n",faces.size());
  vector<coord3d> coordinates(d.N);
  for(node_t u=0;u<d.N;u++){
    const face_t& f = faces[u];
    coord3d avg;
    for(int i=0;i<f.size();i++) avg += points[f[i]];
    coordinates[u] = avg/double(f.size());
  }
  
  return Polyhedron(d,points);
}

double Polyhedron::C20_points[20][3] = {{-1.376381920471174,0,0.2628655560595668},{1.376381920471174,0,-0.2628655560595668},{-0.4253254041760200,-1.309016994374947,0.2628655560595668},{-0.4253254041760200,1.309016994374947,0.2628655560595668},{1.113516364411607,-0.8090169943749474,0.2628655560595668},{1.113516364411607,0.8090169943749474,0.2628655560595668},{-0.2628655560595668,-0.8090169943749474,1.113516364411607},{-0.2628655560595668,0.8090169943749474,1.113516364411607},{-0.6881909602355868,-0.5000000000000000,-1.113516364411607},{-0.6881909602355868,0.5000000000000000,-1.113516364411607},{0.6881909602355868,-0.5000000000000000,1.113516364411607},{0.6881909602355868,0.5000000000000000,1.113516364411607},{0.8506508083520399,0,-1.113516364411607},{-1.113516364411607,-0.8090169943749474,-0.2628655560595668},{-1.113516364411607,0.8090169943749474,-0.2628655560595668},{-0.8506508083520399,0,1.113516364411607},{0.2628655560595668,-0.8090169943749474,-1.113516364411607},{0.2628655560595668,0.8090169943749474,-1.113516364411607},{0.4253254041760200,-1.309016994374947,-0.2628655560595668},{0.4253254041760200,1.309016994374947,-0.2628655560595668}};
