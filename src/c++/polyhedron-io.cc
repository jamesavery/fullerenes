#include "fullerenes/polyhedron.hh"

//////////////////////////// FORMAT MULTIPLEXING ////////////////////////////
vector<string> Polyhedron::formats{{"ascii","planarcode","xyz","mol2","mathematica","latex","cc1","turbomole","gaussian","wavefront","spiral"}};
vector<string> Polyhedron::format_alias{{"txt","ply","xyz","mol2","m","tex","cc1","turbomole","com","obj","rspi"}};
vector<string> Polyhedron::input_formats{{"xyz","mol2"}}; // TODO: "ascii","planarcode","obj"
vector<string> Polyhedron::output_formats{{"ascii","xyz","mol2","cc1","turbomole","gaussian","spiral"}};

int Polyhedron::format_id(string name)
{
  for(int i=0;i<formats.size();i++) if(name == formats[i]) return i;
  for(int i=0;i<formats.size();i++) if(name == format_alias[i]) return i;  
  return -1;
}

Polyhedron Polyhedron::from_file(FILE *file, string format)
{
  switch(format_id(format)){
    // case ASCII:
    //   return from_ascii(file);
    //  case PLANARCODE:
    // PlanarGraph G = PlanarGraph::from_planarcode(file);
    // COMPUTE EMBEDDING
    // return P(G,points);
  case XYZ:
    return from_xyz(file);
  case MOL2:
    return from_mol2(file);
  default: {
    PlanarGraph G = PlanarGraph::from_file(file,format);

    if(!G.N){    
      cerr << "Input format is '"<<format<<"'; must be one of: " << input_formats << " or " << PlanarGraph::input_formats << "\n";
      abort();
    } else {
      G.layout2d = G.tutte_layout();
      Polyhedron P(G,G.zero_order_geometry(),6);
      P.optimize();

      return P;
    }
  }
  }
}


bool Polyhedron::to_file(const Polyhedron &G, FILE *file, string format)
{
  switch(format_id(format)){
  case ASCII:
    return Polyhedron::to_ascii(G,file);
  case XYZ:
    return Polyhedron::to_xyz(G,file);
  case MOL2:
    return Polyhedron::to_mol2(G,file);
  case CC1:
    return Polyhedron::to_cc1(G,file);
  case TURBOMOLE:
    return Polyhedron::to_turbomole(G,file);
  case GAUSSIAN:
    return Polyhedron::to_gaussian(G,file);
  case WAVEFRONT_OBJ:
    return Polyhedron::to_wavefront_obj(G,file);
  case SPIRAL:
    return PlanarGraph::to_spiral(G,file);    
  default:
    cerr << "Output format is '"<<format<<"'  must be one of: " << output_formats << "\n";
    return false;
  }
}


Polyhedron Polyhedron::from_file(string filename)
{
  FILE *file = fopen(filename.c_str(),"rb");
  string extension = filename_extension(filename);
  Polyhedron G = from_file(file,extension);
  fclose(file);
  return G;
}

bool Polyhedron::to_file(const Polyhedron &G, string filename)
{
  FILE *file = fopen(filename.c_str(),"wb");
  string extension = filename_extension(filename);  
  to_file(G,file,extension);
  fclose(file);
  return true;			// TODO: Check success
}


////////////////////////////// OUTPUT ROUTINES //////////////////////////////
bool Polyhedron::to_ascii(const Polyhedron &P, FILE *file)  {
  string s = LIST_OPEN + to_string(P.neighbours) + "," + to_string(P.points) + "," + to_string(P.faces) + LIST_CLOSE;
  fputs(s.c_str(),file);
  return ferror(file) == 0;
}
  

bool Polyhedron::to_turbomole(const Polyhedron &P, FILE *file)  {
  const double aa2bohr = 1.889716164632;
  fprintf(file,"$coord\r\n");
  for(int i=0; i<P.N; ++i){
    const coord3d p = P.points[i] * aa2bohr;
    fprintf(file,"%f %f %f  c\r\n", p[0],p[1],p[2]);
  }
  fprintf(file,"$end\r\n");

  return true;			// TODO: Check file status
}

bool Polyhedron::to_gaussian(const Polyhedron &P, FILE *file, string header)  {
  
  
  if(header.empty()){
    // TODO: Make general! Autodetect fullerene, triangulation, cubic, etc.
    // TODO: Move to to_file, add as parameter.
    auto naming_scheme = spiral_nomenclature::FULLERENE;
    auto construction_scheme = spiral_nomenclature::CUBIC;
    spiral_nomenclature name(P,naming_scheme,construction_scheme);
    
    header ="\nGeometry for C" + to_string(P.N) + " fullerene " + name.to_string() + "\n\n0 1\n";
  }

  fprintf(file,"%s",header.c_str());

  // Atom section
  for(node_t u=0; u<P.N; u++){
    const coord3d p = P.points[u];
    fprintf(file," C %f %f %f\n", p[0],p[1],p[2]);
  }
  fprintf(file,"\n");

  // Connectivity section
  for(node_t u=0; u<P.N;u++){
    auto nu = P.neighbours[u];
    for(auto v: nu) fprintf(file,"%d %d B\n",u+1, v+1);
  }

  return true;			// TODO: Check file status
}

bool Polyhedron::to_xyz(const Polyhedron &P, FILE *file) {
  fprintf(file,"%d\r\n",P.N);
  fprintf(file,"# Created by libgraph from Fullerene (http://tinyurl.com/fullerenes)\r\n");
  for(node_t u=0; u<P.N; ++u){
    const coord3d p = P.points[u];
    fprintf(file,"C  %f  %f  %f\r\n",p[0],p[1],p[2]);
  }
  return true;
}

bool Polyhedron::to_wavefront_obj(const Polyhedron &P, FILE *file)
{
  fprintf(file,"# Vertices:\n");    
  for(auto p: P.points)
    fprintf(file,"v %f %f %f\n",p[0],p[1],p[2]);

  for(auto f: P.faces){
    fprintf(file,"f ");
    for(auto v: f) fprintf(file,"%d ",v);
    fprintf(file,"\n");
  }
  // fprintf(file,"# Pentagons:\n"
  // 	       "g pentagons\n");
  // for(auto f: P.faces)
  //   if(f.size()==5)
  //     fprintf(file,"f %d %d %d %d %d\n",f[0]+1,f[1]+1,f[2]+1,f[3]+1,f[4]+1);

  // fprintf(file,"# Hexagons:\n"
  // 	       "g hexagons\n");
  // for(auto f: P.faces)
  //   if(f.size()==6)
  //     fprintf(file,"f %d %d %d %d %d %d\n",f[0]+1,f[1]+1,f[2]+1,f[3]+1,f[4]+1,f[5]+1);  
}

bool Polyhedron::to_mol2(const Polyhedron &P, FILE *file)
{
  size_t Nedges = P.count_edges();
  fprintf(file,
	  "# Created by libgraph from Fullerene (http://tinyurl.com/fullerenes)\r\n"
	  "@<TRIPOS>MOLECULE\r\n"
	  "Fullerene\r\n"
	  "\t%d\t%ld\t0\t0\t0\r\n"
	  "SMALL\r\n"
	  "NO_CHARGES\r\n\r\n",P.N,Nedges);

  fprintf(file,"@<TRIPOS>ATOM\r\n");
  
  for(node_t u=0; u < P.N; u++){
    const coord3d p = P.points[u];
    fprintf(file,"%d\t C%d\t %f\t %f\t %f\t C\t 1\t Unk\t 0\r\n",u+1,u,p[0],p[1],p[2]);
  }

  fprintf(file,"@<TRIPOS>BOND\r\n");
  int i = 1;
  for(node_t u=0;u<P.N;u++){
    for(node_t v: P.neighbours[u])
      if(v>=u)
	fprintf(file,"%d\t %d\t %d\t un\r\n",i++,u+1,v+1);
  }

  return true;
}

bool Polyhedron::to_cc1(const Polyhedron &P, FILE *file) 
{
  const int weird_constant = 2;

  fprintf(file,"%d\r\n",P.N);

  for(node_t u=0; u < P.N; u++){
    const coord3d p         = P.points[u];
    
    fprintf(file,"C\t %d\t %f\t %f\t %f\t %d", u+1,p[0],p[1],p[2],weird_constant);
    for(node_t v: P.neighbours[u])
      fprintf(file,"\t%d",v);
    fprintf(file,"\r\n");
  }

  return true;
}

// TODO: Decide on consistent I/O interface
string Polyhedron::to_latex(bool show_dual, bool number_vertices, bool include_latex_header) const 
{
  ostringstream s;
  s.precision(2);
  s << fixed;

  vector<edge_t> edge_set = undirected_edges();

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
  for(int i=0;i<edge_set.size();i++){
    edge_t e = edge_set[i];
    s << "{v"<<e.first<<"/v"<<e.second<<"}";
    if(i+1<edge_set.size()) s << ", ";
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
    PlanarGraph dual(dual_graph(face_max));        // TODO: This breaks for everything else than fullerenes
    s << "\\foreach \\place/\\name/\\lbl in {";
    for(node_t u=0;u<dual.N;u++){
      const coord2d& xs(dual.layout2d[u]);
      s << "{(" << xs.first << "," << xs.second << ")/v" << u << "/$" << u << "$}" << (u+1<dual.N? ", ":"}\n\t");
    }    
    s << "\\node[dualvertex] (\\name) at \\place {"<<(number_vertices?"\\lbl":"")<<"};\n";
    s << "\\foreach \\u/\\v in {";
    vector<edge_t> dual_edges = dual.undirected_edges();
    for(int i=0;i<dual_edges.size();i++){
      edge_t e = dual_edges[i];
      s << "{v"<<e.first<<"/v"<<e.second<<"}";
      if(i+1<dual_edges.size()) s << ", ";
    }
    s << "}\n\t\\draw[dualedge] (\\u) -- (\\v);\n";
  }

  s<<"\\end{tikzpicture}\n";
  if(include_latex_header)
    s << "\\end{document}\n";

  return s.str();
}

// TODO: Decide on consistent I/O interface
string Polyhedron::to_povray(double w_cm, double h_cm, 
                   int line_colour, int vertex_colour, int face_colour,
                   double line_width, double vertex_diameter, double face_opacity) const 
{
  //  coord3d whd(width_height_depth()); // TODO: Removed width/height -- much better to use real coordinates and handle layout in host pov file.

  ostringstream s;
  s << "#declare facecolour=color rgb <"<<((face_colour>>16)&0xff)/256.<<","<<((face_colour>>8)&0xff)/256.<<","<<(face_colour&0xff)/256.<<">;\n";
  s << "#declare faceopacity="<<face_opacity<<";\n\n";

  s << PlanarGraph(*this).to_povray(w_cm,h_cm,line_colour,vertex_colour,line_width,vertex_diameter);
  s << "#declare layout3D=array["<<N<<"][3]" << points <<";\n\n";

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
  vector<int>     triface;
  vector<coord3d> centroid_points(points.begin(),points.end());
  vector<coord3d> trinormals(tris.size()), facenormals(faces.size()), vertexnormals(points.size()+faces.size());

  for(int i=0;i<faces.size();i++)
    centroid_points.push_back(faces[i].centroid(points));

  for(int i=0;i<tris.size();i++){
    coord3d n(Tri3D(centroid_points[tris[i][0]],centroid_points[tris[i][1]],centroid_points[tris[i][2]]).n);
    trinormals[i] = n/n.norm();
    for(int j=0;j<3;j++) vertexnormals[tris[i][j]] += trinormals[i];
  }

  for(int i=0;i<N;i++)
    vertexnormals[i] /= vertexnormals[i].norm();

  // Calculate volume
  double V=0;
  for(size_t i=0;i<tris.size();i++){
    const face_t& t(tris[i]);
    Tri3D T(centroid_points[t[0]],centroid_points[t[1]],centroid_points[t[2]]);

    V += ((T.a).dot(T.n))*T.area()/T.n.norm();
  }
  s << "#declare volume="<<fabs(V)<<";\n";

  if(V<0)                        // Calculated normals are pointing inwards!
    for(int i=0;i<tris.size();i++) trinormals[i] *= -1;

  for(int i=0;i<faces.size();i++) {
    coord3d normal;
    if(faces[i].size()>3){
      for(int j=0;j<faces[i].size();j++){
        triface.push_back(i);
        normal += trinormals[triface.size()-1];
      } 
      facenormals[i] = normal/normal.norm();
    } else {
      triface.push_back(i);
      facenormals[i] = trinormals[i];
    }
  }



  s << "#declare Ntris = "<<tris.size()<<";\n";
  s << "#declare tris = array["<<tris.size()<<"][3]" << tris << ";\n\n";
  s << "#declare triface = array["<<triface.size()<<"]" << triface << ";\n\n";
    
  s << "#declare cpoints=array["<<centroid_points.size()<<"][3]" << centroid_points << ";\n\n"; 
  s << "#declare vertexnormals =array["<<vertexnormals.size()<<"][3]" << vertexnormals << ";\n\n";
  s << "#declare trinormals =array["<<tris.size()<<"][3]" << trinormals << ";\n\n";
  s << "#declare facenormals=array["<<faces.size()<<"][3]" << facenormals << ";\n\n";

  //  s << "#include \"drawpolyhedron.pov\"\n\n";
  return s.str();
}

////////////////////////////// INPUT ROUTINES //////////////////////////////
Polyhedron Polyhedron::from_xyz(FILE *file)
{
  int N;
  string Nstring, comment, element,line;
  vector<coord3d> points;

  getline(file,Nstring);
  getline(file,comment);

  N = strtol(Nstring.c_str(),0,0);

  for(int i=0; i < N && getline(file,line); i++){
    stringstream l(line);
    coord3d x;

    l >> element;
    for(int j=0;j<3 && l.good(); j++)
      l >> x[j];

    points.push_back(x);
  }

  assert(points.size() == N);
  return Polyhedron(points);
}

// Read in .mol2 files. 
// NB: Doesn't support full format. Can only read .mol2 files that we've written ourselves!
Polyhedron Polyhedron::from_mol2(FILE *file)
{
  string 
    header_marker = "@<TRIPOS>MOLECULE",
    point_marker  = "@<TRIPOS>ATOM",
    edge_marker   = "@<TRIPOS>BOND";

  int N, Nedges;
  vector<coord3d> points;
  string line;

  // Fast forward to metadata section
  while(getline(file,line) && line.compare(0,header_marker.size(),header_marker)) ;
  getline(file,line);  
  //  assert(!line.compare(0,9,"Fullerene")); // TODO: Fail gracefully if we didn't create the file.

  getline(file,line);
  stringstream l(line);
  l >> N;
  l >> Nedges;

  Graph G(N);
  
  //  cerr << "line="<<line<<"; N="<<N<<", Nedges="<<Nedges<<endl;
  
  // Fast forward to coordinate section
  while(getline(file,line) && line.compare(0,point_marker.size(),point_marker)) ;

  bool file_ok = true;
  for(int i=0;i<N && file_ok;i++){
    getline(file,line);
    file_ok = (ferror(file) == 0);
    stringstream l(line);
    string vid,element;
    coord3d x;

    if(file_ok) l >> vid;
    if(file_ok) l >> element;
    for(int j=0;j<3 && l.good(); j++) l >> x[j];
    points.push_back(x);
    //    cerr << i << " of " << N << ": Read line "<< line;
    //    cerr << "Point " << x << endl;
  }
  assert(points.size() == N);         // TODO: Fail gracefully if file format error.


  // Fast forward to edge section
  while(getline(file,line) && line.compare(0,edge_marker.size(),edge_marker)) ;  

  int i=0;
  for(;i<Nedges && file_ok;i++){
    getline(file,line);
    file_ok = (ferror(file) == 0);
    stringstream l(line);
    int eid, u[2];

    l >> eid;
    for(int j=0;j<2 && l.good(); j++) l >> u[j];
    G.insert_edge(edge_t(u[0]-1,u[1]-1));

    //    cerr << "Edge " << i << " of " << Nedges << ": Read line "<< line <<endl;
  }

  Polyhedron P;  
  P.N = G.N;
  P.neighbours = G.neighbours;
  P.points = points;
  P.layout2d = P.tutte_layout();
  P.orient_neighbours();
  P.faces = P.compute_faces();
  //  cout << "faces = " << P.faces << "\n";
  
  return P;
}
