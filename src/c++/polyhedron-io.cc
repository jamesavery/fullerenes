#include "fullerenes/polyhedron.hh"
#include "fullerenes/layout2d.hh"

#include <algorithm>
#include <array>
#include <cstdint>
#include <cstring>
#include <map>
#include <span>
#include <unordered_map>

//////////////////////////// FORMAT MULTIPLEXING ////////////////////////////
vector<string> Polyhedron::formats{"ascii","planarcode","xyz","mol2","mathematica","latex","cc1","turbomole","gaussian","wavefront","spiral","ply"};
vector<string> Polyhedron::format_alias{"txt","pc","xyz","mol2","m","tex","cc1","turbomole","com","obj","rspi","ply"};
vector<string> Polyhedron::input_formats{"xyz","mol2","ply"}; // TODO: "ascii","planarcode","obj"
vector<string> Polyhedron::output_formats{"ascii","xyz","mol2","cc1","turbomole","gaussian","spiral","ply"};

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
  case PLY:
    return from_ply(file);
  default: {
    // Reject a genuinely-unknown format here: PlanarGraph::from_file still
    // abort()s on its own default, which would bypass the catchable mesh_io_error
    // contract this reader promises.
    if(PlanarGraph::format_id(format) < 0){
      ostringstream msg;
      msg << "Polyhedron::from_file: unknown format '"<<format<<"'; must be one of: "
          << input_formats << " or " << PlanarGraph::input_formats;
      throw mesh_io_error(mesh_io_error::Code::UnknownFormat, msg.str());
    }
    PlanarGraph G = PlanarGraph::from_file(file,format);

    if(!G.N){
      ostringstream msg;
      msg << "Polyhedron::from_file: format '"<<format<<"' parsed an empty graph";
      throw mesh_io_error(mesh_io_error::Code::EmptyGraph, msg.str());
    } else {
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
  case PLY:
    return Polyhedron::to_ply(G,file);
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
  string s = LIST_OPEN + to_string(static_cast<const neighbours_t&>(P)) + "," + to_string(P.points) + "," + to_string(P.faces()) + LIST_CLOSE;
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
    auto nu = P.nbrs(u);
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

  for(auto f: P.faces()){
    fprintf(file,"f ");
    for(auto v: f) fprintf(file,"%d ",v);
    fprintf(file,"\n");
  }

  return true;
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
    for(node_t v: P.nbrs(u))
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
    for(node_t v: P.nbrs(u))
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
    vector<coord2d> dual_layout = dual.tutte_layout();
    s << "\\foreach \\place/\\name/\\lbl in {";
    for(node_t u=0;u<dual.N;u++){
      const coord2d& xs(dual_layout[u]);
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

  {
    vector<coord2d> pov_layout = tutte_layout();
    s << layout2d::to_povray(*this, pov_layout, w_cm,h_cm,line_colour,vertex_colour,line_width,vertex_diameter);
  }
  s << "#declare layout3D=array["<<N<<"][3]" << points <<";\n\n";

  auto fs = faces();
  s << "#declare faces   =array["<<fs.size()<<"]["<<(face_max+1)<<"]{";
  for(int i=0;i<int(fs.size());i++) {
    const face_t& f(fs[i]);
    s << "{";
    for(int j=0;j<int(f.size());j++) s << f[j] << ",";
    for(int j=f.size();j<face_max;j++) s << "-1,";
    s << "-1}" << (i+1<int(fs.size())? ",":"}\n\n");
  }
  s << "#declare facelength=array["<<fs.size()<<"]{";for(int i=0;i<int(fs.size());i++) s<< fs[i].size() << (i+1<int(fs.size())?",":"}\n\n");


  vector<tri_t>   tris(centroid_triangulation(fs));
  vector<int>     triface;
  vector<coord3d> centroid_points(points.begin(),points.end());
  vector<coord3d> trinormals(tris.size()), facenormals(fs.size()), vertexnormals(points.size()+fs.size());

  for(int i=0;i<int(fs.size());i++)
    centroid_points.push_back(fs[i].centroid(points));

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

  for(int i=0;i<int(fs.size());i++) {
    coord3d normal;
    if(fs[i].size()>3){
      for(int j=0;j<int(fs[i].size());j++){
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
  s << "#declare facenormals=array["<<fs.size()<<"][3]" << facenormals << ";\n\n";

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
  static_cast<Owned<PolyhedronView<double>>&>(P) = static_cast<const GraphView&>(G);
  P.owned_points = std::move(points);
  P.repoint();
  {
    vector<coord2d> layout = P.tutte_layout();
    layout2d::orient_neighbours(P, layout);
  }
  // faces are now computed on demand via P.faces()
  //  cout << "faces = " << P.faces << "\n";
  
  return P;
}

size_t file_size(const char *filename) {
  struct stat st;

  if (stat(filename, &st) == 0)
    return st.st_size;

  return -1;
}
size_t file_size(const string filename) { return file_size(filename.c_str()); }


////////////////////////////// PLY MESH I/O //////////////////////////////
// PLY (Stanford polygon) read/write for arbitrary polygon meshes. A PLY face is
// a per-face (count, indices) record, so triangles and n-gons serialise the same
// way; the file's consistent face winding IS the embedding orientation, so the
// oriented graph is built directly (no orient_neighbours / 2D layout). The raw
// read_ply/write_ply + the graph-from-faces builder are file-static; the public
// surface is Polyhedron::{from,to}_ply (and Deltahedron's in deltahedron-io.cc).
namespace {

using Code = mesh_io_error::Code;   // local alias so the throw sites below read cleanly

// A PLY scalar type, resolved from its header name once (parse_ply_scalar_type)
// so the per-value read is an integer switch, not a string compare. `unknown`
// is a recognised-as-unrecognised tag: ascii ignores the type (so it stays
// harmless), binary rejects it at read time -- preserving the prior behaviour.
enum class ply_scalar { i8, u8, i16, u16, i32, u32, f32, f64, unknown };

ply_scalar parse_ply_scalar_type(const string& t) {
  if (t=="char"   || t=="int8")    return ply_scalar::i8;
  if (t=="uchar"  || t=="uint8")   return ply_scalar::u8;
  if (t=="short"  || t=="int16")   return ply_scalar::i16;
  if (t=="ushort" || t=="uint16")  return ply_scalar::u16;
  if (t=="int"    || t=="int32")   return ply_scalar::i32;
  if (t=="uint"   || t=="uint32")  return ply_scalar::u32;
  if (t=="float"  || t=="float32") return ply_scalar::f32;
  if (t=="double" || t=="float64") return ply_scalar::f64;
  return ply_scalar::unknown;
}

int ply_scalar_bytes(ply_scalar t) {
  switch (t) {
    case ply_scalar::i8:  case ply_scalar::u8:                     return 1;
    case ply_scalar::i16: case ply_scalar::u16:                    return 2;
    case ply_scalar::i32: case ply_scalar::u32: case ply_scalar::f32: return 4;
    case ply_scalar::f64:                                         return 8;
    case ply_scalar::unknown:                                     return 0;
  }
  return 0;   // unreachable; all enumerators handled above
}

// The parsed PLY header: body format, element counts, the vertex property list
// (resolved type, name) with the indices of x/y/z within it, and the face-list
// element types. All types are resolved to ply_scalar once, at header parse.
struct ply_header {
  bool binary;
  int  nv, nf;
  vector<pair<ply_scalar,string>> vprops;      // (resolved type, name) in file order
  int  xi, yi, zi;                             // positions of x,y,z within vprops
  ply_scalar face_count_type, face_index_type;
};

// One scalar read for both body formats: binary decodes the resolved type at its
// true width and signedness; ascii reads a whitespace-separated number (type
// irrelevant). An unknown type in a binary body is a hard error, never a silent skip.
double read_ply_scalar(FILE* file, ply_scalar type, bool binary) {
  if (!binary) {
    double v; if (fscanf(file, "%lf", &v) != 1) throw mesh_io_error(Code::MalformedFile, "read_ply: truncated body");
    return v;
  }
  int sz = ply_scalar_bytes(type);
  if (sz == 0) throw mesh_io_error(Code::UnsupportedFormat, "read_ply: unknown binary property type");
  unsigned char buf[8];
  if (fread(buf,1,sz,file) != (size_t)sz) throw mesh_io_error(Code::MalformedFile, "read_ply: truncated body");
  switch (type) {
    case ply_scalar::f32: { float    v; memcpy(&v,buf,4); return v; }
    case ply_scalar::f64: { double   v; memcpy(&v,buf,8); return v; }
    case ply_scalar::i8:  { int8_t   v; memcpy(&v,buf,1); return v; }
    case ply_scalar::u8:  { uint8_t  v; memcpy(&v,buf,1); return v; }
    case ply_scalar::i16: { int16_t  v; memcpy(&v,buf,2); return v; }
    case ply_scalar::u16: { uint16_t v; memcpy(&v,buf,2); return v; }
    case ply_scalar::i32: { int32_t  v; memcpy(&v,buf,4); return v; }
    case ply_scalar::u32: { uint32_t v; memcpy(&v,buf,4); return v; }
    case ply_scalar::unknown: ;   // unreachable: sz==0 already threw above
  }
  return 0;   // unreachable; all readable enumerators return above
}

// Read header lines up to end_header. Throws on an unsupported body format, a
// negative element count, or a vertex element lacking x/y/z.
ply_header parse_ply_header(FILE* file) {
  if (!file) throw mesh_io_error(Code::NullFile, "read_ply: null file handle");

  string fmt;
  ply_header h{}; h.face_count_type = ply_scalar::u8; h.face_index_type = ply_scalar::i32;
  bool in_vertex = false, in_face = false;
  char line[256];
  while (fgets(line, sizeof line, file)) {
    char a[64]={0}, b[64]={0}, c[64]={0}, d[64]={0};
    int n = sscanf(line, "%63s %63s %63s %63s", a, b, c, d);
    string t0 = a;
    if (t0 == "end_header") break;
    if (t0 == "format") fmt = b;
    else if (t0 == "element" && strcmp(b,"vertex")==0) { h.nv = atoi(c); in_vertex=true; in_face=false; }
    else if (t0 == "element" && strcmp(b,"face")==0)   { h.nf = atoi(c); in_vertex=false; in_face=true; }
    else if (t0 == "property" && in_vertex && n>=3) h.vprops.push_back({parse_ply_scalar_type(b), c});
    else if (t0 == "property" && in_face && strcmp(b,"list")==0 && n>=4) {
      h.face_count_type = parse_ply_scalar_type(c); h.face_index_type = parse_ply_scalar_type(d);
    }
  }

  h.binary = (fmt == "binary_little_endian");
  if (fmt != "ascii" && !h.binary) throw mesh_io_error(Code::UnsupportedFormat, "read_ply: unsupported format '" + fmt + "'");
  if (h.nv < 0 || h.nf < 0) throw mesh_io_error(Code::MalformedFile, "read_ply: negative element count");
  h.xi = h.yi = h.zi = -1;
  for (int i=0;i<(int)h.vprops.size();i++) {
    if (h.vprops[i].second=="x") h.xi=i; else if (h.vprops[i].second=="y") h.yi=i; else if (h.vprops[i].second=="z") h.zi=i;
  }
  if (h.xi<0||h.yi<0||h.zi<0) throw mesh_io_error(Code::UnsupportedFormat, "read_ply: missing x/y/z vertex properties");
  return h;
}

// Read h.nv vertices, keeping x/y/z and skipping any other scalar vertex properties.
void read_ply_vertices(FILE* file, const ply_header& h, vector<coord3d>& points) {
  // A bogus vertex count would otherwise throw std::length_error/bad_alloc here
  // (not a mesh_io_error); convert it to the file-failure channel.
  try { points.assign(h.nv, coord3d()); }
  catch (const std::exception&) { throw mesh_io_error(Code::MalformedFile, "read_ply: vertex count too large to allocate"); }
  for (int i=0;i<h.nv;i++) {
    double x=0,y=0,z=0;
    for (int k=0;k<(int)h.vprops.size();k++) {
      double val = read_ply_scalar(file, h.vprops[k].first, h.binary);
      if (k==h.xi) x=val; else if (k==h.yi) y=val; else if (k==h.zi) z=val;
    }
    points[i] = coord3d(x,y,z);
  }
}

// Read h.nf polygon faces. No nf pre-reserve: a too-large nf simply exhausts the
// body and trips read_ply_scalar's truncated-body guard.
void read_ply_faces(FILE* file, const ply_header& h, vector<face_t>& faces) {
  faces.clear();
  for (int i=0;i<h.nf;i++) {
    int cnt = (int)read_ply_scalar(file, h.face_count_type, h.binary);
    if (cnt < 3 || cnt > h.nv) throw mesh_io_error(Code::InvalidTopology, "read_ply: face vertex count out of range [3, nv]");
    face_t f(cnt);
    for (int k=0;k<cnt;k++) {
      int idx = (int)read_ply_scalar(file, h.face_index_type, h.binary);
      if (idx < 0 || idx >= h.nv)            // never trust a raw index: an out-of-range
        throw mesh_io_error(Code::InvalidTopology, "read_ply: face vertex index out of range");  // one is an OOB write downstream
      f[k] = idx;
    }
    faces.push_back(std::move(f));
  }
}

// Parse a PLY (ascii or binary_little_endian) into its x,y,z vertices and polygon
// face list. Throws mesh_io_error on an unsupported header, missing x/y/z, or a
// truncated body.
void read_ply(FILE* file, vector<coord3d>& points, vector<face_t>& faces) {
  ply_header h = parse_ply_header(file);
  read_ply_vertices(file, h, points);
  read_ply_faces(file, h, faces);
}

// Oriented neighbour lists (the rotation system) from a consistently-wound
// polygon face list over vertices 0..N-1: for each boundary arc x -> v -> y the
// neighbour of v after x is y; chaining around each vertex yields its fan. The
// chain is then REVERSED to match the library's rotation convention: face tracing
// follows next_on_face(u,v) = prev(v,u), so reproducing an input face ...u->v->w...
// requires w to sit immediately BEFORE u in nb[v] (prev(v,u)=w), i.e. the reverse
// of the after-chain order [u, after{u,v}=w, ...]. (A symmetric mesh like K4 hides
// the direction; a general mesh does not -- a wrong rotation makes compute_faces
// trace non-closing faces and spin forever.) Satisfies the orientation invariant
// by construction (no orient_neighbours); returned as owned lists since
// neighbours_t is a non-owning view.
vector<vector<node_t>> oriented_neighbours_from_faces(int N, const vector<face_t>& faces) {
  unordered_map<arc_t,int> after;             // arc (x->v) -> next vertex y; O(1) lookup
  vector<int> start(N, -1), incident(N, 0);   // incident[v] = #faces (= arcs) at v
  for (const auto& f : faces) {
    int k = (int)f.size();
    for (int i=0;i<k;i++) {
      int x = f[i], v = f[(i+1)%k], y = f[(i+2)%k];
      if (!after.insert({{x,v}, y}).second)   // a directed arc in two faces == non-manifold
        throw mesh_io_error(Code::InvalidTopology, "from_ply: directed edge shared by two faces (non-manifold or inconsistent winding)");
      incident[v]++;
      if (start[v] < 0) start[v] = x;
    }
  }
  vector<vector<node_t>> nb(N);
  for (int v=0; v<N; v++) {
    if (start[v] < 0) continue;               // unreferenced vertex; caught by check_degrees
    int x = start[v];
    bool closed = false;
    for (int step=0; step <= incident[v]; step++) {
      nb[v].push_back(x);
      auto it = after.find({x, v});
      if (it == after.end()) break;           // boundary edge: the fan has an open end
      x = it->second;
      if (x == start[v]) { closed = true; break; }
    }
    // The rotation must close having visited every incident face -- else the
    // vertex is a boundary (hole) or a non-manifold pinch (two fans meet at v).
    if (!closed || (int)nb[v].size() != incident[v])
      throw mesh_io_error(Code::InvalidTopology, "from_ply: vertex fan does not close (boundary hole or non-manifold pinch)");
    std::reverse(nb[v].begin(), nb[v].end());   // after-chain order -> library prev-convention
  }
  return nb;
}

// Weld vertices that share an EXACT position -- a duplicate-vertex defect that a
// PLY can carry (two vertex records at the same coordinates, woven into the
// triangulation by a ZERO-LENGTH edge and its two zero-area faces). Left in, the
// zero-length edge makes the intrinsic metric non-Euclidean (a + 0 <= a on the
// degenerate face), which every metric consumer must then reject. Weld remaps
// each face to the canonical (first-seen) index at its position, drops faces
// that collapse to a repeated vertex, and compacts the point list; the rotation
// builder below then validates that the collapsed mesh is still a closed
// manifold (it throws otherwise -- a genuinely broken mesh, not a stray
// duplicate). A clean mesh has no coincident vertices, so this is a no-op
// (returns 0) and the point/face arrays are untouched. Returns #vertices merged.
int weld_coincident_vertices(vector<coord3d>& points, vector<face_t>& faces) {
  std::map<std::array<double, 3>, int> canon;    // position -> canonical (lowest) index
  vector<int> to_canon(points.size());
  int merged = 0;
  for (int i = 0; i < (int)points.size(); i++) {
    const std::array<double, 3> key{points[i][0], points[i][1], points[i][2]};
    const auto [it, fresh] = canon.emplace(key, i);
    to_canon[i] = it->second;
    if (!fresh) merged++;
  }
  if (merged == 0) return 0;                      // clean mesh: leave everything as read

  // Relabel faces to canonical indices; drop any that now repeat a vertex (the
  // zero-area faces of a collapsed edge).
  vector<face_t> kept;
  kept.reserve(faces.size());
  for (const face_t& f : faces) {
    face_t g;
    g.reserve(f.size());
    for (int v : f) g.push_back(to_canon[v]);
    bool degen = false;
    for (size_t a = 0; a + 1 < g.size() && !degen; a++)
      for (size_t b = a + 1; b < g.size(); b++)
        if (g[a] == g[b]) { degen = true; break; }
    if (!degen) kept.push_back(std::move(g));
  }

  // Compact: renumber the still-referenced vertices to a dense [0, M) in
  // first-appearance order, rebuilding the point list to match.
  vector<int> new_index(points.size(), -1);
  vector<coord3d> new_points;
  for (face_t& g : kept)
    for (int& v : g)
      if (new_index[v] < 0) {
        new_index[v] = (int)new_points.size();
        new_points.push_back(points[v]);
      }
  for (face_t& g : kept)
    for (int& v : g) v = new_index[v];

  points = std::move(new_points);
  faces  = std::move(kept);
  return merged;
}

// Flatten owned adjacency into the row-major CSR an RSRAdjacencyView wraps.
// `flat`/`deg` must outlive the returned view.
neighbours_t csr_view(const vector<vector<node_t>>& nb, vector<node_t>& flat, vector<uint8_t>& deg) {
  const int N = (int)nb.size();
  int dmax = 1;
  for (const auto& l : nb) dmax = std::max(dmax, (int)l.size());
  flat.assign((size_t)N*dmax, 0);
  deg.assign(N, 0);
  for (int v=0; v<N; v++) {
    deg[v] = (uint8_t)nb[v].size();
    for (int k=0;k<(int)nb[v].size();k++) flat[(size_t)v*dmax + k] = nb[v][k];
  }
  return neighbours_t(N, dmax, std::span<node_t>(flat), std::span<uint8_t>(deg));
}

// A valid closed mesh has every vertex degree in [3, 255]: 3 is the minimum on a
// closed surface, 255 the uint8_t degree ceiling. Rejecting degree < 3 catches a
// stray unreferenced vertex (degree 0) before it produces a malformed graph.
void check_degrees(const vector<vector<node_t>>& nb, const char* who) {
  for (const auto& row : nb) {
    if (row.size() < 3)   throw mesh_io_error(Code::InvalidTopology, string(who) + ": vertex of degree < 3 (isolated or non-manifold)");
    if (row.size() > 255) throw mesh_io_error(Code::InvalidTopology, string(who) + ": vertex degree exceeds 255");
  }
}

// Signed enclosed volume * 6 (fan triangulation): its sign is the orientation
// handedness -- positive for outward-facing (CCW-on-outside) faces.
double signed_volume6(const Polyhedron& P) {
  double V = 0;
  for (const auto& f : P.faces())
    for (size_t i=1; i+1<f.size(); i++)
      V += P.points[f[0]].dot(P.points[f[i]].cross(P.points[f[i+1]]));
  return V;
}

void write_ply(FILE* file, std::span<const coord3d> points, const vector<face_t>& faces, bool binary) {
  if (!file) throw mesh_io_error(Code::NullFile, "write_ply: null file handle");
  for (const auto& f : faces)
    if (f.size() > 255) throw mesh_io_error(Code::FaceTooLarge, "write_ply: face with more than 255 vertices");
  fprintf(file, "ply\nformat %s 1.0\n", binary ? "binary_little_endian" : "ascii");
  fprintf(file, "element vertex %zu\n", points.size());
  fprintf(file, "property float x\nproperty float y\nproperty float z\n");
  fprintf(file, "element face %zu\n", faces.size());
  fprintf(file, "property list uchar int vertex_indices\nend_header\n");
  if (binary) {
    for (const auto& p : points) { float xyz[3]={(float)p[0],(float)p[1],(float)p[2]}; fwrite(xyz,4,3,file); }
    for (const auto& f : faces) {
      unsigned char n = (unsigned char)f.size(); fwrite(&n,1,1,file);
      for (int idx : f) { int32_t v = idx; fwrite(&v,4,1,file); }
    }
  } else {
    for (const auto& p : points) fprintf(file, "%.10g %.10g %.10g\n", p[0], p[1], p[2]);
    for (const auto& f : faces) {
      fprintf(file, "%zu", f.size());
      for (int idx : f) fprintf(file, " %d", idx);
      fprintf(file, "\n");
    }
  }
}

} // namespace

Polyhedron Polyhedron::from_ply(FILE *file)
{
  vector<coord3d> points;
  vector<face_t>  faces;
  read_ply(file, points, faces);
  // Repair duplicate coincident vertices (zero-length edges) before building the
  // rotation system -- a no-op on a clean mesh (see weld_coincident_vertices).
  weld_coincident_vertices(points, faces);

  vector<vector<node_t>> nb = oriented_neighbours_from_faces((int)points.size(), faces);
  check_degrees(nb, "Polyhedron::from_ply");
  vector<node_t> flat; vector<uint8_t> deg;
  neighbours_t view = csr_view(nb, flat, deg);
  PlanarGraphView G(view.N, view.dmax, view.neighbours, view.deg);
  if (!G.is_consistently_oriented())
    throw mesh_io_error(mesh_io_error::Code::InvalidTopology, "Polyhedron::from_ply: input faces are not consistently wound");

  Polyhedron P(G, points);
  if (signed_volume6(P) < 0) P.flip_all_orientations();   // normalise to outward (CCW-on-outside)
  return P;
}

bool Polyhedron::to_ply(const Polyhedron &P, FILE *file, bool binary)
{
  write_ply(file, P.points, P.faces(), binary);
  return ferror(file) == 0;
}
