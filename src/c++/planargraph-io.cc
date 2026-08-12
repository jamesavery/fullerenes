#include "fullerenes/planargraph.hh"
#include "fullerenes/polyhedron.hh"
#include <stdio.h>
#include <cctype>
#include <stdexcept>

//////////////////////////// FORMAT MULTIPLEXING ////////////////////////////
vector<string> PlanarGraph::formats{"ascii","planarcode","xyz","mol2","mathematica","latex","spiral"};
vector<string> PlanarGraph::input_formats{"planarcode","xyz","mol2","spiral"}; // TODO: Add ASCII
vector<string> PlanarGraph::output_formats{"ascii","planarcode","spiral"}; // TODO: Add LaTeX, Mathematica
 
int PlanarGraph::format_id(string name)
{
  cout << "input_formats = " << input_formats << endl;
  for(int i=0;i<formats.size();i++) if(name == formats[i]) return i; else cout << "Format is not " << formats[i] << endl;
  return -1;
}

PlanarGraph PlanarGraph::from_file(FILE *file, string format, int index)
{
  switch(format_id(format)){
  case PLANARCODE:
    return from_planarcode(file,index);
  case XYZ:
    return Polyhedron::from_xyz(file);
  case MOL2:
    return Polyhedron::from_mol2(file);
  case SPIRAL:
    return from_spiral(file, index); 
  default:
    cerr << "Input format is '" << format << "'; must be one of: " << input_formats << "\n";
    abort();
  }
}


bool PlanarGraph::to_file(const PlanarGraph &G, FILE *file, string format)
{
  switch(format_id(format)){
  case ASCII:
    return PlanarGraph::to_ascii(G,file);
  case PLANARCODE:
    return PlanarGraph::to_planarcode(G,file);
  default:
    cerr << "Output format is '" << format << "'; must be one of: " << output_formats << "\n";
    return false;
  }
}

PlanarGraph PlanarGraph::from_file(string filename, int index)
{
  FILE *file = fopen(filename.c_str(),"rb");
  string extension = filename_extension(filename);
  PlanarGraph G = from_file(file,extension);
  fclose(file);
  return G;
}

bool PlanarGraph::to_file(const PlanarGraph &G, string filename)
{
  FILE *file = fopen(filename.c_str(),"wb");
  string extension = filename_extension(filename);  
  to_file(G,file,extension);
  fclose(file);
  return true;			// TODO: Check success
}



////////////////////////////// OUTPUT ROUTINES //////////////////////////////
bool PlanarGraph::to_spiral(const PlanarGraph &P, FILE *file)  {
  // TODO: Make general! Autodetect fullerene, triangulation, cubic, etc.
  // TODO: Move to to_file, add as parameter.
  auto naming_scheme = spiral_nomenclature::FULLERENE;
  auto construction_scheme = spiral_nomenclature::CUBIC;
  spiral_nomenclature name(P,naming_scheme,construction_scheme);  
  string s = name.to_string()+"\n";
  fputs(s.c_str(),file);
  return ferror(file) == 0;
}


// TODO: Where does this belong?
// Assumes file is at position of a graph start
PlanarGraph PlanarGraph::read_hog_planarcode(FILE *file)
{
  // Read the number N of vertices per graph.
  int number_length=1, N=0;
  auto read_int = [&]() -> int {
    int x = fgetc(file);
    if(number_length==2) x |= (fgetc(file) << 8);
    return x;
  };

  N = read_int();
  if(feof(file)) return Graph();
  if(N == 0){ number_length=2; N = read_int(); }
  if(feof(file) || N <= 0) return Graph();

  // First pass: collect adjacency lists to find max degree.
  vector<vector<node_t>> adj_lists(N);
  for(node_t u=0; u<N && !feof(file); ++u){
    int v=0;
    do{
      v = read_int();
      if(v!=0) adj_lists[u].push_back(v-1); // In oriented order
    } while(v!=0 && !feof(file));
  }
  int dmax = GRAPH_DMAX;
  for(int u=0; u<N; u++) dmax = std::max(dmax, (int)adj_lists[u].size());

  Graph adj(N, dmax);
  for(int u=0; u<N; u++)
    for(node_t v : adj_lists[u]) adj.push_back(u, v);
  Graph g(adj);
  
  // Check graph. TODO: does this belong here?
  // for(node_t u=0;u<N;u++){
  //   for(auto v: g.neighbours[u]){
  //     bool found_vu = false;
      
  //     for(node_t w: g.neighbours[v])
  //       if(w == u) found_vu = true;
  //     if(!found_vu){
  //       fprintf(stderr,"Graph is not symmetric: (u,v) = (%d,%d) has\n",u,v);
  //       cerr << "neighbours["<<u<<"] = " << g.neighbours[u] <<";\n";
  //       cerr << "neighbours["<<v<<"] = " << g.neighbours[v] <<";\n";
  //       abort();
  //     }
  //   }
  //}

  return g;
}


bool PlanarGraph::read_hog_metadata(FILE *file, size_t &graph_count, size_t &graph_size)
{
  const int header_size = 15;

  // Get file size
  size_t file_pos  = ftell(file);
  fseek(file, 0, SEEK_END);
  size_t file_size = ftell(file);

  //find number of vertices per graph
  //this only works for files with graphs of equal size
  fseek(file, header_size, SEEK_SET);

  // ASSUMPTION: All graphs in the file have equal record size (same N, same degree sequence).
  // This holds for single-N planarcode files (e.g. all C60 isomers from buckygen/plantri).
  // Random access via from_planarcode(file, index) will give wrong results for mixed-N files.
  Graph first(read_hog_planarcode(file));
  graph_size = ftell(file)-header_size;

  //check if selected graphnumber is valid
  graph_count = (file_size - header_size ) / graph_size;
  fseek(file,file_pos,SEEK_SET); // Back to where we started

  return (ferror(file) == 0);
}


// TODO: Read only a range
vector<PlanarGraph> PlanarGraph::read_hog_planarcodes(FILE *file) {
  const int header_size = 15;
  vector<PlanarGraph> graph_list;

  //the actual parsing of the selected graph:
  //go to the beginning of the selected graph
  fseek(file,  header_size, SEEK_SET);

  //  int i = 1;
  while(!feof(file)){
    //    cerr << "Reading graph " << (i++) << ".\n";
    Graph g = read_hog_planarcode(file);
    //    cerr << "Got graph on " << g.N << " vertices.\n";
    if(g.N != 0){
      graph_list.push_back(g);
    }
  }
    
  return graph_list;
}


// Write House of Graphs planarcode
bool PlanarGraph::to_planarcode(const PlanarGraph &G, FILE *file)
{
  auto write_int = [&](uint16_t x){ fputc(x&0xff,file); if(G.N>255) fputc((x>>8)&0xff,file); };
    
  // TODO: The header should only be written once per file, not per graph.
  // Callers writing multiple graphs should write the header separately.
  fputs(">>planar_code<<",file);
  if(G.N>255) fputc(0,file);

  write_int(G.N);

  for(uint16_t u=0;u<G.N;u++){
    for(uint16_t v: G.nbrs(u))
      write_int(v+1);		// planar_code is 1-indexed; 0 is the terminator
    write_int(0);
  }

  return ferror(file) == 0;
}

bool PlanarGraph::to_ascii(const PlanarGraph &G, FILE *file)
{
  // Neighbour list is unique representation of graph that preserves orientation.
  // N is length of list.
  string s = to_string(static_cast<const neighbours_t&>(G));
  fputs(s.c_str(),file);

  return ferror(file) == 0;
}

bool PlanarGraph::to_mathematica(const PlanarGraph &G, FILE *file)
{
  ostringstream s;
  s << G << "\n";
  fputs(s.str().c_str(),file);

  return ferror(file) == 0;
}

////////////////////////////// INPUT ROUTINES //////////////////////////////

PlanarGraph PlanarGraph::from_ascii(FILE *file)
{
  // Inverse of to_ascii. Reads one neighbour-list record of the form
  //
  //   [[a,b,c,...],[a,b,c,...],...]
  //
  // from the current stream position. Whitespace is tolerated between
  // tokens (the writer emits none, but a human-edited file may have it).
  // Reads up to and including the matching outer ']'; subsequent
  // characters (newline, separator, EOF) are left in the stream so a
  // caller may iterate multiple records.
  //
  // NOT AN ORIENTATION BOUNDARY -- a TRUSTED FORMAT, like planarcode.  The
  // record IS the rotation system: to_ascii writes the rows in their stored
  // cyclic order precisely so that "neighbour list is a unique representation of
  // the graph that preserves orientation" (its comment above), and this reads
  // them back verbatim.  Re-sorting them here would be the opposite of reading
  // the file -- it would discard the one thing the format exists to carry, and
  // (nothing having coordinates) could only put back the mirror embedding.
  //
  // Nor does it VERIFY the genus, though a real predicate now exists.  This
  // format's second job is diagnostic: the claude-projects validator harness
  // dumps the graph of a failing isomer through to_ascii so the failure can be
  // reloaded and studied, and a dump whose bug IS its rotation system must
  // still load.  What is checked below is what a record can be wrong about
  // regardless of any topological claim -- an index naming no vertex, or an arc
  // with no reverse, neither of which is a graph at all.
  //
  // @post   the returned rows are, in order, the rows in the file
  // @pre    oriented: the file's rows are a consistent orientation.  A caller
  //         that does not trust its source checks with oriented_surface() /
  //         require_oriented_surface(); this function will not check for it.
  // @throws std::runtime_error on EOF before the record, on unterminated input,
  //         on a parse failure, on an out-of-range vertex index, or on
  //         asymmetric adjacency

  // 1) Consume the record into a buffer by counting bracket balance.
  std::string buf;
  int balance = 0;
  bool started = false;
  for (int c; (c = fgetc(file)) != EOF; ) {
    if (!started) {
      if (std::isspace(c)) continue;
      if (c != '[') {
        throw std::runtime_error(
          std::string("PlanarGraph::from_ascii: expected '[', got '") +
          char(c) + "'");
      }
      started = true;
      balance = 1;
      buf.push_back(char(c));
      continue;
    }
    buf.push_back(char(c));
    if      (c == '[') balance++;
    else if (c == ']') { if (--balance == 0) break; }
  }
  if (!started)     throw std::runtime_error("PlanarGraph::from_ascii: EOF before record");
  if (balance != 0) throw std::runtime_error("PlanarGraph::from_ascii: unterminated record");

  // 2) Parse the buffered string.
  std::vector<std::vector<node_t>> adj;
  size_t i = 1;  // skip outer '['
  auto skip_ws = [&](){ while (i < buf.size() && std::isspace((unsigned char)buf[i])) i++; };

  skip_ws();
  if (i < buf.size() && buf[i] != ']') {
    while (true) {
      skip_ws();
      if (i >= buf.size() || buf[i] != '[')
        throw std::runtime_error("PlanarGraph::from_ascii: expected '[' for row");
      i++;
      std::vector<node_t> row;
      skip_ws();
      if (i < buf.size() && buf[i] != ']') {
        while (true) {
          skip_ws();
          size_t j = i;
          if (j < buf.size() && (buf[j] == '-' || buf[j] == '+')) j++;
          if (j >= buf.size() || !std::isdigit((unsigned char)buf[j]))
            throw std::runtime_error("PlanarGraph::from_ascii: expected integer in row");
          while (j < buf.size() && std::isdigit((unsigned char)buf[j])) j++;
          row.push_back(node_t(std::stol(buf.substr(i, j - i))));
          i = j;
          skip_ws();
          if (i < buf.size() && buf[i] == ',') { i++; continue; }
          if (i < buf.size() && buf[i] == ']') break;
          throw std::runtime_error("PlanarGraph::from_ascii: expected ',' or ']' in row");
        }
      }
      i++;  // consume row's ']'
      adj.push_back(std::move(row));
      skip_ws();
      if (i < buf.size() && buf[i] == ',') { i++; continue; }
      if (i < buf.size() && buf[i] == ']') break;
      throw std::runtime_error("PlanarGraph::from_ascii: expected ',' or ']' between rows");
    }
  }

  // 3) The two things the record can be wrong about on its own terms.  The range
  // check comes first because an out-of-range index is an out-of-bounds read the
  // moment anything walks the rows -- including the symmetry check below.
  const int N = int(adj.size());
  for(int u=0;u<N;u++)
    for(node_t v: adj[u])
      if(v < 0 || v >= N)
        throw std::runtime_error("PlanarGraph::from_ascii: row " + std::to_string(u)
                                 + " names vertex " + std::to_string(int(v))
                                 + ", which is outside [0," + std::to_string(N) + ")");

  PlanarGraph G{Graph(Spanify::OwnedDenseGraph<node_t>(adj))};
  if(!G.adjacency_is_symmetric())
    throw std::runtime_error("PlanarGraph::from_ascii: adjacency is not symmetric -- some arc "
                             "u->v in the record has no reverse v->u, so it is not a graph");
  return G;
}

PlanarGraph PlanarGraph::from_spiral(FILE *f, const size_t index)  {
  // TODO: Make general! Autodetect fullerene, triangulation, cubic, etc.
  // TODO: Move to to_file, add as parameter.
  auto naming_scheme = spiral_nomenclature::FULLERENE;
  auto construction_scheme = spiral_nomenclature::CUBIC;
  string line;
  bool read_ok;
  printf("from_spiral()\n");
  for(int i=0;i<=index;i++){
    read_ok = getline(f,line);
    if(!read_ok){
      perror("getline");
      return PlanarGraph{};
    } // else {
    //   printf("result: %s, line: %s\n",read_ok,line);
    // }
  }
  
  spiral_nomenclature name(line);  
  return FullereneGraph(name).halma_fullerene(1);
}


// Parse House of Graphs planarcode (not restricted to cubic graphs)
PlanarGraph PlanarGraph::from_planarcode(FILE* file, const size_t index){
  const int header_size = 15;

  size_t graph_count = 0, graph_size = 0;
  read_hog_metadata(file,graph_count,graph_size);

  size_t address = header_size + graph_size * index;
  //check if selected graphnumber is valid
  if(graph_count-1 < index){
    cerr << "You asked for the " << index+1 << (index==0?"st":(index==1?"nd":"th"))<<" graph, but there "
	 <<(graph_count==1?"is":"are")<<" only"
	 << graph_count << " stored in this file." << std::endl;
    abort();
  }
  //Find the beginning of the selected graph and read it
  fseek(file, address, SEEK_SET);
  return read_hog_planarcode(file);
}

////////////////////////////// AUXILIARY //////////////////////////////

string filename_extension(const string& filename)
{
  size_t i = filename.rfind(".");
  bool found = i != string::npos;
  if(found) 
    return filename.substr(i+1,filename.size());
  else 
    return "";
}

bool getline(FILE *file, string& str){
  char *line_ptr     = 0;
  size_t line_length = 0;
  ssize_t error = getline(&line_ptr,&line_length,file);
  //  fprintf(stderr,"error = %ld, line_length = %ld, line = %s\n",error,line_length,line_ptr);
  //  str = string(line_ptr,line_length);
  str = string(line_ptr);
  free(line_ptr);
  return error > 0;
}


