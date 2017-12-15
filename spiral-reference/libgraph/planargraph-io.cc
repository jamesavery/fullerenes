#include "planargraph.hh"
#include "polyhedron.hh"
#include <stdio.h>

//////////////////////////// FORMAT MULTIPLEXING ////////////////////////////
vector<string> PlanarGraph::formats{{"ascii","planarcode","xyz","mol2","mathematica","latex"}};
vector<string> PlanarGraph::input_formats{{"planarcode","xyz","mol2"}}; // TODO: Add ASCII
vector<string> PlanarGraph::output_formats{{"ascii","planarcode"}}; // TODO: Add LaTeX, Mathematica

int PlanarGraph::format_id(string name)
{
  for(int i=0;i<formats.size();i++) if(name == formats[i]) return i;
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
}



////////////////////////////// OUTPUT ROUTINES //////////////////////////////

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
  if(N == 0){ number_length=2; N = read_int(); }
  
  Graph g(N,true);
  for(node_t u=0; u<N && !feof(file); ++u){
    int v=0;
    do{
      v = read_int();
      if(v!=0) g.neighbours[u].push_back(v-1); // In oriented order
    } while(v!=0 && !feof(file));
  }
  
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
  
  //Assumes planarcode files containing only graphs of equal size
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
    
  fputs(">>planar_code<<",file);
  if(G.N>255) fputc(0,file);

  write_int(G.N);

  for(uint16_t u=0;u<G.N;u++){
    for(uint16_t v: G.neighbours[u])
      write_int(v);
    write_int(0);
  }

  return ferror(file) == 0;
}

bool PlanarGraph::to_ascii(const PlanarGraph &G, FILE *file)
{
  // Neighbour list is unique representation of graph that preserves orientation.
  // N is length of list.
  string s = to_string(G.neighbours);
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

