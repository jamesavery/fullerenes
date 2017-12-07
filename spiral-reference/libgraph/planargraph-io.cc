#include "planargraph.hh"

vector<string> PlanarGraph::input_formats_txt{{"ascii","planarcode","xyz","mol2"}};
vector<string> PlanarGraph::output_formats_txt{{"ascii","mathematica","planarcode"}};


int output_format_id(string name)
{
  for(int i=0;i<output_formats.size();i++) if(name == output_formats[i]) return i;
  return -1;
}


static bool PlanarGraph::to_file(const PlanarGraph &G, string path)
{
  switch(output_format_id(output_format)){
  case ASCII:
    PlanarGraph::to_ascii(G,stdout);
    break;
  case MATHEMATICA:
    // TODO: to_mathematica
    cout << G << "\n";
    break;
  case PLANARCODE:
    PlanarGraph::to_planarcode(G,stdout);
    break;
  case XYZ:
  case MOL2:
  default:
    break
  }

}

// TODO: Where does this belong?
// Assumes file is at position of a graph start
PlanarGraph PlanarGraph::read_hog_planarcode(FILE *planarcode_file)
{
  // Read the number N of vertices per graph.
  int number_length=1, N=0;
  fread(reinterpret_cast<unsigned char*>(&N), 1, 1, planarcode_file);
  if(N == 0){
    fread(reinterpret_cast<unsigned char*>(&N), 2, 1, planarcode_file);
    number_length=2;
  }
  
  Graph g(N,true);
  for(node_t u=0; u<N && !feof(planarcode_file); ++u){
    int v=0;
    do{
      int n_read = fread(reinterpret_cast<char*>(&v), number_length, 1, planarcode_file);
      if(n_read != 1 && !feof(planarcode_file)){
        perror("Error reading HoG PlanarCode file: ");
        abort();
      }
      if(v!=0) g.neighbours[u].push_back(v-1); // In oriented order
    } while(v!=0 && !feof(planarcode_file));
  }
  // Check graph
  for(node_t u=0;u<N;u++){
    for(auto v: g.neighbours[u]){
      bool found_vu = false;
      
      for(node_t w: g.neighbours[v])
        if(w == u) found_vu = true;
      if(!found_vu){
        fprintf(stderr,"Graph is not symmetric: (u,v) = (%d,%d) has\n",u,v);
        cerr << "neighbours["<<u<<"] = " << g.neighbours[u] <<";\n";
        cerr << "neighbours["<<v<<"] = " << g.neighbours[v] <<";\n";
        abort();
      }
    }
  }

  return g;
}


// TODO: Read only a range
vector<PlanarGraph> PlanarGraph::read_hog_planarcodes(FILE *planarcode_file) {
  const int header_size = 15;
  vector<PlanarGraph> graph_list;

  //the actual parsing of the selected graph:
  //go to the beginning of the selected graph
  fseek(planarcode_file,  header_size, SEEK_SET);

  //  int i = 1;
  while(!feof(planarcode_file)){
    //    cerr << "Reading graph " << (i++) << ".\n";
    Graph g = read_hog_planarcode(planarcode_file);
    //    cerr << "Got graph on " << g.N << " vertices.\n";
    if(g.N != 0){
      graph_list.push_back(g);
    }
  }
    
  return graph_list;
}

// Parse House of Graphs planarcode (not restricted to cubic graphs)
PlanarGraph PlanarGraph::from_planarcode(FILE* file, const size_t index){
  const int header_size = 15;

  // Get file size
  fseek(file, 0, SEEK_END);
  size_t file_size = ftell(file);

  //find number of vertices per graph
  //this only works for files with graphs of equal size
  fseek(file, header_size, SEEK_SET);

  //Assumes planarcode files containing only graphs of equal size
  size_t step = 0;
  if(index != 0){
    Graph first(read_hog_planarcode(file));
    step = ftell(file);
  }

  size_t address = header_size + step * index;

  //check if selected graphnumber is valid
  unsigned int graphs_per_file = (file_size - header_size ) /step;
  if(graphs_per_file -1 < index){
    cerr << "You asked for the " << index+1 << "th graph, but there are only " << graphs_per_file << " stored in this file." << std::endl;
    abort();
  }

  //Find the beginning of the selected graph and read it
  fseek(file, address+1, SEEK_SET);
  return read_hog_planarcode(file);
}
