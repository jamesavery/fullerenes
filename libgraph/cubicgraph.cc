#include "cubicgraph.hh"

#include <cstdlib>
#include <iostream>
#include <vector>
#include <fstream>


CubicGraph::CubicGraph(FILE *file){
    char line[0x300];
    while(!feof(file)){
      node_t n, ns[3];
      double x, y;
      char *p = fgets(line,0x2ff,file);
      if(!p){
	if(feof(file)) continue;
	else {
	  fprintf(stderr,"File read error.\n");
	  abort();
	}
      }

      int count = sscanf(line,"%d %lf %lf %d %d %d",&n,&x,&y,ns,ns+1,ns+2);

       if(count == 6){
	 // Change index start from 1 to 0
	 n--;
	 for(int i=0;i<3;i++) ns[i]--;
	 if(n>=neighbours.size()){
	   neighbours.resize(n+1);
	   layout2d.resize(n+1);
	 }
	 neighbours[n] = vector<node_t>(ns,ns+3);
	 layout2d[n] = coord2d(x,y);
       } else {			// Try
	 char dummy[0x300];
	 count = sscanf(line,"%d %s %s %d %d %d",&n,dummy,dummy,ns,ns+1,ns+2);
	 if(count == 6){
	   n--;
	   for(int i=0;i<3;i++) ns[i]--;
	   if(n>=neighbours.size())
	     neighbours.resize(n+1);
	   
	   neighbours[n] = vector<node_t>(ns,ns+3);
	 } else {
	   //	   fprintf(stderr,"Skipped line: %s\n",line);
	 }
       }
    }
    update_from_neighbours();
  }

// parse house of graphs
CubicGraph::CubicGraph(const unsigned int index, FILE* file){
  const int header_size = 15;	
	
  // Get file size
  fseek(file, 0, SEEK_END);
  size_t file_size = ftell(file);

  //find number of vertices per graph
  //this only works for files with graphs of the equal size
  fseek(file, header_size, SEEK_SET);

  // Read the number N of vertices per graph.
  fread(reinterpret_cast<char*>(&N), 1, 1, file);
  if(N == 0){
    fread(reinterpret_cast<char*>(&N), 2, 1, file);
  }

  neighbours.resize(N);

  //only for files with graphs of the equal size
  unsigned int step;
  if(N<=255)
    {step = N * 4 + 1;}
  else
    {step = N * 8 + 3;}
  size_t address = header_size + step * index;

  //check if selected graphnumber is valid
  unsigned int graphs_per_file = (file_size - header_size ) /step;
  if(graphs_per_file -1 < index){
    cerr << "You asked for the " << index+1 << "th fullerene, but there are only " << graphs_per_file << " stored in this file." << std::endl;
    abort();
  }

  //the actual parsing of the selected graph:
  //go to the beginning of the selected graph
  fseek(file, address+1, SEEK_SET);

  if(N<=255){
    for(node_t u=0; u<N; ++u){
      for(int neighbour=0; neighbour<3; ++neighbour){
	    unsigned char v;
        fread(reinterpret_cast<char*>(&v), 1, 1, file);
	    neighbours[u].push_back(v-1);
      }
      //skip one byte because there is no interesting content
      fseek(file, 1, SEEK_CUR);
    }
  } else{
    fseek(file, 2, SEEK_CUR);//because three bytes are not read
    for(node_t u=0; u<N; ++u){
      for(int neighbour=0; neighbour<3; ++neighbour){
	    unsigned short v;
        fread(reinterpret_cast<char*>(&v), 2, 1, file);
	    neighbours[u].push_back(v-1);
      }
      fseek(file, 2, SEEK_CUR);
    }
  }

  update_from_neighbours();
}

