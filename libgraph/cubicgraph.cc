#include "cubicgraph.hh"

#include <cstdlib>
#include <iostream>
#include <vector>
#include <fstream>


CubicGraph::CubicGraph(FILE *file = stdin){
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

CubicGraph::CubicGraph(unsigned int *index, FILE* file){
  const int header_size = 15;	
	
  //std::cout<< "entering constuctor" << std::endl;
  // Get file size
  //file->seekg(0,ios::end);
  fseek(file, 0, SEEK_END);
  //size_t file_size = file->tellg();
  size_t file_size = ftell(file);
  //std::cout<< "size has been read" << std::endl;

  //find number of vertices per graph
  //only for files with graphs of the equal size
  //file->seekg(header_size, ifstream::beg);
  fseek(file, header_size, SEEK_SET);
  // Read the number N of vertices per graph.
  //file->read(reinterpret_cast<char*>(&N), 1);
  //if(N == 0){
  //  file->read(reinterpret_cast<char*>(&N),2);
  //}
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
  size_t address = header_size + step * *index;

  //check if selected graph is valid
  unsigned int graphs_per_file = (file_size - header_size ) /step;
  //std::cout << graphs_per_file << std::endl;
  if(graphs_per_file -1 < *index){
    cerr << "You asked for the " << *index+1 << "th fullerene, but there are only " << graphs_per_file << " stored in this file.\n";
    abort();}

  //the actual parsing of the selected graph
  //file->seekg(address+1, ifstream::beg);//because the size is known
  fseek(file, address+1, SEEK_SET);

  if(N<=255){
    for(node_t u=0; u<N; ++u){
      for(int neighbour=0; neighbour<3; ++neighbour){
	unsigned char v;

	//fprintf(stderr, "Adding edge (%d,%d)\n",u,v);
	//file->read(reinterpret_cast<char*>(&v), 1);
    fread(reinterpret_cast<char*>(&v), 1, 1, file);
	neighbours[u].push_back(v-1);
      }
      //file->seekg(1,ifstream::cur);
      fseek(file, 1, SEEK_CUR);
    }
  } else{
    //file->seekg(2, std::ifstream::cur);//because three bytes are not read
    fseek(file, 2, SEEK_CUR);
    for(node_t u=0; u<N; ++u){
      for(int neighbour=0; neighbour<3; ++neighbour){
	unsigned short v;
	//file->read(reinterpret_cast<char*>(&v),2);
    fread(reinterpret_cast<char*>(&v), 2, 1, file);
	neighbours[u].push_back(v-1);
      }
      //file->seekg(2,ifstream::cur);
      fseek(file, 2, SEEK_CUR);
    }
  }

  update_from_neighbours();
}


