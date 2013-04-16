#include <stdio.h>
#include <stdlib.h>
#include <fstream>
#include <vector>
#include "libgraph/cubicgraph.hh"
#include "libgraph/fullerenegraph.hh"
#include "libgraph/polyhedron.hh"

vector<coord3d> readxyz(const string filename)
{
  ifstream xyzfile(filename.c_str());
  string line,comment,dummy;
  int N, linenumber = 0,k=0;
  vector<coord3d> coordinates;

  if(xyzfile.fail()){
    cerr << "Couldn't open file " << filename << " for reading.\n";
    abort();
  }
  
  while(getline(xyzfile,line)){  
    stringstream l(line);
    ++linenumber;

    if(linenumber == 1) {
      l >> N;
      coordinates.resize(N);
    }
    if(linenumber == 2) comment = line;
    if(linenumber > 2) {
      l >> dummy;
      for(int i=0;i<3;i++)
	l >> coordinates[k].x[i];
      if(l.fail()){ cerr << "Malformed line: " << line << endl;  abort(); }      
      k++;
    }
    
  }
  xyzfile.close();
  return coordinates;
}

Graph connect_neighbours(double threshold, const vector<coord3d>& xs)
{
  set<edge_t> edges;
  for(node_t u=0;u<xs.size();u++)
    for(node_t v=u+1;v<xs.size();v++)
      if((xs[u]-xs[v]).norm() <= threshold)
	edges.insert(edge_t(u,v));

  return Graph(edges);
}

int main(int ac, char **av)
{
  if(ac<2){
    cerr << "Syntax: " << av[0] << " <filename.xyz>\n";
    return -1;
  }
  string filename(av[1]);
  string outname(filename.substr(0,filename.rfind("."))+".pov");
  cerr << "Reading from: " << filename << endl
       << "Writing to:   " << outname  << endl;

  vector<coord3d> coordinates(readxyz(filename));
  FullereneGraph G(connect_neighbours(1.6,coordinates));
  G.layout2d = G.tutte_layout();
  vector<face_t> faces(G.compute_faces_flat(6,true));
  Polyhedron P(G,coordinates,6,faces);

  double C20diameter = 4.08, diameter = P.diameter();
  double width = 2.5*diameter/C20diameter;

  cout << P.volume() << endl;
  ofstream outfile(outname.c_str());
  int red = min(int(200.0*(exp(4.0-diameter)/40.0)),0xd0);
  int blue = 0;//min(int(256.0*(1-exp((diameter-40.0)))),0x90);
  int green = min(int(256.0*(exp(diameter-40))),0xa0);

  int vertex_colour = (red<<16) | (green << 8) | blue;
  printf("vertex colour: %x%x%x = %x\n",red,green,blue,vertex_colour);

  outfile << P.to_povray(width,width,0x6a5acd,vertex_colour);
  outfile.close();

  return 0;
}


