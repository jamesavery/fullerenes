#include <stdio.h>
#include <stdlib.h>
#include <fstream>
#include <vector>
#include "libgraph/cubicgraph.hh"
#include "libgraph/fullerenegraph.hh"
#include "libgraph/polyhedron.hh"

int main(int ac, char **av)
{
  if(ac<2){
    cerr << "Syntax: " << av[0] << " <filename.mol2>\n";
    return -1;
  }

  for(int i=1;i<ac;i++){
    string filename(av[i]);
    string outname(filename.substr(0,filename.rfind("."))+".pov");
    //  cerr << "Reading from: " << filename << endl
    //       << "Writing to:   " << outname  << endl;

    Polyhedron P(filename);

    double min_diameter = 4, max_diameter = 30, diameter = P.diameter();
    double width = 2.5*diameter/min_diameter;

    //  printf("{%d, %f, %f}\n",G.N,P.volume(),P.diameter());

    ofstream outfile(outname.c_str());
    int red   = int(256.0*min(diameter/max_diameter-0.1,0.8));
    int green = int(256.0*min(min_diameter/diameter-0.1,0.85));
    int blue = 0;

    int vertex_colour = (red<<16) | (green << 8) | blue;

    outfile << P.to_povray(width,width,0x6a5acd,vertex_colour);
    outfile.close();
  }

  return 0;
}


