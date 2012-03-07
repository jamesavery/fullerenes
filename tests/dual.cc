#include "libgraph/cubicgraph.hh"

int main()
{
  CubicGraph g(stdin);
  g.layout2d = g.tutte_layout();
  Graph dual(g.dual_graph(6,g.layout2d));
  dual.name = "d";

  fprintf(stderr,"Computing faces.\n");
  map<unsigned int, set<face_t> > facemap(dual.compute_faces(3));

  set<int> degrees;
  for(unsigned int i=0;i<dual.N;i++)
    degrees.insert(dual.neighbours[i].size());

  fprintf(stderr,"%d vertices of degrees: ",dual.N);
  for(set<int>::const_iterator d(degrees.begin());d!=degrees.end();d++)
    fprintf(stderr,"%d ",*d);
  fprintf(stderr,"\n");

  fprintf(stderr,"Faces:\n");
  for(unsigned int i=2;i<10;i++){
    set<face_t> faces = facemap[i];
    if(!faces.empty()){
      fprintf(stderr,"%d %d-gons\n",int(faces.size()),i);
      if(i > 3)
	cerr << *faces.begin() << endl;
    }
  }
  
  cout << g << endl;
  cout << dual << endl;

  return 0;
}
