#include "libgraph/cubicgraph.hh"
#include "libgraph/polyhedron.hh"

int main()
{
  CubicGraph g(stdin);

  fprintf(stderr,"%d vertices and %d = %d edges\n",g.N,int(g.edge_set.size()),g.N*3/2);

  vector<node_t> cycle(g.shortest_cycle(0,g.neighbours[0][0],6));
  fprintf(stderr,"Outer face cycle: ");
  for(int i=0;i<cycle.size();i++)
    fprintf(stderr,"%d ",cycle[i]+1);
  fprintf(stderr,"\n");

  vector<unsigned int> vertex_depth(g.multiple_source_shortest_paths(cycle,vector<bool>(g.N*(g.N-1)/2),vector<bool>(g.N)));
  fprintf(stderr,"Depths:\n");
  for(int i=0;i<g.N;i++) fprintf(stderr,"\t%d : %d (%d,%d,%d)\n",i+1,vertex_depth[i],g.neighbours[i][0]+1,g.neighbours[i][1]+1,g.neighbours[i][2]+1);   

  g.layout2d = g.tutte_layout();
  map<unsigned int, set<face_t> > facemap(g.compute_faces_oriented(g.layout2d));

  fprintf(stderr,"Faces:\n");
  for(unsigned int i=2;i<10;i++){
    set<face_t> faces = facemap[i];
    if(!faces.empty()){
      fprintf(stderr,"%d %d-gons\n",int(faces.size()),i);
      if(i > 3)
  	cerr << *faces.begin() << endl;
    }
  }

  g.spherical_layout = g.spherical_projection(g.layout2d);

  cout << g << endl;


  pair< set<face_t>,set<face_t> > faces(g.compute_faces56());
  const set<face_t>& pentagons(faces.first), hexagons(faces.second);
  fprintf(stderr,"%d pentagons and %d hexagons.\n",int(pentagons.size()),int(hexagons.size()));

  set<face_t> all_faces(pentagons.begin(),pentagons.end());
  copy(hexagons.begin(),hexagons.end(),inserter(all_faces,all_faces.end()));
  printf("faces = {\n");
  for(set<face_t>::const_iterator f(all_faces.begin());f!=all_faces.end();){
    printf("{");
    for(unsigned int i=0;i<f->size();i++) printf("%d%s",(*f)[i]+1,i+1<f->size()?", ":"}");
    if(++f != all_faces.end())
      printf(",\n");
    else 
      printf("\n};\n");
  }

  return 0;
}
