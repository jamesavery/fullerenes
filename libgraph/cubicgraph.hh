#ifndef GRAPH3_HH
# define GRAPH3_HH

#include "planargraph.hh"
#include <iostream>
// TODO: Assumes planarity. Should perhaps split into cubic class and planar class?
struct CubicGraph : public PlanarGraph {

  CubicGraph(const PlanarGraph& g) : PlanarGraph(g) {
    for(node_t u=0;u<N;u++)
      if(neighbours[u].size() != 3){
	fprintf(stderr,"Graph not cubic: deg(%d) = %d\n",u,int(neighbours[u].size()));
	abort();
      }
  }
  CubicGraph(const Graph& g, const vector<coord2d>& layout) : PlanarGraph(g,layout) {}

  CubicGraph(FILE *file = stdin) {
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
	   //	   fprintf(stderr,"Skipped line %d\n",l++);
	 }
       }
    }
    N = neighbours.size();
    edges.resize(N*(N-1)/2,false);

    for(node_t n=0;n<N;n++){
      const vector<node_t>& ns(neighbours[n]);
      edge_t e1(n,ns[0]), e2(n,ns[1]), e3(n,ns[2]);
      edge_set.insert(e1);
      edge_set.insert(e2);
      edge_set.insert(e3);
      edges[e1.index()] = true;
      edges[e2.index()] = true;
      edges[e3.index()] = true;
    }
  }



  friend ostream& operator<<(ostream& s, const CubicGraph& g);
};

#endif
