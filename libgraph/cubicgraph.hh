#ifndef GRAPH3_HH
# define GRAPH3_HH

#include "graph.hh"
#include <iostream>
// TODO: Assumes planarity. Should perhaps split into cubic class and planar class?
struct CubicGraph : public Graph {

  vector<coord2d> spherical_layout;

  CubicGraph(const Graph& g) : Graph(g) {
    for(node_t u=0;u<N;u++)
      if(neighbours[u].size() != 3){
	fprintf(stderr,"Graph not cubic: deg(%d) = %d\n",u,int(neighbours[u].size()));
	abort();
      }
  }

  CubicGraph(const unsigned int N, const vector<node_t>& neighbours) : Graph(N) {
    assert(neighbours.size() == 3*N);
    for(int i=0;i<N;i++)
      for(int j=0;j<3;j++)
	edge_set.insert(edge_t(i,neighbours[3*i+j]));

    update_auxiliaries();
  }

  CubicGraph(FILE *file = stdin) {
    int l = 0;
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



  vector<coord2d> tutte_layout(const node_t s=0, const node_t t=0, const node_t r=0) const;
  vector<coord2d> spherical_projection(const vector< coord2d >& layout2d) const;

  friend ostream& operator<<(ostream& s, const CubicGraph& g);
};

#endif
