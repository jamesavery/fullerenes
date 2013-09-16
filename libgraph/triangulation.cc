#include "triangulation.hh"

pair<node_t,node_t> Triangulation::adjacent_tris(const edge_t& e) const
{
  const node_t &u(e.first), &v(e.second);
  const vector<node_t>& nv(neighbours[v]);

  pair<node_t,node_t> tris;

  for(int i=0, t=0;i<nv.size();i++){
    const node_t& w(nv[i]);
    const vector<node_t>& nw(neighbours[w]);
    for(int j=0;j<nw.size();j++)
      if(nw[j] == u) { 
	if(++t == 1)  tris.first  = w;
	else if(t==2) tris.second = w;
      } else {
	  fprintf(stderr,"Triangulation is not orientable, edge %d--%d part of more than two faces.\n",u,v);
	  abort();
	}
  }
  return tris;
}

#define place_dedge(u,v) {   \
  dedge_t uv(u,v), vu(v,u);  \
  dedge_done[uv] = true;     \
  workset.erase(uv);         \
  if(dedge_done[vu] != true) \
    workset.insert(vu);      \
  }

// Computes the faces of a triangulation, returning the faces in a consistent ordering
vector<tri_t> Triangulation::compute_faces() const
{
  vector<tri_t> faces;
  set<dedge_t> workset;
  map<dedge_t,bool> dedge_done;
  node_t u(0), v(neighbours[0][0]), w;

  pair<node_t,node_t> tris(adjacent_tris(edge_t(u,v)));
  w = tris.first;

  // Place the first triangle, its orientation defines the global orientation of the graph
  place_dedge(u,v);
  place_dedge(v,w);
  place_dedge(w,u);

  while(!workset.empty()){
    dedge_t uv(*workset.begin());
    node_t u(uv.first), v(uv.second), w;

    pair<node_t,node_t> tris(uv);
    if(dedge_done[dedge_t(v,tris.first)])      // Already placed triangle defined by v->w0
      if(dedge_done[dedge_t(v,tris.second)]){  // Already placed triangle defined by v->w1
	fprintf(stderr,"Triangulation is not planar: Edge %d-%d is part of more than two faces.\n",u,v);
	abort();
      } else w = tris.second;
    else w = tris.first;

    place_dedge(u,v);
    place_dedge(v,w);
    place_dedge(w,u);

    faces.push_back(tri_t(u,v,w));
  }
  return faces;
}

vector<tri_t> Triangulation::compute_faces_oriented() const 
{
  vector<tri_t> faces;
  map<dedge_t,bool> dedge_done;

  for(node_t u=0;u<N;u++){
    const vector<node_t>& nu(neighbours[u]);
    for(int i=0;i<nu.size();i++){
      const node_t &v(nu[i]);
      dedge_t uv(u,v);
      if(!dedge_done[uv]){
	node_t w = nextCW.find(uv)->second;
	faces.push_back(tri_t(u,v,w));
	dedge_done[uv] = true;
	dedge_done[dedge_t(v,w)] = true;
	dedge_done[dedge_t(w,u)] = true;
      }
    }
  }
  return faces;
}

void Triangulation::orient_neighbours()
{
  vector<tri_t> faces(compute_faces());

  for(int i=0;i<faces.size();i++){
    const tri_t& t(faces[i]);
    nextCW[dedge_t(t[0],t[1])] = t[2];
    nextCW[dedge_t(t[1],t[2])] = t[0];
    nextCW[dedge_t(t[2],t[0])] = t[1];
  }

  for(node_t u=0;u<N;u++){
    int d = neighbours[u].size();

    node_t v = neighbours[u][0]; 
    for(int i=1;i<d;i++){
      node_t w = nextCW[dedge_t(u,v)];
      neighbours[u][i] = w;
      v = w;
    }
  }
}

// Takes full spiral string, e.g. 5666436634665 
Triangulation::Triangulation(const vector<int>& spiral_string)
{
  N = spiral_string.size();
  // Lukas: Can you implement a fast spiral windup?
  
}


Triangulation::Triangulation(const vector<int>& spiral_string, const jumplist_t jumps)
{
  N = spiral_string.size();
  // Lukas: Can you implement a fast general spiral windup?
}
