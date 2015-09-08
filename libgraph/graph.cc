#include "graph.hh"

// Returns true if edge existed prior to call, false if not
bool Graph::remove_edge(const edge_t& e)
{
  node_t u = e.first, v = e.second;
  vector<node_t> &nu(neighbours[u]), &nv(neighbours[v]);

  bool value = false;

  for(int i=0;i<nu.size();i++) if(nu[i] == v){ nu.erase(nu.begin()+i); value = true; break; }
  for(int i=0;i<nv.size();i++) if(nv[i] == u){ nv.erase(nv.begin()+i); value = true; break; }

  return value;
}

// Returns true if edge existed prior to call, false if not
bool Graph::insert_edge(const dedge_t& e, const node_t suc_uv, const node_t suc_vu)
{
  if(edge_exists(e)) return true;	// insert_edge must be idempotent

  const node_t u = e.first, v = e.second;

  assert(u>=0 && v>=0);
  vector<node_t> &nu(neighbours[u]), &nv(neighbours[v]);

  size_t oldsize[2] = {nu.size(),nv.size()};
  
  vector<node_t>::iterator pos_uv = suc_uv<0? nu.end() : find(nu.begin(),nu.end(),suc_uv);
  vector<node_t>::iterator pos_vu = suc_vu<0? nv.end() : find(nv.begin(),nv.end(),suc_vu);

  nu.insert(pos_uv,v);
  if(u!=v) nv.insert(pos_vu,u);

  assert(nu.size() == oldsize[0]+1 && nv.size() == oldsize[1]+1);

  return false;
}

bool Graph::edge_exists(const edge_t& e) const
{
  const vector<node_t> &nu(neighbours[e.first]);
  return find(nu.begin(),nu.end(),e.second) != nu.end();
}

node_t Graph::next(const node_t& u, const node_t& v) const
{
  const vector<node_t>& nu(neighbours[u]);
  for(int j=0;j<nu.size(); j++) if(nu[j] == v) return nu[(j+1)%nu.size()];

  return -1;            // u-v is not an edge in a triangulation
}

node_t Graph::prev(const node_t& u, const node_t& v) const
{
  const vector<node_t>& nu(neighbours[u]);
  for(int j=0;j<nu.size(); j++) if(nu[j] == v) return nu[(j-1+nu.size())%nu.size()];

  return -1;            // u-v is not an edge in a triangulation
}


bool Graph::is_consistently_oriented() const 
{
  map<dedge_t,bool> seen_dedge;

  set<dedge_t> work;
  for(node_t u=0;u<N;u++)
    for(auto v: neighbours[u])
      work.insert({u,v});

  while(!work.empty()){
    const dedge_t e = *work.begin();
    node_t u(e.first), v(e.second);

    // Process CW face starting in u->v
    const node_t u0 = u;
    work.erase(dedge_t(u,v));
    while(v != u0){
      node_t w = next(u,v);	// u--v--w is CW-most corner
      u = v;
      v = w;
      if(work.find(dedge_t(u,v)) == work.end()){ // We have already processed arc u->v
	//	cerr << "Directed edge " << dedge_t(u,v) << " is part of two faces.\n";
	return false;
      }

      work.erase(dedge_t(u,v));
    }
  }
  // Every directed edge is part of exactly one face <-> orientation is consistent
  return true;
}

bool Graph::adjacency_is_symmetric() const
{
  for(node_t u=0;u<N;u++){
    const vector<node_t> &nu = neighbours[u];
    for(int i=0;i<nu.size();i++){
      const vector<node_t> &nv = neighbours[nu[i]];

      bool symmetric = false;
      for(int j=0;j<nv.size();j++) if(nv[j] == u) symmetric = true;
      if(!symmetric) return false;
    }
  }
  return true;
}

 // TODO: Should make two functions: one that takes subgraph (empty is trivially connected) and one that works on full graph.
bool Graph::is_connected(const set<node_t> &subgraph) const 
{
  if(!subgraph.empty()){
    node_t s = *subgraph.begin();
    const vector<int> dist(shortest_paths(s));

    for(set<node_t>::const_iterator u(subgraph.begin()); u!=subgraph.end();u++) 
      if(dist[*u] == INT_MAX) return false;

  } else {
    node_t s = 0; // Pick a node that is part of an edge
    for(;neighbours[s].empty();s++) ;
    assert(s < N);

    const vector<int> dist(shortest_paths(s));

    for(int i=0;i<dist.size();i++) 
      if(dist[i] == INT_MAX) return false;
  }

  return true;
}

list< list<node_t> > Graph::connected_components() const 
{
  vector<bool> done(N);

  list<list<node_t> > components;
  for(node_t u=0;u<N;u++)
    if(!done[u]){
      list<node_t> component;
  
      done[u] = true;
      component.push_back(u);
      list<node_t> q(neighbours[u].begin(),neighbours[u].end());
      while(!q.empty()){
	node_t v = q.back(); q.pop_back(); done[v] = true;
	component.push_back(v);
	for(int i=0;i<neighbours[v].size();i++) if(!done[neighbours[v][i]]) q.push_back(neighbours[v][i]);
      }
      components.push_back(component);
    }
  return components;
}

vector<int> Graph::shortest_path(const node_t& source, const node_t& dest, const vector<int>& dist) const
{
  // Fill in shortest paths -- move to own function.
  node_t vi = source;
  vector<node_t> path;
  //  fprintf(stderr,"Shortest path %d -> %d\n",source,dest);
  for(int i=0;i<N;i++)
    //    fprintf(stderr,"%d: %d\n",i,dist[i]);

  do {
    path.push_back(vi);
    const vector<node_t>& ns(neighbours[vi]);
    // Find next vertex in shortest path
    int kmin = 0;
    for(int k=1, d=dist[ns[0]];k<ns.size();k++) {
      //      fprintf(stderr,"node %d has distance %d\n",ns[k],dist[ns[k]]);
      if(dist[ns[k]] < d){ 
	d = dist[ns[k]];
	kmin = k;
      }
    }
    //    fprintf(stderr,"Choosing neighbour %d = %d\n",kmin,ns[kmin]);
    vi = ns[kmin];
  } while(vi != dest);
  //  cerr << face_t(path) << endl;
  return path;
}

vector<int> Graph::shortest_paths(const node_t& source, const vector<bool>& used_edges, 
					   const vector<bool>& used_nodes, const unsigned int max_depth) const
{
  vector<int> distances(N,INT_MAX);
  list<node_t> queue;    

  distances[source] = 0;
  queue.push_back(source);

  while(!queue.empty()){
    node_t v = queue.front(); queue.pop_front();
      
    const vector<node_t> &ns(neighbours[v]);
    for(unsigned int i=0;i<ns.size();i++){
      const edge_t edge(v,ns[i]);
      if(!used_nodes[ns[i]] && !used_edges[edge.index()] && distances[ns[i]] == INT_MAX){
	distances[ns[i]] = distances[v] + 1;
	if(distances[ns[i]] < max_depth) queue.push_back(ns[i]);
      }
    }
  }
  return distances;
}

vector<int> Graph::shortest_paths(const node_t& source, const unsigned int max_depth) const
{
  vector<int> distances(N,INT_MAX);
  list<node_t> queue;    

  distances[source] = 0;
  queue.push_back(source);

  while(!queue.empty()){
    node_t v = queue.front(); queue.pop_front();
      
    const vector<node_t> &ns(neighbours[v]);
    for(unsigned int i=0;i<ns.size();i++){
      const edge_t edge(v,ns[i]);
      if(distances[ns[i]] == INT_MAX){
	distances[ns[i]] = distances[v] + 1;
	if(distances[ns[i]] < max_depth) queue.push_back(ns[i]);
      }
    }
  }
  return distances;
}

// Returns NxN matrix of shortest distances (or INT_MAX if not connected)
vector<int> Graph::all_pairs_shortest_paths(const unsigned int max_depth) const
{
  vector<int> distances(N*N);
  vector<bool> dummy_edges(N*(N-1)/2), dummy_nodes(N);

  for(node_t u=0;u<N;u++){
    const vector<int> row(shortest_paths(u,dummy_edges,dummy_nodes,max_depth));
    memcpy(&distances[u*N],&row[0],N*sizeof(unsigned int));
  }

  return distances;
}

vector<node_t> Graph::shortest_cycle(const node_t& s, const int max_depth) const 
{
  const vector<node_t> ns(neighbours[s]);
  assert(ns.size() > 0);

  face_t cycle;
  int Lmax = 0;
  for(int i=0;i<ns.size();i++){
    face_t c(shortest_cycle(s,ns[i],max_depth));
    if(c.size() >= Lmax){ Lmax = c.size(); cycle = c; }
  }
  return cycle;
}

vector<node_t> Graph::shortest_cycle(const node_t& s, const node_t& t, const int max_depth) const 
{ 
  vector<bool> used_edges(N*(N-1)/2);
  vector<bool> used_nodes(N);

  // Find shortest path t->s in G excluding edge (s,t)
  used_edges[edge_t(s,t).index()] = true;
  vector<int> distances(shortest_paths(t,used_edges,used_nodes,max_depth));

  // If distances[s] is uninitialized, we have failed to reach s and there is no cycle <= max_depth.
  if(distances[s] == INT_MAX) return vector<node_t>();

  // Then reconstruct the cycle by moving backwards from s to t
  vector<node_t> cycle(distances[s]+1);
  node_t u = s, v = -1;
  for(unsigned int i=0;i<distances[s]+1;i++){
    const vector<node_t> &ns(neighbours[u]);
    unsigned int dmin = INT_MAX;
    for(unsigned int j=0;j<ns.size();j++)
      if(distances[ns[j]] < dmin && (edge_t(u,ns[j]) != edge_t(s,t))){
	dmin = distances[ns[j]];
	v = ns[j];
      } 
    u = v;
    cycle[distances[s]-i] = u;
  }
  cycle[0] = s;
  return cycle;
}

vector<node_t> Graph::shortest_cycle(const node_t& s, const node_t& t, const node_t& r, const int max_depth) const 
{ 
  //  fprintf(stderr,"3: shortest_cycle(%d,%d,%d,max_depth=%d)\n",s,t,r,max_depth);
  assert(s >= 0 && t >= 0 && r >= 0);
  if(max_depth == 3){		// Triangles need special handling
    vector<node_t> cycle(3);
    node_t p=-1,q=-1;
    cycle[0] = s;
    for(int i=0;i<neighbours[s].size();i++) {
      if(neighbours[s][i] == t){ p = t; q = r; break; }
      if(neighbours[s][i] == r){ p = r; q = t; break; }
    }
    bool foundq = false, founds = false;
    if(p<0) return vector<node_t>();
    for(int i=0;i<neighbours[p].size();i++) if(neighbours[p][i] == q) foundq = true;
    for(int i=0;i<neighbours[q].size();i++) if(neighbours[q][i] == s) founds = true;
    
    if(!foundq || !founds) return vector<node_t>();
    cycle[1] = p;
    cycle[2] = q;
    return cycle;
  }

  vector<bool> used_edges(N*(N-1)/2);
  vector<bool> used_nodes(N);
  //  fprintf(stderr,"shortest_cycle(%d,%d,%d), depth = %d\n",s,t,r,max_depth);
  // Find shortest path r->s in G excluding edge (s,t), (t,r)
  used_edges[edge_t(s,t).index()] = true;
  used_edges[edge_t(t,r).index()] = true;

  vector<int> distances(shortest_paths(r,used_edges,used_nodes,max_depth));

  // If distances[s] is uninitialized, we have failed to reach s and there is no cycle <= max_depth.
  if(distances[s] == INT_MAX) return vector<node_t>();

  // Else reconstruct the cycle by moving backwards from s to t
  vector<node_t> cycle(distances[s]+2);
  node_t u = s, v = -1;
  for(unsigned int i=0;i<distances[s]+1;i++){
    const vector<node_t> &ns(neighbours[u]);
    unsigned int dmin = INT_MAX;
    for(unsigned int j=0;j<ns.size();j++)
      if(distances[ns[j]] < dmin && (edge_t(u,ns[j]) != edge_t(s,t))){
	dmin = distances[ns[j]];
	v = ns[j];
      } 
    u = v;
    cycle[distances[s]+1-i] = u;
  }
  cycle[0] = s;
  cycle[1] = t;
  return cycle;
}


vector<int> Graph::multiple_source_shortest_paths(const vector<node_t>& sources, const vector<bool>& used_edges, 
							   const vector<bool>& used_nodes, const unsigned int max_depth) const
{
  vector<int> distances(N,INT_MAX);
  list<node_t> queue;
    
  for(unsigned int i=0;i<sources.size();i++){
    distances[sources[i]] = 0;
    queue.push_back(sources[i]);
  }

  while(!queue.empty()){
    node_t v = queue.front(); queue.pop_front();
      
    const vector<node_t> &ns(neighbours[v]);
    for(unsigned int i=0;i<ns.size();i++){
      const edge_t edge(v,ns[i]);
      if(!used_nodes[ns[i]] && !used_edges[edge.index()] && distances[ns[i]] == INT_MAX){
	distances[ns[i]] = distances[v] + 1;
	if(distances[ns[i]] < max_depth) queue.push_back(ns[i]);
      }
    }
  }
  return distances;
}

vector<int> Graph::multiple_source_shortest_paths(const vector<node_t>& sources, const unsigned int max_depth) const
{
  return multiple_source_shortest_paths(sources,vector<bool>(N*(N-1)/2),vector<bool>(N),max_depth);
}


int Graph::max_degree() const
{
  int max_d = 0;
  for(node_t u=0;u<N;u++) if(neighbours[u].size() > max_d) max_d = neighbours[u].size();
  return max_d;
}

int Graph::degree(const node_t& u) const { 
  return neighbours[u].size(); 
}

void Graph::update_from_edgeset(const set<edge_t>& edge_set) 
{
  // Instantiate auxiliary data strutures: sparse adjacency matrix and edge existence map.
  map<node_t,set<node_t> > ns;
  //  fprintf(stderr,"Initializing edge map.\n");

  // Update node count
  N = 0;
  for(set<edge_t>::const_iterator e(edge_set.begin()); e!= edge_set.end(); e++){
    N = max(N,max(e->first,e->second)+1);
  }

  neighbours.resize(N);

  for(set<edge_t>::const_iterator e(edge_set.begin()); e!= edge_set.end(); e++){
    ns[e->first].insert(e->second);
    ns[e->second].insert(e->first);
  }

  //  fprintf(stderr,"Initializing adjacencies\n");
  for(int u=0;u<N;u++)
    neighbours[u] = vector<node_t>(ns[u].begin(),ns[u].end());

}

set<edge_t> Graph::undirected_edges() const {
  set<edge_t> edges;
  for(node_t u=0;u<N;u++)
    for(int i=0;i<neighbours[u].size();i++)
      edges.insert(edge_t(u,neighbours[u][i]));
  return edges;
}

set<dedge_t> Graph::directed_edges() const {
  set<dedge_t> edges;
  for(node_t u=0;u<N;u++)
    for(int i=0;i<neighbours[u].size();i++)
      edges.insert(dedge_t(u,neighbours[u][i]));
  return edges;
}


ostream& operator<<(ostream& s, const Graph& g) 
{
  set<edge_t> edge_set = g.undirected_edges();

  s << "Graph[Range["<<(g.N)<<"],\n\tUndirectedEdge@@#&/@{";
  for(set<edge_t>::const_iterator e(edge_set.begin()); e!=edge_set.end(); ){    
    s << "{" << (e->first+1) << "," << (e->second+1) << "}";
    if(++e != edge_set.end())
      s << ", ";
  } 
  s << "}]";

  return s;
}

