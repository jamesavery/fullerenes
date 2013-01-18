#include "graph.hh"

vector<node_t> Graph::shortest_path(const node_t& source, const node_t& dest, const vector<unsigned int>& dist) const
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

vector<unsigned int> Graph::shortest_paths(const node_t& source, const vector<bool>& used_edges, 
					   const vector<bool>& used_nodes, const unsigned int max_depth) const
{
  vector<unsigned int> distances(N,INT_MAX);
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

vector<unsigned int> Graph::shortest_paths(const node_t& source, const unsigned int max_depth) const
{
  vector<unsigned int> distances(N,INT_MAX);
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
vector<unsigned int> Graph::all_pairs_shortest_paths(const unsigned int max_depth) const
{
  vector<unsigned int> distances(N*N);
  vector<bool> dummy_edges(N*(N-1)/2), dummy_nodes(N);

  for(node_t u=0;u<N;u++){
    const vector<unsigned int> row(shortest_paths(u,dummy_edges,dummy_nodes,max_depth));
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
  vector<unsigned int> distances(shortest_paths(t,used_edges,used_nodes,max_depth));

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

vector<node_t> remove_node(const vector<node_t>& ns, const node_t& v)
{
  vector<node_t> r(ns.size()-1);
  for(int i=0,j=0;i<ns.size();i++)
    if(ns[i] != v) r[j++] = ns[i];
  return r;
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

  vector<unsigned int> distances(shortest_paths(r,used_edges,used_nodes,max_depth));

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


vector<unsigned int> Graph::multiple_source_shortest_paths(const vector<node_t>& sources, const vector<bool>& used_edges, 
							   const vector<bool>& used_nodes, const unsigned int max_depth) const
{
  vector<unsigned int> distances(N,INT_MAX);
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

vector<unsigned int> Graph::multiple_source_shortest_paths(const vector<node_t>& sources, const unsigned int max_depth) const
{
  return multiple_source_shortest_paths(sources,vector<bool>(N*(N-1)/2),vector<bool>(N),max_depth);
}





void Graph::update_auxiliaries() 
{
  // Instantiate auxiliary data strutures: sparse adjacency matrix and edge existence map.
  map<node_t,set<node_t> > ns;
  //  fprintf(stderr,"Initializing edge map.\n");

  // Update node count
  N = 0;
  for(set<edge_t>::const_iterator e(edge_set.begin()); e!= edge_set.end(); e++){
    N = max(N,max(e->first,e->second)+1);
  }

  edges.resize(N*(N-1)/2);
  neighbours.resize(N);

  for(set<edge_t>::const_iterator e(edge_set.begin()); e!= edge_set.end(); e++){
    edges[e->index()] = true;
    ns[e->first].insert(e->second);
    ns[e->second].insert(e->first);
  }

  //  fprintf(stderr,"Initializing adjacencies\n");
  for(int u=0;u<N;u++)
    neighbours[u] = vector<node_t>(ns[u].begin(),ns[u].end());

}

void Graph::update_from_neighbours() 
{
  edge_set.clear();
  for(node_t u=0;u<neighbours.size();u++)
    for(unsigned int i=0;i<neighbours[u].size();i++)
      edge_set.insert(edge_t(u,neighbours[u][i]));
  update_auxiliaries();
}


ostream& operator<<(ostream& s, const Graph& g) 
{
  s << g.name<< "Graph[Range["<<(g.N)<<"],\n\tUndirectedEdge@@#&/@{";
  for(set<edge_t>::const_iterator e(g.edge_set.begin()); e!=g.edge_set.end(); ){    
    s << "{" << (e->first+1) << "," << (e->second+1) << "}";
    if(++e != g.edge_set.end())
      s << ", ";
    else
      s << "}";
  } 
  s << "]";

  return s;
}

