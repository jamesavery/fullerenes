#include "fullerenes/graph.hh"

char LIST_OPEN='[';
char LIST_CLOSE=']';

// Returns true if edge existed prior to call, false if not
bool Graph::remove_edge(const edge_t& e)
{
  auto [u,v] = e;		
  bool success_uv = neighbours.remove_arc(u,v);
  bool success_vu = neighbours.remove_arc(v,u);

  if(success_uv ^ success_vu){
    cerr << "Corrupt graph: only half of edge ("<<u<<","<<v<<") could be removed.\n";
    abort();
  }
  
  return success_uv && success_vu;
}

// Returns true if edge existed prior to call, false if not
// insert v right before suc_uv in the list of neighbours of u
// insert u right before suc_vu in the list of neighbours of v
// NB: arc_t type is intended, as order matters due to successors.
bool Graph::insert_edge(const arc_t& e, const node_t suc_uv, const node_t suc_vu)
{
  auto [u,v] = e;		
  bool success_uv = neighbours.insert_arc(u,v, suc_uv);
  bool success_vu = neighbours.insert_arc(v,u, suc_vu);

  if(success_uv ^ success_vu){
    cerr << "Corrupt graph: only half of edge ("<<u<<","<<v<<") could be inserted.\n";
    abort();
  }
  
  return success_uv && success_vu;
}

bool Graph::edge_exists(const edge_t& e) const
{
  auto [u,v] = e;
  return neighbours.arc_exists(u,v) && neighbours.arc_exists(v,u);
}

// // remove all vertices without edges from graph
// void Graph::remove_isolated_vertices(){
//   int new_id[N];

//   int u_new = 0;
//   for(int u=0; u<N; u++)
//     if(!neighbours[u].empty())
//       new_id[u] = u_new++;

//   int N_new = u_new;
//   Graph g(N_new);
//   // cerr << "n new: " << N_new << endl;
//   for(int u=0; u<N; u++)
//     for(int v: neighbours[u])
//       g.neighbours[new_id[u]].push_back(new_id[v]);

//   *this = g;
// }

// // completely remove all vertices in sv from the graph
// void Graph::remove_vertices(set<int> &sv){
//   const int N_naught(N);
//   for(int u: sv){
//     while(neighbours[u].size()){
//       const int v = neighbours[u][0];
//       remove_edge({u,v});
//     }
//   }

//   remove_isolated_vertices();

//   // let's see if the graph remained in a sane state
//   // cerr << "N: " << N << endl;
//   if(N_naught != sv.size() + N)
//     cerr << "removed more vertices than intended" << endl;
//   assert(is_connected());
// }

int  Graph::arc_ix(node_t u, node_t v) const
{
  return neighbours.arc_index(u,v); // TODO: Graph should just inherit from View::sparsity<node_t>.
}

// Successor to v in oriented neigbhours of u

node_t Graph::next(node_t u, node_t v) const {
  const auto& ru = neighbours[u];
  int  j   = arc_ix(u,v);
  if(j>=0) return ru[(j+1)%ru.size()];
  else return -1;		// u-v is not an edge in this graph
}

// Predecessor to v in oriented neigbhours of u
node_t Graph::prev(node_t u, node_t v) const {
  const auto& ru = neighbours[u];
  int  j   = arc_ix(u,v);
  if(j>=0) return ru[(j+ru.size()-1)%ru.size()];
  else return -1;		// u-v is not an edge in this graph
}  

// Successor to v in face containing directed edge u->v
node_t Graph::next_on_face(node_t u, node_t v) const
{
  return prev(v,u);
}

// Predecessor to v in face containing directed edge u->v
node_t Graph::prev_on_face(node_t u, node_t v) const
{
  return next(v,u);
}

bool Graph::is_consistently_oriented() const 
{
  map<arc_t,bool> seen_arc;

  set<arc_t> work;
  for(node_t u=0;u<N;u++)
    for(auto v: neighbours[u])
      work.insert({u,v});

  while(!work.empty()){
    const arc_t e = *work.begin();
    node_t u(e.first), v(e.second);

    // Process CW face starting in u->v
    const node_t u0 = u;
    work.erase(arc_t(u,v));
    while(v != u0){
      node_t w = next(v,u);	// u--v--w is CW-most corner
      u = v;
      v = w;
      if(work.find(arc_t(u,v)) == work.end()){ // We have already processed arc u->v
	//	cerr << "Directed edge " << arc_t(u,v) << " is part of two faces.\n";
	return false;
      }

      work.erase(arc_t(u,v));
    }
  }
  // Every directed edge is part of exactly one face <-> orientation is consistent
  return true;
}

// TODO: Doesn't need to be planar and oriented, but is easier to write if it is. Make it work in general.
bool Graph::has_separating_triangles() const
{
  assert(is_oriented);

  for(node_t u=0;u<N;u++){
    const auto &nu = neighbours[u];
    
    for(int i=0;i<nu.size();i++){
      node_t t = nu[i];
      node_t v = prev(u,t), w = next(u,t); // edges: u--t, u--v, u--w
      if(edge_exists({t,w}) && edge_exists({t,v}) && edge_exists({v,w})) return true;
    }
  }
  return false;
}


bool Graph::adjacency_is_symmetric() const
{
  for(node_t u=0;u<N;u++){
    const auto &nu = neighbours[u];
    for(int i=0;i<nu.size();i++){
      const auto &nv = neighbours[nu[i]];

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
  vector<int> dist(N);
  if(!subgraph.empty()){
    node_t s = *subgraph.begin();
    single_source_shortest_paths(s,&dist[0]);

    for(node_t u: subgraph)
      if(dist[u] == INT_MAX) return false;

  } else {
    node_t s = 0; // Pick a node that is part of an edge
    for(;neighbours[s].empty();s++) ;
    assert(s < N);

    single_source_shortest_paths(s,&dist[0]);

    for(int i=0;i<dist.size();i++) 
      if(dist[i] == INT_MAX) return false;
  }

  return true;
}

#include <queue>
vector<vector<node_t>> Graph::connected_components() const 
{
  vector<bool> done(N);

  vector<vector<node_t>> components;
  for(node_t u=0;u<N;u++)
    if(!done[u]){
      vector<node_t> component;
  
      done[u] = true;
      component.push_back(u);
      queue<node_t> Q;
      for(auto v: neighbours[u]) Q.push(v);

      while(!Q.empty()){
	node_t v = Q.front(); Q.pop();
	if(!done[v]){
	  done[v] = true;
	  component.push_back(v);
	  
	  for(int i=0;i<neighbours[v].size();i++)
	    if(!done[neighbours[v][i]]) Q.push(neighbours[v][i]);
	}
      }
      sort(component.begin(), component.end());
      components.push_back(component);
    }
  return components;
}


void Graph::single_source_shortest_paths(node_t source, int *distances, size_t max_depth) const
{
  Deque<node_t> queue(N);
  for(int u=0;u<N;u++) distances[u] = INT_MAX;
  
  distances[source] = 0;
  queue.push_back(source);

  while(!queue.empty()){
    node_t u = queue.pop_front();
    
    for(node_t v: neighbours[u]){
      if(distances[v] == INT_MAX){ // Node is not previously visited
	distances[v] = distances[u] + 1; 
	if(distances[v] < max_depth) queue.push_back(v);
      }
    }
  }
}

// Returns NxN matrix of shortest distances (or INT_MAX if not connected)
// N^2: allocating d
// N*(N-1)/2 steps
matrix<int> Graph::all_pairs_shortest_paths(const unsigned int max_depth) const
{
  matrix<int>   d(N,N,INT_MAX);
  Deque<node_t> queue(N);

  for(node_t u=0;u<N;u++){
    queue.push_back(u);    	// Enqueue source u
    d(u,u) = 0;

    while(!queue.empty()){
      node_t v = queue.pop_front();

      // Process children w of node v
      for(node_t w: neighbours[v]){
	if(d(u,w) == INT_MAX){ // Node is not previously visited
	  int distance = d(u,v)+1;
	  d(u,w)  = distance;
	  //	  d(w,u)  = distance;
	  if(distance < max_depth) queue.push_back(w);
	}
      }
    }
  }

  return d;
}

// Returns MxM matrix of shortest distances between vertices in V
// M^2 memory, O(MN) operations
matrix<int> Graph::all_pairs_shortest_paths(const vector<node_t> &V,
					    const unsigned int max_depth) const
{
  size_t M = V.size();
  matrix<int> D(M,M,INT_MAX);
  vector<int> d(N);

  Deque<node_t> queue(N);

  for(int i=0;i<M;i++){
    for(int j=0;j<N;j++) d[j] = INT_MAX; // Mark all nodes as unvisited    

    node_t source = V[i];
    d[source] = 0;
    queue.push_back(source);

    while(!queue.empty()){
      node_t u = queue.pop_front();

      // Process children of node u
      for(node_t v: neighbours[u]){
	if(d[v] == INT_MAX){ 	// Node is not previously visited
	  d[v] = d[u]+1;
	  if(d[v] < max_depth) queue.push_back(v);
	}
      }
    }
    // Queue is empty, now pick out nodes from V
    for(int j=0;j<M;j++)
      D(i,j) = d[V[j]];
  }

  return D;
}


vector<node_t> Graph::shortest_cycle(node_t s, const int max_depth) const 
{
  face_t cycle;
  int Lmin = INT_MAX;
  for(node_t t: neighbours[s]){
    face_t c(shortest_cycle({s,t},max_depth));
    if(c.size() < Lmin){ Lmin = c.size(); cycle = c; }
  }
  return cycle;
}

// Find shortest cycle s->t->r->...->s, prefix = {s,t,r,...}
// This is the same as: Find shortest path t-*->s in G excluding edges {(s,t),(t,r),...}
// This is the same as:
//   1. Compute d(t,*) when excluding the edge (s,t)
//   2. If d(t,s) == INT_MAX, there is no such cycle <= max_depth
//   3. Otherwise, back-trace as
vector<node_t> Graph::shortest_cycle(const vector<node_t>& prefix, const int max_depth) const 
{
  // Is this a valid start?
  for(int i=0;i+1<prefix.size();i++) assert(edge_exists({prefix[i],prefix[i+1]}));

  node_t s = prefix[0], t = prefix[1];  
  
  if(max_depth == 3){ // Triangles need special handling
    switch(prefix.size()){
    case 2:
      // t must have a neighbor r that neighbours s
      for(node_t r: neighbours[t])
	if(edge_exists({r,s})) return {s,t,r};
      return {};
    case 3:
      return prefix;		// {s,t,r} given and is a valid triangle
    default:
      cerr << "shortest_cycle(): Prefix " << prefix << " is not an appropriate triangle start.\n";
      abort();
    }
  }

  // Now we can assume max_depth >= 4
  vector<int> distances(N);
  Graph G(*this);		// TODO: Get rid of this copy
  for(int i=0;i+1<prefix.size();i++)
    G.remove_edge({prefix[i],prefix[i+1]});

  G.single_source_shortest_paths(t,&distances[0],max_depth);
  // If distances[s] is uninitialized, we have failed to reach s and there is no cycle <= max_depth.
  if(distances[s] == INT_MAX) return {};

  // Otherwise reconstruct the cycle by moving backwards from s to t
  vector<node_t> cycle(distances[s]+1);
  node_t u = s, v = -1;
  for(unsigned int i=0;i<distances[s]+1;i++){
    unsigned int dmin = INT_MAX;
    for(node_t w: neighbours[u])
      if(distances[w] < dmin && (edge_t(u,w) != edge_t(s,t))){
  	dmin = distances[w];
  	v = w;
      } 
    u = v;
    cycle[distances[s]-i] = u;
  }
  cycle[0] = s;
  return cycle;
}


vector<int> Graph::multiple_source_shortest_paths(const vector<node_t>& sources, const unsigned int max_depth) const
{
  vector<int>   distances(N,INT_MAX);
  Deque<node_t> queue(N);
    
  for(node_t s: sources){
    distances[s] = 0;
    queue.push_back(s);
  }

  while(!queue.empty()){
    node_t v = queue.pop_front();
      
    for(node_t w: neighbours[v]){
      const edge_t edge(v,w);
      if(distances[w] == INT_MAX){ // Node not previously visited
	distances[w] = distances[v] + 1;
	if(distances[w] < max_depth) queue.push_back(w);
      }
    }
  }
  return distances;
}


int Graph::max_degree() const
{
  int max_d = 0;
  for(node_t u=0;u<N;u++) if(neighbours[u].size() > max_d) max_d = neighbours[u].size();
  return max_d;
}

int Graph::degree(node_t u) const { 
  return neighbours[u].size(); 
}

// void Graph::update_from_edgeset(const set<edge_t>& edge_set) 
// {
//   // Instantiate auxiliary data strutures: sparse adjacency matrix and edge existence map.
//   map<node_t,set<node_t> > ns;
//   //  fprintf(stderr,"Initializing edge map.\n");

//   // Update node count
//   N = 0;
//   for(set<edge_t>::const_iterator e(edge_set.begin()); e!= edge_set.end(); e++){
//     N = max(N,max(e->first,e->second)+1);
//   }

//   neighbours.resize(N);

//   for(set<edge_t>::const_iterator e(edge_set.begin()); e!= edge_set.end(); e++){
//     ns[e->first].insert(e->second);
//     ns[e->second].insert(e->first);
//   }

//   //  fprintf(stderr,"Initializing adjacencies\n");
//   for(int u=0;u<N;u++)
//     neighbours[u] = vector<node_t>(ns[u].begin(),ns[u].end());

// }

vector<edge_t> Graph::undirected_edges() const {
  set<edge_t> edges;
  for(node_t u=0;u<N;u++)
    for(int i=0;i<neighbours[u].size();i++)
      edges.insert(edge_t(u,neighbours[u][i]));
  return vector<edge_t>(edges.begin(),edges.end());
}

vector<arc_t> Graph::directed_edges() const {
  set<arc_t> edges;
  for(node_t u=0;u<N;u++)
    for(int i=0;i<neighbours[u].size();i++)
      edges.insert(arc_t(u,neighbours[u][i]));
  return vector<arc_t>(edges.begin(),edges.end());
}

size_t Graph::count_edges() const {
  // Don't use edge_set -- it's slow
  size_t twoE = 0;
  for(node_t u=0;u<N;u++)
    twoE += neighbours[u].size();

  return twoE/2;
}

ostream& operator<<(ostream& s, const Graph& g) 
{
  vector<edge_t> edges = g.undirected_edges();

  s << "Graph[Range["<<(g.N)<<"],\n\tUndirectedEdge@@#&/@{";
  for(size_t i=0;i<edges.size();i++){    
    s << "{" << (edges[i].first+1) << "," << (edges[i].second+1) << "}" << (i+1<edges.size()? ", ":"");
  } 
  s << "}]";

  return s;
}

