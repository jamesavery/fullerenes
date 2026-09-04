#include "fullerenes/graph.hh"
#include <stdexcept>

char LIST_OPEN='[';
char LIST_CLOSE=']';

// Returns true if edge existed prior to call, false if not
bool GraphView::remove_edge(const edge_t& e)
{
  node_t u = e.first, v = e.second;
  bool value = false;

  int pos_uv = find(u, v);
  if(pos_uv >= 0){ erase_at(u, pos_uv); value = true; }
  int pos_vu = find(v, u);
  if(pos_vu >= 0){ erase_at(v, pos_vu); value = true; }

  return value;
}

// Returns true if edge existed prior to call, false if not
// insert v right before suc_uv in the list of neighbours of u
// insert u right before suc_vu in the list of neighbours of v
bool GraphView::insert_edge(const arc_t& e, const node_t suc_uv, const node_t suc_vu)
{
  if(edge_exists(e)) return true;	// insert_edge must be idempotent

  const node_t u = e.first, v = e.second;

  assert(u>=0 && v>=0);
  int oldsize_u = degree(u), oldsize_v = degree(v);

  if(suc_uv < 0) push_back(u, v);
  else {
    int pos = find(u, suc_uv);
    insert_at(u, v, pos >= 0 ? pos : degree(u));
  }

  if(u != v){
    if(suc_vu < 0) push_back(v, u);
    else {
      int pos = find(v, suc_vu);
      insert_at(v, u, pos >= 0 ? pos : degree(v));
    }
  }

  assert(degree(u) == oldsize_u+1 && degree(v) == oldsize_v+1);

  return false;
}

bool GraphView::edge_exists(const edge_t& e) const
{
  return find(e.first, e.second) >= 0;
}

// remove all vertices without edges from graph (requires owned storage)
void Graph::remove_isolated_vertices(){
  vector<int> new_id(N);

  int u_new = 0;
  for(int u=0; u<N; u++)
    if(!(*this)[u].empty())
      new_id[u] = u_new++;

  int N_new = u_new;
  Graph g(N_new, dmax);
  // cerr << "n new: " << N_new << endl;
  for(int u=0; u<N; u++)
    for(int v: (*this)[u])
      g.push_back(new_id[u], new_id[v]);

  *this = g;
}

// completely remove all vertices in sv from the graph (requires owned storage)
void Graph::remove_vertices(set<int> &sv){
  const int N_naught(N);
  for(int u: sv){
    while((*this)[u].size()){
      const int v = (*this)[u][0];
      remove_edge({u,v});
    }
  }

  remove_isolated_vertices();

  // let's see if the graph remained in a sane state
  // cerr << "N: " << N << endl;
  if(N_naught != sv.size() + N)
    cerr << "removed more vertices than intended" << endl;
  assert(is_connected());
}

void Graph::apply_permutation(const Permutation& pi)
{
  assert(owns_memory());
  assert((int)pi.size() == N);
  std::vector<node_t> new_neighbours(N * dmax, node_t(-1));
  std::vector<uint8_t> new_deg(N, 0);
  for (int u_old = 0; u_old < N; ++u_old) {
    const int u_new = pi[u_old];
    new_deg[u_new] = owned_deg[u_old];
    for (int i = 0; i < owned_deg[u_old]; ++i) {
      const int t_old = owned_neighbours[u_old * dmax + i];
      new_neighbours[u_new * dmax + i] = pi[t_old];
    }
  }
  owned_neighbours = std::move(new_neighbours);
  owned_deg        = std::move(new_deg);
  repoint();
}

// arc_ix / next / prev / next_on_face / prev_on_face are inline in
// graphview.hh so device code can call them.

// ---------------------------------------------------------------------------
// Orientation: which surface the rotation system embeds the graph in.  The
// mathematics and the contracts are at the declarations in graphview.hh.
// ---------------------------------------------------------------------------

// The faces of that embedding -- the orbits of phi(u->v) = v -> next(v,u) --
// handed to `visit` one at a time, as the vertex each orbit was entered at.
//
// Each orbit is traversed exactly ONCE, following phi until the walk returns to
// its STARTING ARC.  Stopping at the starting VERTEX instead (as this walk's
// predecessor did) cuts an orbit into one segment per visit to that vertex --
// and a face through a cut vertex legitimately visits it more than once -- so a
// counter over such walks counts one face several times.
//
// Per-arc state is a flat vector<uint8_t> indexed by the dense arc id (arcid):
// no per-arc tree node, no O(log E) insert/find/erase.  It is also the
// termination guard: under a permutation an orbit is a simple cycle, so
// stepping onto an already-visited arc other than the start proves phi is not
// one -- which is what parallel arcs and self-loops do to it, since find(v,u)
// resolves every arc between the same pair to the SAME slot.  Without that test
// such a walk enters a cycle not containing a0 and never returns.
//
// @pre  symmetric: G.adjacency_is_symmetric()   (else phi is not total)
// @pre  simple:    no parallel arcs, no self-loops  (else phi is not injective)
// @throws std::logic_error when either @pre is violated
// @time O(E_dmax)
template<typename Visit>
static void for_each_face(const GraphView& G, Visit visit)
{
  vector<uint8_t> visited(size_t(G.N)*G.dmax, 0);

  for(node_t u0=0;u0<G.N;u0++)
    for(int i0=0;i0<G.deg[u0];i0++){
      const size_t a0 = G.arcid(u0,i0);
      if(visited[a0]) continue;

      visit(u0);
      size_t a = a0;
      do {
        visited[a] = 1;
        const node_t u = G.arc_of(a).first, v = G.neighbours[a];
        const int j = G.find(v,u);                   // v's slot holding u
        if(j < 0)
          throw std::logic_error("for_each_face: adjacency is not symmetric");
        a = G.arcid(v, (j+1)%G.deg[v]);              // phi(u->v) = v -> next(v,u)
        if(a != a0 && visited[a])
          throw std::logic_error("for_each_face: the arc successor is not a permutation "
                                 "-- the graph has parallel arcs or self-loops and this "
                                 "predicate is only valid for simple graphs");
      } while(a != a0);
    }
}

// Euler's verdict on one connected component, once its faces are counted:
// chi = N - E + F = 2 - 2g, so the realised genus is (E - N + 2 - F)/2 -- an
// integer, because every rotation system embeds in an orientable surface.
static OrientedSurface euler_verdict(int N, int E, int F, int genus)
{
  const int g = (E - N + 2 - F)/2;
  return {g == genus ? OrientedSurface::Code::Ok : OrientedSurface::Code::GenusMismatch, F, g};
}

OrientedSurface GraphView::oriented_surface(int genus) const
{
  const int E = int(count_edges());
  if(N == 0 || E == 0)          return {};                                       // Degenerate
  if(!adjacency_is_symmetric()) return {OrientedSurface::Code::AsymmetricAdjacency};

  int F = 0;
  for_each_face(*this, [&](node_t){ F++; });
  return euler_verdict(N, E, F, genus);
}

bool GraphView::is_consistently_oriented(int genus) const
{
  return oriented_surface(genus).code == OrientedSurface::Code::Ok;
}

vector<OrientedSurface> GraphView::component_surfaces(const vector<int>& genus) const
{
  const vector<vector<node_t>> components = connected_components();
  const int C = int(components.size());
  vector<OrientedSurface> surfaces(C);                       // Degenerate until measured

  if(!adjacency_is_symmetric()){
    // A missing reverse arc breaks the ARC SET, not one component's topology:
    // phi is then not a permutation, so no component has faces to compare
    // against a claim and none gets a topological verdict.
    for(OrientedSurface& s: surfaces) s.code = OrientedSurface::Code::AsymmetricAdjacency;
    return surfaces;
  }

  vector<int> component_of(N,-1);
  for(int c=0;c<C;c++) for(node_t u: components[c]) component_of[u] = c;

  // E_i off the partition (N_i is the partition), and -- because a phi-orbit
  // never leaves its component -- every face off ONE arc walk, attributed
  // through the vertex it was entered at.  No subgraph is built.
  vector<int> twoEc(C,0), Fc(C,0);
  for(node_t u=0;u<N;u++) twoEc[component_of[u]] += deg[u];
  for_each_face(*this, [&](node_t u){ Fc[component_of[u]]++; });

  // Independent per component -- the tallies above are read-only from here --
  // so this is a par::for_each / omp parallel for the day one is worth its
  // overhead.  (par:: lives in claude-projects/parallel-primitives, which the
  // library must not depend on; plain OpenMP is the in-library form.)
  for(int c=0;c<C;c++){
    const int N_c = int(components[c].size()), E_c = twoEc[c]/2;
    surfaces[c] = E_c == 0 ? OrientedSurface{}
                           : euler_verdict(N_c, E_c, Fc[c],
                                           c < int(genus.size()) ? genus[c] : 0);
  }
  return surfaces;
}

// What the caller was handed instead of the surface it asked for.  One sentence
// per OrientedSurface code, written once because every operation that refuses an
// unoriented graph refuses it for one of these three reasons -- and written with
// NUMBERS, because "not oriented" tells whoever reads the failure nothing they
// can act on.
static string surface_complaint(const OrientedSurface& S, int genus)
{
  switch(S.code){
  case OrientedSurface::Code::Degenerate:
    return "there is nothing to embed (no vertices, or no edges)";
  case OrientedSurface::Code::AsymmetricAdjacency:
    return "adjacency is not symmetric -- some arc u->v has no reverse v->u, so it "
           "has no faces at all";
  default:
    return "its rotation system has " + to_string(S.faces) + " faces, hence genus "
         + to_string(S.genus) + ", not " + to_string(genus);
  }
}

void require_oriented_surface(const GraphView& G, const char* operation, int genus)
{
  const OrientedSurface S = G.oriented_surface(genus);
  if(S.code == OrientedSurface::Code::Ok) return;

  throw unoriented_surface_error(string(operation) + ": needs a consistently oriented genus-"
                                 + to_string(genus) + " surface, but on these "
                                 + to_string(G.N) + " vertices and "
                                 + to_string(G.count_edges()) + " edges "
                                 + surface_complaint(S, genus),
                                 S, genus);
}

// TODO: Doesn't need to be planar and oriented, but is easier to write if it is. Make it work in general.
bool GraphView::has_separating_triangles() const
{
  require_oriented_surface(*this, "GraphView::has_separating_triangles");

  for(node_t u=0;u<N;u++){
    auto nu = (*this)[u];

    for(int i=0;i<nu.size();i++){
      node_t t = nu[i];
      node_t v = prev(u,t), w = next(u,t); // edges: u--t, u--v, u--w
      if(edge_exists({t,w}) && edge_exists({t,v}) && edge_exists({v,w})) return true;
    }
  }
  return false;
}


bool GraphView::adjacency_is_symmetric() const
{
  for(node_t u=0;u<N;u++){
    auto nu = (*this)[u];
    for(int i=0;i<nu.size();i++){
      auto nv = (*this)[nu[i]];

      bool symmetric = false;
      for(int j=0;j<nv.size();j++) if(nv[j] == u) symmetric = true;
      if(!symmetric) return false;
    }
  }
  return true;
}

// TODO: Should make two functions: one that takes subgraph (empty is trivially connected) and one that works on full graph.
bool GraphView::is_connected(const set<node_t> &subgraph) const
{
  vector<int> dist(N);
  if(!subgraph.empty()){
    node_t s = *subgraph.begin();
    single_source_shortest_paths(s,&dist[0]);

    for(node_t u: subgraph)
      if(dist[u] == INT_MAX) return false;

  } else {
    node_t s = 0; // Pick a node that is part of an edge
    for(;(*this)[s].empty();s++) ;
    assert(s < N);

    single_source_shortest_paths(s,&dist[0]);

    for(int i=0;i<dist.size();i++)
      if(dist[i] == INT_MAX) return false;
  }

  return true;
}

#include <queue>
vector<vector<node_t>> GraphView::connected_components() const
{
  vector<bool> done(N);

  vector<vector<node_t>> components;
  for(node_t u=0;u<N;u++)
    if(!done[u]){
      vector<node_t> component;

      done[u] = true;
      component.push_back(u);
      queue<node_t> Q;
      for(auto v: (*this)[u]) Q.push(v);

      while(!Q.empty()){
	node_t v = Q.front(); Q.pop();
	if(!done[v]){
	  done[v] = true;
	  component.push_back(v);

	  for(int i=0;i<(*this)[v].size();i++)
	    if(!done[(*this)[v][i]]) Q.push((*this)[v][i]);
	}
      }
      sort(component.begin(), component.end());
      components.push_back(component);
    }
  return components;
}


void GraphView::single_source_shortest_paths(node_t source, int *distances, size_t max_depth) const
{
  vector<node_t> queue_buf(N);
  Deque<node_t> queue(queue_buf);
  for(int u=0;u<N;u++) distances[u] = INT_MAX;

  distances[source] = 0;
  queue.push_back(source);

  while(!queue.empty()){
    node_t u = queue.pop_front();

    for(node_t v: (*this)[u]){
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
matrix<int> GraphView::all_pairs_shortest_paths(const unsigned int max_depth) const
{
  matrix<int>   d(N,N,INT_MAX);
  vector<node_t> queue_buf(N);
  Deque<node_t> queue(queue_buf);

  for(node_t u=0;u<N;u++){
    queue.push_back(u);    	// Enqueue source u
    d(u,u) = 0;

    while(!queue.empty()){
      node_t v = queue.pop_front();

      // Process children w of node v
      for(node_t w: (*this)[v]){
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
matrix<int> GraphView::all_pairs_shortest_paths(const vector<node_t> &V,
					    const unsigned int max_depth) const
{
  size_t M = V.size();
  matrix<int> D(M,M,INT_MAX);
  vector<int> d(N);

  vector<node_t> queue_buf(N);
  Deque<node_t> queue(queue_buf);

  for(int i=0;i<M;i++){
    for(int j=0;j<N;j++) d[j] = INT_MAX; // Mark all nodes as unvisited

    node_t source = V[i];
    d[source] = 0;
    queue.push_back(source);

    while(!queue.empty()){
      node_t u = queue.pop_front();

      // Process children of node u
      for(node_t v: (*this)[u]){
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


vector<node_t> GraphView::shortest_cycle(node_t s, const int max_depth) const
{
  face_t cycle;
  int Lmin = INT_MAX;
  for(node_t t: (*this)[s]){
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
vector<node_t> GraphView::shortest_cycle(const vector<node_t>& prefix, const int max_depth) const
{
  // Is this a valid start?
  for(int i=0;i+1<prefix.size();i++) assert(edge_exists({prefix[i],prefix[i+1]}));

  node_t s = prefix[0], t = prefix[1];

  if(max_depth == 3){ // Triangles need special handling
    switch(prefix.size()){
    case 2:
      // t must have a neighbor r that neighbours s
      for(node_t r: (*this)[t])
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
    for(node_t w: (*this)[u])
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


vector<int> GraphView::multiple_source_shortest_paths(const vector<node_t>& sources, const unsigned int max_depth) const
{
  vector<int>   distances(N,INT_MAX);
  vector<node_t> queue_buf(N);
  Deque<node_t> queue(queue_buf);

  for(node_t s: sources){
    distances[s] = 0;
    queue.push_back(s);
  }

  while(!queue.empty()){
    node_t v = queue.pop_front();

    for(node_t w: (*this)[v]){
      const edge_t edge(v,w);
      if(distances[w] == INT_MAX){ // Node not previously visited
	distances[w] = distances[v] + 1;
	if(distances[w] < max_depth) queue.push_back(w);
      }
    }
  }
  return distances;
}


void GraphView::flip_all_orientations()
{
  for(node_t u=0;u<N;u++) reverse((*this)[u].begin(), (*this)[u].end());
}

int GraphView::max_degree() const
{
  int max_d = 0;
  for(node_t u=0;u<N;u++) if(degree(u) > max_d) max_d = degree(u);
  return max_d;
}


vector<edge_t> GraphView::undirected_edges() const {
  set<edge_t> edges;
  for(node_t u=0;u<N;u++)
    for(int i=0;i<(*this)[u].size();i++)
      edges.insert(edge_t(u,(*this)[u][i]));
  return vector<edge_t>(edges.begin(),edges.end());
}

vector<arc_t> GraphView::directed_edges() const {
  set<arc_t> edges;
  for(node_t u=0;u<N;u++)
    for(int i=0;i<(*this)[u].size();i++)
      edges.insert(arc_t(u,(*this)[u][i]));
  return vector<arc_t>(edges.begin(),edges.end());
}

size_t GraphView::count_edges() const {
  // Don't use edge_set -- it's slow
  size_t twoE = 0;
  for(node_t u=0;u<N;u++)
    twoE += (*this)[u].size();

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
