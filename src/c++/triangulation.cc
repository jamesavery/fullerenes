#include <algorithm>
#include <sstream>
#include <stdexcept>
#include "fullerenes/triangulation.hh"
#include "fullerenes/unfold.hh"
#include "fullerenes/buckygen-wrapper.hh"

pair<node_t,node_t> TriangulationView::adjacent_tris(const arc_t& e) const
{
  node_t u  = e.first, v = e.second;
  node_t w1 = next_on_face(u,v), w2 = next_on_face(v,u);
  return make_pair(w1,w2);
}

vector<tri_t> TriangulationView::compute_faces_oriented() const
{
  // Per-arc "done" flags as a flat vector<uint8_t> indexed by the dense arc id
  // (arcid) -- no hashing, no allocation.
  vector<uint8_t> arc_done(size_t(N)*dmax, 0);
  vector<tri_t> triangles;
  triangles.reserve(2*(N-2));        // Most common case is cubic dual, but we no longer know it for sure.

  for(node_t u=0;u<N;u++){
    auto nu = nbrs(u);
    for(int i=0;i<(int)nu.size();i++){
      const size_t s_uv = arcid(u,i);          // directed edge u->v
      if(arc_done[s_uv]) continue;
      const node_t v = nu[i];
      const int jv = find(v,u);                // w = next_on_face(u,v) = prev(v,u)
      if(jv == -1){                            // reverse arc v->u missing: not an oriented triangulation
        ostringstream msg;
        msg << "compute_faces_oriented: arc " << u << "->" << v << " has no reverse arc "
            << v << "->" << u << " (graph is not a consistently oriented triangulation)";
        throw logic_error(msg.str());
      }
      const int    jw   = (jv - 1 + deg[v]) % deg[v];
      const size_t s_vw = arcid(v,jw);         // arc v->w
      const node_t w    = neighbours[s_vw];
      const int    ju   = find(w,u);           // arc w->u, the triangle's third edge
      const bool   has_wu = ju >= 0;           // a real edge in a valid triangulation; guarded for robustness
      const size_t s_wu = has_wu ? arcid(w,ju) : 0;
      // Every directed edge is part of exactly one triangle, so arc_done[s_uv]
      // alone would suffice; the v->w / w->u checks mirror the original guard.
      if(!arc_done[s_vw] && (!has_wu || !arc_done[s_wu])){
        triangles.push_back(tri_t(u,v,w));
        arc_done[s_uv] = 1;
        arc_done[s_vw] = 1;
        if(has_wu) arc_done[s_wu] = 1;
      }
    }
  }
  return triangles;
}

/*
// Build precalculated lookup tables
//   1. cubic arc -> dual arc   (transverse arc a,i -> u,j)     : Nx3     node_t,uint8_t 
//   2. dual arc  -> cubic arc  (transverse arc u,j -> a,i)     : NfxFmax node_t,uint8_t 
//
// TODO: * Resolve function overlap with cubic_faces() and translate_arcs. This produces faster data structures.
//       * Gather precalculated tables into one data structure (separate from graph itself)
void Triangulation::compute_lookup_tables(const PlanarGraph&            cubic_graph,
					  vector<pair<node_t,uint8_t>>& carc_to_darc,
					  vector<pair<node_t,uint8_t>>& darc_to_carc,
					  vector<uint8_t>&              degrees,
					  int Fmax = 6
					  ) const
{
  size_t Nf = N;		// To remind ourselves that this is the number of faces
  size_t Nc = 2*(N-2);		// This function assumes the triangulation is the dual of a cubic graph 

  //TODO: Convert asserts to throws with useful information
  assert(triangles.size() == Nc); // Make sure triangles array is initialized before calling
  assert(is_consistently_oriented());	  // We don't bother with non-oriented surfaces

  if(Fmax<1) Fmax = max_degree(); // Calculate fmax
  
  carc_to_darc.resize(Nc*3);
  darc_to_carc.resize(Nf*Fmax);
  darc_to_cnode.resize(Nf*Fmax);
  degrees.resize(Nf);
  
  IDCounter tid(triangles);	  // Look up cubic node id by triangle  
  
  for(node_t u=0;u<Nf;u++){
    const vector<node_t>& nu(neighbours[u]);
    degrees[u] = nu.size();
    
    for(int j=0;j<nu.size();j++){
      const node_t& v  = nu[j];    // Process triangulation arc f->g
        node_t wa = next_on_face(u,v); // Third triangle node in CW direction
        node_t wb = next_on_face(v,u); // Third triangle node in CCW direction

	// If either hl or hr fails, triangulation data is corrupted.
	if((wa<0) || (wb<0)){
	  printf("next_on_face(%d,%d) = prev(%d,%d) = -1\n",u,v,v,u);
	  printf("next_on_face(%d,%d) = prev(%d,%d) = -1\n",v,u,u,v);
	  printf("\tneighbours[%d] = ",u); for(auto x: (*this)[u]) cout << x << ' '; cout << endl;
	  printf("\tneighbours[%d] = ",v); for(auto x: (*this)[v]) cout << x << ' '; cout << endl;	  
	  assert((wa<0) || (wb<0));
	}

	tri_t t1 = {u,g,CW}, t2 = {g,u,CCW};

	node_t a = tid(t1.sorted()), b = tid(t2.sorted);
	
	darc_to_cnode[u*Fmax + j] = u; // Triangle/cubic node defined by {u,g,CW}

	// The transverse arc to u->g in the cubic graph is u->v.
	// Now we want to find the index i such that v = cubic_graph.neighbours[u][i].
	int i = cubic_graph.arc_ix({u,v});
	darc_to_carc[u*Fmax + j] = {a,i};
	carc_to_darc[a*3    + i] = {u,j};
      }
    }
  }
}
*/

// TODO: arc_t -> arc_t everywhere
//       arc   -> arc   everywhere
//       cubic nodes: a,b,c,...
//       dual  nodes: u,v,w,...
unordered_map<arc_t,arc_t> TriangulationView::arc_translation() const
{
  // TODO: Common metadata, calculate once
  IDCounter<tri_t> tri_numbers;
  auto tris = triangles();
  unordered_map<arc_t,arc_t> arc_translate(tris.size()*3);

  if(tris.size() != (N-2)*2){
    cout << "triangles = " << tris << endl;
    assert(tris.size() == (N-2)*2);
  }
  
  // Dual arcs
  for(node_t u=0;u<N;u++)
    for(node_t v: nbrs(u)){
      node_t wa = next_on_face(u,v), wb = next_on_face(v,u);
      tri_t  Ta = {u,v,wa}, Tb = {v,u,wb};
      node_t  a = tri_numbers(Ta.sorted()), b = tri_numbers(Tb.sorted());

      arc_translate[{u,v}] = {a,b};
    }
  return arc_translate;
}


// TODO: Factor out tri_numbers to do only once?
PlanarGraph TriangulationView::dual_graph() const
{
  IDCounter<tri_t> tri_numbers;

  auto tris = triangles();
  //assert(tris.size() == (N-2)*2);
  if(tris.size() != (N-2)*2){
    cout << "triangles = " << tris << endl;
    assert(tris.size() == (N-2)*2);
  }
  for(int i=0;i<tris.size();i++) tri_numbers.insert(tris[i].sorted());

  Graph A(tris.size(), vector<node_t>(3));

  for(node_t U=0;U<tris.size();U++){
    const tri_t& t(tris[U]);

    for(int i=0;i<3;i++){
      const node_t& u(t[i]), v(t[(i+1)%3]);
      node_t w(prev(u,v)); // TODO: CCW for buckygen -- will this give problems elsewhere?
      // TODO: Should this not be prev(v,u)? Or next(v,u)?

      A[U][i] = tri_numbers(tri_t(u,v,w).sorted());

      if(A[U][i] < 0){
          cerr << "Triangle " << tri_t(u,v,w).sorted() << " (opposite " << t << ") not found!\n";
          abort();
      }
    }
  }
  //  cerr << "A = " << A << endl;
  Graph G(A);
  return PlanarGraph(G);
};


vector<face_t> TriangulationView::cubic_faces() const
{
  vector<face_t> dfaces(N);

  IDCounter<tri_t> tri_numbers;
  auto tris = triangles();
  for(int i=0;i<tris.size();i++) tri_numbers.insert(tris[i].sorted());

  for(node_t u=0;u<N;u++){
    auto nu = nbrs(u);
    face_t f(nu.size());
    for(int i=0;i<nu.size();i++){
      node_t v=nu[i], w = nu[(i+1)%nu.size()]; // next_on_face(u,v);
      f[i] = tri_numbers(tri_t(u,v,w).sorted());
    }
    dfaces[u] = f;
  }
  return dfaces;
}


// ── Layer 1: SpiralBoundary ─────────────────────────────────────────
// The boundary is a deque representing the frontier between placed and
// unplaced faces in a triangulation under construction or extraction.
// Each entry is (node, open_valencies), where open_valencies counts
// how many more edges the node needs before it is fully enclosed.
//
// When a new node is placed, it connects to both boundary endpoints
// (front and back). This may trigger a chain reaction: if an endpoint
// becomes fully saturated (open_valencies hits 0), it is absorbed into
// the interior, exposing the next node — which also receives a connection
// from the new node, and may itself saturate, and so on.
//
// This chain reaction is called:
//   cascading (during windup) — pops saturated endpoints, inserting edges
//   draining  (during unwind) — pops saturated endpoints, read-only
//
// Two symmetric public operations:
//   wind(k, ins)  — windup: attach node k, insert edges via ins callback.
//   unwind()      — unwind: absorb saturated endpoints after peeling apex.
// Both return the total number of connections consumed, which is
// subtracted from the node's degree to get its open valency on the boundary.
struct SpiralBoundary : Deque<pair<node_t,int>> {
  using Deque::Deque;

  // Drain saturated endpoints from both ends after the apex connects
  // to the boundary. Each drain decrements the endpoint's open valency
  // (one connection to the apex); if it hits 0, the node is fully
  // enclosed — pop it, and the new endpoint also loses one valency
  // (it too is adjacent to the apex). Continue until an endpoint with
  // remaining valency is reached. Returns total connections consumed (≥ 2).
  int unwind() { return drain(FRONT) + drain(BACK); }

  void rotate() { push_back(front()); pop_front(); }

  // Attach node k to the boundary during windup construction.
  // First connects k to both endpoints (front and back), inserting
  // the edges into the graph. Then cascades: if either endpoint is
  // now saturated, pop it and connect k to the next endpoint, repeating
  // until both ends have remaining open valency.
  // Returns total connections made (≥ 2).
  template<class InsertEdge>
  int wind(node_t k, InsertEdge&& ins, bool best_effort = false) {
    connect(FRONT, k, other(FRONT).first, ins);
    connect(BACK,  k, other(BACK).first,  ins);
    int pu = 2;
    pu += cascade(FRONT, k, best_effort, ins);
    pu += cascade(BACK,  k, best_effort, ins);
    return pu;
  }

  // Close the final node during windup: it connects to every remaining
  // boundary node, completing the surface.
  template<class InsertEdge>
  void closeLast(node_t last, InsertEdge&& ins) {
    int n = size();
    for(int i = 0; i < n; i++)
      ins({last, front(i).first}, front((i+1)%n).first, front((i-1+n)%n).first);
  }

private:
  // Decrement this end's open valency (one connection consumed by the apex).
  // If it hits 0, the node is fully enclosed: pop it off the boundary,
  // and the new endpoint also loses a valency (the apex triangle spans it
  // too). Repeat until the endpoint has remaining open valency.
  // Returns the number of connections absorbed from this end (≥ 1).
  int drain(End d) {
    --end(d).second;
    int k = 1;
    while(end(d).second == 0) { pop(d); --end(d).second; k++; }
    return k;
  }

  // Insert edge {k → endpoint} into the planar embedding with oriented
  // successor hints, and decrement the endpoint's open valency.
  // The two successor slots (suc_uv, suc_vu) get the cascade_from node
  // and a structural neighbour; which slot gets which depends on direction.
  template<class InsertEdge>
  void connect(End d, node_t k, node_t cascade_from, InsertEdge& ins) {
    node_t target = end(d).first;
    if(d == FRONT)
      ins({k, target}, other(d).first, cascade_from);
    else
      ins({k, target}, cascade_from, end(d, 1).first);
    --end(d).second;
  }

  // After connect saturated the endpoint (open == 0), it is fully
  // enclosed: pop it and connect k to the new endpoint. Repeat while
  // the new endpoint is also saturated. This is the windup analog of
  // drain — same chain reaction, but with actual edge insertion.
  // The best_effort guard prevents closing the last two boundary nodes
  // prematurely when building partial triangulations.
  template<class InsertEdge>
  int cascade(End d, node_t k, bool best_effort, InsertEdge& ins) {
    int n = 0;
    while(end(d).second == 0) {
      if(best_effort && size() == 2 && other(d).second == 0) break;
      node_t old = end(d).first;
      pop(d);
      connect(d, k, old, ins);
      n++;
    }
    return n;
  }
};


// Windup: construct an oriented triangulation from a face-degree sequence.
// Mirrors Haskell reference: windupGeneralSpiral spiral jumps = init2 >> foldl' stepK >> closeLast
Triangulation::Triangulation(const vector<int>& spiral_string, const jumplist_t& j, const bool best_effort):
  base_t(int(spiral_string.size()))
{
  jumplist_t jumps = j;

  int max_boundary = 3 * static_cast<int>(ceil(sqrt(N))) + 12;
  vector<pair<node_t,int>> boundary_buf(max_boundary);
  SpiralBoundary B(boundary_buf);

  // NB: Capture via GraphView& reference, not Triangulation* this —
  // calling insert_edge through the derived pointer produces
  // incorrect neighbour ordering in the planar embedding.
  GraphView& g = *this;
  auto ins = [&g](const arc_t& e, node_t su, node_t sv){ g.insert_edge(e, su, sv); };

  // ── Initialize: place first two nodes ────────────────────────────────
  B.push_back({0, spiral_string[0] - 1});
  B.push_back({1, spiral_string[1] - 1});
  insert_edge({0, 1});

  // ── Main fold: stepK for each remaining node (Haskell reference: foldl' stepK) ─
  int N_final = best_effort ? N : N-1;
  for(int k = 2; k < N_final; k++){
    if(!jumps.empty() && k == jumps.front().first){
      for(int i = jumps.front().second; i > 0; --i) B.rotate();
      jumps.erase(jumps.begin());
    }

    int open = spiral_string[k] - B.wind(k, ins, best_effort);
    B.push_back({k, open});

    if(open <= 0){
      if(best_effort){
        if(open != 0){ cerr << "Windup: open=" << open << ", expected 0\n"; abort(); }
        for(int i = 0; i < B.size(); i++)
          if(B.front(i).second != 0){ cerr << "Windup: boundary not fully saturated\n"; abort(); }
        break;
      }
      cerr << "Cage closed but faces left (or otherwise invalid spiral)\n";
      abort();
    }
  }

  // ── Validate terminal state and close last face ──────────────────────
  if(!best_effort){
    if(B.size() != spiral_string.back()){
      cerr << "Cage not closed: " << B.size() << " boundary nodes, expected " << spiral_string.back() << "\n";
      cerr << "Incomplete triangulation = " << *this << "\n";
      abort();
    }
    for(int i = 0; i < B.size(); i++)
      if(B.front(i).second != 1){
        cerr << "Cage not closed: boundary node has " << B.front(i).second << " open valencies\n";
        abort();
      }
    B.closeLast(N-1, ins);
  } else {
    remove_isolated_vertices();
  }

}


Triangulation TriangulationView::GCtransform(unsigned k, unsigned l) const
{
  if(l==0) return halma_transform(k-1);

  Unfolding u(*this);
  Unfolding gcu(u*Eisenstein(k,l));
  Folding gcf(gcu);
  Triangulation t(gcf.fold());
   return t;
}

Triangulation TriangulationView::halma_transform(int m, vector<map<edge_t,node_t>>* face_grids) const {
  if(m<0) return Triangulation(*this);

  map<arc_t,vector<node_t>> arc_nodes;
  node_t v_new = N;

  // Create m new vertices for each edge, stored in both directions
  vector<edge_t> dual_edges = undirected_edges();
  for(const auto &e: dual_edges){
    vector<node_t> nodes;
    for(unsigned int i=0;i<m;i++) nodes.push_back(v_new++);
    arc_nodes[{e.first, e.second}] = nodes;
    arc_nodes[{e.second, e.first}] = vector<node_t>(nodes.rbegin(), nodes.rend());
  }

  // Build per-face halma grids and create interior nodes.
  auto tris = triangles();
  int n_faces = tris.size();
  vector<map<edge_t,node_t>> face_grid(n_faces);

  for(int i = 0; i < n_faces; i++){
    auto& grid = face_grid[i];
    const face_t& T = tris[i];
    const vector<node_t>& ns0(arc_nodes[{T[0],T[1]}]);
    const vector<node_t>& ns1(arc_nodes[{T[1],T[2]}]);
    const vector<node_t>& ns2(arc_nodes[{T[0],T[2]}]);

    grid[edge_t(0,0)]     = T[0];
    grid[edge_t(m+1,0)]   = T[1];
    grid[edge_t(m+1,m+1)] = T[2];
    for(int j=0;j<m;j++){
      grid[edge_t(0,j+1)]   = ns0[j];
      grid[edge_t(j+1,m+1)] = ns1[j];
      grid[edge_t(j+1,j+1)] = ns2[j];
    }
    for(int j=1;j<m;j++)
      for(int k=j+1;k<=m;k++)
        grid[edge_t(j,k)] = v_new++;
  }

  if(face_grids) *face_grids = face_grid;

  node_t N_new = v_new;

  // Build next_on_face map from micro-triangles.
  // Each original face is subdivided into up-triangles and down-triangles.
  // Both types inherit the original face's consistent winding.
  // nof[{u,v}] = w means: in the face (u,v,w), w is the CW successor
  // of v in u's neighbour list.
  map<arc_t, node_t> nof;
  vector<node_t> first_nb(N_new, -1);

  for(int i = 0; i < n_faces; i++){
    const auto& grid = face_grid[i];

    // Up-triangles: face (v, left, dn)
    for(int j=0;j<=m;j++)
      for(int k=j+1;k<=m+1;k++){
        node_t v    = grid.at(edge_t(j,k));
        node_t left = grid.at(edge_t(j,k-1));
        node_t dn   = grid.at(edge_t(j+1,k));

        nof[{v,left}]  = dn;
        nof[{left,dn}] = v;
        nof[{dn,v}]    = left;

        if(first_nb[v]    < 0) first_nb[v]    = left;
        if(first_nb[left] < 0) first_nb[left] = dn;
        if(first_nb[dn]   < 0) first_nb[dn]   = v;
      }

    // Down-triangles: face (v, dn, dnr)
    for(int j=0;j<m;j++)
      for(int k=j+1;k<=m;k++){
        node_t v   = grid.at(edge_t(j,k));
        node_t dn  = grid.at(edge_t(j+1,k));
        node_t dnr = grid.at(edge_t(j+1,k+1));

        nof[{v,dn}]   = dnr;
        nof[{dn,dnr}] = v;
        nof[{dnr,v}]  = dn;
      }
  }

  // Build CW neighbour lists by chaining next_on_face.
  Graph neighbours(N_new, GRAPH_DMAX);
  for(node_t u = 0; u < N_new; u++){
    node_t v = first_nb[u];
    node_t w = v;
    do {
      neighbours.push_back(u, w);
      w = nof.at({u,w});
    } while(w != v);
  }

  return Triangulation(neighbours);
}


// *********************************************************************
//                 SPIRAL STUFF
// *********************************************************************


bool TriangulationView::get_spiral(const node_t f1, const node_t f2, const node_t f3,
                               vector<int> &spiral, jumplist_t& jumps, vector<node_t>& permutation,
                               const bool general) const {
  return get_spiral_implementation(f1,f2,f3,spiral,jumps,permutation,general);
}


// ── Layer 2: RemainingGraph ─────────────────────────────────────────
// Tracks which nodes and edges are still present during spiral extraction
// using a per-node bitmask into the original neighbour list (max degree 16).
//
// Removing node v clears its bit in each neighbour's active mask,
// effectively deleting all edges incident to v. The last remaining
// neighbour is cached in last_nbr for the terminal validation step.
struct RemainingGraph {
  const TriangulationView& tri;
  vector<uint16_t> active;
  int count;
  node_t last_nbr = -1;

  explicit RemainingGraph(const TriangulationView& t)
    : tri(t), count(t.N), active(t.N)
  {
    assert(t.max_degree() <= 16);
    for(int u = 0; u < count; u++)
      active[u] = (1u << t.degree(u)) - 1;
  }

  // Remove v from the graph: clear v's bit in each active neighbour's
  // mask (deleting all incident edges) and zero v's own mask.
  // Caches an arbitrary active neighbour before clearing — this becomes
  // the final node when only one remains.
  void remove(node_t v) {
    last_nbr = tri.nbrs(v)[__builtin_ctz(active[v])];
    count--;
    for(uint16_t m = active[v]; m; m &= m-1) {
      int    i = __builtin_ctz(m);
      node_t w = tri.nbrs(v)[i];
      active[w] &= ~(1u << tri.arc_ix(w, v));
    }
    active[v] = 0;
  }

  // Test whether v is a cut vertex: would removing v disconnect its
  // active neighbours? Since neighbours are in cyclic (planar) order,
  // the active ones form a contiguous arc iff every consecutive pair
  // shares an edge. If any gap exists, the path is broken and v is
  // a cut vertex.
  // Used by generalPeel to skip cut vertices (which would strand faces).
  bool is_cut_vertex(node_t v) const {
    uint16_t m = active[v];
    int n = __builtin_popcount(m);
    if(n < 2) return false;

    node_t nbrs[16];
    int k = 0;
    for(uint16_t t = m; t; t &= t-1)
      nbrs[k++] = tri.nbrs(v)[__builtin_ctz(t)];

    int edges = 0;
    for(int i = 0; i < k; i++) {
      node_t u = nbrs[i], w = nbrs[(i+1) % k];
      int j = tri.arc_ix(u, w);
      if(j >= 0 && (active[u] & (1u << j))) edges++;
    }
    return edges < k - 1;  // connected arc of k nodes needs exactly k-1 edges
  }
};


// ── Layer 3: SpiralState ────────────────────────────────────────────
// Composes RemainingGraph + SpiralBoundary + jump bookkeeping for
// spiral extraction (unwind direction).
//
// The "apex" is the unplaced node completing the triangle formed by
// the two boundary endpoints (back, front). Each peel step identifies
// the apex, removes it from the remaining graph, unwinds the boundary,
// and records the node's degree in the spiral output.
//
// regularPeel: the apex is always unique and peelable.
// generalPeel: if the apex is a cut vertex (removing it would disconnect
//   the remaining graph), rotate the boundary and record a jump.
//   Try successive rotations until a non-cut apex is found.
struct SpiralState {
  RemainingGraph R;
  vector<pair<node_t,int>> buf;
  SpiralBoundary B;
  jumplist_t& jumps;
  bool CCW;
  int jump_count = 0;
  int step = 3;

  SpiralState(const TriangulationView& t, node_t f1, node_t f2, node_t f3,
              jumplist_t& jumps_out)
    : R(t)
    , buf(3 * (int)ceil(sqrt(t.N)) + 12)
    , B(buf)
    , jumps(jumps_out)
    , CCW(t.next(f1,f2) == f3)
  {}

  bool valid(node_t f1, node_t f2, node_t f3) const {
    return CCW || R.tri.prev(f1,f2) == f3;
  }

  int deg(node_t v) const { return R.tri.degree(v); }

  // The apex is the unplaced node completing the triangle between the
  // boundary endpoints. Which side of the (back→front) arc it lies on
  // depends on the winding orientation (CCW vs CW).
  node_t apex() const {
    node_t u = B.back().first, w = B.front().first;
    return CCW ? R.tri.prev(u,w) : R.tri.next(u,w);
  }

  // Place node f at position pos in the spiral. Remove it from the
  // remaining graph and add it to the boundary with deg(f)-2 open
  // valencies (2 are consumed connecting to the existing boundary endpoints).
  void place(int pos, node_t f, vector<int>& spiral, vector<node_t>& perm) {
    spiral[pos] = deg(f);
    perm[pos]   = f;
    R.remove(f);
    B.push_back({f, deg(f) - 2});
  }

  node_t peel(bool general) {
    node_t v = general ? generalPeel() : regularPeel();
    step++;
    return v;
  }

private:
  node_t regularPeel() {
    node_t v = apex();
    if(v == -1 || R.active[v] == 0) return -1; // Invalid or already-removed apex

    R.remove(v);  // Remove the apex vertex v from the remaining graph R...
    int open = deg(v) - B.unwind(); // ...and place it on the boundary B with its open valency
    if(open <= 0) return -1;   
    B.push_back({v, open});
    return v;
  }

  node_t generalPeel() {
    for(int rot = 0; rot < B.size(); rot++) {
      node_t v = apex();
      if(v == -1) return -1;

      // Rotate past removed nodes and cut vertices (matches Haskell generalPeelStep).
      if(R.active[v] == 0 || R.is_cut_vertex(v)){
	B.rotate(); jump_count++;
	continue;
      }

      // Commit accumulated jump if any, then peel v from the boundary.
      if(jump_count > 0) { jumps.push_back({step, jump_count}); jump_count = 0; }

      R.remove(v);
      int open = deg(v) - B.unwind();
      if(open < 1) return -1;
      B.push_back({v, open});
      return v;
    }
    return -1;
  }
};


// Extract a (general) spiral from a starting triple.  Mirrors the Haskell
// decomposition: orient → init3 → fold peel → validateTerminal.
bool TriangulationView::get_spiral_implementation(const node_t f1, const node_t f2, const node_t f3, vector<int> &spiral,
                                              jumplist_t& jumps, vector<node_t> &permutation,
                                              const bool general, const vector<int>& S0, const jumplist_t &J0) const {
  spiral.assign(N, 0);  permutation.assign(N, 0);  jumps.clear();

  SpiralState S(*this, f1, f2, f3, jumps);
  if(!S.valid(f1, f2, f3)) return false;

  S.place(0, f1, spiral, permutation);
  S.place(1, f2, spiral, permutation);
  S.place(2, f3, spiral, permutation);

  if(!S0.empty() && (spiral[0] != S0[0] || spiral[1] != S0[1] || spiral[2] != S0[2])) return false;

  for(int i = 3; i < N-1; ++i) {
    node_t v = S.peel(general);
    if(v == -1) return false;
    spiral[i]      = S.deg(v);
    permutation[i] = v;
    if(!S0.empty() && spiral[i] != S0[i]) return false;
  }

  // ── Terminal: exactly one node remains. The boundary must form its
  // complete star — one open valency per boundary node, and boundary
  // size equals the last node's degree.
  if(S.R.count != 1) return false;
  int last_valency = S.deg(S.R.last_nbr);
  if(S.B.size() != last_valency) return false;
  for(int i = 0; i < S.B.size(); i++)
    if(S.B.front(i).second != 1) return false;

  spiral[N-1]      = last_valency;
  permutation[N-1] = S.R.last_nbr;
  return true;
}

void TriangulationView::get_all_spirals(vector<vector<int>>& spirals, vector<jumplist_t>& jumps,
                     vector<vector<node_t>>& permutations,
                     const bool only_special, const bool general) const
{
  vector<node_t> node_starts;

  vector<int> spiral(N);
  jumplist_t  jump;
  vector<node_t> permutation;

  // Prefer special nodes. TODO: Automatic renumber in order of degrees.
  for(node_t u=0;u<N;u++) if(degree(u) != 6) node_starts.push_back(u);
  for(node_t u=0;u<N;u++)
    if(!only_special && degree(u) == 6) node_starts.push_back(u);

  for(int i=0; i<node_starts.size(); i++){ // Looks like O(N^3), is O(N)
    const node_t u=node_starts[i];
    auto nu = nbrs(u);

    for(int j=0;j<nu.size();j++){
      node_t v=nu[j], w[2];
      w[0] = prev(u,v);
      w[1] = next(u,v);

      for(int k=0;k<2;k++){
        if(get_spiral(u,v,w[k],spiral,jump,permutation,general)){
          spirals.push_back(spiral);
          jumps.push_back(jump);
          permutations.push_back(permutation);
        }
      }
    }
  }
}


bool TriangulationView::get_spiral(vector<int>& spiral, jumplist_t& jumps, const bool only_rarest_special, const bool general, const bool CW_only) const
{
  vector<vector<node_t>> permutations;
  bool success = get_spiral(spiral,jumps,permutations,only_rarest_special,general,CW_only);
  return success;
}

general_spiral TriangulationView::get_general_spiral(const bool only_rarest_special, const bool CW_only) const
{
  general_spiral gs;
  bool success = get_spiral(gs.spiral_code,gs.jumps,only_rarest_special,true,CW_only);
  assert(success); 		// General spirals should *always* succeed
  return gs;
}

// Canonical spiral search: find the lexicographically smallest (general) spiral.
// Mirrors Haskell: canonicalGeneralSpiral g = tryRegular triples <|> tryGeneral triples
//   where triples = startingTriples g
bool TriangulationView::get_spiral(vector<int> &spiral, jumplist_t &jumps, vector<vector<node_t>>& permutations, const bool only_rarest_special, const bool general, const bool CW_only) const
{
  permutations.clear();

  // ── Starting nodes: rarest non-hexagonal degree (Haskell: rarestSpecial) ──
  vector<node_t> node_starts;
  if(only_rarest_special){
    int max_deg = 0;
    for(node_t u = 0; u < N; u++) max_deg = max(max_deg, (int)degree(u));
    vector<int> count(max_deg + 1, 0);
    for(node_t u = 0; u < N; u++) count[degree(u)]++;

    int rarest_deg = 0, rarest_count = INT_MAX;
    for(int d = 0; d <= max_deg; d++)
      if(d != 6 && count[d] > 0 && count[d] < rarest_count)
        { rarest_deg = d; rarest_count = count[d]; }

    for(node_t u = 0; u < N; u++)
      if((int)degree(u) == rarest_deg) node_starts.push_back(u);
  } else {
    for(node_t u = 0; u < N; u++) node_starts.push_back(u);
  }

  // ── Temporaries and best-tracking ─────────────────────────────────────────
  vector<int> spiral_tmp(N), permutation_tmp(N);
  jumplist_t jumps_tmp;
  spiral = vector<int>(1, INT_MAX);                    // sentinel: any real spiral wins
  jumps  = jumplist_t(100, make_pair(0, 0));           // sentinel: any real jumps win

  // ── Try all starting triples, track lexicographic minimum ─────────────────
  // Haskell reference: mapMaybe (regularSpiral g) (startingTriples g), then minimum
  auto tryAllTriples = [&](bool use_general) -> bool {
    bool found = false;
    for(node_t u : node_starts) {
      for(node_t v : nbrs(u)) {
        node_t w[2];
        int n_tries;
        if(CW_only) { w[0] = prev(u,v); n_tries = 1; }          // CW only: prev(f1,f2)==f3
        else        { w[0] = prev(v,u); w[1] = next(v,u); n_tries = 2; }

        for(int k = 0; k < n_tries; k++){
          if(!get_spiral_implementation(u,v,w[k], spiral_tmp,jumps_tmp,permutation_tmp,
                                        use_general))
          {
            if(use_general){ fprintf(stderr, "General spiral failed -- this should never happen!\n"); abort(); }
            continue;
          }
          found = true;
          general_spiral gs_tmp{jumps_tmp,spiral_tmp}, gs_best{jumps,spiral};
          if(gs_tmp < gs_best){
            jumps  = jumps_tmp;
            spiral = spiral_tmp;
            permutations.clear();
          }
          if(!(gs_best < gs_tmp))
            permutations.push_back(permutation_tmp);
        }
      }
    }
    return found;
  };

  // ── Two-phase search: regular first, general only if needed ───────────────
  bool found = tryAllTriples(false);
  if(general && !found)
    found = tryAllTriples(true);

  return spiral.size() == N;
}


void TriangulationView::symmetry_information(int N_generators, Graph& coxeter_diagram, vector<int>& coxeter_labels) const
{
  vector<vector<int>> spirals;
  vector<vector<node_t>> permutations;
  vector<jumplist_t>    jumps;

  get_all_spirals(spirals,jumps,permutations,true,true);

  // // Now get the spiral corresponding to this triangulation
  // // TODO: define S,J,P
  // get_spiral(0,1,2,S,J,P,true);
  // vector< int > group_actions; // Referring to the 'permutations' array.
  // for(int i=0;i<spirals.size();i++)
  //   if(J == jumps[i] && S == spirals[i]) group_actions.push_back(i);

  // Since finite point groups are Coxeter groups, they are generated
  // by a set of involutions. But... not every involution is a generator! What do?

  // Stuff for symmetry:
  //- Multiplication table can be represented as a directed graph.
  //-- Given any spanning tree, the leafs generate the group
  //-- Point groups are generated by involutions. Can we build up a spanning tree ending in involutions only?

}

vector<int> draw_path(int major, int minor)
{
  if(minor == 0) return {1,major};

  int slope = major/minor, slope_remainder = major%minor, slope_accumulator = 0;

  vector<int> paths(minor+1,0), runs(minor);

  for(int i=0; i<minor; i++){
    slope_accumulator += slope_remainder;

    paths[i+1] = paths[i] + slope + (slope_accumulator != 0);

    if((i+1<minor) && (slope_accumulator >= minor || slope_remainder == 0)){
      paths[i+1]++;
      slope_accumulator %= minor;
    }

    runs[i]    = paths[i+1]-paths[i];
  }

  //  cout << make_pair(major,minor) << " runlengths is " << runs << endl;

  return runs;
}




// Given start node u0 and adjacent face F_i, lay down triangles along the the straight
// line to Eisenstein number (a,b), and report what the final node is.
//
// Assumes a,b >= 1.
// TODO: Add special cases for (a,0) and (b,0) to make more general.
// TODO: Better name.
node_t TriangulationView::end_of_the_line(node_t u0, int i, int a, int b) const
{
  // Axis-aligned walks: a-th vertex along axis i is the iterated graph-edge
  // continuation through the oriented neighbour list.  The Eisenstein
  // "straight line" along (a, 0) is, at each step v coming from prev,
  //   axis_step(prev, v) = next(next(next(prev, v), v), v),
  // which is the same triple-`next` the rolling-square pumps; spelled out
  // here it's simply "walk a graph edges along the axis at u0".
  auto axis_step = [&](node_t prev, node_t v) {
    return next(next(next(prev, v), v), v);
  };
  if (a == 0 && b == 0) return u0;
  if (b == 0) {
    node_t prev_v = u0, v = nbrs(u0)[i];
    for (int k = 1; k < a; k++) { node_t nv = axis_step(prev_v, v); prev_v = v; v = nv; }
    return v;
  }
  if (a == 0) {
    // (0, b) at axis i is (b, 0) at axis (i+1) mod deg(u0).
    int j = (i + 1) % nbrs(u0).size();
    node_t prev_v = u0, v = nbrs(u0)[j];
    for (int k = 1; k < b; k++) { node_t nv = axis_step(prev_v, v); prev_v = v; v = nv; }
    return v;
  }

  // Interior walks (a, b >= 1): track the 4-vertex rolling square as a
  // cursor and step it according to the run-length encoding of the line.
  node_t q,r,s,t;

  auto go_north = [&](){
    const node_t S(s), T(t);
    q = S; r = T; s = next(S,T); t = next(s,r);
  };

  auto go_east = [&](){
    const node_t R(r), T(t);
    q = R; s = T; r = next(s,q); t = next(s,r);
  };

  // Initial square
  q = u0;                  // (0, 0)
  r = nbrs(u0)[i];         // (1, 0)
  s = next(q, r);          // (0, 1)
  t = next(s, r);          // (1, 1)

  vector<int> runlengths = draw_path(max(a,b), min(a,b));

  for(int i=0;i<runlengths.size();i++){
    int L = runlengths[i];

    if(a>=b){                        // a is major axis
      for(int j=0;j<L-1;j++)    go_east();
      if(i+1<runlengths.size()) go_north();
    } else {                        // b is major axis
      for(int j=0;j<L-1;j++)    go_north();
      if(i+1<runlengths.size()) go_east();
    }
  }

  return t;                        // End node is upper right corner.
}

// Returns:
// For each run-length, returns the following:
// t0--t1--t2--t3  -> [r0,t0, r1,t1, r2,t2, r3,t3] 
// | \ |  \ | \ |
// r0--r1--r2--r3
// The complete result is a vector of these vertex sequences.
// They are connected as follows when a>b
//        tt0--tt1-- ...
//          | \ |
// t0--t1--t2--t3--- ...
// | \ |  \ | \ |
// r0--r1--r2--r3
// I.e., always sharing an edge to the north, and when a<b
//     .. ..
//     | \ |
// s---t---tt
// | \ | \ | 
// s---t---tt
// | \ |
// s---t
// i.e., always sharing an edge to the East.
//   s---t
//  / \ /
// q---r
vector<vector<node_t>> TriangulationView::quads_of_the_line(node_t u0, int i, int a, int b,
                                                            vector<vector<Eisenstein>>* coords_out) const
{
  node_t q,r,s,t;                // Current square

  auto go_north = [&](){
    const node_t S(s), T(t); // From old square
    q = S; r = T; s = next(S,T); t = next(s,r);
  };

  auto go_east = [&](){
    const node_t R(r), T(t); // From old square
    q = R; s = T; r = next(s,q); t = next(s,r);
  };

  // Square one
  q = u0;                         // (0,0)
  r = nbrs(u0)[i];          // (1,0)
  s = next(q,r);                  // (0,1)
  t = next(s,r);                  // (1,1)

  // Otherwise, draw the line
  vector<int>            runlengths = draw_path(max(a,b), min(a,b));
  vector<vector<node_t>> quad_runs(runlengths.size());
  
  for(int i=0;i<runlengths.size();i++){
    int L = runlengths[i];
    vector<node_t> quad_run(2*L+2);
    
    if(a>=b){                        // a is major axis
      quad_run[0] = q; quad_run[1] = s;
      quad_run[2] = r; quad_run[3] = t;
      for(int j=0;j<L-1;j++)    {
	go_east();
	quad_run[2*(j+2)+0] = r;
	quad_run[2*(j+2)+1] = t;
      }

      quad_runs[i] = quad_run;      
      go_north();

    } else {                        // b is major axis
      quad_run[0] = q; quad_run[2] = s; 
      quad_run[1] = r; quad_run[3] = t;
      for(int j=0;j<L-1;j++) {
	go_north();
	quad_run[2*(j+2)+0] = s;
	quad_run[2*(j+2)+1] = t;	
      }

      quad_runs[i] = quad_run;      
      go_east();
    }
  }

  //  {
    //    node_t t = end_of_the_line(u0,i,a,b);
    //    printf("quads_of_the_line(%d,%d,(%d,%d)) -> %d: ",u0,i,a,b,t);
    //    cout << quad_runs << endl;
  //  }
  // Fill the per-vertex Eisenstein lattice coords in lockstep with the vertex
  // quads above (origin u0 = (0,0), axis i along (1,0); go_east -> +x, go_north
  // -> +y; corners q=(x,y) r=(x+1,y) s=(x,y+1) t=(x+1,y+1)). Mirrors the loop
  // above exactly; the caller's cross-check (terminal corner == (a,b),
  // vid == end_of_the_line) guards against any drift. Valid for a,b >= 1
  // (the axis-aligned b==0 case is the plain graph walk, handled by callers).
  if(coords_out){
    coords_out->assign(runlengths.size(), vector<Eisenstein>());
    int cx=0, cy=0;
    for(size_t ri=0; ri<runlengths.size(); ri++){
      int L = runlengths[ri];
      vector<Eisenstein> cr(2*L+2);
      if(a>=b){                                    // a is major axis (go_east in-run)
        cr[0]=Eisenstein(cx,cy);   cr[1]=Eisenstein(cx,cy+1);
        cr[2]=Eisenstein(cx+1,cy); cr[3]=Eisenstein(cx+1,cy+1);
        for(int j=0;j<L-1;j++){ cx++;
          cr[2*(j+2)+0]=Eisenstein(cx+1,cy); cr[2*(j+2)+1]=Eisenstein(cx+1,cy+1); }
        (*coords_out)[ri]=cr; cy++;               // go_north between runs
      } else {                                     // b is major axis (go_north in-run)
        cr[0]=Eisenstein(cx,cy);   cr[2]=Eisenstein(cx,cy+1);
        cr[1]=Eisenstein(cx+1,cy); cr[3]=Eisenstein(cx+1,cy+1);
        for(int j=0;j<L-1;j++){ cy++;
          cr[2*(j+2)+0]=Eisenstein(cx,cy+1); cr[2*(j+2)+1]=Eisenstein(cx+1,cy+1); }
        (*coords_out)[ri]=cr; cx++;               // go_east between runs
      }
    }
  }
  return quad_runs;                        // End node is upper right corner.
}

matrix<int> TriangulationView::pentagon_distance_mtx() const {
  vector<int> pentagon_indices(12);
  for(int u=0, i=0;u<N;u++) if(degree(u) == 5) pentagon_indices[i++] = u;
  return all_pairs_shortest_paths(pentagon_indices);
}


// TODO: Do we need to do Dijkstra on sqrt(H) after all?
matrix<TriangulationView::simple_geodesic>
TriangulationView::simple_geodesics(vector<node_t> nodes,
				bool calculate_self_geodesics,
				node_t trace_u,
				vector<geodesic_step>* trace_out,
				int* M_out) const
{
  if(nodes.empty()){
    nodes.resize(N);
    for(int i=0;i<N;i++) nodes[i] = i;
  }

  vector<int> nodes_inverse(N,-1);
  for(node_t U=0;U<nodes.size();U++){
    node_t     u = nodes[U];
    nodes_inverse[u] = U;
  }

  // Initialize H to squared graph distances (upper bound to squared surface
  // distance: 3/4 d_g^2 <= d_surface^2 <= d_g^2).
  matrix<int>             H(nodes.size(),nodes.size(),all_pairs_shortest_paths(nodes));
  matrix<simple_geodesic> G(nodes.size(),nodes.size());

  vector<int> M(nodes.size(),0);	// M[u] = max_v(d_g(u,v)) is upper bound to surface distance from u
  for(node_t U=0; U<nodes.size();U++)
    for(node_t V=0;V<nodes.size();V++)
      M[U] = max(M[U], H(U,V));

  for(int i=0;i<H.size();i++) H[i] *= H[i];     // Work with squared distances

  if(calculate_self_geodesics) for(node_t U=0;U<nodes.size();U++){
      H(U,U) = INT_MAX;     // Self-geodesics: search for the shortest non-trivial loop
      M[U]  *= 2;           // ...by walking up to twice the diameter
    }

  // Report the traced source's search radius (after any self-geodesic doubling).
  if(trace_u>=0 && M_out && nodes_inverse[trace_u]!=-1) *M_out = M[nodes_inverse[trace_u]];

  // Loop bounds inclusive (a <= M[U], a^2+ab+b^2 <= M[U]^2) so the
  // axis-aligned (a,0) geodesic of length a == d_g(u,v) is recorded
  // when it matches H(U,V) = d_g(u,v)^2. Tie-breaking: last walk wins.
  for(node_t u: nodes){
    for(int i=0;i<degree(u);i++){
      node_t U  = nodes_inverse[u];

      for(int a=1; a<=M[U]; a++){
	for(int b=0; a*a + a*b + b*b <= M[U]*M[U]; b++){
	  const node_t v = end_of_the_line(u,i,a,b);

	  if(nodes_inverse[v] != -1){ // Endpoint v is in nodes
	    node_t V = nodes_inverse[v];
	    int d_sqr = a*a + a*b + b*b;
	    const bool improved = (d_sqr <= H(U,V));
	    if(trace_out && u==trace_u)        // record the probe before the update
	      trace_out->push_back({i,a,b,v,d_sqr,H(U,V),improved});
	    if(improved){
	      H(U,V) = d_sqr;
	      G(U,V) = simple_geodesic(a,b,i);
	    }
	  }
	}
      }
    }
  }
  return G;
}

TriangulationView::geodesic
TriangulationView::compose_simple_geodesics(const vector<int>& path,
                                        const matrix<simple_geodesic>& simple)
{
  geodesic G;
  if (path.size() < 2) return G;
  G.segments.reserve(path.size() - 1);
  for (size_t k = 0; k + 1 < path.size(); k++)
    G.segments.push_back(simple(path[k], path[k+1]));
  return G;
}

// surface_geodesics: simple-geodesics + APSP-with-paths + path
// reconstruction + per-pair composition.
matrix<TriangulationView::geodesic>
TriangulationView::surface_geodesics(vector<node_t> nodes,
                                 bool calculate_self_geodesics) const
{
  matrix<geodesic> G(0, 0, geodesic());
  surface_distances(nodes, calculate_self_geodesics, &G);
  return G;
}

// ── File-local helpers for the level-synchronous simple_square_surface_distances ──
//
// The Eisenstein ray sweep is reordered so rays are grouped by squared length
// d2 = a^2+ab+b^2 (the Loeschian numbers) and processed in increasing-d2
// "levels". This is bit-identical to the old source-major (u,i,a,b) loop but
// (a) hoists each rep's run-length tape out of the per-(source,sector) fan and
// (b) lets a dynamic d2 > max(H) cut and active-source compaction skip work.
// It is also parallel-ready: within one level every write to H stores the SAME
// d2, so the per-level sweep is a race-free parallel-for and the levels are the
// only synchronisation points. The `PARALLEL:` markers below show where the
// parallel-primitives ctx (par::seq / par::omp / SYCL) drops in once that API
// settles -- until then this runs serially with no atomics or OpenMP.
namespace {

// Loeschian level table: reps {(a,b) : a>=1, b>=0, a^2+ab+b^2 < M_max^2},
// sorted by norm2() and grouped into constant-d2 levels.  The (a,0) axis reps
// are included so axis self-closures reach the self-diagonal (see the loop
// comment below); off-diagonal they idempotently re-min the graph distance.
struct EisensteinLevels {
  vector<Eisenstein> reps;     // sorted by norm2(), tie-break (a,b); level-contiguous
  vector<int>        first;    // level l = reps[first[l], first[l+1)); size()+1 entries
  int size()  const { return (int)first.size() - 1; }
  int d2(int l) const { return reps[first[l]].norm2(); }
  std::span<const Eisenstein> level(int l) const {
    return {reps.data() + first[l], size_t(first[l+1] - first[l])};
  }
};

EisensteinLevels loeschian_levels(int M_max)
{
  EisensteinLevels L;
  for(int a=1; a*a < M_max*M_max; a++)                          // (a,0): smallest is a^2
    for(int b=0; Eisenstein(a,b).norm2() < M_max*M_max; b++)
      L.reps.push_back(Eisenstein(a,b));

  std::sort(L.reps.begin(), L.reps.end(), [](Eisenstein g, Eisenstein h){
    return std::tuple(g.norm2(), g) < std::tuple(h.norm2(), h);
  });
  for(int r=0; r<(int)L.reps.size(); r++)
    if(r==0 || L.reps[r].norm2() != L.reps[r-1].norm2()) L.first.push_back(r);
  L.first.push_back(L.reps.size());
  return L;
}

// Concatenated per-rep turn tapes for one level: rep r's tape is
// tape[tape_first[r], tape_first[r+1)).  Axis reps (b==0) get an empty segment
// -- they are walked directly, and draw_path is undefined for minor==0.
void build_tapes(std::span<const Eisenstein> reps, vector<int>& tape, vector<int>& tape_first)
{
  tape.clear();
  tape_first.assign(1, 0);
  for(auto [a,b] : reps){
    if(b >= 1){
      vector<int> t = draw_path(max(a,b), min(a,b));            // reuse the file-local DDA
      tape.insert(tape.end(), t.begin(), t.end());
    }
    tape_first.push_back(tape.size());
  }
}

// Tape-driven straight-line walk: the rolling-square pump of end_of_the_line,
// but driven by a precomputed tape (interior b>=1) or walked directly (axis
// b==0).  Hoisting the tape is the point; for a one-off walk use the member
// end_of_the_line instead.
node_t walk_line(const TriangulationView& G, node_t u0, int i,
                 Eisenstein g, std::span<const int> tape)
{
  if(g.second == 0){                                            // axis ray (a,0)
    auto axis_step = [&](node_t prev, node_t v){ return G.next(G.next(G.next(prev,v),v),v); };
    node_t prev = u0, v = G[u0][i];
    for(int k=1; k<g.first; k++){ node_t nv = axis_step(prev,v); prev = v; v = nv; }
    return v;
  }
  node_t q,r,s,t;
  auto go_north = [&]{ const node_t S(s),T(t); q=S; r=T; s=G.next(S,T); t=G.next(s,r); };
  auto go_east  = [&]{ const node_t R(r),T(t); q=R; s=T; r=G.next(s,q); t=G.next(s,r); };
  q = u0; r = G[u0][i]; s = G.next(q,r); t = G.next(s,r);       // initial square
  const bool a_major = g.first >= g.second;
  for(size_t k=0; k<tape.size(); k++){
    for(int j=0;j<tape[k]-1;j++) a_major ? go_east() : go_north();
    if(k+1 < tape.size())        a_major ? go_north() : go_east();
  }
  return t;                                                     // upper-right corner
}

} // anonymous namespace

// Squared surface distances via straight Eisenstein rays, refined by the
// APSP pass in surface_distances().
//
// A simple geodesic is a straight ray with no cone (deg != 6 vertex) in its
// interior. The segment u->(a,b) meets an interior lattice vertex iff
// g = gcd(a,b) > 1, at the points (j*a/g, j*b/g), j=1..g-1; a coprime (a,b)
// therefore crosses only flat edges/faces and is exactly a simple geodesic of
// length sqrt(a^2+ab+b^2), while a non-coprime ray is the g-fold repeat of its
// primitive and is a *composition* of simple geodesics, not a simple one.
//
// This routine shoots every (a,b), not just the coprimes. That is exact here
// because every cone has curvature kappa = (pi/3)(6 - deg) >= 0 (deg <= 6): a
// ray developed through a non-negative cone records sqrt(a^2+ab+b^2), which is
// an UPPER bound on the true distance to the landed vertex -- developing
// ignores the removed wedge, and re-gluing the wedge is 1-Lipschitz (distances
// only shrink). So the non-simple rays can only over-estimate; the entrywise
// min discards them, and the APSP pass recomposes the exact distance from the
// primitive cone-free segments. (For a negative cone, deg >= 7, a through-cone
// ray could under-estimate; only the coprime, cone-free enumeration would then
// be valid -- see coprime-simple-geodesics.tex, and the cone-free fan-sweep in
// geodesics/refactor-debt.md.)
matrix<int> TriangulationView::simple_square_surface_distances(vector<node_t> nodes,
							   bool calculate_self_geodesics) const
{
  if(nodes.empty()){ nodes.resize(N); for(int i=0;i<N;i++) nodes[i] = i;  } // If no node list is given, calculate all distances
  const int n = nodes.size(), max_deg = max_degree();

  vector<int> nodes_inverse(N,-1);
  for(int i=0;i<n;i++) nodes_inverse[nodes[i]] = i;

  // Initialize H to graph distances, an upper bound to surface distances:
  // 3/4 d_g^2 <= d_surface^2 <= d_g^2.  Off-diagonal axis-aligned rays (a,0)
  // have length == graph distance, so they are already captured here.
  matrix<int> H(n, n, all_pairs_shortest_paths(nodes));

  // M[U] = max_V d_g(u,v): per-source search radius (upper bound to surface dist).
  vector<int> M(n,0);
  for(node_t U=0; U<n; U++) for(node_t V=0; V<n; V++) M[U] = max(M[U], H(U,V));

  for(int i=0;i<H.size();i++) H[i] *= H[i];                     // work in squared distances

  if(calculate_self_geodesics) for(node_t U=0;U<n;U++){
      M[U] *= 2;             // a self-geodesic loops there and back: look twice as far
    }

  vector<int> M2(n);                                           // per-source search radii, squared
  for(node_t U=0; U<n; U++) M2[U] = M[U]*M[U];

  if(calculate_self_geodesics) for(node_t U=0;U<n;U++)
      H(U,U) = M2[U];        // Diagonal sentinel (2*M[U])^2: "no self-geodesic strictly
                             // inside the search radius".  Finite (unlike INT_MAX) so the
                             // max(H)-based cut below stays usable.  The self-diagonal is
                             // NOT covered by the graph-SP init (it would be 0 then this
                             // sentinel), so axis self-closures end_of_the_line(u,i,a,0)==u
                             // must be shot explicitly -- hence the (a,0) reps in the table.

  // Sources sorted by search radius: as d2 grows a contiguous prefix retires
  // (M2[U] <= d2 => U can shoot no further), so each level sweeps only the
  // still-active suffix order[lo..n) instead of pruning every item.
  vector<int> order(n);
  for(int U=0;U<n;U++) order[U] = U;
  std::sort(order.begin(), order.end(), [&](int a,int b){ return M2[a] < M2[b]; });

  const int M_max = M.empty() ? 0 : *std::max_element(M.begin(), M.end());
  const EisensteinLevels levels = loeschian_levels(M_max);
  vector<int> tape, tape_first;

  // Levels in increasing d2.  Each level shoots {(U,i,g) : g.norm2()==d2(l) <
  // M2[U]} over the active sources; the old bound a < Mu is implied by d2 < Mu^2,
  // so the union over levels is exactly the old candidate set.
  constexpr int bound_refresh_period = 16;     // n^2 max-scan period (guard, not a knob)
  int bound = INT_MAX, lo = 0;
  for(int l=0; l<levels.size(); l++){
    const int d2 = levels.d2(l);
    while(lo < n && M2[order[lo]] <= d2) lo++;                 // retire sources past their radius
    const int n_active = n - lo;
    if(n_active == 0) break;                                   // no source can shoot further

    if(l % bound_refresh_period == 0){
      // Dynamic global cut: once d2 exceeds every current entry, no later
      // (longer) ray can improve any entry.
      int hmax = 0;
      // PARALLEL: hmax = ctx.reduce(n*n, 0, par::max_op{}, [&](int k){ return H[k]; });
      for(int k=0;k<n*n;k++) hmax = max(hmax, H[k]);
      bound = hmax;
    }
    if(d2 > bound) break;

    const std::span<const Eisenstein> reps = levels.level(l);
    const int R = reps.size();
    build_tapes(reps, tape, tape_first);                       // hoisted once per level

    // One work item = one ray (active source, sector i, rep): walk, then min.
    // rep innermost so consecutive items share a tape.  Within a level every
    // write stores the same d2, so this loop is a race-free parallel-for.
    const int n_items = n_active * max_deg * R;
    // PARALLEL: ctx.for_each(n_items, shoot); ctx.barrier();  (shoot = loop body below,
    //           with H(U,V)=min(...) becoming a relaxed atomic_ref store of the equal d2)
    for(int flat=0; flat<n_items; flat++){
      const int rep = flat % R, i = (flat / R) % max_deg;
      const node_t U = order[lo + flat / (R * max_deg)];
      const node_t u = nodes[U];
      if(i >= degree(u)) continue;
      const node_t v = walk_line(*this, u, i, reps[rep],
                                 std::span<const int>(tape).subspan(tape_first[rep],
                                                                    tape_first[rep+1]-tape_first[rep]));
      const int V = nodes_inverse[v];
      if(V != -1) H(U,V) = min(H(U,V), d2);
    }
  }
  return H;
}



matrix<double> TriangulationView::surface_distances(vector<node_t> nodes,
						bool calculate_self_geodesics,
						matrix<geodesic>* geodesics_out) const
{
  if (geodesics_out == nullptr) {
    // Fast path: distances only, via the existing min-plus APSP.
    matrix<double> H(simple_square_surface_distances(nodes,calculate_self_geodesics));
    H = H.sqrt_elementwise();
    return H.APSP(false).square_elementwise();
  }

  // Pipeline: simple geodesics -> sqrt of simple distances -> path-tracking
  // APSP -> per-pair path reconstruction + composition.
  const matrix<simple_geodesic> simple = simple_geodesics(nodes, calculate_self_geodesics);
  const int n = simple.m;

  matrix<double> H(n, n, 0.0);
  for (int U = 0; U < n; U++)
    for (int V = 0; V < n; V++)
      H(U, V) = sqrt((double)simple(U, V).g.norm2());
  if (!calculate_self_geodesics)
    for (int U = 0; U < n; U++) H(U, U) = 0.0;

  APSPResult<double> apsp = H.APSP_with_paths();

  *geodesics_out = matrix<geodesic>(n, n, geodesic());
  for (int U = 0; U < n; U++)
    for (int V = 0; V < n; V++)
      (*geodesics_out)(U, V) =
        compose_simple_geodesics(reconstruct_path(apsp.preds, U, V), simple);

  return apsp.dist.square_elementwise();
}

// --- Thurston star unfolding (@anchor tri-star-unfold in graphview.hh;
//     construction + proofs: claude-projects/visualization/STAR-UNFOLDING.md).

namespace {

// w^k as an Eisenstein unit: CCW rotation by k * 60 degrees; rotation is
// Eisenstein multiplication.
Eisenstein w_pow(int k) {
  static const Eisenstein w[6] = {{1,0},{0,1},{-1,1},{-1,0},{0,-1},{1,-1}};
  return w[((k % 6) + 6) % 6];
}

// Lattice shoelace of a closed labelled (a,b)-polygon: counts covered unit
// lattice triangles, each contributing +-1 ((0,0)(1,0)(0,1) -> +1).
long long lattice_shoelace(const vector<pair<Eisenstein, node_t>>& poly) {
  long long s = 0;
  for (size_t i = 0; i < poly.size(); ++i) {
    const Eisenstein& p = poly[i].first;
    const Eisenstein& q = poly[(i + 1) % poly.size()].first;
    s += (long long)p.first * q.second - (long long)q.first * p.second;
  }
  return s;
}

// One cut source -> cone: the fan data (axis, g) kept raw for the exact
// angular comparator, the developed fan direction delta = w^axis * g, the
// Loeschian squared length, and the leaf's degree (cone angle / (pi/3)).
struct StarCut {
  node_t     cone;
  int        axis, deg;
  Eisenstein g, delta;
  long long  d2;
};

// Angular order within the developed fan [0, deg(source)*60): by sector,
// then by the in-sector cross product.  Total and exact; equal directions
// are the CollinearCuts degeneracy, detected before sorting.
bool fan_less(const StarCut& x, const StarCut& y) {
  if (x.axis != y.axis) return x.axis < y.axis;
  return (long long)x.g.first * y.g.second - (long long)y.g.first * x.g.second > 0;
}

}  // namespace

TriangulationView::star_unfolding
TriangulationView::star_unfold(node_t source, const vector<node_t>& cones,
                               const matrix<simple_geodesic>& simple,
                               const matrix<double>& dist) const
{
  using SU = star_unfolding;

  const int n = (int)cones.size();
  int       S = -1;
  for (int i = 0; i < n; ++i)
    if (cones[i] == source) S = i;
  if (S < 0) return SU::error(SU::Code::SourceNotCone, std::to_string(source));

  SU out;
  out.source = source;

  vector<StarCut> cuts;
  for (int I = 0; I < n; ++I) {
    if (I == S) continue;
    const simple_geodesic& sg = simple(S, I);
    const long long d2 = sg.g.norm2();
    if (d2 <= 0)
      return SU::error(SU::Code::NoSimpleGeodesic, std::to_string(cones[I]));
    if ((double)d2 > dist(S, I) + 1e-6) out.cuts_globally_shortest = false;
    cuts.push_back({cones[I], sg.axis, (int)degree(cones[I]),
                    sg.g, w_pow(sg.axis) * sg.g, d2});
  }
  for (const StarCut& x : cuts)              // collinear pre-scan (k is tiny)
    for (const StarCut& y : cuts)
      if (&x != &y && x.axis == y.axis && !fan_less(x, y) && !fan_less(y, x))
        return SU::error(SU::Code::CollinearCuts,
                         std::to_string(x.cone) + "," + std::to_string(y.cone));
  sort(cuts.begin(), cuts.end(), fan_less);

  // The walk (STAR-UNFOLDING.md: LEAF TURN + FAN FRAME, telescoped):
  //   v_i = s^i + w^(D_(i-1)) delta_i,   s^(i+1) = v_i - w^(D_i) delta_i.
  // Emitted directly as the standard CW labelled outline.
  Eisenstein s(0, 0);
  int        D = 0;
  for (const StarCut& c : cuts) {
    const Eisenstein v = s + w_pow(D) * c.delta;
    out.outline.push_back({s, source});
    out.outline.push_back({v, c.cone});
    D = (D + c.deg) % 6;
    s = v - w_pow(D) * c.delta;
  }

  if (!(s == Eisenstein(0, 0)))
    throw std::runtime_error("star_unfold: boundary does not close at (" +
                             std::to_string(s.first) + "," +
                             std::to_string(s.second) + ")");
  const long long n_tri = 2LL * N - 4;
  if (lattice_shoelace(out.outline) != -n_tri)
    throw std::runtime_error("star_unfold: area check failed: shoelace " +
                             std::to_string(lattice_shoelace(out.outline)) +
                             " != -" + std::to_string(n_tri));
  return out;
}

TriangulationView::star_unfolding
TriangulationView::star_unfold(node_t source) const
{
  vector<node_t> cones;
  for (node_t u = 0; u < (node_t)N; ++u)
    if (degree(u) != 6) cones.push_back(u);
  return star_unfold(source, cones, simple_geodesics(cones),
                     surface_distances(cones));
}

Triangulation TriangulationView::sort_nodes() const
{
  vector< pair<int,int> > degrees(N);

  for(int u=0;u<N;u++) degrees[u] = make_pair(degree(u), u);

  sort(degrees.begin(), degrees.end());

  vector<int> newname(N), oldname(N);

  int deg_val, u_old;
  for(node_t u_new=0;u_new<N;u_new++){
    tie(deg_val,u_old) = degrees[u_new];
    newname[u_old] = u_new;
    oldname[u_new] = u_old;
  }

  Graph new_neighbours(N, GRAPH_DMAX);
  for(int u=0;u<N;u++)
    for(int i=0;i<degree(u);i++)
      new_neighbours.push_back(newname[u], newname[nbrs(u)[i]]);

  return Triangulation(new_neighbours);
}

Permutation TriangulationView::sort_flat_last() const
{
  // Sort vertices so cone points (degree != 6) come first in original order,
  // then flat (degree-6) vertices in original order.  This is the correct
  // ordering for iDT vertex removal, which removes from the back.
  return argsort([this](node_t a, node_t b){
    const int ka = (degree(a) == 6) ? 1 : 0;
    const int kb = (degree(b) == 6) ? 1 : 0;
    if (ka != kb) return ka < kb;
    return a < b;
  });
}

spiral_nomenclature FullereneDualView::name(bool rarest_start) const
{
  return spiral_nomenclature(dual_graph(), spiral_nomenclature::FULLERENE,
			     spiral_nomenclature::CUBIC,
			     rarest_start);  
}

// call for one general spiral and extract the pentagon indices
bool FullereneDualView::get_rspi(const node_t f1, const node_t f2, const node_t f3, vector<int>& rspi, jumplist_t& jumps, const bool general) const
{
  rspi.resize(12);
  jumps.clear();
  vector<int> spiral;
  vector<node_t> permutation;
  if(!get_spiral(f1, f2, f3, spiral, jumps, permutation, general)) return false;

  for(int i=0,j=0;i<spiral.size();i++)
    if(spiral[i] == 5)
      rspi[j++] = i;

  return true;
}

// call for the canonical general spiral and extract the pentagon indices
bool FullereneDualView::get_rspi(vector<int>& rspi, jumplist_t& jumps, const bool general, const bool pentagon_start) const
{
  rspi.resize(12);
  jumps.clear();
  vector<int> spiral;
  if(!get_spiral(spiral, jumps, pentagon_start, general)) return false;

  for(int i=0,j=0;i<spiral.size();i++)
    if(spiral[i] == 5)
      rspi[j++] = i;

  return true;
}

general_spiral FullereneDualView::get_rspi(const bool rarest_start) const
{
  general_spiral S;
  get_rspi(S.spiral_code,S.jumps,true,rarest_start);
  return S;
}

// permutation of vertex numbers (ie, replace v by vertex_numbers[v], to get numbered vertices)
// where permutations are as returned by T.get_spiral(J,S,perm)
// locants are vertices that should have small vertex numbers (as far as permitted by symmetry equivalent canonical spirals)
vector<node_t> TriangulationView::vertex_numbers(vector<vector<node_t>> &permutations, const vector<node_t> &locants) const{
  assert(is_triangulation());
  vector<node_t> vertex_numbers(N);
  vector<node_t> vertex_numbers_inv(N,INT_MAX);
  for(int p=0; p<permutations.size(); p++){
    const vector<node_t> &vertex_numbers_tmp=permutations[p];
    // invert
    vector<node_t> vertex_numbers_inv_tmp(N);
    for(int i=0; i<vertex_numbers_tmp.size(); i++) vertex_numbers_inv_tmp[vertex_numbers_tmp[i]] = i;
    // copy to vertex_numbers_inv?
    if(locants.size()==0){
      vertex_numbers_inv = vertex_numbers_inv_tmp;
      break;
    }
    // compare two vectors, but only at chosen positions
    for(int l=0; l<locants.size(); l++){
      if(vertex_numbers_inv_tmp[locants[l]] > vertex_numbers_inv[locants[l]]) break;
      if(vertex_numbers_inv_tmp[locants[l]] < vertex_numbers_inv[locants[l]]){
        vertex_numbers_inv = vertex_numbers_inv_tmp;
        break;
      }
    }
  }
  //invert
  for(int i=0; i<vertex_numbers.size(); i++) vertex_numbers[vertex_numbers_inv[i]] = i;
  return vertex_numbers; 
}

// takes a triangulation, and returns a dual of the inverse leapfrog
// this is easy because we just remove a set of faces
// the resulting planar graph is oriented because the input is oriented und we only remove vertices
PlanarGraph TriangulationView::inverse_leapfrog_dual() const
{
  assert(is_consistently_oriented());
  PlanarGraph PG(*this);
  set<int> face_vertices, to_do_set;

  // find all vertices with degree < 6 (one could additionally find the vertices with odd degree)
  for(int v=0; v<size(); v++){
    if(degree(v) < 6){
      face_vertices.insert(v);
      to_do_set.insert(v);
    }
  }

  // ... find all face vertices
  while(to_do_set.size() != 0){
    const int u = *(to_do_set.begin());
    to_do_set.erase(to_do_set.begin());

    for(int v: nbrs(u)){
      const int w = next(u,v);
      const int s = face_vertices.size();
      const int x = next(w,v);
      face_vertices.insert(x);
      if(face_vertices.size() != s)
        to_do_set.insert(x);
    }
  }
  // cerr << "faces to remove:" << endl;
  // cerr << face_vertices << endl;

  // check number of face_vertices
  // FIXME, optional

  // remove all face_vertices
  PG.remove_vertices(face_vertices);
  return PG;
}




vector<general_spiral> FullereneDual::isomer_search(const Triangulation::predicate_t& predicate, size_t N, size_t verbose_level,
					    bool IPR, bool only_nontrivial_symmetry, size_t N_chunks, size_t chunk_index)
{
  vector<general_spiral> spirals;
  size_t i=0;
  Graph g;

  // cout << "N_chunks = " << N_chunks << "\n"
  //      << "chunk_index = " << chunk_index << "\n\n";
  
  BuckyGen::buckygen_queue Q = BuckyGen::start(N,IPR,only_nontrivial_symmetry,
					       chunk_index,N_chunks);

  while(BuckyGen::next_fullerene(Q,g)){
    FullereneDual G(g);
    vector<int> rspi;
    jumplist_t jumps;

    i++;
    if(verbose_level>0 && i%verbose_level == 0) fprintf(stderr,"Reached isomer %ld\n",i);

    if(predicate(G)){
      spirals.push_back(FullereneDual(G).get_rspi());
      spiral_nomenclature name(G, spiral_nomenclature::FULLERENE,spiral_nomenclature::TRIANGULATION);
      //      spirals.push_back(name.to_string());
      cerr << name.to_string() << endl;
    }
  }
  BuckyGen::stop(Q);

  if(verbose_level>0) fprintf(stderr,"Generated %ld %d-vertex graphs.\n",i,int(g.size()));
  
  return spirals; 
}
