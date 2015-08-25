#include "triangulation.hh"
#include <algorithm>

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
    //    printf("%d/%d: %d->%d->%d\n",i,t,u,v,w);
    t++;
    if(t == 1)  tris.first  = w;
    else if(t==2) tris.second = w;
    else {
      fprintf(stderr,"Triangulation is not orientable, edge %d--%d part of more than two faces.\n",u,v);
      cerr << "neighbours = " << neighbours << ";\n";
      abort();
    }
      }
  }
  return tris;
}

vector<tri_t> Triangulation::compute_faces() const
// Does not assume graph is oriented
// Produces oriented triangles
{
  set<tri_t> triangle_set;

  for(node_t u=0;u<N;u++)
    for(int i=0;i<neighbours[u].size();i++){
      node_t v = neighbours[u][i];
      pair<node_t,node_t> ws(adjacent_tris(edge_t(u,v)));

      triangle_set.insert(tri_t(u,v,ws.first ).sorted());
      triangle_set.insert(tri_t(u,v,ws.second).sorted());
    }


  vector<tri_t> triangles(triangle_set.begin(),triangle_set.end());
  orient_triangulation(triangles);
  return triangles;
}


void Triangulation::orient_neighbours()
{
  map<dedge_t,node_t> next;

  for(int i=0;i<triangles.size();i++){
    const tri_t& t(triangles[i]);
    next[dedge_t(t[0],t[1])] = t[2];
    next[dedge_t(t[1],t[2])] = t[0];
    next[dedge_t(t[2],t[0])] = t[1];
  }
  for(node_t u=0;u<N;u++){
    int d = neighbours[u].size();

    node_t v = neighbours[u][0];
    for(int i=1;i<d;i++){
      node_t w = next[dedge_t(u,v)];
      neighbours[u][i] = w;
      v = w;
    }
  }
}

node_t Triangulation::nextCW(node_t u, node_t v) const
{
  return next(u,v);
}

node_t Triangulation::nextCCW(node_t u, node_t v) const
{
  return prev(u,v);
}


vector<tri_t> Triangulation::compute_faces_oriented() const
{
  vector<tri_t> triangles(2*(N-2)); // Assumes cubic
  map<dedge_t,bool> dedge_done;    // Change to vector<bool> dedge_done[N*3] to increase speed.

  int k=0;
  for(node_t u=0;u<N;u++){
    const vector<node_t>& nu(neighbours[u]);
    for(int i=0;i<nu.size();i++){
      const node_t& v(nu[i]);    // Process directed edge u->v
      const dedge_t uv(u,v);

      if(!dedge_done[uv]){
        node_t w = nextCW(uv);

        if(!dedge_done[dedge_t(v,w)] && !dedge_done[dedge_t(w,u)]){

          triangles[k++] = tri_t(u,v,w);

          dedge_done[dedge_t(u,v)] = true;
          dedge_done[dedge_t(v,w)] = true;
          dedge_done[dedge_t(w,u)] = true;
        }
      }
    }
  }
  return triangles;
}



PlanarGraph Triangulation::dual_graph() const
{
  IDCounter<tri_t> tri_numbers;

  for(int i=0;i<triangles.size();i++) tri_numbers.insert(triangles[i].sorted());

  neighbours_t A(triangles.size(),vector<node_t>(3));

  for(node_t U=0;U<triangles.size();U++){
    const tri_t& t(triangles[U]);

    for(int i=0;i<3;i++){
      const node_t& u(t[i]), v(t[(i+1)%3]);
      node_t w(nextCCW(u,v)); // TODO: CCW for buckygen -- will this give problems elsewhere?


      A[U][i] = tri_numbers(tri_t(u,v,w).sorted());

      //      printf("A[%d][%d] = %d (%s)\n",U,i,tri_numbers(tri_t(u,v,w).sorted()),to_string(tri_t(u,v,w).sorted()).c_str());

      if(A[U][i] < 0){
	cerr << "Triangle " << tri_t(u,v,w).sorted() << " (opposite " << t << ") not found!\n";
	abort();
      }
    }
  }
  return PlanarGraph(Graph(A));
};

vector<face_t> Triangulation::dual_faces() const
{
  vector<face_t> dfaces(N);

  IDCounter<tri_t> tri_numbers;
  for(int i=0;i<triangles.size();i++) tri_numbers.insert(triangles[i].sorted());

  for(node_t u=0;u<N;u++){
    const vector<node_t> &nu(neighbours[u]);
    face_t f(nu.size());
    for(int i=0;i<nu.size();i++){
      node_t v=nu[i], w = nextCW(dedge_t(u,v));
      f[i] = tri_numbers(tri_t(u,v,w).sorted());
    }
    dfaces[u] = f;
  }
  return dfaces;
}


// Takes full spiral string, e.g. 566764366348665
// where the degrees are between 3 and 8 (or anything larger, really)
Triangulation::Triangulation(const vector<int>& spiral_string, const jumplist_t& j): PlanarGraph()
{
  jumplist_t jumps = j; // we need a local copy to remove elements
  N = spiral_string.size();

  set<edge_t> edge_set;

  // open_valencies is a list with one entry per node that has been added to
  // the spiral but is not fully saturated yet.  The entry contains the number
  // of the node and the number of open valencies
  list<pair<node_t,int> > open_valencies;

  // set up first two nodes
  open_valencies.push_back(make_pair(0,spiral_string[0]-1));
  open_valencies.push_back(make_pair(1,spiral_string[1]-1));
  edge_set.insert(edge_t(0, 1));

  //iterate over atoms
  //k=0, k=1 have been done already
  // omit the last one because it requires special treatment
  for (int k=2; k<N-1; ++k){
    int pre_used_valencies=0;
//    cout << "k: " << k << endl;

    // should a cyclic shift be applied before adding the next atom?
    if(jumps.size() != 0 && k == jumps.front().first){
      // perform cyclic shift on open_valencies
      for(int i = jumps.front().second; i>0; --i){ // 0 is no jump
        open_valencies.push_back(open_valencies.front());
        open_valencies.pop_front();
      }
      jumps.pop_front();
    }

    // connect k and <last>
    auto connect_backward = [&](){
      edge_set.insert(edge_t(k, open_valencies.back().first));
      --open_valencies.back().second;
      ++pre_used_valencies;
    };

    // connect k and <first>
    auto connect_forward = [&](){
      edge_set.insert(edge_t(k, open_valencies.front().first));
      --open_valencies.front().second;
      ++pre_used_valencies;
    };

    connect_backward();
    connect_forward();

    // do the remaining connect forwards
    while(open_valencies.front().second==0){
      open_valencies.pop_front();
      connect_forward();
    }

    // do the remaining connect backwards
    while(open_valencies.back().second==0){
      open_valencies.pop_back();
      connect_backward();
    }

    if(spiral_string[k] - pre_used_valencies < 1){//the current atom is saturated (which may only happen for the last one)
      cout << "Cage closed but faces left (or otherwise invalid spiral)" << endl;
      abort();
    }

    // add node to spiral
    open_valencies.push_back(make_pair(k,spiral_string[k]-pre_used_valencies));

  } // iterate over atoms

  // make sure we left the spiral in a sane state
  // open_valencies must be either spiral.back() times '1' at this stage
  if(open_valencies.size() != spiral_string.back()){
    cout << "Cage not closed but no faces left (or otherwise invalid spiral), wrong number of faces left" << endl;
    cout << "Incomplete graph g = " << Graph(edge_set) << "\n";
    abort();
  }
  for(list<pair<node_t,int> >::iterator it = open_valencies.begin(); it!=open_valencies.end(); ++it){
    if(it->second!=1){
      cout << "Cage not closed but no faces left (or otherwise invalid spiral), more than one valency left for at least one face" << endl;
    abort();
    }
  }

  // add remaining edges, we don't care about the valency list at this stage
  for(int i=0; i<spiral_string.back(); ++i){
    edge_set.insert(edge_t(N-1, open_valencies.front().first));
    open_valencies.pop_front();
  }

  *this = Triangulation(PlanarGraph(edge_set));
}


//FIXME remove / replace layout assertion ?
Triangulation Triangulation::GCtransform(const unsigned k, const unsigned l) const
{
  assert(layout2d.size() == N);
  Unfolding u(*this,true);
  Unfolding gcu(u*Eisenstein(k,l));
  Folding gcf(gcu);
  Triangulation t(gcf.fold());
  t.layout2d = t.tutte_layout();
  return t;
}



// *********************************************************************
//                 SPIRAL STUFF
// *********************************************************************


bool Triangulation::get_spiral(const node_t f1, const node_t f2, const node_t f3,
			       vector<int> &spiral, jumplist_t& jumps, vector<node_t>& permutation,
			       const bool general) const {
  bool normal_spiral = get_spiral_implementation(f1,f2,f3,spiral,jumps,permutation,false);

  return normal_spiral || (general && get_spiral_implementation(f1,f2,f3,spiral,jumps,permutation,true));
}

void remove_node(const node_t u, Graph &remaining_graph, set<node_t> &remaining_nodes, vector<node_t> &deleted_neighbours){
  remaining_nodes.erase(u);	// O(log(N)) with big coefficient - is set<node_t> the best data structure to use?
  vector<node_t>& nu(remaining_graph.neighbours[u]);

  //remove u from all neighbour lists and erase all neighbours from the u-list
  for(int i=0;i<nu.size();i++){	// O(1) since neighbour count is bounded by max degree
    const node_t& v = nu[i];
    vector<node_t>& nv(remaining_graph.neighbours[v]);

    for(int j=0;j<nv.size();j++){ // O(1) since neighbour count is bounded by max degree
      if(nv[j] == u){
        nv[j] = nv[nv.size()-1];//shift the last entry to the deleted pos
        nv.pop_back();//delete the last
        break;
      }
    }
  }
  deleted_neighbours = nu;
  nu.clear();
}

// jumps start to count at 0
// perform a general spiral search and return the spiral and the jump positions + their length
// TODO: Add jumps to S0.
// TODO: Make GSpiral data type
// TODO: if S0 is given, no need to test for connectedness at every step - jump positions are predetermined.
// TODO: return GSpiral
bool Triangulation::get_spiral_implementation(const node_t f1, const node_t f2, const node_t f3, vector<int> &spiral,
        				      jumplist_t& jumps, vector<node_t> &permutation,
        				      const bool general, const vector<int>& S0) const {
  //this routine expects empty containers pentagon_indices and jumps.  we make sure they *are* empty
  spiral.clear();
  jumps.clear();
  permutation.clear();

  // TODO: Static allocation everywhere. Simplify function.
  permutation.resize(N);
  spiral.resize(N);

  Graph remaining_graph(neighbours); // remaining_graph consists of all nodes that haven't been added to the result yet
  set<node_t> remaining_nodes; // all the nodes that haven't been added yet, not ordered and starting at 0
  vector<int> valencies(N, 0); // valencies is the N-tuple consisting of the valencies for each node

  // open_valencies is a list with one entry per node that has been added to
  // the spiral but is not fully saturated yet.  The entry contains the number
  // of the node and the number of open valencies
  list<pair<node_t,int> > open_valencies;

  // a backup of the neighbours of the current node ... required in case of a jump
  vector<int> deleted_neighbours_bak;

  int pre_used_valencies=0; // number of valencies removed from the last vertex *before* it is added to the spiral
  int jump_state=0; // the number of cyclic shifts of length 1 in the current series.

  //init of the valency-list and the set of nodes in the remaining graph
  for(node_t u=0; u<remaining_graph.N; u++){
    valencies[u] = remaining_graph.neighbours[u].size();
    remaining_nodes.insert(u);
  }

  bool
    CW  = nextCW(dedge_t(f1,f2)) == f3,
    CCW = nextCCW(dedge_t(f1,f2)) == f3;

  // check if starting nodes are a valid spiral start.
  // TODO: Input f1,f2 and bool CW instead to only allow valid spiral starts.
    if(!(CW || CCW)){
    // TODO: Set a global verbosity level; Usually we don't want to look at all this stuff.
    // cerr << "The requested nodes are not connected." << endl;
    return false;
  }

  auto connect_forward = [&](){
    --open_valencies.front().second;
    ++pre_used_valencies;
  };
  auto connect_backward = [&](){
    --open_valencies.back().second;
    ++pre_used_valencies;
  };

  // add the first three (defining) nodes
  // first node
  spiral[0] = valencies[f1];
  permutation[0] = f1;
  remove_node(f1, remaining_graph, remaining_nodes, deleted_neighbours_bak);
  open_valencies.push_back(make_pair(f1,valencies[f1]));

  // second node
  spiral[1] = valencies[f2];
  permutation[1] = f2;
  remove_node(f2, remaining_graph, remaining_nodes, deleted_neighbours_bak);
  connect_backward();
  open_valencies.push_back(make_pair(f2,valencies[f2]-1));

  // third node
  spiral[2] = valencies[f3];
  permutation[2] = f3;
  remove_node(f3, remaining_graph, remaining_nodes, deleted_neighbours_bak);
  connect_backward();
  connect_forward();
  open_valencies.push_back(make_pair(f3,valencies[f3]-2));

  if(!S0.empty() && (spiral[0] != S0[0] || spiral[1] != S0[1] || spiral[2] != S0[2])) return false;

  // iterate over all nodes (of the initial graph) but not by their respective number
  // starting at 3 because we added 3 already
  for(int i=3; i<N-1; ++i){
    pre_used_valencies=0;
    //    list<pair<int,int> > open_valencies_bak(open_valencies); // Makes the whole thing O(N^2), but is never used!

    // find *the* node in *this (not the remaining_graph), that is connected to open_valencies.back() und open_valencies.front()
    // we can't search in the remaining_graph because there are some edges deleted already
    node_t u = open_valencies.back().first, w = open_valencies.front().first;
    node_t v = CCW? nextCW(dedge_t(u,w)) : nextCCW(dedge_t(u,w));
    //    cout << "u->v: " << u << " -> " << v << endl;
    if(v == -1) return false; // non-general spiral failed

    //remove all edges of which *j is part from the remaining graph
    remove_node(v, remaining_graph, remaining_nodes, deleted_neighbours_bak);

    if(general){
      bool is_connected = remaining_graph.is_connected(remaining_nodes);

      if(!is_connected){//further cyclic rotation required
        //revert the last operations
        //	cout << "Reverting deleted node " << v << " to neighbours " << deleted_neighbours_bak << endl;
        remaining_graph.neighbours[v] = deleted_neighbours_bak;
        deleted_neighbours_bak.clear();
        for(vector<node_t>::iterator it = remaining_graph.neighbours[v].begin(); it != remaining_graph.neighbours[v].end(); ++it){
          remaining_graph.neighbours[*it].push_back(v);
        }

        //perform cyclic shift on open_valencies
        open_valencies.push_back(open_valencies.front());
        open_valencies.pop_front();
        //	cout << "open valencies = " << open_valencies << endl;
        //there was no atom added, so 'i' must not be incremented
        --i;
        ++jump_state;
        continue;
      } else if(jump_state!=0){//end of cyclic rotation
        //	cout << "//end of cyclic rotation\n";
        jumps.push_back(make_pair(i,jump_state));
        jump_state=0;
      }
    }

    connect_forward();
    while (open_valencies.front().second==0){
      open_valencies.pop_front();
      connect_forward();
    }

    connect_backward();
    while (open_valencies.back().second==0){
      open_valencies.pop_back();
      connect_backward();
    }

    spiral[i] = valencies[v];
    permutation[i] = v;
    if(!S0.empty() && spiral[i] != S0[i]) return false;

    if(valencies[v]-pre_used_valencies < 1) return false;
    open_valencies.push_back(make_pair(v,valencies[v]-pre_used_valencies));
    if(open_valencies.back().second < 1) return false; //i.e., the spiral is stuck. This can only happen if the spiral missed a jump
  }

  // make sure we left the loop in a sane state
  // this probably requires some proper error handling: throw and catch and so on ...
  if(remaining_nodes.size() != 1){
    // cerr << "more than one node left ... exiting." << endl;
    return false;
  }
  const int last_valency = valencies[*remaining_nodes.begin()];
  if(open_valencies.size() != last_valency){
    // cerr << "wrong number of nodes with open valencies: " << open_valencies.size() << " ... exiting." << endl;
    return false;
  }
  for(list<pair<node_t,int> >::const_iterator it=open_valencies.begin(); it!=open_valencies.end(); ++it){
    if(it->second != 1){
      // cerr << "number of open valencies is not 1 (but it should be) ... exiting." << endl;
      return false;
    }
  }

  spiral[N-1] = last_valency;
  permutation[N-1] = *remaining_nodes.begin();
  return true;
}

void Triangulation::get_all_spirals(vector< vector<int> >& spirals, vector<jumplist_t>& jumps, // TODO: Should only need to supply jumps when general=true
                     vector< vector<int> >& permutations,
                     const bool only_special, const bool general) const
{
  vector<node_t> node_starts;

  vector<int> spiral(N);
  jumplist_t  jump;
  vector<int> permutation;

  // Prefer special nodes. TODO: Automatic renumber in order of degrees.
  for(node_t u=0;u<N;u++) if(neighbours[u].size() != 6) node_starts.push_back(u);
  for(node_t u=0;u<N;u++)
    if(!only_special && neighbours[u].size() == 6) node_starts.push_back(u);

  for(int i=0; i<node_starts.size(); i++){ // Looks like O(N^3), is O(N)
    const node_t u=node_starts[i];
    const vector<node_t>& nu(neighbours[u]);

    for(int j=0;j<nu.size();j++){
      node_t v=nu[j], w[2];
      w[0] = nextCW(dedge_t(u,v));
      w[1] = nextCCW(dedge_t(u,v));

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

// perform the canonical general spiral search and the spiral and the jump positions + their length
// special_only is a switch to search for spirals starting at non-hexagons only
bool Triangulation::get_spiral(vector<int> &spiral, jumplist_t &jumps, const bool canonical, const bool only_special, const bool general) const
{
  vector<node_t> node_starts;

  for(node_t u=0;u<N;u++)
    if(neighbours[u].size() != 6) node_starts.push_back(u);

  for(node_t u=0;u<N;u++)
    if(!only_special && neighbours[u].size() == 6) node_starts.push_back(u);

  vector<int> spiral_tmp,permutation_tmp;
  jumplist_t jumps_tmp;
  spiral = vector<int>(1,INT_MAX); // so it gets overwritten
  jumps = jumplist_t(100,make_pair(0,0)); // so it gets overwritten

  for(int i=0; i<node_starts.size(); i++){
    const node_t u=node_starts[i];
    const vector<node_t>& nu(neighbours[u]);

    for(int j=0;j<nu.size();j++){
      node_t v=nu[j], w[2];
      w[0] = nextCW(dedge_t(u,v));
      w[1] = nextCCW(dedge_t(u,v));

      for(int k=0;k<2;k++){        // Looks like O(N^3), is O(N) (or O(1) if only_special is set)
        if(!get_spiral(u,v,w[k],spiral_tmp,jumps_tmp,permutation_tmp,general))
          continue;

        // + If we don't need the canonical spiral, just return the first one that works
        if(!canonical){
          jumps  = jumps_tmp;
          spiral = spiral_tmp;
          return true;
        }

        // store the shortest / lexicographically smallest (general) spiral
        if(jumps_tmp.size() < jumps.size() ||
           (jumps_tmp.size() == jumps.size() && lexicographical_compare(jumps_tmp.begin(), jumps_tmp.end(), jumps.begin(), jumps.end())) ||
           (jumps_tmp.size() == jumps.size() && jumps_tmp == jumps &&
            lexicographical_compare(spiral_tmp.begin(), spiral_tmp.end(), spiral.begin(), spiral.end()))){
          jumps = jumps_tmp;
          spiral = spiral_tmp;
        }
      }
    }
  }

  if(spiral.size()!=N) return false;

  return true;
}

void Triangulation::symmetry_information(int N_generators, Graph& coxeter_diagram, vector<int>& coxeter_labels) const
{
  vector< vector<int> > spirals, permutations;
  vector<jumplist_t>    jumps;

  get_all_spirals(spirals,jumps,permutations,true,false);

  if(spirals.empty())
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
node_t Triangulation::end_of_the_line(node_t u0, int i, int a, int b) const
{
  node_t q,r,s,t;		// Current square

  auto go_north = [&](){
    const node_t S(s), T(t); // From old square
    q = S; r = T; s = next(S,T); t = next(s,r);
  };

  auto go_east = [&](){
    const node_t R(r), T(t); // From old square
    q = R; s = T; r = next(s,q); t = next(s,r);
  };

  // Square one
  q = u0; 			// (0,0)
  r = neighbours[u0][i];	// (1,0)
  s = next(q,r);	// (0,1)
  t = next(s,r);	// (1,1)

  if(a==1 && b==1) return t;

  vector<int> runlengths = draw_path(max(a,b), min(a,b));

  for(int i=0;i<runlengths.size();i++){
    int L = runlengths[i];

    if(a>=b){			// a is major axis
      for(int j=0;j<L-1;j++)    go_east();
      if(i+1<runlengths.size()) go_north();
    } else {			// b is major axis
      for(int j=0;j<L-1;j++)    go_north();

      if(i+1<runlengths.size()) go_east();
    }
  }

  return t;			// End node is upper right corner.
}

matrix<int> Triangulation::convex_square_surface_distances() const
{
  matrix<int> H = matrix<int>(N,N,all_pairs_shortest_paths());
  int M = *max_element(H.begin(),H.end());      // M is upper bound to path length

  for(int i=0;i<H.size();i++) H[i] *= H[i]; // Work with square distances, so that all distances are integers.

  for(node_t u=0;u<N;u++)
    for(int i=0;i<neighbours[u].size();i++){

      // Note: All Eisenstein numbers of the form (a,0) or (0,b) yield same lengths
      //       as graph distance, and are hence covered by initial step. So start from 1.
      //       M is upper bound for distance, so only need to do a^2+ab+b^2 strictly less than M.
      for(int a=1; a<M;a++)
        for(int b=1; a*a + a*b + b*b < M*M; b++){
          // Check: if(gcd(a,b) != 1) continue.
          const node_t v = end_of_the_line(u,i,a,b);

          // printf("min(H(%d,%d),|(%d,%d)|^2)  = min(%d,%d)\n",
          // 	 u,v,a,b,H(u,v), a*a+a*b+b*b);
          H(u,v) = min(H(u,v), a*a + a*b + b*b);
        }
    }
  return H;
}

matrix<double> Triangulation::surface_distances() const
{
  matrix<double> H(convex_square_surface_distances());
  for(int i=0;i<N*N;i++) H[i] = sqrt(H[i]);

  bool nonconvex = false;
  for(node_t u=0;u<N;u++) if(neighbours[u].size() > 6) nonconvex = true;

  if(nonconvex) return H.APSP();
  else return H;
}

Triangulation Triangulation::sort_nodes() const
{
  vector< pair<int,int> > degrees(N);

  for(int u=0;u<N;u++) degrees[u] = make_pair(neighbours[u].size(), u);

  sort(degrees.begin(), degrees.end());

  cout << "degrees = " << degrees << endl;

  vector<int> newname(N);
  for(int u=0;u<N;u++) newname[degrees[u].second] = u;

  neighbours_t new_neighbours(N);
  for(int u=0;u<N;u++)
    for(int i=0;i<neighbours[u].size();i++)
      new_neighbours[newname[u]].push_back(newname[neighbours[u][i]]);

  return Triangulation(new_neighbours);
}

// call for the canonical general spiral and extract the pentagon indices
bool FullereneDual::get_rspi(const node_t f1, const node_t f2, const node_t f3, vector<int>& rspi, jumplist_t& jumps, const bool general) const
{
  rspi.clear();
  jumps.clear();
  vector<int> spiral;
  vector<node_t> permutation;
  if(!get_spiral(f1, f2, f3, spiral, jumps, permutation, general)) return false;

  int i=0;
  vector<int>::const_iterator it=spiral.begin();
  for(; it!=spiral.end(); ++it,++i){
    if(*it==5){
      rspi.push_back(i);
    }
  }
  //assert(rspi.size()==12);
  return true;
}

// call for the canonical general spiral and extract the pentagon indices
bool FullereneDual::get_rspi(vector<int>& rspi, jumplist_t& jumps, bool canonical, bool general) const
{
  rspi.clear();
  jumps.clear();
  vector<int> spiral;
  bool special_only = false;
  if(!get_spiral(spiral, jumps, canonical, special_only, general)) return false;

  int i=0;
  vector<int>::const_iterator it=spiral.begin();
  for(; it!=spiral.end(); ++it,++i){
    if(*it==5){
      rspi.push_back(i);
    }
  }
  //assert(rspi.size()==12);
  return true;
}

