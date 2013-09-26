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
	//	printf("%d/%d: %d->%d->%d\n",i,t,u,v,w);
	t++;
	if(t == 1)  tris.first  = w;
	else if(t==2) tris.second = w;
	else {
	  fprintf(stderr,"Triangulation is not orientable, edge %d--%d part of more than two faces.\n",u,v);
	  cerr << "edges = " << edge_set << ";\n";
	  abort();
	}
      }
  }
  return tris;
}

vector<tri_t> Triangulation::compute_faces() const // non-oriented triangles
{
  set<tri_t> triangles;

  int i=0;
  for(set<edge_t>::const_iterator e(edge_set.begin()); e!=edge_set.end(); e++,i++){
    pair<node_t,node_t> t(adjacent_tris(*e));
    triangles.insert(tri_t(e->first,e->second,t.first).sorted());
    triangles.insert(tri_t(e->first,e->second,t.second).sorted());
  }
  return vector<tri_t>(triangles.begin(),triangles.end());
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

node_t Triangulation::nextCW(const dedge_t& uv) const
{
  const node_t &u(uv.first), &v(uv.second);
  const vector<node_t>& nv(neighbours[v]);

  for(int j=0;j<nv.size(); j++) if(nv[j] == u) return nv[(j+1)%nv.size()];

  return -1;			// u-v is not an edge in a triangulation
}

node_t Triangulation::nextCCW(const dedge_t& uv) const
{
  const node_t &u(uv.first), &v(uv.second);
  const vector<node_t>& nv(neighbours[v]);

  for(int j=0;j<nv.size(); j++) if(nv[j] == u) return nv[(j-1+nv.size())%nv.size()];

  return -1;			// u-v is not an edge in a triangulation
}


vector<tri_t> Triangulation::compute_faces_oriented() const 
{
  vector<tri_t> triangles(edge_set.size()-N+2);
  map<dedge_t,bool> dedge_done;	// Change to vector<bool> dedge_done[N*3] to increase speed.

  int k=0;
  for(node_t u=0;u<N;u++){
    const vector<node_t>& nu(neighbours[u]);
    for(int i=0;i<nu.size();i++){ 
      const node_t& v(nu[i]);	// Process directed edge u->v
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
      node_t w(nextCCW(dedge_t(u,v)));      

      A[U][i] = tri_numbers(tri_t(u,v,w).sorted());

      if(A[U][i] < 0){
	cerr << "Triangle " << tri_t(u,v,w).sorted() << " (opposite " << t << ") not found!\n";
	abort();
      }
    }
  }
  return PlanarGraph(Graph(A));
};


extern void wg_connect_backward(set<edge_t> &edge_set, list<pair<node_t, int> > &ov);
extern void wg_connect_forward(set<edge_t> &edge_set, list<pair<node_t, int> > &ov);

// // debug only (do not remove, please [lukas])
// void pdp(list<pair<int,int> > &open_valencies){
//   for(list<pair<int,int> >::iterator it(open_valencies.begin()); it!= open_valencies.end(); ++it){
//      cout << it->first << ": " << it->second << endl;
//   }
// }


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
  list<pair<int,int> > open_valencies;

  // set up first two nodes
  open_valencies.push_back(make_pair(0,spiral_string[0]));
  open_valencies.push_back(make_pair(1,spiral_string[1]));
  //connect first two faces
  wg_connect_backward(edge_set, open_valencies);

  //iterate over atoms
  //k=0, k=1 have been done already
  // omet the last one because it requires special treatment
  for (int k=2; k<N-1; ++k){
//    cout << "k: " << k << endl;

    if(jumps.size() != 0 && k == jumps.front().first){
      // perform cyclic shift on open_valencies
      for(int i = jumps.front().second; i>0; --i){ // 0 is no jump
        open_valencies.push_back(open_valencies.front());
        open_valencies.pop_front();
      }
      jumps.pop_front();
    }

    // add node to spiral
    open_valencies.push_back(make_pair(k,spiral_string[k]));

    // connect k to k-1
    wg_connect_backward(edge_set, open_valencies);

    // connect k to k-2, etc
    wg_connect_forward(edge_set, open_valencies);

    // do the remaining connect forwards
    while(open_valencies.front().second==0){
      open_valencies.pop_front();
      wg_connect_forward(edge_set, open_valencies);
    }
    // do the remaining connect backwards //not neat but the most simple way to emulate 'while second_last->second==0) ...'
    while(true){
      list<pair<int,int> >::iterator second_last(open_valencies.end());
      --second_last;
      --second_last;
      if(second_last->second==0){
        open_valencies.erase(second_last);
        wg_connect_backward(edge_set, open_valencies);
      } else break;
      
    }
//    pdp(open_valencies);

    if (open_valencies.back().second == 0){//the current atom is saturated (which may only happen for the last one)
      cout << "Cage closed but faces left (or otherwise invalid spiral)" << endl;
      abort();
    }
   
  }//iterate over atoms

  // make sure we left the spiral in a sane state
  // open_valencies must be either spiral.back() times '1' at this stage
  if(open_valencies.size() != spiral_string.back()){
    cout << "Cage not closed but no faces left (or otherwise invalid spiral), wrong number of faces left" << endl;
    abort();
  }
  for(list<pair<int,int> >::iterator it = open_valencies.begin(); it!=open_valencies.end(); ++it){
    if(it->second!=1){
      cout << "Cage not closed but no faces left (or otherwise invalid spiral), more than one valency left for at least one face" << endl;
    abort();
    }
  }

  // add last node to spiral // the name of the last node is N -1 (because it's the last one)
  open_valencies.push_back(make_pair(N-1,spiral_string[N-1]));

  for(int i=0; i<spiral_string.back(); ++i){
    wg_connect_forward(edge_set, open_valencies);
    open_valencies.pop_front();
  }

  *this = Triangulation(PlanarGraph(edge_set));
//  FIXME is there any update_from ... required?
}


// *********************************************************************
//			     SPIRAL STUFF
// *********************************************************************
// gpi is for 'get pentagon indices'
inline void gpi_connect_forward(list<pair<int,int> > &open_valencies){
  --open_valencies.back().second;
  --open_valencies.front().second;
}

inline void gpi_connect_backward(list<pair<int,int> > &open_valencies){
  list<pair<int,int> >::iterator second_last(open_valencies.end());
  second_last--;
  second_last--;

  --open_valencies.back().second;
  --second_last->second;//decrement the last but one entry
}

inline void gpi_remove_node(const int i, PlanarGraph &remaining_graph, set<int> &remaining_nodes, vector<int> &deleted_neighbours){
  remaining_nodes.erase(i);
  //remove i from all neighbour lists and erase all neighbours from the i-list
  for(vector<int>::iterator it = remaining_graph.neighbours[i].begin(); it != remaining_graph.neighbours[i].end(); ++it){
    remaining_graph.neighbours[*it].erase(find(remaining_graph.neighbours[*it].begin(),remaining_graph.neighbours[*it].end(),i));
  }
  deleted_neighbours = remaining_graph.neighbours[i];
  remaining_graph.neighbours[i].clear();
}

// jumps start to count at 0
// perform a general spiral search and return the spiral and the jump positions + their length
void Triangulation::get_spiral(const node_t f1, const node_t f2, const node_t f3, vector<int> &spiral, jumplist_t& jumps, bool general) const {

  //this routine expects empty containers pentagon_indices and jumps.  we make sure they *are* empty
  spiral.clear();
  jumps.clear();

  // remaining_graph is the graph that consists of all nodes that haven't been
  // added to the graph yet
  Triangulation remaining_graph(*this);
  // all the nodes that haven't been added yet, not ordered and starting at 0
  set<node_t> remaining_nodes;

  // valencies is a list of length N and contains the valencies of each node
  vector<node_t> valencies(N, 0);  
  // open_valencies is a list with one entry per node that has been added to
  // the spiral but is not fully saturated yet.  The entry contains the number
  // of the node and the number of open valencies
  list<pair<node_t,int> > open_valencies;
  // a backup of the neighbours of he current node ... required in case of a
  // jump
  vector<int> deleted_neighbours_bak;

  //the current jumping state
  int x=0;

  //init of the valency-list and the set of nodes in the remaining graph
  for(int i=0; i!=remaining_graph.N; ++i){
    valencies[i] = remaining_graph.neighbours[i].size();
    //cout << i << ": " << valencies[i]<< endl;
    remaining_nodes.insert(i);
  }

  //check if starting nodes share a face
  if(edge_set.find(edge_t(f1,f2)) == edge_set.end() ||
     edge_set.find(edge_t(f1,f3)) == edge_set.end() ||
     edge_set.find(edge_t(f2,f3)) == edge_set.end()){
    cerr << "The requested nodes are not connected.  Aborting ..." << endl;
    abort();
  }

  // add the first three (defining) nodes
  // first node
  spiral.push_back(valencies[f1]);
  gpi_remove_node(f1, remaining_graph, remaining_nodes, deleted_neighbours_bak);
  open_valencies.push_back(make_pair(f1,valencies[f1]));

  // second node
  spiral.push_back(valencies[f2]);
  gpi_remove_node(f2, remaining_graph, remaining_nodes, deleted_neighbours_bak);
  open_valencies.push_back(make_pair(f2,valencies[f2]));
  gpi_connect_backward(open_valencies);

  // third node
  spiral.push_back(valencies[f3]);
  gpi_remove_node(f3, remaining_graph, remaining_nodes, deleted_neighbours_bak);
  open_valencies.push_back(make_pair(f3,valencies[f3]));
  gpi_connect_backward(open_valencies);
  gpi_connect_forward(open_valencies);

  // iterate over all nodes (of the initial graph) but not by their respective number
  // starting at 3 because we added 3 already
  for(int i=3; i<N-1; ++i){

    list<pair<int,int> > open_valencies_bak(open_valencies);

    // find *the* node in *this (not the remaining_graph), that is connected to open_valencies.back() und open_valencies.front()
    // we can't search in the remaining_graph because there are some edges deleted already
    set<int>::iterator j=remaining_nodes.begin();
    node_t u = open_valencies.back().first, w = open_valencies.front().first;
    for( ; j!=remaining_nodes.end(); ++j){
      if(edge_set.find(edge_t(u,*j)) != edge_set.end() &&
         edge_set.find(edge_t(w,*j)) != edge_set.end()) break;
    }
    // there is allways a node to be added next
    // even in the non-general spiral fails should be caught in the connectedness test
    assert(j!=remaining_nodes.end());

    spiral.push_back(valencies[*j]);
    open_valencies.push_back(make_pair(*j,valencies[*j]));
    gpi_connect_backward(open_valencies);
    gpi_connect_forward(open_valencies);

    // there are three positions in open_valencies that can be 0---one shouldn't happen, the other two cases require interaction.
    while(open_valencies.front().second==0){
      open_valencies.pop_front();
      gpi_connect_forward(open_valencies);
    }
    while(true){
      list<pair<int,int> >::iterator second_last(open_valencies.end());
      second_last--;
      second_last--;
      
      if(second_last->second==0){
        open_valencies.erase(second_last);
        gpi_connect_backward(open_valencies);
      }
      else break;
    }
//    assert(open_valencies.back().second!=0);//i.e., the spiral is stuck. This can only happen if the spiral missed a jump

    node_t v = *j;
    //remove all edges of which *j is part from the remaining graph
    gpi_remove_node(v, remaining_graph, remaining_nodes, deleted_neighbours_bak);

    bool is_connected = remaining_graph.is_connected(remaining_nodes);
    if(!general && !is_connected){//failing spiral
      spiral.front() = INT_MAX; // as an error code tht behaves correctly with respect to lexicographical sorting
      return;
    }
    else if(general && !is_connected){//further cyclic rotation required
      //revert the last operations
      remaining_nodes.insert(v);
      spiral.pop_back();
      open_valencies = open_valencies_bak;
      remaining_graph.neighbours[v] = deleted_neighbours_bak;
      for(vector<node_t>::iterator it = remaining_graph.neighbours[v].begin(); it != remaining_graph.neighbours[v].end(); ++it){
        remaining_graph.neighbours[*it].push_back(v);
      }
      //perform cyclic shift on open_valencies
      open_valencies.push_back(open_valencies.front());
      open_valencies.pop_front();
      //there was no atom added, so 'i' must not be incremented
      --i;
      ++x;
    }
    else if(general && x!=0 && is_connected){//end of cyclic rotation
      jumps.push_back(make_pair(i,x));
      x=0;
    }
  }

  // make sure we left the loop in a sane state
  // this probably requires some proper error handling thow and catch and so on ...
  assert(remaining_nodes.size() == 1);
  const int last_valency = valencies[*remaining_nodes.begin()];
  assert(open_valencies.size() == last_valency);
  for(list<pair<int,int> >::const_iterator it=open_valencies.begin(); it!=open_valencies.end(); ++it){
    assert(it->second == 1);
  }
  spiral.push_back(last_valency);
 
}


// perform the canonical general spiral search and the spiral and the jump positions + their length
void Triangulation::get_canonical_spiral(vector<int> &spiral, jumplist_t &jumps, bool general) const {

//  vector<int> pentagon_indices_tmp;
  vector<int> spiral_tmp;
  jumplist_t jumps_tmp;
  spiral = vector<int>(1,INT_MAX); // so it gets overwritten
  jumps = jumplist_t(100,make_pair(0,0));
  
  vector<face_t> faces(compute_faces_flat(3));

//  cout << "generating all spirals ";

  for(int i=0; i<faces.size(); i++){
    int permutations[6][3] = {{0,1,2},{0,2,1},{1,0,2},{1,2,0},{2,0,1},{2,1,0}};
    const face_t& f = faces[i];
    for(int j=0; j<6; j++){

      int f1 = f[permutations[j][0]], f2 = f[permutations[j][1]], f3 = f[permutations[j][2]];

      get_spiral(f1, f2, f3, spiral_tmp, jumps_tmp, general);

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

//  cout << "got spiral, size: " << spiral.size() << endl;
//  cout << "got spiral: " << spiral << endl;

//  cout << "got jumps, size: " << jumps.size() << endl;
//  cout << "got jumps: " << jumps << endl;

}


void FullereneDual::get_canonical_fullerene_rspi(vector<int>& rspi, jumplist_t& jumps, bool general) const {

  rspi.clear();
  jumps.clear();
  vector<int> spiral;
  get_canonical_spiral(spiral, jumps, general);

  int i=0;
  vector<int>::const_iterator it=spiral.begin();
  for(; it!=spiral.end(); ++it,++i){
    if(*it==5){
      rspi.push_back(i);
    }
  }
}


