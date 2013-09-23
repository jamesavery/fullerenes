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
  for(set<edge_t>::const_iterator e(edge_set.begin()); e!=edge_set.end(); e++){
    pair<node_t,node_t> t(adjacent_tris(*e));
    triangles.insert(tri_t(e->first,e->second,t.first));
    triangles.insert(tri_t(e->first,e->second,t.second));
  }
  return vector<tri_t>(triangles.begin(),triangles.end());
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
  orient_triangulation(faces);

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


void wg_connect_backward(set<edge_t> &edge_set, list<pair<node_t, int> > &ov)
{
  list< pair<node_t,int> >::iterator second_last(ov.end());
  --second_last;
  --second_last;

  edge_set.insert(edge_t(ov.back().first, second_last->first));
  --ov.back().second;
  --(second_last->second);//decrement the last but one entry
}

void wg_connect_forward(set<edge_t> &edge_set, list<pair<node_t, int> > &ov)
{
  edge_set.insert(edge_t(ov.back().first, ov.front().first));
  --ov.back().second;
  --ov.front().second;
}

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
}

PlanarGraph Triangulation::dual_graph() const {
  //  assert(is_consistently_oriented());
  IDCounter<tri_t> tri_numbers;
  vector<tri_t> triangles(compute_faces());

  for(int i=0;i<triangles.size();i++) tri_numbers.insert(triangles[i].sorted());

  neighbours_t A(triangles.size(),vector<node_t>(3));

  for(node_t U=0;U<triangles.size();U++){
    const face_t& t(triangles[U]);
    
    for(int i=0;i<3;i++){
      const node_t& u(t[i]), v(t[(i+1)%3]);
      node_t w(nextCW.find(dedge_t(u,v))->second);      

      A[U][i] = tri_numbers(tri_t(u,v,w).sorted());

      if(A[U][i] < 0){
	cerr << "Triangle " << tri_t(u,v,w).sorted() << " (opposite " << t << ") not found!\n";
	abort();
      }
    }
  }
  return PlanarGraph(Graph(A));
};
