#include "graph.hh"

vector<unsigned int> Graph::shortest_paths(const node_t& source, const vector<bool>& used_edges, 
					   const vector<bool>& used_nodes, const unsigned int max_depth) const
{
    //    if(used_nodes.size() == 0) used_nodes = vector<bool>(N,false);
    //    if(used_edges.size() == 0) used_edges = vector<bool>(N*(N-1)/2,false);
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
  return shortest_paths(source,vector<bool>(N*(N-1)/2),vector<bool>(N),max_depth);
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

vector<node_t> Graph::shortest_cycle(const node_t& s, const node_t& t, const node_t& r, const int max_depth) const 
{ 
  vector<bool> used_edges(N*(N-1)/2);
  vector<bool> used_nodes(N);

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

Graph Graph::dual_graph(unsigned int Fmax, const vector<coord2d> layout) const {
  Graph dual;
  unsigned int Nfaces = edge_set.size()-N+2;
  dual.N = Nfaces;
  dual.neighbours.resize(Nfaces);
  dual.edges.resize(Nfaces*(Nfaces-1)/2);

  const vector<face_t> allfaces(compute_faces_flat(Fmax,layout));

  if(Nfaces != allfaces.size()){
    fprintf(stderr,"%d != %d faces: Graph is not polyhedral.\n",Nfaces,int(allfaces.size()));
    cout << "errgraph = " << *this << endl;
  }

  // Construct mapping e -> faces containing e (these are mutually adjacent)
  map< edge_t, set<int> > facenodes;
  for(unsigned int i=0;i<allfaces.size(); i++){
    const face_t& face(allfaces[i]);
    //  cerr << "Face "<<i<<": " << face << endl;
    for(unsigned int j=0;j<face.size();j++)
      facenodes[edge_t(face[j],face[(j+1)%face.size()])].insert(i);
  }

  for(map<edge_t,set<int> >::const_iterator fs(facenodes.begin());fs!=facenodes.end();fs++){
    const edge_t&   e(fs->first);
    const set<int>& connects(fs->second);
    if(connects.size() != 2)
      fprintf(stderr,"Edge (%d,%d) connects %d faces: Graph is not planar.\n",e.first,e.second,int(connects.size()));
  }
  
  // Insert edge between each pair of faces that share an edge
  for(set<edge_t>::const_iterator e(edge_set.begin()); e!= edge_set.end(); e++){
    const set<int>& adjacent_faces(facenodes[*e]);
    for(set<int>::const_iterator f(adjacent_faces.begin()); f!= adjacent_faces.end(); f++){
      set<int>::const_iterator g(f);
      for(++g; g!= adjacent_faces.end(); g++)
	dual.edge_set.insert(edge_t(*f,*g));
    }
  }
  fprintf(stderr,"%d nodes, and %d edges in dual graph.\n",int(dual.N), int(dual.edge_set.size()));

  dual.update_auxiliaries();

  // If original graph was planar with 2D layout, there's a corresponding layout for the dual grap
  if(layout.size() == N){
    dual.layout2d = vector<coord2d>(Nfaces);
    for(unsigned int i=0;i<Nfaces;i++){
      face_t face(allfaces[i]);
      coord2d center = 0;
      for(unsigned int j=0;j<face.size();j++) center += layout[face[j]];
      dual.layout2d[i] = center / face.size();
    }
  }

  // for(unsigned int u=0;u<Nfaces;u++){
  //   printf("%d: ",u);
  //   for(int i=0;i<dual.neighbours[u].size();i++)
  //     printf("%d ",dual.neighbours[u][i]);
  //   printf("\n");
  // }  
  return dual;
}


Graph::facemap_t Graph::compute_faces(unsigned int Nmax, const vector<coord2d> layout) const 
{
  facemap_t facemap;

  // TODO: This is a much better and faster method, but needs to be debugged.
  //  if(layout.size() == N) return compute_faces_oriented(layout);

  for(set<edge_t>::const_iterator e(edge_set.begin()); e!= edge_set.end(); e++){
    const node_t s = e->first, t = e->second;

    const vector<node_t>& ns(neighbours[t]);
    for(unsigned int i=0;i<ns.size();i++)
      if(ns[i] != s) {
	const node_t u = ns[i];

	face_t face(shortest_cycle(s,t,u,Nmax));  
	if(face.size() <= Nmax){
	  facemap[face.size()].insert(face);
	} //else {
	  //	  fprintf(stderr,"Erroneous face starting at (%d -> %d -> %d) found: ",s,t,u); 
	  //	  cerr << face << endl;
	  
	//	}
      }
  }
  return facemap;
}


vector<face_t> Graph::compute_faces_flat(unsigned int Nmax, const vector<coord2d> layout) const 
{
  vector<face_t> result;
  facemap_t facemap(compute_faces(Nmax,layout));
  for(facemap_t::const_iterator fs(facemap.begin()); fs != facemap.end(); fs++)
    copy(fs->second.begin(),fs->second.end(),inserter(result,result.end()));

  return result;
}


void Graph::update_auxiliaries() 
{
  // Instantiate auxiliary data strutures: sparse adjacency matrix and edge existence map.
  map<node_t,set<node_t> > ns;
  //  fprintf(stderr,"Initializing edge map.\n");

  // Update node count
  for(set<edge_t>::const_iterator e(edge_set.begin()); e!= edge_set.end(); e++){
    N = max(N,max(e->first,e->second)+1);
  }

  edges.resize(N*(N-1)/2);
  neighbours.resize(N);
  node_t nmax = 0;
  for(set<edge_t>::const_iterator e(edge_set.begin()); e!= edge_set.end(); e++){
    edges[e->index()] = true;
    ns[e->first].insert(e->second);
    ns[e->second].insert(e->first);
  }

  //  fprintf(stderr,"Initializing adjacencies\n");
  for(unsigned int u=0;u<N;u++)
    neighbours[u] = vector<node_t>(ns[u].begin(),ns[u].end());

}

ostream& operator<<(ostream& s, const Graph& g) 
{
  s << g.name<< "Graph[Range[0,"<<(g.N-1)<<"],\n\tUndirectedEdge@@#&/@{";
  for(set<edge_t>::const_iterator e(g.edge_set.begin()); e!=g.edge_set.end(); ){    
    s << "{" << e->first << "," << e->second << "}";
    if(++e != g.edge_set.end())
      s << ", ";
    else
      s << "}";
  }

  if(g.layout2d.size() == g.N){
    s << g.name << ",\n\tVertexCoordinates->{";
    for(unsigned int i=0;i<g.N;i++){
      coord2d xy(g.layout2d[i]);
      s << xy << (i+1<g.N?", ":"}");
    }
  }
  s << "\n]";

  return s;
}
