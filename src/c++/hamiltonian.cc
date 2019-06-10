#include "graph.hh"

int Graph::hamiltonian_count() const 
  {
    node_t start_node = 0;
    vector<bool> used_edges(N*(N-1)/2);
    vector<bool> used_nodes(N);
    vector<node_t> path;
    used_nodes[start_node] = true;
    path.push_back(start_node);

    vector<int> distances(shortest_paths(start_node,used_edges,used_nodes,N));
    
    fprintf(stderr,"distances:\n");
    for(unsigned int i=0;i<N;i++) fprintf(stderr,"\t%d : %d\n",i+1,distances[i]);
    return hamiltonian_count(start_node,used_edges,used_nodes,path,distances);
  }

int Graph::hamiltonian_count(const node_t& current_node, vector<bool>& used_edges, vector<bool>& used_nodes, vector<node_t>& path, const vector<int>& distances) const {
    if(path.size() == N){
      if(edge_set.find(edge_t(current_node,path[0])) != edge_set.end()){ // Hamiltonian circuit
	//	fprintf(stderr,"Hamiltonian: ");
	for(unsigned int i=0;i<path.size();i++) fprintf(stderr,"%d ",path[i]+1); fprintf(stderr,"%d\n",path[0]+1);
	return 1;
      } else {
	// fprintf(stderr,"Non-Hamiltonian: "); for(unsigned int i=0;i<path.size();i++) fprintf(stderr,"%d ",path[i]+1); printf("\n");
	return 0;
      }
    }
  
    int count = 0;
    const vector<node_t> &ns(neighbours[current_node]);

    // Exhaustive depth first search
    for(int i=0;i<3;i++){
      const edge_t edge(current_node,ns[i]);
      if(!used_nodes[ns[i]] && !used_edges[edge.index()]  && distances[ns[i]] <= (N-path.size()) ){
	used_edges[edge.index()] = true;
	used_nodes[ns[i]] = true;
	path.push_back(ns[i]);

	// Periodically update distance vector
	// if(path.size() < N-15)
	//   count += hamiltonian_count(ns[i],used_edges,used_nodes,path, shortest_paths(0,used_edges,used_nodes,N-path.size()));
	// else
	count += hamiltonian_count(ns[i], used_edges, used_nodes, path,distances);
	
	path.pop_back();
	used_nodes[ns[i]] = false;
	used_edges[edge.index()] = false;
      }
    }
    return count;
  }
