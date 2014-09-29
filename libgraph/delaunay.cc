#include "delaunay.hh"

vector<dedge_t> FulleroidDelaunay::triangulate_hole(const vector<node_t>& hole0)
{
  vector<node_t> hole(hole0);
  vector<dedge_t> new_arcs;
  
  int i=0;
  while(hole.size()>3){
    // Is it safe to connect h[i] to h[i+2]?
    node_t u = hole[i], v = hole[(i+1)%hole.size()], w = hole[(i+2)%hole.size()];
    dedge_t uw(u,v);
    
    if(!edge_exists(uw)){
	// Edge doesn't already exist - but does it lead to a consistent triangle?
	insert_edge(uw,v,v);
	
	// New triangle uvw has already-Delaunay neighbour triangles uvs and vwt. Is uvw screwy?
	node_t s = nextCCW(v,u), t = nextCCW(w,v);

	Quad q1(u,s,v,w), q2(v,t,w,u);
	
	if(!is_consistent(q1) || !is_consistent(q2)){ // New triangle is screwy -- discard edge
	  remove_edge(dedge_t(u,w));
	} else { 			// New triangle is OK. Delaunayify!
	  hole.erase(hole.begin()+((i+1)%hole.size()));       // Remove v from hole
	  
	  // We have possibly removed an element before the current one - find current element again.
	  i = hole.begin() - find(hole.begin(),hole.end(),u); 
	  new_arcs.push_back(dedge_t(u,w));
	}
      }
    i = (i+1)%hole.size();
  }
  return new_arcs;
}


vector<dedge_t> FulleroidDelaunay::delaunayify_hole(const vector<dedge_t>& edges)
{
  vector<dedge_t> new_edges(edges);

  int flips = 0;
  bool done = false;
  while(!done){			// Naive |edges|^2 algorithm
    done = true;
    for(int i=0;i<edges.size();i++){
      node_t A = edges[i].first, C = edges[i].second;
      node_t B = nextCCW(C,A), D = nextCCW(A,C);

      Quad q(A,B,C,D);

      if(!is_delaunay(q) && is_consistent(q.flipped())){ // Do a flip!
	remove_edge(dedge_t(A,C));
	insert_edge(dedge_t(B,D),A,C);
	
	new_edges.erase(find(new_edges.begin(),new_edges.end(), dedge_t(A,C)));
	new_edges.push_back(dedge_t(B,D));
	
	flips++;
	done = false;
      }

    }
  }

  return new_edges;
}

void FulleroidDelaunay::remove_last_vertex()
{
  vector<node_t> hole(neighbours.back());

  node_t h = N-1;
  for(int i=0;i<hole.size();i++) remove_edge(edge_t(h,hole[i]));
  neighbours.pop_back();

  vector<dedge_t> triangle_edges = triangulate_hole(hole);
  vector<dedge_t> delaunay_edges = delaunayify_hole(triangle_edges);
}

void FulleroidDelaunay::remove_flat_vertices()
{
  // TODO: Assume that vertices are sorted such that hexagons appear
  // at the end. This needs a fix for fulleroids with negative
  // curvature.
  while(neighbours.back().size() == 6) remove_last_vertex();
  // TODO: Perform a final Delaunayification.
}
