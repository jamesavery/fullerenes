#include "delaunay.hh"

double FulleroidDelaunay::angle(node_t A, node_t B, node_t C) const {
  double 
    a = distances(B,A),
    b = distances(B,C),
    c = distances(A,C);
  // printf("dist(%d,%d) = %g\n"
  // 	 "dist(%d,%d) = %g\n"
  // 	 "dist(%d,%d) = %g\n",
  // 	 B,A,distances(B,A),
  // 	 B,C,distances(B,C),
  // 	 A,C,distances(A,C));
  //  printf("(a,b,c) = (%g,%g,%g) ~> acos(%g) = %g\n",a,b,c,(a*a+b*b-c*c)/(2*a*b),acos((a*a+b*b-c*c)/(2*a*b)));
  return acos((a*a+b*b-c*c)/(2*a*b));
}


double FulleroidDelaunay::angle(const Quad& Q, int i, int subangle) const {
  int A(Q.v[(i+3-(subangle==2))%4]), B(Q.v[i]), C(Q.v[(i+1+(subangle==1))%4]);
  // printf("angle({%d,%d,%d,%d},%d(%d),%d) = angle(%d,%d,%d)\n",
  // 	 Q.v[0],Q.v[1],Q.v[2],Q.v[3],i,Q.v[i],subangle,
  // 	 A,B,C);

  return angle(A,B,C);
}

bool FulleroidDelaunay::is_consistent(const Quad& Q,int i) const { 
  //  printf("%g =?= %g+%g\n",angle(Q,i,0),angle(Q,i,1),angle(Q,i,2));
  //  return true;
  return fabs(angle(Q,i,0) - angle(Q,i,1) - angle(Q,i,2)) < epsilon; 
}

bool FulleroidDelaunay::is_consistent(const Quad& Q)       const { 
  // printf("{%d,%d,%d,%d}: %d-%d-%d-%d\n",Q.v[0],Q.v[1],Q.v[2],Q.v[3],
  // 	 is_consistent(Q,0), is_consistent(Q,1), is_consistent(Q,2), is_consistent(Q,3));
    return is_consistent(Q,0) && is_consistent(Q,1) && is_consistent(Q,2) && is_consistent(Q,3); 
}


vector<dedge_t> FulleroidDelaunay::triangulate_hole(const vector<node_t>& hole0)
{
  vector<node_t> hole(hole0);
  vector<dedge_t> new_arcs;
  
  int i=0;
  int steps = 0;
  while(hole.size()>3){
    // Is it safe to connect h[i] to h[i+2]?
    node_t u = hole[i], v = hole[(i+1)%hole.size()], w = hole[(i+2)%hole.size()];
    dedge_t uw(u,w);
    
    cout << "Trying edge " << uw << ": ";
    if(!edge_exists(uw)){
	// Edge doesn't already exist - but does it lead to a consistent triangle?
	insert_edge(uw,v,v);
	
	// New triangle uvw has already-Delaunay neighbour triangles uvs and vwt. Is uvw screwy?
	node_t s = nextCCW(v,u), t = nextCCW(w,v);

	Quad q1(u,s,v,w), q2(v,t,w,u);
	
	if(!is_consistent(q1) || !is_consistent(q2)){ // New triangle is screwy -- discard edge
	  cout << "leads to inconsistent neighbourhood " << vector<int>(q1.v,q1.v+4) << "("<< is_consistent(q1) << ")"
	       << " or " << vector<int>(q2.v,q2.v+4)  << "("<< is_consistent(q2) << ")\n";

	  remove_edge(dedge_t(u,w));
	} else { 			// New triangle is OK. Delaunayify!
	  hole.erase(hole.begin()+((i+1)%hole.size()));       // Remove v from hole
	  
	  // We have possibly removed an element before the current one - find current element again.
	  i = hole.begin() - find(hole.begin(),hole.end(),u); 
	  new_arcs.push_back(dedge_t(u,w));
	}
    } else {
      cout << "edge exists.\n";
    }
    i = (i+1)%hole.size();
    if(++steps > 10){
      cout << "Got stuck; graph is:\n"
	   <<"g = " << *this << ";\n"
	   <<"hole = " << hole << ";\n"
	   <<"newarcs = " << new_arcs <<";\n";
      abort();
    }
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
