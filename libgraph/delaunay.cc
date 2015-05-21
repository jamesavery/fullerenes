#include "delaunay.hh"
#include "debug.hh"

ostream& operator<<(ostream& s, const FulleroidDelaunay::Quad& q) 
{
  s << "{" << q.v[0] << ", " << q.v[1] << ", " << q.v[2] << ", " << q.v[3] << "}";
  return s;
}

double FulleroidDelaunay::angle(node_t A, node_t B, node_t C) const {
// returns angle at B
  double 
    a = distances(B,A),
    b = distances(B,C),
    c = distances(A,C);
  double cos_theta = (a*a+b*b-c*c)/(2*a*b);

  Debug("Delaunay",Debug::INFO3) 
    << "dist(" << make_pair(B,A) << ") = " << distances(B,A) << "; "
    << "dist(" << make_pair(B,C) << ") = " << distances(B,A) << "; "
    << "dist(" << make_pair(A,C) << ") = " << distances(B,A) << "\n"
    << "(a,b,c) = " << vector<double>({a,b,c}) << " ~> acos("<<cos_theta<<") = " << acos(cos_theta) << "\n\n";

  return acos(cos_theta);
}


double FulleroidDelaunay::angle_d6y(node_t A, node_t B, node_t C) const {
// returns angle at B
  double 
    a = edge_lengths_d6y(B,A),
    b = edge_lengths_d6y(B,C),
    c = edge_lengths_d6y(A,C);
  double cos_theta = (a*a+b*b-c*c)/(2*a*b);

  Debug("Delaunay",Debug::INFO3) 
    << "dist(" << make_pair(B,A) << ") = " << distances(B,A) << "; "
    << "dist(" << make_pair(B,C) << ") = " << distances(B,A) << "; "
    << "dist(" << make_pair(A,C) << ") = " << distances(B,A) << "\n";

    assert(a>epsilon && b>epsilon && c>epsilon);
  if( cos_theta >= 1.0 ){ // catch 0 degree case + numerical inaccuracy
    Debug("Delaunay",Debug::INFO2)  
      << "(a,b,c) = " << vector<double>({a,b,c}) 
      << " ~> acos("<<cos_theta<<") = " << acos(cos_theta) << " (exception used)\n\n";

    return 0;
  }
  if( cos_theta <= -1.0 ){ // catch 180 degree case + numerical inaccuracy
    Debug("Delaunay",Debug::INFO2)  
      << "(a,b,c) = " << vector<double>({a,b,c}) 
      << " ~> acos("<<cos_theta<<") = " << acos(cos_theta) << " (exception used)\n\n";

    return M_PI;
  }

    Debug("Delaunay",Debug::INFO3)  
      << "(a,b,c) = " << vector<double>({a,b,c}) 
      << " ~> acos("<<cos_theta<<") = " << acos(cos_theta) << " (exception used)\n\n";

  return acos(cos_theta);
}


double FulleroidDelaunay::angle(const Quad& Q, int i, int subangle) const {
  int A(Q.v[(i+3-(subangle==2))%4]), B(Q.v[i]), C(Q.v[(i+1+(subangle==1))%4]);
  // printf("angle({%d,%d,%d,%d},%d(%d),%d) = angle(%d,%d,%d)\n",
  //          Q.v[0],Q.v[1],Q.v[2],Q.v[3],i,Q.v[i],subangle,
  //          A,B,C);

  return angle(A,B,C);
}


double FulleroidDelaunay::angle_d6y(const Quad& Q, int i, int subangle) const {
  int A(Q.v[(i+3-(subangle==2))%4]), B(Q.v[i]), C(Q.v[(i+1+(subangle==1))%4]);
  // printf("angle({%d,%d,%d,%d},%d(%d),%d) = angle(%d,%d,%d)\n",
  //          Q.v[0],Q.v[1],Q.v[2],Q.v[3],i,Q.v[i],subangle,
  //          A,B,C);

  return angle_d6y(A,B,C);
}

bool FulleroidDelaunay::is_consistent(const Quad& Q,int i) const { 
  //  printf("%g =?= %g+%g\n",angle(Q,i,0),angle(Q,i,1),angle(Q,i,2));
  //  return true;
  return fabs(angle(Q,i,0) - angle(Q,i,1) - angle(Q,i,2)) < epsilon; 
}

bool FulleroidDelaunay::is_consistent(const Quad& Q)       const { 
  // printf("{%d,%d,%d,%d}: %d-%d-%d-%d\n",Q.v[0],Q.v[1],Q.v[2],Q.v[3],
  //          is_consistent(Q,0), is_consistent(Q,1), is_consistent(Q,2), is_consistent(Q,3));
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
     
     cout << "Trying edge " << make_pair(u+1,w+1) << ": ";
     if(!edge_exists(uw)){
       // Edge doesn't already exist - but does it lead to a consistent triangle?
       insert_edge(uw,next(u,v),v);  // u: ..v,w,x,.. w: ...,u,v,... ; 
         
       // New triangle uvw has already-Delaunay neighbour triangles uvs and vwt. Is uvw screwy?
       node_t s = nextCW(v,u), t = nextCW(w,v);
       Quad q1(u,s,v,w), q2(v,t,w,u);
         
       if(!is_consistent(q1) || !is_consistent(q2)){ // New triangle is screwy -- discard edge
         cout << "leads to inconsistent neighbourhood " << (vector<int>(q1.v,q1.v+4)+1) << "("<< is_consistent(q1) << ")"
              << " or " << (vector<int>(q2.v,q2.v+4)+1)  << "("<< is_consistent(q2) << ")\n";
 
         remove_edge(dedge_t(u,w));
       } else {         		// New triangle is OK. Delaunayify!
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
            <<"hole = " << hole+1 << ";\n"
            <<"newarcs = " << new_arcs <<";\n";
       abort();
     }
   }
  return new_arcs;
}


// vector<dedge_t> FulleroidDelaunay::delaunayify_hole(const vector<dedge_t>& edges)
// {
//   vector<dedge_t> new_edges(edges);
// 
//   int flips = 0;
//   bool done = false;
//   while(!done){ // Naive |edges|^2 algorithm
//     done = true;
//     for(int i=0;i<edges.size();i++){
//       node_t A = edges[i].first, C = edges[i].second;
//       node_t B = nextCW(C,A), D = nextCW(A,C);
// 
//       Quad q(A,B,C,D);
// 
//       if(!is_delaunay(q) && is_consistent(q.flipped())){ // Do a flip!
//         remove_edge(dedge_t(A,C));
//         insert_edge(dedge_t(B,D),A,C);
//         
//         new_edges.erase(find(new_edges.begin(),new_edges.end(), dedge_t(A,C)));
//         new_edges.push_back(dedge_t(B,D));
//         
//         flips++;
//         done = false;
//       }
// 
//     }
//   }
// 
//   return new_edges;
// }

void FulleroidDelaunay::delaunayify_hole_2(const vector<edge_t>& edges)
{
  vector<edge_t> new_edges(edges);

  int flips = 0;
  bool done = false;
  while(!done){ // Naive |edges|^2 algorithm
    done = true;
    for(int i=0;i<new_edges.size();i++){
      cout << "edges considered for dealaunayification: " << new_edges << ". current index: " << i << endl;

      node_t A = new_edges[i].first, C = new_edges[i].second;
      node_t B = nextCW(C,A), D = nextCW(A,C);

      Quad q(A,B,C,D);

      if(!is_delaunay_d6y(q)){ // Do a flip!
        cout << q << " is not delaunay -- flipping." << endl;

        // TODO: implement void flip()

        double ab = edge_lengths_d6y(A,B);
        double ad = edge_lengths_d6y(D,A);
        double alpha1 = angle_d6y(D,A,C);
        double alpha2 = angle_d6y(C,A,B);
        double bd = sqrt( ad*ad + ab*ab - 2.0*ad*ab*cos(alpha1+alpha2) );
        cout << ab << ", " << ad << ", " << alpha1 << ", " << alpha2 << ", " << bd << endl;

        remove_edge_d6y(edge_t(A,C));
        if(B > D) swap(A,C);

        insert_edge_d6y(edge_t(B,D),A,C,bd);
        
        new_edges.erase(new_edges.begin()+i);
        new_edges.push_back(edge_t(B,D));

        flips++;
        done = false;
        cout << "flip done" << endl;
      }
      else{ cout << q << " is delaunay, all good." << endl; }
    }
  }
}

void FulleroidDelaunay::align_hole(vector<node_t>& hole) const
{
  bool done = false;
  while(!done){ // Rotate hole until hole[0] is connected only to hole[-1] and hole[1].
    const vector<node_t>& n0(neighbours[hole[0]]);    
    done = true;

    for(int i=2; i< hole.size()-1; i++)
      if(find(n0.begin(), n0.end(), hole[i]) != n0.end()){
	hole.push_back(hole[0]);
	hole.erase(hole.begin());
	done = false;
	break;
      }
  }
}

vector<double> FulleroidDelaunay::hole_angles(const node_t&v, const vector<node_t>& nv) const 
{
  vector<double> angles(nv.size());
  for(int i=0;i<nv.size();i++) angles[i] = angle_d6y(nv[i],v,nv[(i+1)%nv.size()]);
  return angles;
}

// Distance from hole[0] to every other element in the hole
vector<double> FulleroidDelaunay::new_distances_d6y(const node_t& v, const vector<node_t>& hole) const
{
  const size_t n = hole.size();
  vector<double> distances(n);
  
  double accumulated_angle = 0;//angle_d6y(hole[0], v, hole[1]);//angles[0]
  double d0 = edge_lengths_d6y(v, hole[0]);    

  distances[0]   = 0;
  for (int i=1; i<n; i++){
    double di    = edge_lengths_d6y(v, hole[i+1]);
    double angle = angle_d6y(hole[i], v, hole[(i+1)%n]); // angles[i-1];
    accumulated_angle += angle;
    distances[i] = sqrt( d0*d0 + di*di - 2.0*d0*di*cos(accumulated_angle) ) ;
    cerr << distances[i] << endl;
  }  
  return distances;
}
  

vector<edge_t> FulleroidDelaunay::triangulate_hole_d6y(const vector<node_t>& hole, const vector<double>& new_distances) 
{
  vector<edge_t> triangle_edges;
  cerr << "triangulate hole " << hole << endl;

  for (int i=2; i< hole.size()-1; i++){
    cerr << "hole[" << i << "]: " << hole[i] << endl;

    node_t a=hole[0], b=hole[hole.size()-1], c=hole[i], d=hole[i-1];
    if(hole[0] > hole[i]) swap(b,d);

    cerr << vector<int>({a,b,c,d}) << endl;

    insert_edge_d6y(edge_t(a,c),b,d,new_distances[i-2]);
    triangle_edges.push_back(edge_t(a,c));

    cerr << neighbours << endl;
  }
  return triangle_edges;
}

void FulleroidDelaunay::remove_flat_vertex(node_t v)
{
  cerr << "begin remove flat vertex" << endl;
  vector<node_t> hole(neighbours[v]);
  cerr << "hole: " << hole << endl;

  // check if hole[0] is already connected to any of the other hole-nodes in
  // which case we have to start the fan-connecting from somewhere else
  align_hole(hole);

  // get distances of hole[0] to all hole[2]--hole[n-2] within the hole and before removing the vertex
  vector<double> new_distances = new_distances_d6y(v,hole);
  cerr << "new distances: " << new_distances << endl;

  // remove vertices from graph and distance mtx
  for(int i=0; i<hole.size(); i++)
    remove_edge_d6y(edge_t(v,hole[i]));

  neighbours.pop_back();
  N--;

  //triangulate hole
  vector<edge_t> triangle_edges = triangulate_hole_d6y(hole,new_distances);
  
  // delaunayify hole (not sure if it's enough to only delaunayify the hole)
  delaunayify_hole_2(triangle_edges);

  cerr << "g=" << *this << endl;
}

void FulleroidDelaunay::remove_flat_vertices()
{

  // TODO: Assume that vertices are sorted such that hexagons appear at the
  // end, i.e., vertices 12 -- N-1 are hexagons. This needs a fix for
  // fulleroids with negative curvature.

  cerr << "neighbours: " << neighbours << endl;
  cerr << "neighbours-size: " << neighbours.size() << endl;


  while(neighbours.size() > 12){
    cerr << "-----" << endl;
    assert(edge_lengths_d6y_are_symmetric());
    cerr << "edge_lengths_d6y: " << edge_lengths_d6y << endl;
    node_t remove = neighbours.size()-1;
    cerr << "remove node-" << remove << endl;
    remove_flat_vertex(remove);
    cerr << "removed node-" << remove << endl;
    cerr << "neighbours: " << neighbours << endl;
    cerr << "-----" << endl;
  }

  cerr << "--- done ---" << endl;

  cerr << edge_lengths_d6y << endl;

  set<edge_t> set_ue = undirected_edges();
  vector<edge_t> vec_ue;
  for(set<edge_t>::iterator it = set_ue.begin(), to = set_ue.end(); it!=to; it++){
    vec_ue.push_back(*it);
  }
  delaunayify_hole_2(vec_ue);

  cerr << edge_lengths_d6y << endl;
  cerr << "--- done ---" << endl;
  
  // TODO: Perform a final Delaunayification.

}

