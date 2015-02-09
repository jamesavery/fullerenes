#include "delaunay.hh"


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
  // printf("dist(%d,%d) = %g\n"
  //          "dist(%d,%d) = %g\n"
  //          "dist(%d,%d) = %g\n",
  //          B,A,distances(B,A),
  //          B,C,distances(B,C),
  //          A,C,distances(A,C));
  //  printf("(a,b,c) = (%g,%g,%g) ~> acos(%g) = %g\n",a,b,c,(a*a+b*b-c*c)/(2*a*b),acos((a*a+b*b-c*c)/(2*a*b)));
  return acos((a*a+b*b-c*c)/(2*a*b));
}


double FulleroidDelaunay::angle_d6y(node_t A, node_t B, node_t C) const {
// returns angle at B
  double 
    a = edge_lengths_d6y(B,A),
    b = edge_lengths_d6y(B,C),
    c = edge_lengths_d6y(A,C);
    printf("dist(%d,%d) = %g\n"
           "dist(%d,%d) = %g\n"
           "dist(%d,%d) = %g\n",
           B,A,edge_lengths_d6y(B,A),
           B,C,edge_lengths_d6y(B,C),
           A,C,edge_lengths_d6y(A,C));
    assert(a>epsilon && b>epsilon && c>epsilon);
  if( (a*a+b*b-c*c)/(2.0*a*b) >= 1.0 ){ // catch 0 degree case + numerical inaccuracy
    printf("(a,b,c) = (%g,%g,%g) ~> acos(%g) = %g (exception used)\n",a,b,c,(a*a+b*b-c*c)/(2*a*b),0.0);
    return 0;
  }
  if( (a*a+b*b-c*c)/(2.0*a*b) <= -1.0 ){ // catch 180 degree case + numerical inaccuracy
    printf("(a,b,c) = (%g,%g,%g) ~> acos(%g) = %g (exception used)\n",a,b,c,(a*a+b*b-c*c)/(2*a*b),M_PI);
    return M_PI;
  }
  printf("(a,b,c) = (%g,%g,%g) ~> acos(%g) = %g\n",a,b,c,(a*a+b*b-c*c)/(2*a*b),acos((a*a+b*b-c*c)/(2*a*b)));
  return acos((a*a+b*b-c*c)/(2.0*a*b));
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
        if(B > D){ // FIXME write elegantly
          node_t buf=A; A=C; C=buf;
        }
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

void FulleroidDelaunay::remove_flat_vertex(node_t v)
{
  cout << "begin remove flat vertex" << endl;
  vector<node_t> hole(neighbours[v]);
  cout << "hole: " << hole << endl;

  
  // check if hole[0] is already connected to any of the other hole-nodes in
  // which case we have to start the fan-connecting from somewhere else
  //TODO what about the second one

  for(int i=2; i< hole.size()-1; i++){
    if(find(neighbours[hole[0]].begin(), neighbours[hole[0]].end(), hole[i]) != neighbours[hole[0]].end()){
      hole.push_back(hole[0]);
      hole.erase(hole.begin());
    }
  }

  // get distances of hole[0] to all hole[2]--hole[n-2] within the hole and before removing the vertex
  vector<double> new_distances;
  double accumulated_angle = angle_d6y(hole[0], v, hole[1]);
  cout << "angle: " << accumulated_angle << endl;
  double d0 = edge_lengths_d6y(v, hole[0]);
  for (int i=2; i<hole.size()-1; i++){
    double di = edge_lengths_d6y(v, hole[i]);
    accumulated_angle += angle_d6y(hole[i-1], v, hole[i]);
    cout << "angle: " << accumulated_angle << endl;
    new_distances.push_back( sqrt( d0*d0 + di*di - 2.0*d0*di*cos(accumulated_angle) ) );
  }
  cout << "new distances: " << new_distances << endl;

  // remove vertices from graph and distance mtx
  for(int i=0; i<hole.size(); i++){
    remove_edge_d6y(edge_t(v,hole[i]));
  }
  neighbours.pop_back();
  N--;

  //triangulate hole
  // vector<dedge_t> triangle_edges = triangulate_hole(hole);
  vector<edge_t> triangle_edges;
  cout << "triangulate hole " << hole << endl;
  for (int i=2; i< hole.size()-1; i++){
    cout << "hole[" << i << "]: " << hole[i] << endl;
    node_t a=hole[0], b=hole[hole.size()-1], c=hole[i], d=hole[i-1];
    if(hole[0] > hole[i]){ // FIXME write elegantly
      node_t buf = b; b=d; d=buf;
    }
    cout << a << ", " << b << ", " << c << ", " << d << endl;
    insert_edge_d6y(edge_t(a,c),b,d,new_distances[i-2]);
    cout << neighbours << endl;
    triangle_edges.push_back(edge_t(a,c));
  }
  
  // delaunayify hole (not sure if it's enough to only delaunayify the hole)
  delaunayify_hole_2(triangle_edges);

  cout << "g=" << *this << endl;
}

void FulleroidDelaunay::remove_flat_vertices()
{

  // TODO: Assume that vertices are sorted such that hexagons appear at the
  // end, i.e., vertices 12 -- N-1 are hexagons. This needs a fix for
  // fulleroids with negative curvature.

  cout << "neighbours: " << neighbours << endl;
  cout << "neighbours-size: " << neighbours.size() << endl;


  while(neighbours.size() > 12){
    cout << "-----" << endl;
    assert(edge_lengths_d6y_are_symmetric());
    cout << "edge_lengths_d6y: " << edge_lengths_d6y << endl;
    node_t remove = neighbours.size()-1;
    cout << "remove node-" << remove << endl;
    remove_flat_vertex(remove);
    cout << "removed node-" << remove << endl;
    cout << "neighbours: " << neighbours << endl;
    cout << "-----" << endl;
  }

  cout << "--- done ---" << endl;

  cout << edge_lengths_d6y << endl;

  set<edge_t> set_ue = undirected_edges();
  vector<edge_t> vec_ue;
  for(set<edge_t>::iterator it = set_ue.begin(), to = set_ue.end(); it!=to; it++){
    vec_ue.push_back(*it);
  }
  delaunayify_hole_2(vec_ue);

  cout << edge_lengths_d6y << endl;
  cout << "--- done ---" << endl;
  
  // TODO: Perform a final Delaunayification.

}

