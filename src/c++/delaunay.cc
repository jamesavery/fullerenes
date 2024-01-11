#include "fullerenes/delaunay.hh"
#include "fullerenes/debug.hh"

struct cos_sin: public pair<double,double> {

  cos_sin(double a, double b, double c){
    double th = tan_half(a,b,c);
    double ct = cot(a,b,c);
    double cs = (1 - th*th)/(1 + th*th);
    double sn = cs/ct;
   
    first = cs; second = sn;
  }

  cos_sin() : pair<double,double>(1,0) {}
  cos_sin(double cs, double sn) : pair<double,double>(cs,sn) {}
  cos_sin(double angle) : pair<double,double>(cos(angle),sin(angle)) {}

  cos_sin operator+(const cos_sin& b) const {
    return cos_sin(first*b.first - second*b.second, first*b.second+second*b.first);
  }

  cos_sin& operator+=(const cos_sin& b){
    cos_sin result(*this + b);
    return *this = result;
  }
  
  double angle() const { return atan2(first,second); }

  static double tan_half(double a, double b, double c){ 
    return sqrt((a - b + c)*(a + b - c)/((a + b + c)*(-a + b + c))); 
  }

  static double cot(double a, double b, double c){
    double tha = tan_half(a,b,c);
    return (1 - tha*tha)/(2*tha);
  }

};

ostream& operator<<(ostream& s, const FulleroidDelaunay::Quad& q) 
{
  s << "{" << q.v[0] << ", " << q.v[1] << ", " << q.v[2] << ", " << q.v[3] << "}";
  return s;
}


double FulleroidDelaunay::tan_halfangle(node_t vi, node_t vj, node_t vk) const
  { 
    /*   b /vi
     *    / |
     *  vk  | a
     *    \ |
     *   c \vj
     */
    double a = edge_lengths(vi,vj), b = edge_lengths(vk,vi), c = edge_lengths(vk,vj);

    // We must never encounter the need to read an edge length of an edge that isn't there.
    if(!(a != 0 && b != 0 && c !=0 )){
      fprintf(stderr,"Assertion failed in tan_halfangle(%d,%d,%d): a != 0 && b != 0 && c !=0.\n",vi+1,vj+1,vk+1);
      cerr << "(a,b,c) = " << vector<double>{a,b,c} << ";\n";
      abort();
    } 
    
    return sqrt((a-b+c)*(a+b-c)/((a+b+c)*(-a+b+c)));
  }


double FulleroidDelaunay::cot_angle(node_t vi, node_t vj, node_t vk) const 
{
  double tha = tan_halfangle(vi,vj,vk);
  return (1-tha*tha)/(2*tha);
}

double FulleroidDelaunay::add_tan(double tan_a, double tan_b)
{
  return (tan_a + tan_b)/(1-tan_a*tan_b);
}

  
double FulleroidDelaunay::tan_adh(const Quad& Q) const
{
  /*  tan((alpha+delta)/2), where:
   *
   *   / alpha
   *  v0--- e ---
   *   \ delta
   */
  double tha1  = tan_halfangle(Q.v[1],Q.v[2],Q.v[0]), tha2 = tan_halfangle(Q.v[2],Q.v[3],Q.v[0]);
  double t_adh = (tha1+tha2)/(1-tha1*tha2);
  return t_adh;
}

double FulleroidDelaunay::cos_ad(const Quad& Q) const
{
  /*  cos(alpha+delta), where:
   *
   *   / alpha
   *  v0--- e ---
   *   \ delta
   */
  double t_adh = tan_adh(Q);    
  double c_ad = (1-t_adh*t_adh)/(1+t_adh*t_adh);
  return c_ad;
}

double FulleroidDelaunay::flipped_length(const Quad& Q) const
{
  /*  length of side f, where:
   *       
   *  b /  f
   *   /   |
   *  v0-e-- 
   *   \   |
   *  c \  f
   */
  double b = edge_lengths(Q.v[0],Q.v[1]), c = edge_lengths(Q.v[0],Q.v[3]), c_ad = cos_ad(Q);

  // We must never encounter the need to read an edge length of an edge that isn't there.
  assert(b != 0 && c !=0); 

  double f = sqrt(b*b+c*c-2*b*c*c_ad);

  return f;
}

bool FulleroidDelaunay::is_delaunay(const Quad& Q) const {
  double d[4];
  const double eps = 1e-8;
  /*   B 
   * 0/ \1
   * A-e-C
   * 3\ /2
   *   D
   *
   * cot(∠ABC) + cot(∠ADC) >= 0
   */
  double e = edge_lengths(Q.v[0],Q.v[2]);
  for(int i=0;i<4;i++) d[i] = edge_lengths(Q.v[i],Q.v[(i+1)%4]);

  return cos_sin::cot(e,d[0],d[1]) + cos_sin::cot(e,d[2],d[3]) + eps >= 0;
}


bool FulleroidDelaunay::is_consistent(const Quad& Q) const { 
  double d[4];
  for(int i=0;i<4;i++) d[i] = edge_lengths(Q.v[i],Q.v[(i+1)%4]);

  double e = edge_lengths(Q.v[0],Q.v[2]);
  double f = flipped_length(Q);
  
  cos_sin cs_a(f,d[0],d[3]), cs_a1(d[1],d[0],e), cs_a2(d[2],d[3],e);
  cos_sin cs_b(e,d[0],d[1]), cs_b1(d[3],d[0],f), cs_b2(d[2],d[1],f);
  cos_sin cs_c(f,d[1],d[2]), cs_c1(d[0],e,d[1]), cs_c2(d[3],d[2],e);
  cos_sin cs_d(e,d[1],d[3]), cs_d1(d[0],f,d[3]), cs_d2(d[1],d[2],f);
  
  cos_sin angle[4] = {cs_a,cs_b,cs_c,cs_d};
  cos_sin angle1[4] = {cs_a1,cs_b1,cs_c1,cs_d1};
  cos_sin angle2[4] = {cs_a2,cs_b2,cs_c2,cs_d2}; 

  for(int i=0;i<4;i++){
    if(fabs(angle[i].angle()-angle1[i].angle()-angle2[i].angle()) > 1e-5) return false;
  }
  return true;
}


#include <stack>
vector<arc_t> FulleroidDelaunay::delaunayify_hole(const vector<edge_t>& edges)
{
  Debug debug("Delaunay",Debug::INFO1);
  vector<arc_t> new_edges(edges.begin(), edges.end());
  stack<arc_t,vector<arc_t> >  S(new_edges);

  map<edge_t,bool> mark;
  for(auto e: new_edges) mark[edge_t(e)] = true;

  int flips = 0;
  while(!S.empty()){ 
    const arc_t AC = S.top(); S.pop(); mark[edge_t(AC)] = false;
    debug << "Next edge to check is " << AC << ".\n";

    node_t A = AC.first, C = AC.second;
    node_t B = next(C,A), D = next(A,C);

    Quad Q(A,B,C,D);

    if(!is_delaunay(Q)){ // Do a flip!
      debug << "Flipping " << (Q.to_vector()+1) << endl;
      debug << "gg = " << *this << ";\n";
      
      debug << "AC = " << edge_lengths(A,C) << ";\n";

      flip(Q);
    
      if(!mark[{A,B}]){ S.push(arc_t({A,B})); mark[{A,B}] = true; }
      if(!mark[{B,C}]){ S.push(arc_t({B,C})); mark[{B,C}] = true; }
      if(!mark[{C,D}]){ S.push(arc_t({C,D})); mark[{C,D}] = true; }
      if(!mark[{D,A}]){ S.push(arc_t({D,A})); mark[{D,A}] = true; }

      new_edges.erase(find(new_edges.begin(),new_edges.end(), arc_t(A,C)));
      new_edges.push_back(arc_t(B,D));
      
      flips++;
    } else{ debug << (Q.to_vector()+1) << " is delaunay, all good." << endl; }
  }

  return new_edges;
}

void FulleroidDelaunay::flip(const Quad& Q) 
{
  Debug debug("Delaunay",Debug::INFO1);
  node_t A = Q.v[0], B = Q.v[1], C = Q.v[2], D = Q.v[3];
  double f = flipped_length(Q);

  debug << "edge_length["<<(B+1)<<","<<(D+1)<<"] = " << f << ";\n";

  insert_edge(edge_t(B,D),C,A,f); // Check predecessors
  remove_edge(edge_t(A,C));
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
        cout << "rotated hole" << endl;
        break;
      }
  }
}


// Distance from hole[0] to every other element in the hole
vector<double> FulleroidDelaunay::new_distances(const node_t& v, const vector<node_t>& hole) const
{
  const size_t n = hole.size();
  vector<double> distances(n);
  
  cos_sin cossin_acc(1,0);

  double d0 = edge_lengths(v, hole[0]);    

  distances[0]   = 0;
  for (int i=1; i<n; i++){
    double di = edge_lengths(v, hole[i+1]);
    
    double 
      a = edge_lengths(hole[i],hole[(i+1)%n]),
      b = edge_lengths(v,hole[i]),
      c = edge_lengths(v,hole[(i+1)%n]);

    cossin_acc += cos_sin(a,b,c);

    distances[i] = sqrt( d0*d0 + di*di - 2.0*d0*di*cossin_acc.first ) ;
  }  
  return distances;
}
  

vector<edge_t> FulleroidDelaunay::triangulate_hole(const vector<node_t>& hole, const vector<double>& new_distances) 
{
  vector<edge_t> triangle_edges;
  Debug debug("Delaunay",Debug::INFO1);
  
  debug << "(* Triangulate hole *) " << "hole = " << (hole+1) << endl;

  for (int i=2; i< hole.size()-1; i++){
    debug << "hole[" << i << "]: " << (hole[i]+1) << endl;

    node_t a=hole[0], b=hole[hole.size()-1], c=hole[i], d=hole[i-1];
    if(hole[0] > hole[i]) swap(b,d);

    debug << "(a,b,c,d) = " << (vector<int>({a,b,c,d})+1) << endl;

    insert_edge(edge_t(a,c),b,d,new_distances[i]);
    triangle_edges.push_back(edge_t(a,c));

    debug << "neighbours = " << (neighbours+1) << ";\n";
    //    debug << "lengths    = " << edge_lengths << ";\n";
  }
  return triangle_edges;
}

void FulleroidDelaunay::remove_flat_vertex(node_t v)
{
  Debug debug("Delaunay",Debug::INFO1);

  debug << "(*begin remove flat vertex*)" << endl;
  vector<node_t> hole(neighbours[v]);
  debug << "hole=" << (hole+1) << "\n";

  // check if hole[0] is already connected to any of the other hole-nodes in
  // which case we have to start the fan-connecting from somewhere else
  align_hole(hole);

  // get distances of hole[0] to all hole[2]--hole[n-2] within the hole and before removing the vertex
  vector<double> dist = new_distances(v,hole);
  debug << "dist=" << dist << ";\n";

  // remove edges from graph and distance mtx
  debug << "(* Remove edges from vertex "<<(v+1)<<" *)\n";
  for(int i=0; i<hole.size(); i++)
    remove_edge(edge_t(v,hole[i]));

  neighbours.pop_back();
  N--;

  //triangulate hole
  debug << "(* Triangulating hole.*)\n";
  vector<edge_t> triangle_edges = triangulate_hole(hole,dist);
  
  // delaunayify triangulation
  debug << "(* Delaunayifying hole.*)\n";
  vector<edge_t> indel_edges(triangle_edges);
  for(int i=0;i<hole.size();i++)
    indel_edges.push_back({hole[i],hole[(i+1)%hole.size()]});
  delaunayify_hole(indel_edges);
}

void FulleroidDelaunay::remove_flat_vertices()
{
  MathematicaDebug debug("Delaunay",0);
  // Assumes that vertices are sorted such that hexagons appear at the
  // end. Procedure incrementally removes vertices from the back until 
  // reaching a vertex that was not of degree 6 in the initial graph.

  debug << "neighbours=" << neighbours << ";\n";
  debug << "neighbourssize=" << neighbours.size() << ";\n";

  vector<int> original_degrees(N);
  for(node_t v=0;v<N;v++) original_degrees[v] = neighbours[v].size();

  int i=1;
  node_t v = neighbours.size()-1;
  while(original_degrees[v] == 6){
    assert(edge_lengths_are_symmetric());
    debug << "(* removing node " << (v+1) << " *)\n";
    MathematicaDebug::channel_stream["Delaunay"]->flush();
    remove_flat_vertex(v);
    debug << "g["<<(i++)<<"] = " << *this << ";\n";
    v = neighbours.size()-1;
    MathematicaDebug::channel_stream["Delaunay"]->flush();
  }

  debug << "(* --- done --- *)" << endl;
}

