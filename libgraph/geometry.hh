#ifndef GEOMETRY_HH
# define GEOMETRY_HH

#include <string.h>
#include <iostream>
#include <vector>
#include <set>
#include <map>
#include <math.h>
using namespace std;

typedef int node_t;
typedef vector< vector<node_t> > neighbours_t;
typedef vector< bool > edges_t;

template <typename T> ostream& operator<<(ostream& s, const vector<T>& v)
{
  s << "{";
  for(int i=0;i<v.size();i++) s << v[i] << (i+1<v.size()? ",":"");
  s << "}";
  return s;
}

template <typename S, typename T> ostream& operator<<(ostream& s, const pair<S,T>& p)
{
  s << "{" << p.first << "," << p.second << "}";
  return s;
}

// Undirected edge is an unordered pair of nodes
struct edge_t : public pair<node_t,node_t> {
  edge_t(const pair<node_t,node_t>& p) : pair<node_t,node_t>(min(p.first,p.second),max(p.first,p.second)) {}
  edge_t(const node_t& u, const node_t& v): pair<node_t,node_t>(min(u,v),max(u,v)) {}
  edge_t(const int index) {
    node_t u=0;
    for(;u*(u-1)/2<=index;u++) ;
    u--;
    first = u;
    second = index-u*(u-1)/2;
  }
  inline size_t index() const { 
    const node_t v = first, u = second;
    return u*(u-1)/2 + v; 
  }
  friend ostream& operator<<(ostream& s, const edge_t& e){ s<<"{"<<e.first <<","<<e.second<<"}"; return s; }
};

// Directed edge is an ordered pair of nodes
typedef pair<node_t,node_t> dedge_t;

struct coord2d : public pair<double,double> {
  coord2d(const double x=0, const double y=0) : pair<double,double>(x,y) {}
  coord2d operator/(const double s)   const { return coord2d(first/s,second/s); }
  coord2d operator*(const double s)   const { return coord2d(first*s,second*s); }
  coord2d operator*(const coord2d& y)   const { return coord2d(first*y.first,second*y.second); }
  coord2d operator+(const coord2d& y) const { return coord2d(first+y.first,second+y.second); }
  coord2d operator-(const coord2d& y) const { return coord2d(first-y.first,second-y.second); }
  coord2d& operator+=(const coord2d& y){ first += y.first; second += y.second; return *this; }
  coord2d& operator*=(const coord2d& y){ first *= y.first; second *= y.second; return *this; }
  coord2d& operator*=(const double s)  { first*=s;second*=s; return *this;}
  double dot(const coord2d& y) const { return first*y.first+second*y.second; }

  double line_angle(const coord2d& v) const { // CCW between two lines ([-pi;pi] where -*this is +/-pi and *this is 0 -- i.e. pi is "leftmost", -pi is "rightmost")
    double angle = fmod(atan2(first*v.second-v.first*second,first*v.first+second*v.second)+2*M_PI,2*M_PI)-M_PI;
    //    fprintf(stderr,"angle[{%f,%f},{%f,%f}] == %f\n",first,second,v.first,v.second,angle);
    return angle;
  }

  double point_angle(const coord2d& y=0, const bool periodic=false) const { 
    if(periodic) return point_angle_periodic(y);
    else {
      const coord2d vx(*this-y);
      return atan2(vx.second,vx.first);
    }
  }

  // CCW angle of y around point on a periodic surface [0;2pi[ x [0;pi[. 
  double point_angle_periodic(const coord2d& y=0) const {
    const coord2d dy[4] = {coord2d(0,0), coord2d(2*M_PI,0), coord2d(0,M_PI), coord2d(2*M_PI,M_PI)};
    const coord2d& x(*this);

    // Step 1: There are four potential "proper" coordinates for y:  y, y+(2pi,0), y+(0,pi), or y+(2pi,pi)
    //         The appropriate one has the smallest distance to x.
    double rsqrmin = INFINITY;
    int imin = 0;
    for(int i=0;i<4;i++){
      const coord2d y0(y+dy[i]-x);
      const double rsqr = y0.dot(y0);
      if(rsqr < rsqrmin){
	rsqrmin = rsqr;
	imin = i;
      }
    }
    // Step 2: Now the angle is calculated as usual.
    return point_angle(y+dy[imin],false);
  }
  
  double norm() const { return sqrt(first*first+second*second); }

  friend ostream& operator<<(ostream &s, const coord2d& x){ s << fixed << "{" << x.first << "," << x.second << "}"; return s; }
  friend istream& operator>>(istream &s, coord2d& x){ s >> x.first; s >> x.second; return s; }
};


struct coord3d {
  double x[3];

  coord3d(const double y[3]) { memcpy(x,y,3*sizeof(double)); }
  coord3d(const double x_=0, const double y=0, const double z=0) { x[0] = x_; x[1] = y; x[2] = z; }
  coord3d operator/(const double s)   const { return coord3d(*this) /= s; }
  coord3d operator*(const double s)   const { return coord3d(*this) *= s; }
  coord3d operator+(const coord3d& y) const { return coord3d(*this) += y; }
  coord3d operator-(const coord3d& y) const { return coord3d(*this) -= y; }
  coord3d& operator+=(const coord3d& y){ x[0] += y[0]; x[1] += y[1]; x[2] += y[2]; return *this; }
  coord3d& operator-=(const coord3d& y){ x[0] -= y[0]; x[1] -= y[1]; x[2] -= y[2]; return *this; }
  coord3d& operator*=(const coord3d& y){ x[0] *= y[0]; x[1] *= y[1]; x[2] *= y[2]; return *this; }
  coord3d& operator*=(const double& y){ x[0] *= y; x[1] *= y; x[2] *= y; return *this; }
  coord3d& operator/=(const double& y){ x[0] /= y; x[1] /= y; x[2] /= y; return *this; }

  coord3d cross(const coord3d& y) const {
    return coord3d(x[1]*y[2]-x[2]*y[1], x[2]*y[0]-x[0]*y[2], x[0]*y[1] - x[1]*y[0]);
  }
  double dot(const coord3d& y) const { return x[0]*y[0]+x[1]*y[1]+x[2]*y[2]; }
  double norm() const { return sqrt(dot(*this)); }
  coord2d polar_angle(const coord3d& centre = coord3d()) const
  { 
    const coord3d& x(*this - centre);
    const double r = norm();
    return coord2d(acos(x[2]/r),atan2(x[1],x[0]));
  }

  double& operator[](unsigned int i){ return x[i]; }
  double  operator[](unsigned int i) const { return x[i]; }



  friend ostream& operator<<(ostream &s, const coord3d& x){ s << fixed << "{" << x[0] << "," << x[1] << "," << x[2]<< "}"; return s; }
  friend istream& operator>>(istream &s, coord3d& x){ for(int i=0;i<3;i++){ s >> x[i]; } return s; }
};


struct face_t : public vector<node_t> {
  face_t(const size_t size=0) : vector<node_t>(size) {}
  face_t(const vector<node_t>& vertices) : vector<node_t>(vertices) {}
  
  bool operator==(const face_t& B) const { 
    // Two faces are the same if they contain the same vertices
    // (I.e., we disregard the orientation)
    return set<node_t>(begin(),end()) == set<node_t>(B.begin(),B.end());
  }
  bool operator<(const face_t& B) const {
    return set<node_t>(begin(),end()) < set<node_t>(B.begin(),B.end());
  }
  
  coord2d centroid(const vector<coord2d>& layout) const { 
    coord2d c(0);
    for(size_t i=0;i<size();i++) c += layout[(*this)[i]];
    return c/size();
  }
  coord3d centroid(const vector<coord3d>& layout) const { 
    coord3d c(0.0);
    for(size_t i=0;i<size();i++) c += layout[(*this)[i]];
    return c/size();
  }

  // http://www.softsurfer.com/Archive/algorithm_0103/algorithm_0103.htm#Winding%20Number
  double winding_number(const vector<coord2d>& layout, const coord2d& x) const {
    vector<coord2d> Cp(size());
    for(int i=0;i<size();i++)
      Cp[i] = layout[(*this)[i]]-x;

    double wn = 0;
    for(int i=0;i<size();i++){
      coord2d segment = Cp[(i+1)%size()] - Cp[i];
      double theta = atan2(segment.second,segment.first);
      wn += theta;
    }
    wn /= 2*M_PI;
    return wn;
  }

  bool point_inside(const vector<coord2d>& layout, const coord2d& x) const {
    return winding_number(layout,x) != 0;
  }
  bool contains(const node_t v) const { for(int i=0;i<size();i++) if(v == (*this)[i]) return true; return false; }
  
  

  friend ostream& operator<<(ostream &s, const face_t& f){
    s << "{"; for(unsigned int i=0;i<f.size();i++) s << f[i] << (i+1<f.size()?", ":"}"); return s;
  }
};

struct tri_t : public face_t {
  tri_t(const node_t a,const node_t b,const node_t c) : face_t(3) { u(0) = a; u(1) = b; u(2) = c; }
  tri_t(const vector<node_t>& f) : face_t(f) {}
  node_t& u(const unsigned int i)  { return (*this)[i]; }
  const node_t& u(const unsigned int i) const  { return (*this)[i]; }

  coord3d centroid(const vector<coord3d>& points) const { return (points[u(0)]+points[u(1)]+points[u(2)])/3.0; }
  coord2d centroid(const vector<coord2d>& points) const { return (points[u(0)]+points[u(1)]+points[u(2)])/3.0; }
  void flip(){ node_t t = u(1); u(1) = u(2); u(2) = t; }
};

typedef map<unsigned int,set<face_t> > facemap_t;

struct Tri3D {
  coord3d a,b,c,u,v,n;
  typedef pair<coord3d,coord3d> segment_t;

  Tri3D(const coord3d *T) : a(T[0]), b(T[1]), c(T[2]), u(b-a), v(c-a), n(u.cross(v)) {}
  Tri3D(const coord3d& a,const coord3d& b,const coord3d& c) : a(a), b(b), c(c), u(b-a), v(c-a), n(u.cross(v)) {}
  double plane_intersection(const segment_t& segment) const { 
    return n.dot(a-segment.first)/n.dot(segment.second-segment.first);
  }
  Tri3D(const vector<coord3d>& points, const tri_t& t) : a(points[t[0]]),b(points[t[1]]),c(points[t[2]]), u(b-a), v(c-a), n(u.cross(v)) {}

  double distance(const coord3d& x) const {
    return fabs((x-a).dot(n))/n.norm();
  }
  bool intersects(const coord3d& x) const {
    const coord3d w(x-a);
    double d = u.dot(v)*u.dot(v) - u.dot(u)*v.dot(v);
    double s = (u.dot(v)*w.dot(v)-v.dot(v)*w.dot(u))/d;
    double t = (u.dot(v)*w.dot(u)-u.dot(u)*w.dot(v))/d;
    return s>=0 && t>=0 && s+t<=1;
  }
  bool intersects(const segment_t& segment) const {
    double r = plane_intersection(segment);
    if(r<0 || r>1) return false;
    
    const coord3d x(segment.first + (segment.second-segment.first)*r);
    return intersects(x);
  }

  void flip(){ coord3d t = b; b = c; c = t; t = u; u = v; v = u; n = n*-1; }
  void flip_normal(){ n = n*-1.0; }

  bool back_face(const coord3d& p) const {
    const coord3d centre((a+b+c)/3.0);
    const coord3d line(centre-p);

    return line.dot(n) > 0;
  }

  double area() const {
    return (((b-a).cross(c-a)).norm()/2.0);
  }

  coord3d centroid() const { return coord3d((a+b+c)/3.0); }

  friend ostream& operator<<(ostream& s, const Tri3D& T){
    s << "{" << T.a << "," << T.b << "," << T.c << "}";
    return s;
  }
};

struct Tetra3D {
  coord3d a,b,c,d;

  Tetra3D(const coord3d *T): a(T[0]),  b(T[1]), c(T[2]), d(T[3]) {}
  Tetra3D(const coord3d& a, const coord3d& b, const coord3d& c, const coord3d& d): a(a),  b(b), c(c), d(d) {}

  coord3d centroid() const { return (a+b+c+d)/4.0; }

  double volume() const {
    return fabs((a-d).dot((b-d).cross(c-d)))/6.0;
  }
};

struct sort_ccw_point {
  const vector<coord2d>& layout;
  const coord2d& x;
  const bool periodic;
  sort_ccw_point(const vector<coord2d>& layout, const coord2d& x, const bool periodic = false)
    : layout(layout), x(x), periodic(periodic)
  { }
  
  bool operator()(const node_t& s, const node_t& t) const {
    double angs = x.point_angle(layout[s],periodic), 
           angt = x.point_angle(layout[t],periodic);

    // printf("compare: %d:{%g,%g}:%g %d:{%g,%g}:%g\n",
    // 	   s,layout[s].first,layout[s].second,angs,
    // 	   t,layout[t].first,layout[t].second,angt);
    return angs >= angt; 	// TODO: Is the sign here correct?
  }
};


#endif
