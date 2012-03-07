#ifndef GEOMETRY_HH
# define GEOMETRY_HH

#include <string.h>

using namespace std;

typedef unsigned int node_t;
typedef vector< vector<node_t> > neighbours_t;
typedef vector< bool > edges_t;


struct edge_t : public pair<node_t,node_t> {
  edge_t(const node_t& u, const node_t& v): pair<node_t,node_t>(min(u,v),max(u,v)) {}
  edge_t(const unsigned int index) {
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
};

struct face_t : public vector<node_t> {
  face_t() {}
  face_t(const vector<node_t>& vertices) : vector<node_t>(vertices) {}
  bool operator==(const face_t& B) const { 
    // Two faces are the same if they contain the same vertices
    // (I.e., we disregard the orientation)
    return set<node_t>(begin(),end()) == set<node_t>(B.begin(),B.end());
  }
  bool operator<(const face_t& B) const {
    return set<node_t>(begin(),end()) < set<node_t>(B.begin(),B.end());
  }
  friend ostream& operator<<(ostream &s, const face_t& f){
    s << "["; for(unsigned int i=0;i<f.size();i++) s << f[i] << (i+1<f.size()?", ":"]"); return s;
  }
};


struct coord2d : public pair<double,double> {
  coord2d(const double x=0, const double y=0) : pair<double,double>(x,y) {}
  coord2d operator/(const double s)   const { return coord2d(first/s,second/s); }
  coord2d operator*(const double s)   const { return coord2d(first*s,second*s); }
  coord2d operator+(const coord2d& y) const { return coord2d(first+y.first,second+y.second); }
  coord2d operator-(const coord2d& y) const { return coord2d(first-y.first,second-y.second); }
  coord2d& operator+=(const coord2d& y){ first += y.first; second += y.second; return *this; }
  double dot(const coord2d& y) const { return first*y.first+second*y.second; }

  double line_angle(const coord2d& v) const { // CCW between two lines ([-pi;pi] where -*this is +/-pi and *this is 0 -- i.e. pi is "leftmost", -pi is "rightmost")
    double angle = fmod(atan2(first*v.second-v.first*second,first*v.first+second*v.second)+2*M_PI,2*M_PI)-M_PI;
    //    fprintf(stderr,"angle[{%f,%f},{%f,%f}] == %f\n",first,second,v.first,v.second,angle);
    return angle;
  }
  double point_angle(const coord2d& x=0) const { // CCW angle of x around *this ([-pi;pi])
    const coord2d vx(x-*this);
    double angle = atan2(vx.second,vx.first);
    return angle;
  }
  
  double norm() const { return sqrt(first*first+second*second); }

  friend ostream& operator<<(ostream &s, const coord2d& x){ s << fixed << "{" << x.first << "," << x.second << "}"; return s; }
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
  coord3d& operator*=(const double& y){ x[0] *= y; x[1] *= y; x[2] *= y; return *this; }
  coord3d& operator/=(const double& y){ x[0] /= y; x[1] /= y; x[2] /= y; return *this; }

  coord3d cross(const coord3d& y) const {
    return coord3d(x[1]*y[2]-x[2]*y[1], x[2]*y[0]-x[0]*y[2], x[0]*y[1] - x[1]*y[0]);
  }
  double dot(const coord3d& y) const { return x[0]*y[0]+x[1]*y[1]+x[2]*y[2]; }
  double norm() const { return sqrt(dot(*this)); }

  double& operator[](unsigned int i){ return x[i]; }
  double  operator[](unsigned int i) const { return x[i]; }

  friend ostream& operator<<(ostream &s, const coord3d& x){ s << fixed << "{" << x[0] << "," << x[1] << "," << x[2]<< "}"; return s; }
};


struct Tri3D {
  const coord3d a,b,c,u,v,n;
  typedef pair<coord3d,coord3d> segment_t;

  Tri3D(const coord3d *T) : a(T[0]), b(T[1]), c(T[2]), u(b-a), v(c-a), n(u.cross(v)) {}
  double plane_intersection(const segment_t& segment) const { 
    return n.dot(a-segment.first)/n.dot(segment.second-segment.first);
  }
  bool   intersects(const coord3d& x) const {
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

  friend ostream& operator<<(ostream& s, const Tri3D& T){
    s << "{" << T.a << "," << T.b << "," << T.c << "}";
    return s;
  }
};

struct sort_ccw_point {
  const vector<coord2d>& layout;
  const coord2d& x;
  sort_ccw_point(const vector<coord2d>& layout, const coord2d& x) : layout(layout), x(x)
  { }
  
  bool operator()(const node_t& s, const node_t& t) const {
    double angs = x.point_angle(layout[s]), 
           angt = x.point_angle(layout[t]);

    // printf("compare: %d:{%g,%g}:%g %d:{%g,%g}:%g\n",
    // 	   s,layout[s].first,layout[s].second,angs,
    // 	   t,layout[t].first,layout[t].second,angt);
    return angs <= angt; 
  }
};


#endif
