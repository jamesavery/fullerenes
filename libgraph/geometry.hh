#ifndef GEOMETRY_HH
#define GEOMETRY_HH

#include <string.h>
#include <iostream>
#include <vector>
#include <set>
#include <map>
#include <math.h>
#include <sstream>
#include <list>
#include <complex>
#include <algorithm>
#include "auxiliary.hh"

using namespace std;

typedef int node_t;
typedef vector< vector<node_t> > neighbours_t;
typedef vector< bool > edges_t;

struct matrix3d; // required for coord3d.outer(coord3d)

// TODO: geometry.hh is getting huge. Move most of the implementation to geometryc.cc

// Directed edge is an ordered pair of nodes
typedef pair<node_t,node_t> dedge_t;

struct edge_t : public pair<node_t,node_t> {
  edge_t() {}
  edge_t(const pair<node_t,node_t>& p) : pair<node_t,node_t>(min(p.first,p.second),max(p.first,p.second)) {}
  edge_t(const node_t u, const node_t v): pair<node_t,node_t>(min(u,v),max(u,v)) {}
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
};


struct coord2d : public pair<double,double> {
  coord2d(const complex<double>& x) : pair<double,double>(x.real(),x.imag()) {}
  coord2d(const pair<double,double>& x) : pair<double,double>(x) {}
  coord2d(const double x=0, const double y=0) : pair<double,double>(x,y) {}
  coord2d operator/(const double s)   const { return coord2d(first/s,second/s); }
  coord2d operator*(const double s)   const { return coord2d(first*s,second*s); }
  coord2d operator*(const coord2d& y)   const { return coord2d(first*y.first,second*y.second); }
  coord2d operator+(const coord2d& y) const { return coord2d(first+y.first,second+y.second); }
  coord2d operator-(const coord2d& y) const { return coord2d(first-y.first,second-y.second); }
  coord2d& operator+=(const coord2d& y){ first += y.first; second += y.second; return *this; }
  coord2d& operator-=(const coord2d& y){ first -= y.first; second -= y.second; return *this; }
  coord2d& operator*=(const coord2d& y){ first *= y.first; second *= y.second; return *this; }
  coord2d& operator*=(const double s)  { first*=s;second*=s; return *this;}
  coord2d operator-() const {coord2d y(-first,-second); return y;}
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
  // d/dx_i ||x|| = x_i/||x||.
  static coord2d dnorm(const coord2d& x){ return x/x.norm(); }

  static coord2d displacement(const coord2d& x, const coord2d& y, bool layout_is_spherical)
  {
    if(!layout_is_spherical) return x-y;

    int i0=0,j0=0;
    double dmin = INFINITY;

    for(int i=0;i<=1;i++)
      for(int j=0;j<=1;j++){
	double d = (x+coord2d(i*2*M_PI,j*M_PI) - y).norm();
	if(d < dmin){ i0 = i; j0 = j; }
      }
    return x+coord2d(i0*2*M_PI,j0*M_PI) - y;
  }

  friend ostream& operator<<(ostream &s, const coord2d& x){ s << fixed << "{" << x.first << "," << x.second << "}"; return s; }
  friend istream& operator>>(istream &s, coord2d& x){ s >> x.first; s >> x.second; return s; }
};



struct coord3d {
  double x[3];

  coord3d(const double y[3]) { x[0] = y[0]; x[1] = y[1]; x[2] = y[2]; }
  explicit coord3d(const double x_=0, const double y=0, const double z=0) { x[0] = x_; x[1] = y; x[2] = z; }
  coord3d operator/(const double s)   const { return coord3d(*this) /= s; }
  coord3d operator*(const double s)   const { return coord3d(*this) *= s; }
  coord3d operator*(const coord3d& y) const { return coord3d(*this) *= y; }
  coord3d operator+(const coord3d& y) const { return coord3d(*this) += y; }
  coord3d operator-(const coord3d& y) const { return coord3d(*this) -= y; }
  coord3d& operator+=(const coord3d& y){ x[0] += y[0]; x[1] += y[1]; x[2] += y[2]; return *this; }
  coord3d& operator-=(const coord3d& y){ x[0] -= y[0]; x[1] -= y[1]; x[2] -= y[2]; return *this; }
  coord3d& operator*=(const coord3d& y){ x[0] *= y[0]; x[1] *= y[1]; x[2] *= y[2]; return *this; }
  coord3d& operator*=(const double& y){ x[0] *= y; x[1] *= y; x[2] *= y; return *this; }
  coord3d& operator/=(const double& y){ x[0] /= y; x[1] /= y; x[2] /= y; return *this; }
  coord3d operator-() const {coord3d y(-x[0],-x[1],-x[2]); return y;}
  coord3d operator*(const matrix3d& m) const;

  coord3d cross(const coord3d& y) const {
    return coord3d(x[1]*y[2]-x[2]*y[1], x[2]*y[0]-x[0]*y[2], x[0]*y[1]-x[1]*y[0]);
  }
  matrix3d outer(const coord3d& y) const;
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

  static double dist(const coord3d& x, const coord3d& y){ return (x-y).norm(); }
  // d/dx_i ||x|| = x_i/||x||.
  static coord3d dnorm(const coord3d& x){ return x/x.norm(); }
  // d^2/(dx_i dx_j) ||x|| = -x_i x_j/||x||^3 + [i==j]/||x||
  static void ddnorm(const coord3d& x, vector<double> &H)
  {
    const double n = 1.0/x.norm(), n3 = n*n*n;

    for(int i=0;i<3;i++)
      for(int j=0;j<3;j++)
	H[i*3+j] = -x[i]*x[j]*n3 + (i==j? n : 0);
  }

  // calculation of the angle beta at b(0,0,0)
  static double angle(const coord3d& a, const coord3d& c);
  // calculation of the derivative of angle beta at b(0,0,0) according to coordinates a and c with fixed b
  static void dangle(const coord3d& a, const coord3d& c, coord3d& da, coord3d& dc);
  // calculation of the dihedral angle theta at a(0,0,0), b, c and d,  the result is an angle between -\pi and +\pi (in radians)
  static double dihedral(const coord3d& b, const coord3d& c, const coord3d& d);
  // calculation of the derivative of dihedral angle theta at a(0,0,0), b, c and d according to coordinates b, c and d with fixed a
  static void ddihedral(const coord3d& b, const coord3d& c, const coord3d& d, coord3d& db, coord3d& dc, coord3d& dd);

  static double ideal_dihedral(const int lA, const int lB, const int lC, const double ur=1.0, const double us=1.0, const double ut=1.0);


  friend vector<coord3d> &operator-=(vector<coord3d>& xs, const coord3d& y)
  {
    for(int i=0;i<xs.size();i++) xs[i] -= y;
    return xs;
  }

  friend vector<coord3d> &operator*=(vector<coord3d>& xs, const double& y)
  {
    for(int i=0;i<xs.size();i++) xs[i] *= y;
    return xs;
  }

  // NB: Does this belong here?
  static coord3d line_plane_intersect(const coord3d& x0, const coord3d& x1, const coord3d& X0, const coord3d& n)
  {
    const coord3d dx(x1-x0);
    double t = (X0-x0).dot(n)/dx.dot(n);
    return x0+dx*t;
  }

  friend ostream& operator<<(ostream &s, const coord3d& x){ s << fixed << "{" << x[0] << "," << x[1] << "," << x[2]<< "}"; return s; }
  friend istream& operator>>(istream &s, coord3d& x){ for(int i=0;i<3;i++){ s >> x[i]; } return s; }
};



struct matrix3d {
  double values[9];

//  matrix3d()                { memset(values,0,9*sizeof(double)); }
  matrix3d(const double *v) { memcpy(values,v,9*sizeof(double)); }
  explicit matrix3d(const double r=0, const double s=0, const double t=0, const double u=0, const double v=0, const double w=0, const double x=0, const double y=0, const double z=0) {
    values[0]=r; values[1]=s; values[2]=t; values[3]=u; values[4]=v; values[5]=w; values[6]=x; values[7]=y; values[8]=z;
  }

  double& operator()(int i, int j)       { return values[i*3+j]; }
  double  operator()(int i, int j) const { return values[i*3+j]; }
  matrix3d operator+(const matrix3d& y) const { return matrix3d(*this) += y; }
  matrix3d operator-(const matrix3d& y) const { return matrix3d(*this) -= y; }
  matrix3d& operator+=(const matrix3d& y){ for(int i=0;i<3;++i){for(int j=0;j<3;++j){values[3*i+j] += y(i,j);}}; return *this; }
  matrix3d& operator-=(const matrix3d& y){ for(int i=0;i<3;++i){for(int j=0;j<3;++j){values[3*i+j] -= y(i,j);}}; return *this; }
  matrix3d operator*(const double s)   const { return matrix3d(*this) *= s; }
  matrix3d& operator*=(const double& s){ for(int i=0;i<3;++i){for(int j=0;j<3;++j){values[3*i+j] *= s;}}; return *this; }
  matrix3d operator-() const {matrix3d m(-values[0],-values[1],-values[2],-values[3],-values[4],-values[5],-values[6],-values[7],-values[8]); return m;}

  matrix3d transpose() const {
    const matrix3d &M(*this);
    matrix3d Mt;
    for(int i=0;i<3;i++)
      for(int j=0;j<3;j++)
	Mt(i,j) = M(j,i);
    return Mt;
  }

  double norm() const {
    const matrix3d &M(*this);

    double norm = 0;
    for(int i=0;i<3;i++)
      for(int j=0;j<3;j++)
	norm += M(i,j)*M(i,j);

    return sqrt(norm);
  }
  
  matrix3d operator*(const matrix3d& B) const {
    const matrix3d &A(*this);
    matrix3d C;

    for(int i=0;i<3;i++){
      for(int j=0;j<3;j++){
        double sum = 0;
        for(int k=0;k<3;k++){
          sum += A(i,k)*B(k,j);
        }
        C(i,j) = sum;
      }
    }
    return C;
  }

  coord3d operator*(const coord3d& x) const {
    coord3d y;
    for(int j=0;j<3;j++)
      y += coord3d(values[j]*x[j],values[3+j]*x[j],values[6+j]*x[j]);
    return y;
  }
  
  vector<coord3d> operator*(const vector<coord3d>& xs) const {
    const matrix3d &A(*this);
    vector<coord3d> ys(xs.size());

    for(int i=0;i<xs.size();i++) ys[i] = A*xs[i];
    return ys;
  }

  // Eigenvalue solver specialized to symmetric real 3x3 matrices using Viete's 
  // closed form solution to cubic polynomials with three real roots.
  coord3d eigenvalues() const
  {
    const matrix3d &M(*this);
    // Make sure that matrix is symmetric. TODO: FP comparison, not exact.
    assert(M(0,1) == M(1,0) && M(0,2) == M(2,0) && M(1,2) == M(2,1));

    // Coefficients up to symmetry
    double a(M(0,0)), b(M(0,1)), c(M(0,2)), d(M(1,1)), e(M(1,2)), f(M(2,2));
    
    // Coefficients of characteristic polynomial, calculated with Mathematica
    long double 
      A = -1.L,
      B = a+d+f,
      C = b*b + c*c - a*d + e*e - a*f - d*f,
      D = -c*c*d + 2*b*c*e - a*e*e - b*b*f + a*d*f;

    if(D==0){			// Second order equation. TODO: FP comparison
      long double Disc = sqrtl(B*B-4*A*C);
      cout << "D = " << Disc << endl;
      return coord3d(0,(-B-Disc)/(2.L*A),(-B+Disc)/(2.L*A));
    }

    // Depress characteristic polynomial - see http://en.wikipedia.org/wiki/Cubic_equation#Reduction_to_a_depressed_cubic
    long double 
      p  = (3.L*A*C - B*B)/(3.L*A*A),
      q  = (2.L*B*B*B - 9.L*A*B*C + 27.L*A*A*D)/(27.L*A*A*A),
      xc = B/(3.L*A);

    // François Viète's solution to cubic polynomials with three real roots. 
    coord3d t;
    long double K = 2*sqrtl(-p/3.L), 
                theta0 = (1.L/3.L)*acosl((3.L*q)/(2.L*p)*sqrtl(-3.L/p));
    for(int k=0;k<3;k++) t[k] = K*cosl(theta0-k*2.L*M_PI/3.L);

    // lambda = t - B/(3A)
    return t - coord3d(xc,xc,xc);
  }
  
  coord3d eigenvector(const double lambda) const {
// [ a_12 * a_23 - a_13 * (a_22 - r) ]
// [ a_12 * a_13 - a_23 * (a_11 - r) ]
// [(a_11 - r) * (a_22 - r) - a_12^2 ]

    const matrix3d &M(*this);
    coord3d x(M(0,1)*M(1,2) - M(0,2)*(M(1,1)-lambda),
	      M(0,1)*M(0,2) - M(1,2)*(M(0,0)-lambda),
	      (M(0,0)-lambda)*(M(1,1)-lambda)-M(0,1)*M(0,1));
    return x/x.norm();
  }

  pair<coord3d,matrix3d> eigensystem() const {
    coord3d lambda(eigenvalues());
    matrix3d C;

    // Sort eigenvalues by absolute value, smallest first
    if(fabs(lambda[0]) > fabs(lambda[1])) std::swap(lambda[0],lambda[1]);
    if(fabs(lambda[1]) > fabs(lambda[2])) std::swap(lambda[1],lambda[2]);
    if(fabs(lambda[0]) > fabs(lambda[1])) std::swap(lambda[0],lambda[1]);

    // Build eigenvector matrix
    for(int i=0;i<3;i++){
      coord3d c(eigenvector(lambda[i]));
      for(int j=0;j<3;j++) C(i,j) = c[j];
    }
    return make_pair(lambda,C);
  }

  static matrix3d unit_matrix(){
    return matrix3d(1,0,0,0,1,0,0,0,1);
  }
  
  friend ostream& operator<<(ostream& S, const matrix3d &M)
  {
    S << "{"; for(int i=0;i<3;i++) S << vector<double>(&M.values[i*3],&M.values[(i+1)*3]) << (i+1<3?",":"}");
    return S;
  }
};



struct tri_t {
  int x_[3];
  tri_t(const node_t a=0,const node_t b=0,const node_t c=0) { u(0) = a; u(1) = b; u(2) = c; }
  tri_t(const vector<node_t>& f) { for(int i=0;i<3;i++) x_[i] = f[i]; }

  node_t& operator[](const unsigned int i)  { return x_[i]; }
  const node_t& operator[](const unsigned int i) const  { return x_[i]; }
  
  node_t& u(const unsigned int i)  { return x_[i]; }
  const node_t& u(const unsigned int i) const  { return x_[i]; }

  coord3d centroid(const vector<coord3d>& points) const { return (points[u(0)]+points[u(1)]+points[u(2)])/3.0; }
  coord2d centroid(const vector<coord2d>& points) const { return (points[u(0)]+points[u(1)]+points[u(2)])/3.0; }
  void flip(){ node_t t = u(1); u(1) = u(2); u(2) = t; }

  bool operator!=(const tri_t& x) const { return x_[0] != x[0] || x_[1] != x[1] || x_[2] != x[2]; }
  bool operator==(const tri_t& x) const { return x_[0] == x[0] && x_[1] == x[1] && x_[2] == x[2]; }
  bool operator<(const tri_t& x)  const { return x_[0] < x[0] || (x_[0] == x[0] && (x_[1] < x[1] || (x_[1] == x[1] && x_[2] < x[2]))); }

  tri_t sorted() const { 
    tri_t t(*this); 
    if(t[0] > t[1]) std::swap(t[0],t[1]);
    if(t[1] > t[2]) std::swap(t[1],t[2]);
    if(t[0] > t[1]) std::swap(t[0],t[1]);
    return t;
  }

  friend ostream& operator<<(ostream& S, const tri_t& t){ S << vector<int>(t.x_,t.x_+3); return S; }
};


struct face_t : public vector<node_t> {
  face_t(const size_t size=0) : vector<node_t>(size) {}
  face_t(const vector<node_t>& vertices) : vector<node_t>(vertices) {}
  face_t(const tri_t& t) : vector<node_t>(t.x_,t.x_+3) {}
  
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
  // TODO: This seems to fail when compiled with -std=c++0x and -O0!
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

  face_t sorted() const {
    face_t f(*this);
    sort(f.begin(),f.end());
    return f;
  }
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
  const coord2d x;
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
    // ascending sorting of angles, implies CCW orientation, and sort requires 'less' for ascending sorting
    return angs <= angt; 	// TODO: Is the sign here correct? // lnw: changed, i think less/equal is correct
  }
};

#include "eisenstein.hh"


class polygon { // Given in CW order
public:
  vector<Eisenstein> outline;
  vector<Eisenstein> reduced_outline;

  polygon(const vector<Eisenstein>& outline) : outline(outline), reduced_outline(reduce()) {  }
  polygon(const vector<pair<int,int> >& outline) : outline(outline.begin(),outline.end()), reduced_outline(reduce()) {  }

  class scanline {
  public:
    int minY;
    vector< vector<int> > xs;
  };


  static vector<Eisenstein> draw_line(const Eisenstein& x0,const Eisenstein& x1); 
  scanline scanConvert() const;  

  pair<int,int> slope(int i,bool reduced=false) const;
  int turn_direction(int j,bool reduced=false) const;
  bool peak(int j,bool reduced=false) const;
  bool saddle(int j,bool reduced=false) const;

  double winding_number(const Eisenstein& x) const;
  bool point_inside(const Eisenstein& x) const;

  set<Eisenstein> allpoints() const;
  vector<Eisenstein> controlpoints() const;

  friend ostream& operator<<(ostream& S, const polygon& P);
private:
  vector<Eisenstein> reduce() const;
};

#endif
