#pragma once

#include <string.h>
#include <iostream>
#include <vector>
#include <set>
#include <map>
#include <math.h>
#include <sstream>
#include <list>
#include <algorithm>
#include <array>
#include <span>
#include "auxiliary.hh"

using namespace std;

// node_t and neighbours_t defined in auxiliary.hh (included above)
typedef vector< bool > edges_t;

template<typename T> struct coord2; // forward declaration for coord3<T>::polar_angle
struct matrix3d;                   // required for coord3<T>::outer(coord3<T>)

// TODO: geometry.hh is getting huge. Move most of the implementation to geometryc.cc

namespace std {
  template<>
  struct hash<edge_t> {
    std::size_t operator()(const edge_t &p) const {
      std::size_t seed1(0);
      ::hash_combine(seed1, p.first);
      ::hash_combine(seed1, p.second);

      std::size_t seed2(0);
      ::hash_combine(seed2, p.second);
      ::hash_combine(seed2, p.first);

      return std::min(seed1, seed2);
    }
  };
}


template<typename T = double>
struct coord2 : public pair<T,T> {
  coord2(const pair<T,T>& p) : pair<T,T>(p) {}
  coord2(const T x=0, const T y=0) : pair<T,T>(x,y) {}
  coord2 operator/(const T s)        const { return coord2(this->first/s, this->second/s); }
  coord2 operator*(const T s)        const { return coord2(this->first*s, this->second*s); }
  coord2 operator*(const coord2& y)  const { return coord2(this->first*y.first, this->second*y.second); }
  coord2 operator+(const coord2& y)  const { return coord2(this->first+y.first, this->second+y.second); }
  coord2 operator-(const coord2& y)  const { return coord2(this->first-y.first, this->second-y.second); }
  coord2& operator+=(const coord2& y){ this->first += y.first; this->second += y.second; return *this; }
  coord2& operator-=(const coord2& y){ this->first -= y.first; this->second -= y.second; return *this; }
  coord2& operator*=(const coord2& y){ this->first *= y.first; this->second *= y.second; return *this; }
  coord2& operator*=(const T s)      { this->first*=s; this->second*=s; return *this; }
  coord2 operator-() const { return coord2(-this->first, -this->second); }

  T operator()(int i) const { return i? this->second : this->first; }
  T& operator()(int i)      { return i? this->second : this->first; }

  T dot(const coord2& y) const { return this->first*y.first + this->second*y.second; }

  // CCW angle between two lines ([-pi;pi])
  double line_angle(const coord2& v) const {
    double angle = fmod(atan2((double)(this->first*v.second - v.first*this->second),
                              (double)(this->first*v.first  + this->second*v.second))
                        + 2*M_PI, 2*M_PI) - M_PI;
    return angle;
  }

  double point_angle(const coord2& y=0, const bool periodic=false) const {
    if(periodic) return point_angle_periodic(y);
    else {
      const coord2 vx(*this-y);
      return atan2((double)vx.second, (double)vx.first);
    }
  }

  // CCW angle of y around point on a periodic surface [0;2pi[ x [0;pi[.
  double point_angle_periodic(const coord2& y=0) const {
    const coord2 dy[4] = {coord2(0,0), coord2((T)(2*M_PI),0), coord2(0,(T)M_PI), coord2((T)(2*M_PI),(T)M_PI)};
    const coord2& self(*this);

    double rsqrmin = INFINITY;
    int imin = 0;
    for(int i=0;i<4;i++){
      const coord2 y0(y+dy[i]-self);
      const double rsqr = (double)y0.dot(y0);
      if(rsqr < rsqrmin){ rsqrmin = rsqr; imin = i; }
    }
    return point_angle(y+dy[imin],false);
  }

  T    norm()  const { return (T)sqrt((double)norm2()); }
  T    norm2() const { return this->first*this->first + this->second*this->second; }
  static coord2 dnorm(const coord2& x){ return x/x.norm(); }

  static coord2 displacement(const coord2& x, const coord2& y, bool layout_is_spherical)
  {
    if(!layout_is_spherical) return x-y;
    int i0=0, j0=0;
    double dmin = INFINITY;
    for(int i=0;i<=1;i++)
      for(int j=0;j<=1;j++){
        double d = (x + coord2((T)(i*2*M_PI),(T)(j*M_PI)) - y).norm();
        if(d < dmin){ i0 = i; j0 = j; }
      }
    return x + coord2((T)(i0*2*M_PI),(T)(j0*M_PI)) - y;
  }

  friend ostream& operator<<(ostream &s, const coord2<T>& x){
    s << std::fixed << LIST_OPEN << x.first << "," << x.second << LIST_CLOSE; return s; }
  friend istream& operator>>(istream &s, coord2<T>& x){ s >> x.first; s >> x.second; return s; }
};

using coord2d = coord2<double>;



template<typename T = double>
struct coord3 {
  T x[3];

  coord3() { x[0] = 0; x[1] = 0; x[2] = 0; }
  coord3(const T y[3]) { x[0] = y[0]; x[1] = y[1]; x[2] = y[2]; }
  coord3(const T x_, const T y, const T z) { x[0] = x_; x[1] = y; x[2] = z; }
  coord3 operator/(const T s)    const { return coord3(*this) /= s; }
  coord3 operator*(const T s)    const { return coord3(*this) *= s; }
  coord3 operator*(const coord3& y) const { return coord3(*this) *= y; }
  coord3 operator+(const coord3& y) const { return coord3(*this) += y; }
  coord3 operator-(const coord3& y) const { return coord3(*this) -= y; }
  coord3& operator+=(const coord3& y){ x[0] += y[0]; x[1] += y[1]; x[2] += y[2]; return *this; }
  coord3& operator-=(const coord3& y){ x[0] -= y[0]; x[1] -= y[1]; x[2] -= y[2]; return *this; }
  coord3& operator*=(const coord3& y){ x[0] *= y[0]; x[1] *= y[1]; x[2] *= y[2]; return *this; }
  coord3& operator*=(const T& y){ x[0] *= y; x[1] *= y; x[2] *= y; return *this; }
  coord3& operator/=(const T& y){ x[0] /= y; x[1] /= y; x[2] /= y; return *this; }
  coord3 operator-() const { return coord3(-x[0],-x[1],-x[2]); }
  coord3 operator*(const matrix3d& m) const;

  coord3 cross(const coord3& y) const {
    return coord3(x[1]*y[2]-x[2]*y[1], x[2]*y[0]-x[0]*y[2], x[0]*y[1]-x[1]*y[0]);
  }
  matrix3d outer(const coord3& y) const;
  T dot(const coord3& y)  const { return x[0]*y[0]+x[1]*y[1]+x[2]*y[2]; }
  T norm2()               const { return dot(*this); }                       // squared length
  T norm()                const { return (T)sqrt((double)dot(*this)); }
  coord2<T> polar_angle(const coord3& centre = coord3()) const
  {
    const coord3 rel(*this - centre);
    const T r = norm();  // norm of *this (preserved from original)
    return coord2<T>((T)acos((double)rel[2]/r), (T)atan2((double)rel[1],(double)rel[0]));
  }

  T& operator[](unsigned int i)       { return x[i]; }
  T  operator[](unsigned int i) const { return x[i]; }

  static T dist(const coord3& a, const coord3& b){ return (a-b).norm(); }
  static coord3 dnorm(const coord3& v){ return v/v.norm(); }
  // d^2/(dx_i dx_j) ||x|| = -x_i x_j/||x||^3 + [i==j]/||x||  (double precision)
  static void ddnorm(const coord3& v, vector<double>& H)
  {
    const double n = 1.0/v.norm(), n3 = n*n*n;
    for(int i=0;i<3;i++)
      for(int j=0;j<3;j++)
        H[i*3+j] = -(double)v[i]*(double)v[j]*n3 + (i==j? n : 0);
  }

  // calculation of the angle beta at b(0,0,0)
  static double angle(const coord3& a, const coord3& c);
  // calculation of the derivative of angle beta at b(0,0,0)
  static void dangle(const coord3& a, const coord3& c, coord3& da, coord3& dc);
  // dihedral angle theta at a(0,0,0), b, c and d  [-pi; pi]
  static double dihedral(const coord3& b, const coord3& c, const coord3& d);
  // derivative of dihedral angle
  static void ddihedral(const coord3& b, const coord3& c, const coord3& d,
                        coord3& db, coord3& dc, coord3& dd);
  static double ideal_dihedral(const int lA, const int lB, const int lC,
                               const double ur=1.0, const double us=1.0, const double ut=1.0);

  friend vector<coord3<T>>& operator-=(vector<coord3<T>>& xs, const coord3<T>& y)
  {
    for(size_t i=0;i<xs.size();i++) xs[i] -= y;
    return xs;
  }

  friend vector<coord3<T>>& operator*=(vector<coord3<T>>& xs, const T& y)
  {
    for(size_t i=0;i<xs.size();i++) xs[i] *= y;
    return xs;
  }

  static coord3 line_plane_intersect(const coord3& x0, const coord3& x1,
                                     const coord3& X0, const coord3& n)
  {
    const coord3 dx(x1-x0);
    T t = (X0-x0).dot(n)/dx.dot(n);
    return x0+dx*t;
  }

  friend ostream& operator<<(ostream &s, const coord3<T>& v){
    s << std::fixed << LIST_OPEN << v[0] << "," << v[1] << "," << v[2] << LIST_CLOSE; return s; }
  friend istream& operator>>(istream &s, coord3<T>& v){
    for(int i=0;i<3;i++){ s >> v[i]; } return s; }
};

using coord3d = coord3<double>;

struct matrix2d {
  double A[4];
  matrix2d(const double *v) { for(int i=0;i<4;i++) A[i] = v[i]; }
  explicit matrix2d(const double r=0, const double s=0, const double t=0, const double u=0) : A{r,s,t,u} { }
  double& operator()(int i, int j)      { return A[i*2+j]; }
  double operator()(int i, int j) const { return A[i*2+j]; }

  //  coord2d operator*(const coord2d& x){ return {A[0*2+0]*x(0) + A[1*2+0]*x(1), A[0*2+1]*x(0)+A[1*2+1]*x(1)}; }
  coord2d  operator*(const coord2d& x) const { return {A[0*2+0]*x(0)+A[0*2+1]*x(1),
						      A[1*2+0]*x(0)+A[1*2+1]*x(1)}; }
  matrix2d operator*(const matrix2d& B) const {
    return matrix2d(A[0*2+0]*B(0,0) + A[0*2+1]*B(1,0), A[0*2+0]*B(0,1) + A[0*2+1]*B(1,1),
		    A[1*2+0]*B(0,0) + A[1*2+1]*B(1,0), A[1*2+0]*B(0,1) + A[1*2+1]*B(1,1));
  }
  
  static matrix2d rotation(double th){
    double s = sin(th), c = cos(th);
    return matrix2d{c,-s,
	            s, c};
  }
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

    for(size_t i=0;i<xs.size();i++) ys[i] = A*xs[i];
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
      return coord3d(0,(-B-Disc)/(2.L*A),(-B+Disc)/(2.L*A));
    }

    // Depress characteristic polynomial - see http://en.wikipedia.org/wiki/Cubic_equation#Reduction_to_a_depressed_cubic
    long double 
      p  = (3.L*A*C - B*B)/(3.L*A*A),
      q  = (2.L*B*B*B - 9.L*A*B*C + 27.L*A*A*D)/(27.L*A*A*A),
      xc = B/(3.L*A);

    // François Viète's solution to cubic polynomials with three real roots. 
    coord3d t;
    long double cos_arg = (3.L*q)/(2.L*p)*sqrtl(-3.L/p);
    cos_arg = std::max(-1.0L, std::min(1.0L, cos_arg)); // clamp for FP robustness
    long double K = 2*sqrtl(-p/3.L),
                theta0 = (1.L/3.L)*acosl(cos_arg);
    for(int k=0;k<3;k++) t[k] = K*cosl(theta0-k*2.L*M_PI/3.L);

    // lambda = t - B/(3A)
    return t - coord3d(xc,xc,xc);
  }
  

  coord3d eigenvector(const double lambda) const {
    const matrix3d &M(*this);
    coord3d x;
  
    // using the first two eqs
    // [ a_12 * a_23 - a_13 * (a_22 - r) ]
    // [ a_12 * a_13 - a_23 * (a_11 - r) ]
    // [ (a_11 - r) * (a_22 - r) - a_12^2 ]
    x = coord3d( M(0,1)*M(1,2) - M(0,2)*(M(1,1)-lambda),
                 M(0,1)*M(0,2) - M(1,2)*(M(0,0)-lambda),
                 (M(0,0)-lambda)*(M(1,1)-lambda) - M(0,1)*M(0,1) );
    if (x.norm() / (M(0,0) + M(1,1) + M(2,2)) > 1.e-12) // not zero-ish
      return x/x.norm();
  
    // using the first+last eqs
    // [ a_12 * (a_33 - r) - a_13 * a_23 ]
    // [ a_13^2 - (a_11 - r) * (a_33 - r) ]
    // [ a_23 * (a_11 - r) - a_12 * a_13 ]
    x = coord3d( M(0,1)*(M(2,2)-lambda) - M(0,2)*M(1,2),
                 M(0,2)*M(0,2) - (M(0,0)-lambda)*(M(2,2)-lambda),
                 M(1,2)*(M(0,0)-lambda) - M(0,1)*M(0,2) );
    if (x.norm() / (M(0,0) + M(1,1) + M(2,2)) > 1.e-12) // not zero-ish
      return x/x.norm();
  
    // using the last two eqs
    // [ a_23^2 - (a_22 - r) * (a_33 - r) ]
    // [ a_12 * (a_33 - r) - a_13 * a_23 ]
    // [ a_13 * (a_22 - r) - a_12 * a_23 ]
    x = coord3d( M(1,2)*M(1,2) - (M(1,1)-lambda)*(M(2,2)-lambda),
                 M(0,1)*(M(2,2)-lambda) - M(0,2)*M(1,2),
                 M(0,2)*(M(1,1)-lambda) - M(0,1)*M(1,2) );
    if (x.norm() / (M(0,0) + M(1,1) + M(2,2)) > 1.e-12) // not zero-ish
      return x/x.norm();
  
    cerr << "something is very wrong, possibly two degenerate evals" << std::endl;
    return coord3d();
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

  // Skew-symmetric cross product matrix: cross_matrix(a) * b == a.cross(b)
  static matrix3d cross_matrix(const coord3d& a) {
    return matrix3d(   0, -a[2],  a[1],
                    a[2],     0, -a[0],
                   -a[1],  a[0],     0);
  }

  static matrix3d unit_matrix(){
    return matrix3d(1,0,0,0,1,0,0,0,1);
  }
  
  friend ostream& operator<<(ostream& S, const matrix3d &M)
  {
    S << LIST_OPEN; for(int i=0;i<3;i++) S << vector<double>(&M.values[i*3],&M.values[(i+1)*3]) << (i+1<3?',':LIST_CLOSE);
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

  coord3d centroid(std::span<const coord3d> points) const { return (points[u(0)]+points[u(1)]+points[u(2)])/3.0; }
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

// TODO: Hash functions gathered in single file?
namespace std {
  template<> struct hash<tri_t> { // Vectors of integers smaller than 32 bit
    size_t operator()(const tri_t &t) const {
      size_t seed(0);
      hash_combine(seed, t[0]);
      hash_combine(seed, t[1]);
      hash_combine(seed, t[2]);
      return seed;
    }
  };
}


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
  coord3d centroid(std::span<const coord3d> layout) const {
    coord3d c;
    for(size_t i=0;i<size();i++) c += layout[(*this)[i]];
    return c/size();
  }

  // http://www.softsurfer.com/Archive/algorithm_0103/algorithm_0103.htm#Winding%20Number
  // TODO: This seems to fail when compiled with -std=c++0x and -O0!
  double winding_number(const vector<coord2d>& layout, const coord2d& x) const {
    vector<coord2d> Cp(size());
    for(size_t i=0;i<size();i++)
      Cp[i] = layout[(*this)[i]]-x;

    double wn = 0;
    for(size_t i=0;i<size();i++){
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
  bool contains(const node_t v) const { for(size_t i=0;i<size();i++) if(v == (*this)[i]) return true; return false; }

  // Unique representation for sets and maps.
  // To preserve orientation: don't sort, but rotate so as to start at smallest node.
  // NB: Not necessary because of == and < operations?
  face_t normalized() const {
    face_t f(*this);

    // Find smallest node
    node_t i_min = 0;
    for(size_t i=0;i<f.size();i++) if(f[i] < f[i_min]) i_min = i;
    // Rotate so smallest node is first
    rotate(f.begin(),f.begin()+i_min,f.end());
    return f;
  }

  // The directed edge with minimal start vertex gives a unique representation of an oriented face
  arc_t minimal_edge() const {
    face_t f(*this);
    
    arc_t e_min{f[0],f[1]};
    for(size_t i=1;i<f.size();i++) if(f[i] < e_min.first) e_min = {f[i],f[(i+1)%f.size()]};

    return e_min;
  }
};

namespace std {
  template<> struct hash<face_t> { // Vectors of integers smaller than 32 bit
    size_t operator()(face_t const &f) const {
      u32string s(f.begin(),f.end());
      return std::hash<u32string>()(s);      
    }
  };
}


typedef map<unsigned int,set<face_t> > facemap_t; // Never really used -- candidate for retirement

struct Tri3D {
  coord3d a,b,c,u,v,n;
  typedef pair<coord3d,coord3d> segment_t;

  Tri3D(const coord3d *T) : a(T[0]), b(T[1]), c(T[2]), u(b-a), v(c-a), n(u.cross(v)) {}
  Tri3D(const coord3d& a,const coord3d& b,const coord3d& c) : a(a), b(b), c(c), u(b-a), v(c-a), n(u.cross(v)) {}
  double plane_intersection(const segment_t& segment) const { 
    return n.dot(a-segment.first)/n.dot(segment.second-segment.first);
  }
  Tri3D(std::span<const coord3d> points, const tri_t& t) : a(points[t[0]]),b(points[t[1]]),c(points[t[2]]), u(b-a), v(c-a), n(u.cross(v)) {}

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

  void flip(){ coord3d t = b; b = c; c = t; t = u; u = v; v = t; n = n*-1; }
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
    s << LIST_OPEN << T.a << "," << T.b << "," << T.c << LIST_CLOSE;
    return s;
  }
};

// The closest point of a triangle to a query point: its 3D position, barycentric
// coordinates, and squared distance to the query.
struct ClosestPoint { coord3d pos; std::array<double,3> bary; double dist2; };

// Nearest point of triangle t (vertices points[t[0..2]]) to p -- the metric
// projection of p onto the filled triangle: the foot of the perpendicular if it
// lands inside, else the nearest point on an edge or vertex. dist2 = |pos-p|^2,
// and is +infinity for a degenerate (zero-area) face so it is never selected as
// the nearest among several.
ClosestPoint closest_point_on_triangle(const coord3d& p, const tri_t& t,
                                        std::span<const coord3d> points);

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

  polygon operator*(Eisenstein x) const;
  polygon operator+(Eisenstein x) const;
  
  static vector<Eisenstein> draw_line(const Eisenstein& x0,const Eisenstein& x1);
  scanline scanConvert() const;  

  pair<int,int> slope(int i,bool reduced=false) const;
  int turn_direction(int j,bool reduced=false) const;
  bool peak(int j,bool reduced=false) const;
  bool saddle(int j,bool reduced=false) const;

  double winding_number(const Eisenstein& x) const;
  bool point_inside(const Eisenstein& x) const;
  bool point_on_periphery(const Eisenstein& x) const;
  bool point_included(const Eisenstein& x) const; // Inside or on periphery
  
  vector<Eisenstein> allpoints() const;
  vector<Eisenstein> controlpoints() const;

  friend ostream& operator<<(ostream& S, const polygon& P);
private:
  vector<Eisenstein> reduce() const;
};


