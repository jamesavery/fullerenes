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

// The apex of a flat triangle with base a=(0,0), b=(L_ab,0), given its three edge
// lengths (law of cosines); apex placed at y >= 0. The Euclidean analogue of
// place_third_eis (the Z[w] lattice version, eisenstein.hh).
inline coord2d place_apex(double L_ab, double L_bc, double L_ca) {
  const double x  = (L_ab*L_ab + L_ca*L_ca - L_bc*L_bc) / (2.0*L_ab);
  const double y2 = L_ca*L_ca - x*x;
  return coord2d(x, sqrt(y2 > 0 ? y2 : 0.0));
}



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

  // Symmetric part (m_ij + m_ji)/2.  Exact for symmetric input; for input that
  // is symmetric only up to roundoff (e.g. an X^T X product formed under FMA
  // contraction, whose off-diagonals can differ in the last ulp) it is the
  // well-defined symmetric representative.  The eigen* methods below operate on
  // this part, so they never require -- or assert -- exact symmetry of *this.
  matrix3d symmetric_part() const {
    const matrix3d &M(*this);
    matrix3d S;
    for(int i=0;i<3;i++)
      for(int j=0;j<3;j++)
        S(i,j) = 0.5*(M(i,j)+M(j,i));
    return S;
  }

  // Unnormalised null-space direction of a symmetric B = S - lambda*I: the
  // largest-magnitude cross product of two of its rows.  For a simple
  // eigenvalue (B has rank 2) this is the eigenvector direction; for a repeated
  // eigenvalue all row pairs are parallel and the returned norm collapses to ~0.
  static coord3d null_direction(const matrix3d& B) {
    const coord3d r0(B(0,0),B(0,1),B(0,2)), r1(B(1,0),B(1,1),B(1,2)), r2(B(2,0),B(2,1),B(2,2));
    const coord3d cr[3] = { r0.cross(r1), r0.cross(r2), r1.cross(r2) };
    coord3d best = cr[0];
    double best_n2 = best.norm2();
    for(int k=1;k<3;k++){ double n2 = cr[k].norm2(); if(n2 > best_n2){ best_n2 = n2; best = cr[k]; } }
    return best;
  }

  // Eigenvalue solver specialized to symmetric real 3x3 matrices, via the
  // trigonometric closed form applied to the DEVIATORIC part (Smith's method;
  // the same formulation as the device solver symMat3::eigenvalues in
  // sym-mat.cc).  Working on S - q*I with q = Tr(S)/3 avoids the catastrophic
  // cancellation the raw characteristic-polynomial route suffers for clustered
  // spectra: the trace part carries no error into the root separation, and the
  // middle root is recovered exactly from the trace identity.  Operates on the
  // symmetric part (see symmetric_part()), so it is correct for exactly-
  // symmetric input and for input asymmetric only by floating-point roundoff.
  // Returned in DESCENDING order: lambda[0] >= lambda[1] >= lambda[2].
  coord3d eigenvalues() const
  {
    const matrix3d S(symmetric_part());
    const double a(S(0,0)), b(S(0,1)), c(S(0,2)), d(S(1,1)), e(S(1,2)), f(S(2,2));

    const long double q  = (a + (long double)d + f)/3.L;         // Tr(S)/3
    const long double p1 = (long double)b*b + (long double)c*c + (long double)e*e;
    const long double p2 = (a-q)*(a-q) + (d-q)*(d-q) + (f-q)*(f-q) + 2.L*p1;
    if(p2 == 0) return coord3d(q, q, q);   // exactly isotropic: S = q*I

    const long double p = sqrtl(p2/6.L);
    // r = det((S - q*I)/p) / 2 lies in [-1,1] up to roundoff; clamp.
    const long double det_dev = -c*(long double)c*(d-q) + 2.L*(long double)b*c*e
                               - (a-q)*(long double)e*e - (long double)b*b*(f-q)
                               + (a-q)*(d-q)*(f-q);
    long double r = det_dev/(2.L*p*p*p);
    r = std::max(-1.0L, std::min(1.0L, r));

    const long double phi  = acosl(r)/3.L;                // phi in [0, pi/3]
    const long double lam0 = q + 2.L*p*cosl(phi);                    // largest
    const long double lam2 = q + 2.L*p*cosl(phi + 2.L*M_PI/3.L);     // smallest
    const long double lam1 = 3.L*q - lam0 - lam2;                    // trace identity
    return coord3d(lam0, lam1, lam2);
  }
  

  // Unit eigenvector for a SIMPLE eigenvalue lambda of the symmetric part:
  // the normalised largest-norm cross product of two rows of (S - lambda*I).
  // A repeated (degenerate) eigenvalue has a >= 2-dimensional eigenspace in
  // which no single vector is distinguished; this returns the zero vector as a
  // documented sentinel -- use eigensystem() to obtain a full orthonormal
  // eigenbasis, which handles degeneracy directly.
  coord3d eigenvector(const double lambda) const {
    const matrix3d S(symmetric_part());
    matrix3d B(S);
    for(int i=0;i<3;i++) B(i,i) -= lambda;

    const coord3d n(null_direction(B));
    const double nn = n.norm(), scale = S.norm();
    if(nn > 1e-9*(scale>0?scale*scale:1.0)) return n/nn;
    return coord3d();
  }

  // Eigen-decomposition of the symmetric part (see eigenvalues()).  Returns
  // (lambda, C) where lambda holds the eigenvalues sorted by ABSOLUTE VALUE,
  // smallest first (PolyhedronView::principal_axes and other callers rely on
  // this ordering), and row i of C is the unit eigenvector for lambda[i].
  //
  // The eigenvectors are always orthonormal (C is orthogonal).  Degeneracy is
  // handled by extracting the eigenvector of the most isolated eigenvalue first
  // and completing an orthonormal basis in its 2D orthogonal complement via a
  // 2x2 sub-problem (which gracefully covers a two-fold degenerate eigenvalue);
  // the isotropic case -- all eigenvalues equal, e.g. a multiple of the
  // identity or the zero matrix -- returns the standard basis.
  //
  // The closed-form eigenvalues (eigenvalues(), Smith form) seed the
  // eigenvector extraction; the returned eigenvalues are nevertheless
  // recovered from the eigenvectors -- one step of Rayleigh-quotient
  // iteration refines the isolated eigenvalue, and the other two come directly
  // from the 2x2 restriction of S to the complementary plane -- so the result
  // does not inherit even the (small) closed-form rounding.
  pair<coord3d,matrix3d> eigensystem() const {
    const matrix3d S(symmetric_part());
    const coord3d  lam0(eigenvalues());   // closed-form seeds (descending)

    const double scale  = S.norm();
    const double fscale = (scale>0 ? scale : 1.0);

    double  lam[3];
    coord3d vec[3];
    bool    resolved = false;

    // Isotropy is detected from the MATRIX, not from lam0: S is a multiple of I
    // iff its deviatoric part vanishes.  (The closed-form eigenvalues carry a
    // spurious O(sqrt(eps)) spread for a clustered spectrum, so a gap test on
    // them would misclassify an isotropic matrix and then divide by a collapsed
    // null direction.)  Every direction is then an eigenvector.
    const double trace = S(0,0)+S(1,1)+S(2,2);
    matrix3d dev(S); for(int i=0;i<3;i++) dev(i,i) -= trace/3.0;

    if(dev.norm() > 1e-12*fscale){
      // Most isolated eigenvalue: largest minimum gap to the other two.  Its
      // eigenspace is 1-dimensional; the seed only needs to pick it out.
      const double gap[3] = {
        std::min(fabs(lam0[0]-lam0[1]), fabs(lam0[0]-lam0[2])),
        std::min(fabs(lam0[1]-lam0[0]), fabs(lam0[1]-lam0[2])),
        std::min(fabs(lam0[2]-lam0[0]), fabs(lam0[2]-lam0[1]))
      };
      int iso = 0;
      if(gap[1] > gap[iso]) iso = 1;
      if(gap[2] > gap[iso]) iso = 2;

      // Isolated eigenvector: extract, refine the eigenvalue as v.(S v), then
      // sharpen with a second extraction -- but only if the refined eigenvalue
      // has not made (S - li*I) collapse to zero (which happens when li lands
      // exactly on the eigenvalue), in which case the first vi is already exact.
      double li = lam0[iso];
      matrix3d B(S); for(int k=0;k<3;k++) B(k,k) -= li;
      coord3d n(null_direction(B));
      const double n0 = n.norm();
      coord3d vi(n/n0);
      li = vi.dot(S*vi);

      matrix3d B1(S); for(int k=0;k<3;k++) B1(k,k) -= li;
      const coord3d n1(null_direction(B1));
      if(n1.norm() > 1e-8*n0){ vi = n1/n1.norm(); li = vi.dot(S*vi); }

      // Orthonormal basis {u,w} of the plane orthogonal to vi.
      const coord3d ax =
        (fabs(vi[0]) <= fabs(vi[1]) && fabs(vi[0]) <= fabs(vi[2])) ? coord3d(1,0,0)
        : (fabs(vi[1]) <= fabs(vi[2]) ? coord3d(0,1,0) : coord3d(0,0,1));
      coord3d u(vi.cross(ax)); u /= u.norm();
      const coord3d w(vi.cross(u));   // unit: vi and u are orthonormal

      // Diagonalise the symmetric 2x2 restriction [[a2,b2],[b2,d2]] of S to the
      // plane.  Its eigenvalues are the (accurate) remaining eigenvalues; the
      // rotation is the identity when they are degenerate.
      const coord3d Su(S*u), Sw(S*w);
      const double a2 = u.dot(Su), b2 = u.dot(Sw), d2 = w.dot(Sw);
      const double theta = 0.5*atan2(2.0*b2, a2-d2);
      const double cs = cos(theta), sn = sin(theta);

      const int p = (iso+1)%3, q = (iso+2)%3;
      lam[iso] = li;                              vec[iso] = vi;
      lam[p]   = a2*cs*cs + 2.0*b2*cs*sn + d2*sn*sn; vec[p] = u*cs + w*sn;
      lam[q]   = a2*sn*sn - 2.0*b2*cs*sn + d2*cs*cs; vec[q] = u*(-sn) + w*cs;
      resolved = true;
    }

    if(!resolved){
      // Isotropic: standard basis; eigenvalues are the (equal) diagonal.
      const coord3d e[3] = { coord3d(1,0,0), coord3d(0,1,0), coord3d(0,0,1) };
      for(int i=0;i<3;i++){ vec[i] = e[i]; lam[i] = e[i].dot(S*e[i]); }
    }

    // Sort (eigenvalue, eigenvector) pairs by absolute value, smallest first.
    for(int a=0;a<2;a++)
      for(int b=a+1;b<3;b++)
        if(fabs(lam[a]) > fabs(lam[b])){ std::swap(lam[a],lam[b]); std::swap(vec[a],vec[b]); }

    matrix3d C;
    for(int i=0;i<3;i++)
      for(int j=0;j<3;j++) C(i,j) = vec[i][j];
    return make_pair(coord3d(lam[0],lam[1],lam[2]), C);
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

// Convex hull of a 3D point set as a triangle list (indices into pts), outward-
// oriented normals. Incremental construction; empty result for degenerate input
// (< 4 points or all coplanar). Requires >= 4 points not all coplanar for a
// non-empty result; degenerate input returns an empty vector (never throws). Impl
// in deltahedron.cc (the hull the deltahedron pipeline's project_onto_convex_hull
// builds, exported for warm-start convexification in the research sub-projects).
std::vector<std::array<int,3>> convex_hull_tris(std::span<const coord3d> pts);

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

  bool point_inside(const Eisenstein& x) const;
  bool point_on_periphery(const Eisenstein& x) const;
  bool point_included(const Eisenstein& x) const; // Inside or on periphery
  
  vector<Eisenstein> allpoints() const;
  vector<Eisenstein> controlpoints() const;

  friend ostream& operator<<(ostream& S, const polygon& P);
private:
  vector<Eisenstein> reduce() const;
};


