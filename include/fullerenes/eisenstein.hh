#pragma once
#include <assert.h>
#include <stdlib.h>
#include <inttypes.h>
#include <cmath>
#include <algorithm>
#include <vector>

using namespace std;

/*
  Multiplication:
  (a,b) = a*1 + b*w -> (a,b)*(c,d) = (a*1+b*w)(c*1+d*w) = (ac*1 + (bc+ad)*w + b*d*(w-1) = (ac-bd)*1 + (bc+ad+bd)*w
  
  Inverse: (a+b,-b)/(a^2+ab+b^2) = complex_conj(a,b)/norm2(a,b)
  
*/
class Eisenstein: public pair<int,int> {
public:
  Eisenstein(int a=0, int b=0) : pair<int,int>(a,b) {}
  Eisenstein(const pair<int,int>& x) : pair<int,int>(x.first,x.second) {}
  Eisenstein(const pair<double,double>& x) : pair<int,int>(round(x.first-x.second/sqrt(3)), round(2*x.second/sqrt(3)))
  { }
  Eisenstein operator*(const int y) const { return Eisenstein(first*y,second*y); }
  Eisenstein operator/(const int y) const { return Eisenstein(first/y,second/y); }
  Eisenstein operator*(const Eisenstein& y) const { 
    int a = first, b = second, c = y.first, d = y.second;
    return Eisenstein(a*c-b*d, b*c+(a+b)*d); 
  }
  vector<Eisenstein> operator*(const vector<Eisenstein>& v) const {
    vector<Eisenstein> result(v.size());
    for(size_t i=0;i<v.size();i++) result[i] = (*this)*v[i];
    return result;
  }
  Eisenstein operator+(const Eisenstein& y) const { return Eisenstein(first+y.first,second+y.second); }
  Eisenstein operator-(const Eisenstein& y) const { return Eisenstein(first-y.first,second-y.second); }
  Eisenstein operator-()                    const { return Eisenstein(-first,-second); }   
  pair<double,double>    operator-(const pair<double,double>& y)    const { return pair<double,double>(first-y.first,second-y.second); } 
  pair<double,double>    operator+(const pair<double,double>& y)    const { return pair<double,double>(first+y.first,second+y.second); } 
  Eisenstein& operator*=(const Eisenstein& y) { return (*this = *this * y); }
  Eisenstein& operator+=(const Eisenstein& y) { first += y.first; second += y.second; return *this; }
  Eisenstein& operator-=(const Eisenstein& y) { first -= y.first; second -= y.second; return *this; }

  static Eisenstein unit[7];
  
  bool isUnit() const { return norm2() == 1;  }

  // w-reflection: conjugation in Z[w] in its native (1, w) basis,
  Eisenstein eis_conj()  const { return Eisenstein(first,-second);  }
  // i-reflection: complex conjugation in C, restricted to Z[w].
  Eisenstein complex_conj() const { return Eisenstein((first+second), -second); }
  Eisenstein GCtransform(int k, int l) const {  return Eisenstein(k,l) * (*this);  }
  Eisenstein affine(const Eisenstein& x0, const Eisenstein w) const {
    Eisenstein x(*this);
    //    cout << x0 << " + " << w << " * " << *this << " = " << (x0+(w*(*this))) <<endl;
    return x0+w*x;
  }
  //  
  // (-1,1)   \ /  (0,1)
  // (-1,0)  --x-- (1,0)
  // ( 0,-1)  / \  (1,-1)
  // 
  Eisenstein nextCW()    const { return (*this) * Eisenstein(1,-1); }
  Eisenstein nextCCW()   const { return (*this) * Eisenstein(0,1);  }
  Eisenstein transpose() const { return Eisenstein(second,first);   }

  int unit_angle() const {
    assert(norm2() == 1);
    switch(first*10+second){
    case  1*10 + 0: return 0;	
    case  0*10 + 1: return 1;	
    case -1*10 + 1: return 2;	
    case -1*10 + 0: return 3;  	
    case  0*10 - 1: return 4;
    case  1*10 - 1: return 5;	
    default:
      abort();
    }
  }

  int nearest_unit_angle() const {
    int l = first, k = second;
    if(l>=0 && k<=0 && -k < l) return 0;
    if(l>=0 && k>0)            return 1;
    if(l<0  && k>0 && -l <= k) return 2;
    if(l<0  && k>=0 && -l > k) return 3;
    if(l<=0 && k<0)            return 4;
    if(l>0 && k<=0 && -k >= l) return 5;
    return -1;
  }


  pair<double,double> coord() const { 
    return make_pair(first + second/2., sqrt(3/4.)*second);
  }
  int norm2() const { return first*first + first*second + second*second; }
  double norm() const { return sqrt(norm2()); }

  Eisenstein abs() const { return Eisenstein(::abs(first),::abs(second)); }

  Eisenstein div(const Eisenstein& y) const {
    // Naive, possibly non-robust algorithm
    Eisenstein z(*this * y.complex_conj());

    double re = z.first, im = z.second;
    re /= y.norm2(); im /= y.norm2();

    return Eisenstein(round(re),round(im));
  }

  Eisenstein mod(const Eisenstein& y) const {
    Eisenstein q(this->div(y));
    return (*this) - q*y;
  }

  static pair<double,double> average(const vector<Eisenstein>& xs)
  {
    double re = 0, im = 0;
    for(size_t i=0;i<xs.size();i++){ re += xs[i].first; im += xs[i].second; }
    re /= xs.size(); im /= xs.size();
    
    return {re,im};
  }

  static Eisenstein gcd(Eisenstein a, Eisenstein b)  {
    Eisenstein c;
    while ( a != Eisenstein(0,0) ) {
      c = a; a = b.mod(a);  b = c;
    }
    return b;
  }

  static Eisenstein gcd(const vector<Eisenstein>& xs) {
    if(xs.empty()) return Eisenstein(1,0);

    Eisenstein d(xs[0]);
    for(size_t i=1;i<xs.size();i++) d = gcd(xs[i],d);
    return d;
  }


  static int turn(const Eisenstein& a, const Eisenstein& b, const Eisenstein& c) {
    int x = (b.first-a.first)*(c.second-a.second) - (b.second-a.second)*(c.first-a.first);
    return x < 0 ? -1 : (x == 0 ? 0 : 1);
  }

};

namespace std {
  template<> struct hash<Eisenstein> { // Vectors of integers smaller than 32 bit
    size_t operator()(Eisenstein const &f) const {
      uint64_t combined_int = (uint64_t(f.first)<<32) + f.second;
      return std::hash<uint64_t>()(combined_int);
    }
  };
}


// =====================================================================
// Free-function helpers on Z[w].
//
// Things that are about Z[w] but read more naturally as binary or
// query operations than as methods on a single Eisenstein.
// =====================================================================

// d-th sixth-root of unity, with mod-6 wrap on d (handles d < 0 too).
// Eisenstein::unit[i] is the same value when 0 <= i <= 6 (where i=6 is
// a duplicate of i=0); use unit_direction when d may be out of range.
inline Eisenstein unit_direction(int d) {
  return Eisenstein::unit[((d % 6) + 6) % 6];
}

// Integer signed area in shear coords:  u.first*v.second - u.second*v.first.
// Sign matches the Cartesian wedge (Cartesian-wedge = (sqrt(3)/2) * this),
// so trigonometry-free orientation tests are exact.  Companion to
// Eisenstein::turn(a, b, c), which returns the SIGN of wedge(b-a, c-a).
inline long wedge(Eisenstein u, Eisenstein v) {
  return (long)u.first * v.second - (long)u.second * v.first;
}

// Some Eisenstein integer (a, b) with a >= 0, b >= 0 and
// a^2 + a*b + b^2 == N.  Returns the first sector-0 representative
// found by scanning b = 0, 1, ...; aborts if no solution exists.
// Precondition: N >= 0 is a valid Eisenstein norm.
Eisenstein eisenstein_of_norm(int N);

// Enumerate ALL sector-0 Eisenstein reps (a >= 0, b >= 0) of norm N.
// Generic norms return 1 entry; split-prime norms return 2 entries in
// distinct rotation orbits.
std::vector<Eisenstein> sector0_reps_of_norm(int N);


// =====================================================================
// D6 action between norm-equal Eisensteins.
// =====================================================================

// A D6 affine transform of Z[w] onto itself:
//   T(z) = unit * z                       (if reflect == false)
//        = unit * z.complex_conj()        (if reflect == true)
// `unit` is one of the 6 Eisenstein units.
struct D6Affine {
  Eisenstein unit;
  bool       reflect;

  Eisenstein apply(Eisenstein z) const {
    return (reflect ? z.complex_conj() : z) * unit;
  }
};

// Find the D6Affine T with T(z_from) == z_to.  Both inputs must have
// the same norm.  Exactly one branch (rotation, reflection) gives a
// unit by the D6 symmetry of Z[w]; align returns that branch.  Aborts
// on norm mismatch or non-divisibility.
D6Affine align(Eisenstein z_from, Eisenstein z_to);




