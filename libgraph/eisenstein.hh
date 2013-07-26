#ifndef EISENSTEIN_HH
# define EISENSTEIN_HH

#include <stdlib.h>
#include "geometry.hh"
using namespace std;

/*
  Multiplication:
  (a,b) = a*1 + b*w -> (a,b)*(c,d) = (a*1+b*w)(c*1+d*w) = (ac*1 + (bc+ad)*w + b*d*(w-1) = (ac-bd)*1 + (bc+ad+bd)*w
  
  Inverse: (a+b,-b)/(a^2+ab+b^2)
  
*/
class Eisenstein: public pair<int,int> {
public:
  Eisenstein(int a=0, int b=0) : pair<int,int>(a,b) {}
  Eisenstein(const pair<int,int>& x) : pair<int,int>(x.first,x.second) {}
  Eisenstein(const coord2d& x) : pair<int,int>(round(x.first-x.second/sqrt(3)), round(2*x.second/sqrt(3)))
  { }
  Eisenstein operator*(const int y) const { return Eisenstein(first*y,second*y); }
  Eisenstein operator/(const int y) const { return Eisenstein(first/y,second/y); }
  Eisenstein operator*(const Eisenstein& y) const { 
    int a = first, b = second, c = y.first, d = y.second;
    return Eisenstein(a*c-b*d, b*c+(a+b)*d); 
  }
  vector<Eisenstein> operator*(const vector<Eisenstein>& v) const {
    vector<Eisenstein> result(v.size());
    for(int i=0;i<v.size();i++) result[i] = (*this)*v[i];
    return result;
  }
  Eisenstein operator+(const Eisenstein& y) const { return Eisenstein(first+y.first,second+y.second); }
  Eisenstein operator-(const Eisenstein& y) const { return Eisenstein(first-y.first,second-y.second); } 
  coord2d    operator-(const coord2d& y)    const { return coord2d(first-y.first,second-y.second); } 
  coord2d    operator+(const coord2d& y)    const { return coord2d(first+y.first,second+y.second); } 
  Eisenstein& operator+=(const Eisenstein& y) { first += y.first; second += y.second; return *this; }
  Eisenstein& operator-=(const Eisenstein& y) { first -= y.first; second -= y.second; return *this; }

  // invertn(a,b) * (a,b) == norm2() (1,0)
  Eisenstein invertn() const { return Eisenstein((first+second), -second); }
  Eisenstein GCtransform(int k, int l) const {  return Eisenstein(k,l) * (*this);  }
  Eisenstein affine(const Eisenstein& x0, const Eisenstein w) const {
    //    cout << x0 << " + " << w << " * " << *this << " = " << (x0+(w*(*this))) <<endl;
    return x0+(w*(*this));
  }
  //  
  // (-1,1)   \ /  (0,1)
  // (-1,0)  --x-- (1,0)
  // ( 0,-1)  / \  (1,-1)
  // 
  Eisenstein nextCW() const { return (*this) * Eisenstein(1,-1); }
  Eisenstein nextCCW() const { return (*this) * Eisenstein(0,1); }

  coord2d coord() const { return coord2d(1,0)*first + coord2d(1/2.,sqrt(3/4.))*second; }

  int norm2() const { return first*first + first*second + second*second; }
  double norm() const { return sqrt(norm2()); }

  Eisenstein abs() const { return Eisenstein(::abs(first),::abs(second)); }

  Eisenstein div(const Eisenstein& y) const {
    // Naive, possibly non-robust algorithm
    Eisenstein z(*this * y.invertn());
    coord2d zf(z.first,z.second);
    zf *= 1.0/y.norm2();

    return Eisenstein(round(zf.first),round(zf.second));
  }

  Eisenstein mod(const Eisenstein& y) const {
    Eisenstein q(this->div(y));
    return (*this) - q*y;
  }

  static coord2d average(const vector<Eisenstein>& xs)
  {
    coord2d sum(0,0);
    for(int i=0;i<xs.size();i++) sum += coord2d(xs[i].first,xs[i].second);
    return sum / xs.size(); 
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
    for(int i=1;i<xs.size();i++) d = gcd(xs[i],d);
    return d;
  }
};



#endif
