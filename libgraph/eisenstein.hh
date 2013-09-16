#ifndef EISENSTEIN_HH
# define EISENSTEIN_HH

#include <assert.h>
#include <stdlib.h>
#include <algorithm>
#include <complex>
#include <vector>

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
    for(int i=0;i<v.size();i++) result[i] = (*this)*v[i];
    return result;
  }
  Eisenstein operator+(const Eisenstein& y) const { return Eisenstein(first+y.first,second+y.second); }
  Eisenstein operator-(const Eisenstein& y) const { return Eisenstein(first-y.first,second-y.second); } 
  pair<double,double>    operator-(const pair<double,double>& y)    const { return pair<double,double>(first-y.first,second-y.second); } 
  pair<double,double>    operator+(const pair<double,double>& y)    const { return pair<double,double>(first+y.first,second+y.second); } 
  Eisenstein& operator*=(const Eisenstein& y) { return (*this = *this * y); }
  Eisenstein& operator+=(const Eisenstein& y) { first += y.first; second += y.second; return *this; }
  Eisenstein& operator-=(const Eisenstein& y) { first -= y.first; second -= y.second; return *this; }

  bool isUnit() const {
    Eisenstein unit(1,0);
    for(int i=0;i<6;i++,unit = unit.nextCW()) if((*this) == unit) return true;
    return false;
  }

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
  Eisenstein transpose() const { return Eisenstein(second,first); }
  Eisenstein conj() const { return Eisenstein(first,-second); }

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
    Eisenstein z(*this * y.invertn());

    complex<double> zf(z.first,z.second);
    zf /= y.norm2();

    return Eisenstein(round(zf.real()),round(zf.imag()));
  }

  Eisenstein mod(const Eisenstein& y) const {
    Eisenstein q(this->div(y));
    return (*this) - q*y;
  }

  static pair<double,double> average(const vector<Eisenstein>& xs)
  {
    complex<double> avg(0,0);
    for(int i=0;i<xs.size();i++) avg += complex<double>(xs[i].first,xs[i].second);
    avg /= xs.size(); 
    return make_pair(avg.real(),avg.imag());
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


  static int turn(const Eisenstein& a, const Eisenstein& b, const Eisenstein& c) {
    int x = (b.first-a.first)*(c.second-a.second) - (b.second-a.second)*(c.first-a.first);
    return x < 0 ? -1 : (x == 0 ? 0 : 1);
  }

};




#endif
