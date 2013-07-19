#ifndef EISENSTEIN_HH
# define EISENSTEIN_HH

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
  Eisenstein& operator+=(const Eisenstein& y) { first += y.first; second += y.second; return *this; }
  Eisenstein& operator-=(const Eisenstein& y) { first -= y.first; second -= y.second; return *this; }

  // invertn(a,b) * (a,b) == norm2() (1,0)
  Eisenstein invertn() const { return Eisenstein((first+second), -second); }
  Eisenstein GCtransform(int k, int l) const {  return Eisenstein(k,l) * (*this);  }
  //  
  // (-1,1)   \ /  (0,1)
  // (-1,0)  --x-- (1,0)
  // ( 0,-1)  / \  (1,-1)
  // 
  Eisenstein nextCW() const { return (*this) * Eisenstein(1,-1); }
  Eisenstein nextCCW() const { return (*this) * Eisenstein(0,1); }

  Eisenstein nextCWold() const {
    switch(second + 10*first){
    case  10+0: /*( 1 ,0)*/  return Eisenstein(1,-1);
    case  10-1: /*( 1,-1)*/  return Eisenstein(0,-1); 
    case   0-1: /*( 0,-1)*/  return Eisenstein(-1,0);
    case -10+0: /*(-1, 0)*/  return Eisenstein(-1,1);
    case -10+1: /*(-1, 1)*/  return Eisenstein(0,1);
    case   0+1: /*( 0, 1)*/  return Eisenstein(1,0);
    default:
      cerr << "nextCW(): " << *this << " is not a grid direction.\n";
      abort();
    }
 }

  coord2d coord() const { return coord2d(1,0)*first + coord2d(1/2.,sqrt(3/4.))*second; }

  int norm2() const { return first*first + first*second + second*second; }
  double norm() const { return sqrt(norm2()); }

  Eisenstein GCD(const Eisenstein &y) const;

};



#endif
