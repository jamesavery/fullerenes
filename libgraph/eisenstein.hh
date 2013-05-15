#ifndef EISENSTEIN_HH
# define EISENSTEIN_HH

#include "geometry.hh"
using namespace std;

class Eisenstein: public pair<int,int> {
public:
  Eisenstein(int a=0, int b=0) : pair<int,int>(a,b) {}
  Eisenstein(const coord2d& x) : pair<int,int>(round(x.first-x.second/sqrt(3)), round(2*x.second/sqrt(3)))
  { }
  Eisenstein operator*(const Eisenstein& y) const { return Eisenstein(first*y.first,second*y.second); }
  Eisenstein operator+(const Eisenstein& y) const { return Eisenstein(first+y.first,second+y.second); }
  Eisenstein operator-(const Eisenstein& y) const { return Eisenstein(first-y.first,second-y.second); } 
  Eisenstein& operator+=(const Eisenstein& y) { first += y.first; second += y.second; return *this; }
  Eisenstein& operator-=(const Eisenstein& y) { first -= y.first; second -= y.second; return *this; }

  Eisenstein GCtransform(int k, int l) const {
    return Eisenstein(k*first - l*second, l*first + (k+l)*second);
  }
  //  
  // (-1,1)   \ /  (1,1)
  // (-1,0)  --x-- (1,0)
  // (-1,-1)  / \  (1,-1)
  // 
  Eisenstein nextCW() const {
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

  static vector<int> rasterize_line(const Eisenstein &x0, const Eisenstein &x1)
  {
#define ROUND(x) (round(x*1e4)*1e-4)
    const int 
      dx = x1.first - x0.first,
      dy = x1.second - x0.second;

    // Degenerate case
    if(dy == 0) return vector<int>();
    int top  = max(x0.second,x1.second),
      bottom = min(x0.second,x1.second);
    vector<int> result(top-bottom+1);

     // NB: No, wait: What about exact point i + epsilon -> i+1?
    // UNSTABLE! Rewrite with standard integer algorithm
    double slope = dx/double(dy);
    //    printf("slope = %d/%d = %g\n",dx,dy,slope);
    if(sgn(dy) < 0){
      double x = x1.first;
      for(int y=x1.second;y<=x0.second;y++,x+=slope){
	//	printf("\t%g\n",x);
	result[y-bottom] = ceil(ROUND(x));
      }
    } else {
      double x = x0.first;
      for(int y=x0.second;y<=x1.second;y++,x+=slope){
	//	printf("\t%g\n",x);
	result[y-bottom] = floor(ROUND(x));
      }
    }
    return result;
  }

  static vector< set<int> > rasterize_polygonb(const vector<Eisenstein>& outline)
  {
    int top = INT_MIN, bottom = INT_MAX;
    for(int i=0;i<outline.size();i++){
      if(outline[i].second < bottom) bottom = outline[i].second;
      if(outline[i].second > top)    top = outline[i].second;
    }

    vector< set<int> > xvalues(top-bottom+1);

    for(int i=0;i<outline.size();i++){
      const Eisenstein &ux(outline[i]), &vx(outline[(i+1)%outline.size()]);
      const vector<int> line(rasterize_line(ux,vx));
      int line_bottom = min(ux.second,vx.second);

      for(int i=0;i<line.size();i++)
	xvalues[i-bottom+line_bottom].insert(line[i]);
    }
    return xvalues;
  }
  

  coord2d coord() const { return coord2d(1,0)*first + coord2d(1/2.,sqrt(3/4.))*second; }

  int norm2() const { return first*first + first*second + second*second; }
  double norm() const { return sqrt(norm2()); }
  bool operator<(const Eisenstein& y) const { return norm2() < y.norm2(); }

};


class EOp {
public:
  int a[4];
  int denom;

  EOp(int a0,int a1, int a2, int a3, int denom=1) : a{a0,a1,a2,a3}, 
				     denom(denom) {}
  EOp(int a[4], int denom=1) : a{a[0],a[1],a[2],a[3]}, 
				     denom(denom) {}
  
  EOp operator*(const EOp& B) const {
    return EOp(a[0]*B.a[0] + a[1]*B.a[2], a[0]*B.a[1] + a[1]*B.a[3],
	       a[2]*B.a[0] + a[3]*B.a[2], a[2]*B.a[1] + a[3]*B.a[3], 
	       denom*B.denom);
  }
  Eisenstein operator*(const Eisenstein& x) const {
    return Eisenstein((a[0]*x.first + a[1]*x.second)/denom, (a[2]*x.first + a[3]*x.second)/denom);
  }

  vector<Eisenstein> operator*(const vector<Eisenstein>& xs) const {
    vector<Eisenstein> ys(xs.size());
    for(int i=0;i<xs.size();i++) ys[i] = (*this)*xs[i];
    return ys;
  }

  EOp transpose() const { return EOp(a[0],a[2],a[1],a[3],denom); }

  static EOp GC(int k, int l)        { return EOp(k,-l,l,k+l); }
  static EOp GCInverse(int k, int l) { return EOp(k+l,l,-l,k, k*k + k*l + l*l); }

  static EOp Delta(const Eisenstein& d)
  {
    const Eisenstein e(d.nextCW());
    return EOp(d.first, e.first,d.second,e.second);
  }

  friend ostream& operator<<(ostream &s, const EOp &A)
  {
    s << "{{" << A.a[0] << "," << A.a[1] << "},{" << A.a[2] << "," << A.a[3] << "}}";
    if(A.denom != 1) s << "/" << A.denom;
    return s;
  }

};

#endif
