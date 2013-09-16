#include <limits.h>
#include <iostream>
#include <cstdio>
#include <vector>
#include <algorithm>
#include <cassert>
#include <complex>
#include "geometry.hh"
#include "eisenstein.hh"

string pad_string(const string& s, int length, char padchar)
{
  string result(s);
  int pad = length - s.size();
  string padstring;
  for(int i=0;i<pad;i++) padstring += padchar;
  return padstring+result;
}

int gcd(int a, int b)
{
  int c;
  while ( a != 0 ) { c = a; a = b%a;  b = c;  }
  return b;
}

// vector<Eisenstein> draw_line(const Eisenstein &a, const Eisenstein &b)
// {
//   if(a.second == b.second)
// }

vector<Eisenstein> polygon::draw_line(const Eisenstein& x0, const Eisenstein& x1)
{
  vector<Eisenstein> result;

  Eisenstein D(x1-x0);
  int Dx = D.first, Dy = D.second;
  int sx = sgn(Dx), sy = sgn(Dy);
  Eisenstein xy(x0);

  // Invariant: (x-x0)*Dy >= (y-y0)*Dx
  // Or:         sx*i*Dy  >= sy*j*Dx
  // Or:         sx*i*Dy - sy*j*Dx >= 0

  if(sx == 0){
    while(xy.second-sy != x1.second){ result.push_back(xy); xy.second += sy; }
    return result;
  }
  if(sy == 0){
    //      while(xy.first-sx != x1.first){ result.push_back(xy); xy.first += sx; }
    result.push_back(x0);
    result.push_back(x1);
    return result;
  }

  int threshold = 0, t0 = sx*sy>0? 0 : -abs(Dy)+1;
  
  if(abs(Dx) > abs(Dy)){	// X-major
    for(xy = x0; xy.first-sx != x1.first; xy.first += sx, threshold += sx*Dy){
      while(sx*sy*threshold >= t0){
	result.push_back(xy);
	xy.second += sy;
	threshold-=sy*Dx;
      }
    }
  } else {			// Y-major
    for(xy = x0; xy.second-sy != x1.second; xy.second += sy, threshold -= sy*Dx){
      while(sx*sy*threshold < t0){
	xy.first  += sx;
	threshold += sx*Dy;
      }
      result.push_back(xy);
    }
  }

  return result;
}

// Slope of line segment i.
pair<int,int> polygon::slope(int i,bool reduced) const
{
  const vector<Eisenstein>& pts(reduced? reduced_outline : outline);
  Eisenstein x(pts[i]), y(pts[(i+1)%pts.size()]), dx(y-x);
  int d = gcd(dx.first,dx.second);
  return make_pair(dx.first/d, dx.second/d);
}


// Node j is at left-turn (-1), straight (0), or right-turn (1) 
// in CW traversal
int polygon::turn_direction(int j,bool reduced) const {

  const vector<Eisenstein>& pts(reduced? reduced_outline : outline);
  Eisenstein 
    xi(pts[(pts.size()+j-1)%pts.size()]),
    xj(pts[j]),
    xk(pts[(j+1)%pts.size()]),
    dx1(xj-xi), dx2(xk-xj);

  return sgn(dx2.first * dx1.second - dx2.second * dx1.first);
}

// Node j is at peak turn.
bool polygon::peak(int j, bool reduced) const {

  const vector<Eisenstein>& pts(reduced? reduced_outline : outline);
  Eisenstein 
    xi(pts[(pts.size()+j-1)%pts.size()]),
    xj(pts[j]),
    xk(pts[(j+1)%pts.size()]),
    dx1(xj-xi), dx2(xk-xj);

  return sgn(dx1.second)*sgn(dx2.second) < 0;
}

// Exactly one of either segment i--j or j--k is horizontal
bool polygon::saddle(int j, bool reduced) const {
  const vector<Eisenstein>& pts(reduced? reduced_outline : outline);
  Eisenstein 
    xi(pts[(pts.size()+j-1)%pts.size()]),
    xj(pts[j]),
    xk(pts[(j+1)%pts.size()]),
    dx1(xj-xi), dx2(xk-xj);

  return (dx1.second == 0) ^ (dx2.second == 0);
}


vector<Eisenstein> polygon::reduce() const {
  vector<Eisenstein> xs(1,outline[0]);

  pair<int,int> slope0(slope(0));
  for(int i=1;i<outline.size();i++){
    pair<int,int> slopei(slope(i));
    if(slopei != slope0){
      xs.push_back(outline[i]);
      slope0 = slopei;
    }
  }  
  
  return xs;
}

set<Eisenstein> polygon::allpoints() const {
  scanline S = scanConvert();

  set<Eisenstein> points;
  for(int i=0;i<S.xs.size();i++)
    for(int j=0;j<S.xs[i].size()/2;j++){
      int start = S.xs[i][2*j], end = S.xs[i][2*j+1];
      for(int x=start;x<=end;x++)
	points.insert(Eisenstein(x,i+S.minY));
    }

  return points;
}

vector<Eisenstein> polygon::controlpoints() const {
  scanline S = scanConvert();

  vector<Eisenstein> points;
  for(int i=0;i<S.xs.size();i++)
    for(int j=0;j<S.xs[i].size();j++)
      points.push_back(Eisenstein(S.xs[i][j],i+S.minY));

  return points;
}

polygon::scanline polygon::scanConvert() const {
  int minY=INT_MAX, maxY=INT_MIN;
  for(int i=0;i<reduced_outline.size();i++) {
    if(reduced_outline[i].second < minY) minY = reduced_outline[i].second;
    if(reduced_outline[i].second > maxY) maxY = reduced_outline[i].second;
  }
    
  scanline S;
  S.minY = minY;
  S.xs = vector<vector<int> >(maxY-minY+1);

  for(int i=0;i<reduced_outline.size();i++){ 
    vector<Eisenstein> segment(draw_line(reduced_outline[i],reduced_outline[(i+1)%reduced_outline.size()]));

    if(peak(i,true) || (saddle(i,true) && turn_direction(i,true) == -1)) // If peak or left-turning saddle
      {
	int loc = segment[0].second-minY;
	assert(loc >= 0 && loc<S.xs.size());
	S.xs[segment[0].second-minY].push_back(segment[0].first);          // include point twice.
      }

    for(int j=1;j<segment.size();j++){
      const Eisenstein& xy(segment[j]);
      int loc = xy.second-minY;
      assert(loc >= 0 && loc<S.xs.size());
      S.xs[xy.second-minY].push_back(xy.first);
    }
  }

  for(int i=0;i<S.xs.size();i++)
    sort(S.xs[i].begin(),S.xs[i].end());

  return S;
}

// http://www.softsurfer.com/Archive/algorithm_0103/algorithm_0103.htm#Winding%20Number
// TODO: This seems to fail when compiled with -std=c++0x and -O0!
double polygon::winding_number(const Eisenstein& x) const {
  vector<Eisenstein> Cp(reduced_outline.size());
  for(int i=0;i<reduced_outline.size();i++)
    Cp[i] = reduced_outline[i]-x;

  double wn = 0;
  for(int i=0;i<Cp.size();i++){
    Eisenstein segment = Cp[(i+1)%Cp.size()] - Cp[i];
    double theta = atan2(segment.second,segment.first);
    wn += theta;
  }
  wn /= 2*M_PI;
  return wn;
}

bool polygon::point_inside(const Eisenstein& x) const 
{
  return winding_number(x) != 0;
}
  

ostream& operator<<(ostream& S, const polygon& P)
{
  S << make_pair(P.outline,P.reduced_outline);
  return S;
}
