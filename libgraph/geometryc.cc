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

int gcd (int a, int b)
{
  int c;
  while ( a != 0 ) {
     c = a; a = b%a;  b = c;
  }
  return b;
}


pair<int,int> polygon::slope(int i) const
{
  Eisenstein x((*this)[i]), y((*this)[(i+1)%this->size()]), dx(y-x);
  int d = gcd(dx.first,dx.second);
  return make_pair(dx.first/d, dx.second/d);
}

polygon polygon::reduce() const {
  vector<pair<int,int> > xs(1,(*this)[0]);

  pair<int,int> slope0(slope(0));
  for(int i=1;i<this->size();i++){
    pair<int,int> slopei(slope(i));
    if(slopei != slope0){
      xs.push_back((*this)[i]);
      slope0 = slopei;
    }
  }  
  
  return polygon(xs);
}

polygon::scanline polygon::scanConvert() const {
  
  vector<pointinfo> points;
  vector<int> x,y;
  polygon p(reduce());
  transform(p.begin(),p.end(),back_inserter(x), getFirst<int,int>);
  transform(p.begin(),p.end(),back_inserter(y), getSecond<int,int>);

  for(size_t i = 0; i < y.size(); ++i) {
    size_t j = (i+1)%y.size();
    //scanning line i -> j.
    int x0 = x[i], y0 = y[i];
    int x1 = x[j], y1 = y[j];
    if(y0 > y1) { 
      std::swap(x0, x1);
      std::swap(y0, y1);
    }
    int dx = x1-x0, dy = y1-y0;
    int x = x0, y = y0, r = 0, d = dy;
    //the point is (x+r/d, y)
    while(y < y1-1) {
      //update x
      r += dx;
      while(r >= d) { r -= d; ++x; }
      while(r < 0)  { r += d; --x; }
      //update y
      ++y;
      //insert point
      points.push_back(pointinfo(x, y, r == 0, true, 0));
    }
    //the above also works as intended when y0==y1
  }
  vector<int> lineDir(y.size());
  //lineDir[i] holds the direction of the line:
  //from x[i],y[i] to x[i+1],y[i+1] (mod y.size())
  for(size_t i = 0; i < y.size(); ++i) {
    size_t j = (i+1)%y.size();
    lineDir[i] = y[j]-y[i];
    lineDir[i] = lineDir[i] > 0 ? 1 : (lineDir[i] == 0 ? 0 : -1);
  }
  for(size_t k = 0;k < y.size(); ++k) {
    size_t i = (k-2+y.size())%y.size(), j = (k-1+y.size())%y.size(), l = (k+1)%y.size();
    if(lineDir[j] != 0 && lineDir[k] != 0) {
      points.push_back(pointinfo(x[k], y[k], true, lineDir[j]*lineDir[k] == 1, 0));
    } else {
      assert(lineDir[j] != 0 || lineDir[k] != 0);
      if(lineDir[j] == 0) { //lineDir[k] != 0
	assert(x[j] != x[k]);
	assert(lineDir[i] != 0);
	if(x[k] < x[j]) {
	  points.push_back(pointinfo(x[k], y[k], true, lineDir[i]*lineDir[k] == 1, x[j]-x[k]));
	}
      } else { //lineDir[k] == 0, lineDir[j] != 0
	// This fails if polygon is not "reduced", in the sense that there exists three points on the outline
	// such that xy[i]--xy[j] has the same slope as xy[j]--xy[k].
	assert(x[k] != x[l]);
	//	if(lineDir[l] == 0){
	  //	  printf("(i,j,k,l) = (%ld,%ld,%ld,%ld)\n",i,j,k,l);
	  //	  printf("lineDir   = (%d,%d,%d,%d)\n",lineDir[i],lineDir[j],lineDir[k],lineDir[l]);
	  //	  cout << (*this)[i] << "; "<< (*this)[j] << "; "<< (*this)[k] << "; " << (*this)[l] << "; " << (*this)[(k+2)%y.size()] << ";\n";
	  //	}
	assert(lineDir[l] != 0);
	if(x[k] < x[l]) {
	  points.push_back(pointinfo(x[k], y[k], true, lineDir[j]*lineDir[l] == 1, x[l]-x[k]));
	}
      }
    }
  }
  sort(points.begin(), points.end());
  scanline s;
  s.minY = points[0].y;
  int y_cur = s.minY-1;
  size_t i = 0;
  while(i < points.size()) {
    for(int t = 0;t < points[i].y-y_cur;++t) {
      s.xs.push_back(vector<int>());
      s.edge_xs.push_back(vector<int>());
    }
    y_cur = points[i].y;
    size_t j = i+1;
    while(j < points.size() && points[j].y == y_cur) ++j;
    bool inside = false;
    for(;i < j;++i) {
      if(points[i].integral) {
	for(int w = 0;w <= points[i].width; ++w) {
	  s.edge_xs.back().push_back(points[i].x+w);
	}
      }
      //cout << points[i].x << ", " << points[i].y << ", " << points[i].integral<< ", " << points[i].sameDir << ", " << points[i].width << endl;
      if(inside) {
	s.xs.back().push_back(points[i].x-(int)points[i].integral);
      }
      if(points[i].sameDir != inside) {
	s.xs.back().push_back(points[i].x+points[i].width+1);
      }
      if(points[i].sameDir) {
	inside = !inside;
      }
    }
    i=j;
  }
  return s;
}
