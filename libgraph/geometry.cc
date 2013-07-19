#include <iostream>
#include <cstdio>
#include <vector>
#include <algorithm>
#include <cassert>
#include <complex>
#include "geometry.hh"


string pad_string(const string& s, int length, char padchar)
{
  string result(s);
  int pad = length - s.size();
  string padstring;
  for(int i=0;i<pad;i++) padstring += padchar;
  return padstring+result;
}

polygon::scanline polygon::scanConvert() const {
  
  vector<pointinfo> points;
  
  for(size_t i = 0; i < y.size(); ++i) {
    size_t j = (i+1)%y.size();
    //scanning line i -> j.
    int x0 = x[i], y0 = y[i];
    int x1 = x[j], y1 = y[j];
    if(y0 > y1) { 
      swap(x0, x1);
      swap(y0, y1);
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
	assert(x[k] != x[l]);
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
