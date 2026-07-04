#include <limits.h>
#include <limits>
#include <iostream>
#include <cstdio>
#include <vector>
#include <algorithm>
#include <cassert>
#include "fullerenes/geometry.hh"
#include "fullerenes/eisenstein.hh"

// TODO: This stuff does not belong here. Start auxiliary.cc. Also, clean up!

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

vector<Eisenstein> polygon::draw_line(const Eisenstein& xy0, const Eisenstein& xy1)
{
  vector<Eisenstein> result;

  Eisenstein D(xy1-xy0);
  int Dx = D.first, Dy = D.second;
  int sx = sgn(Dx), sy = sgn(Dy);
  Eisenstein xy(xy0);

  // Invariant: (x-xy0)*Dy >= (y-y0)*Dx
  // Or:         sx*i*Dy  >= sy*j*Dx
  // Or:         sx*i*Dy - sy*j*Dx >= 0

  if(sx == 0){
    while(xy.second-sy != xy1.second){ result.push_back(xy); xy.second += sy; }
    return result;
  }
  if(sy == 0){
    while(xy.first-sx != xy1.first){ result.push_back(xy); xy.first += sx; }
    return result;
  }

  int threshold = 0, t0 = sx*sy>0? 0 : -abs(Dy)+1;
  
  if(abs(Dx) > abs(Dy)){	// X-major
    for(xy = xy0; xy.first-sx != xy1.first; xy.first += sx, threshold += sx*Dy){
      while(sx*sy*threshold >= t0){
	result.push_back(xy);
	xy.second += sy;
	threshold-=sy*Dx;
      }
    }
  } else {			// Y-major
    for(xy = xy0; xy.second-sy != xy1.second; xy.second += sy, threshold -= sy*Dx){
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

vector<Eisenstein> polygon::allpoints() const {
  scanline S = scanConvert();

  vector<Eisenstein> points;
  for(int i=0;i<S.xs.size();i++)
    for(int j=0;j<S.xs[i].size()/2;j++){
      int start = S.xs[i][2*j], end = S.xs[i][2*j+1];
      for(int x=start;x<=end;x++)
	points.push_back(Eisenstein(x,i+S.minY));
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

// A rational scanline abscissa num/den, den > 0. Exact comparisons by
// cross-multiplication; outline coordinates are small, so long never overflows.
namespace {
struct rat {
  long num, den;
  rat(long n, long d) : num(d < 0 ? -n : n), den(d < 0 ? -d : d) {}
  bool operator<(const rat& o)  const { return num*o.den <  o.num*den; }
  bool operator<=(const rat& o) const { return num*o.den <= o.num*den; }
  long floor() const { return num >= 0 ? num/den : -((-num + den - 1)/den); }
  long ceil()  const { return -rat(-num, den).floor(); }
};
}  // anonymous: TU-local rational

// STANDARD even-odd scanline polygon fill (edge table + per-row crossing
// pairing; Foley--van Dam sec. 3.6, Abrash's Black Book ch. 38-40 -- the same
// scan scan_triangle implements for triangles), with the standard half-open
// vertex rule: a non-horizontal edge owns scanlines [ymin, ymax).  Valid for
// ANY simple lattice polygon (any edge slope); restores the invariant
//
//     x in allpoints()  <=>  point_included(x)
//
// which the previous version broke by feeding the parity stage draw_line
// PIXELS instead of edge crossings (multiple pixels per row on shallow
// X-major edges corrupt the parity; the C20 star unfolding rotated by (0,1)
// gained two exterior phantom grid points).
//
// Two deviations from the screen-fill version, both forced by the contract
// "all lattice points inside OR ON the polygon":
//   1. crossings are exact rationals (num/den, cross-multiplied compares)
//      instead of an incremental integer DDA -- same numbers, exactness
//      needed, fewer moving parts;
//   2. CLOSED coverage instead of a fill convention (screen fills use
//      top-left-style rules that deliberately drop boundary pixels): spans
//      round ceil(lo)/floor(hi), and the two boundary point classes the
//      parity stage cannot see -- horizontal edges and strict local-max
//      corners -- enter as boundary intervals.  (Local minima fall out of
//      the parity stage as degenerate intervals.)
// Intervals merge at the RATIONAL level: touching intervals share a boundary
// point, while strictly separated ones stay separate even when their integer
// roundings touch -- a sub-unit exterior slit between star arms is never
// bridged.  Oracle-tested in tests/test_seam_step.cc against an independent
// exact inclusion predicate over full bounding boxes.
//
// (Considered and rejected: building this on scan_triangle over the
// ear-clipped triangulation. Integer truncation per triangle loses sub-unit
// continuity at INTERNAL chords -- a contiguous row splits into integer
// spans that no longer reveal interior unit arcs crossing the chord --
// correct as a point set, wrong for span consumers like connect_polygon.)
polygon::scanline polygon::scanConvert() const {
  int minY=INT_MAX, maxY=INT_MIN;
  for(auto xy: reduced_outline) {
    int y = xy.second;
    minY = min(minY,y);
    maxY = max(maxY,y);
  }

  scanline S;
  S.minY = minY;
  S.xs = vector<vector<int> >(maxY-minY+1);

  const int n = (int)reduced_outline.size();
  for(int row=minY; row<=maxY; row++){
    vector<rat> cross;                                  // even-odd crossings, paired after sorting
    vector<pair<rat,rat>> spans;                        // closed rational intervals

    for(int i=0;i<n;i++){
      const Eisenstein A = reduced_outline[i], B = reduced_outline[(i+1)%n];
      if(A.second == B.second){                         // horizontal edge on this row: boundary
        if(A.second == row)
          spans.push_back({rat(min(A.first,B.first),1), rat(max(A.first,B.first),1)});
        continue;
      }
      const int ylo = min(A.second,B.second), yhi = max(A.second,B.second);
      if(row < ylo || row >= yhi) continue;             // half-open ownership [ymin, ymax)
      const long den = B.second - A.second;
      const long num = (long)A.first*den + (long)(row - A.second)*(B.first - A.first);
      cross.push_back(rat(num, den));
    }
    for(int i=0;i<n;i++){                               // strict local-max corners: boundary points
      const Eisenstein P = reduced_outline[i];
      if(P.second != row) continue;
      const int yp = reduced_outline[(i+n-1)%n].second, yn = reduced_outline[(i+1)%n].second;
      if(yp < row && yn < row) spans.push_back({rat(P.first,1), rat(P.first,1)});
    }

    if(cross.size() % 2)
      throw std::logic_error("polygon::scanConvert: odd crossing count (non-simple outline?)");
    sort(cross.begin(), cross.end());
    for(size_t k=0;k+1<cross.size();k+=2) spans.push_back({cross[k], cross[k+1]});

    sort(spans.begin(), spans.end(),
         [](const pair<rat,rat>& a, const pair<rat,rat>& b){ return a.first < b.first; });
    vector<int>& out = S.xs[row-minY];
    for(size_t k=0;k<spans.size();){
      rat lo = spans[k].first, hi = spans[k].second;
      size_t j = k+1;
      while(j<spans.size() && spans[j].first <= hi){    // rational touch/overlap only
        if(hi < spans[j].second) hi = spans[j].second;
        j++;
      }
      const long a0 = lo.ceil(), a1 = hi.floor();
      if(a0 <= a1){ out.push_back((int)a0); out.push_back((int)a1); }
      k = j;
    }
  }

  return S;
}

// Exact integer point-in-polygon, replacing the old floating-point atan2
// winding number (which needed a 1e-10 tolerance and broke the exact-Eisenstein
// invariant). Wedge-based crossing-parity predicate, validated as the reference
// oracle in the fold-unfold seam-step tests. Operates on reduced_outline
// (collinear-free, same polygon).

// On the boundary: exactly collinear with an edge (wedge == 0) and within its
// bounding box.
bool polygon::point_on_periphery(const Eisenstein& x) const
{
  const int n = reduced_outline.size();
  for(int i=0;i<n;i++){
    const Eisenstein s = reduced_outline[i], t = reduced_outline[(i+1)%n];
    if(wedge(t-s, x-s) == 0 &&
       min(s.first, t.first)  <= x.first  && x.first  <= max(s.first, t.first) &&
       min(s.second,t.second) <= x.second && x.second <= max(s.second,t.second))
      return true;
  }
  return false;
}

// Strictly interior: odd crossing parity of the east-going ray, with half-open
// [ymin, ymax) edge ownership so a shared vertex is counted exactly once.
// Boundary points are excluded here -- that is point_on_periphery's job.
bool polygon::point_inside(const Eisenstein& x) const
{
  if(point_on_periphery(x)) return false;
  const int n = reduced_outline.size();
  int parity = 0;
  for(int i=0;i<n;i++){
    const Eisenstein s = reduced_outline[i], t = reduced_outline[(i+1)%n];
    if(s.second == t.second) continue;
    const int ylo = min(s.second, t.second), yhi = max(s.second, t.second);
    if(x.second < ylo || x.second >= yhi) continue;          // half-open [ymin, ymax)
    const long den = t.second - s.second;
    const long num = (long)s.first*den + (long)(x.second - s.second)*(t.first - s.first);
    // the edge crosses the ray strictly right of x iff num/den > x.first (den sign folded in)
    if((den > 0 ? num - (long)x.first*den : (long)x.first*den - num) > 0) parity ^= 1;
  }
  return parity;
}

bool polygon::point_included(const Eisenstein& x) const   // inside or on the boundary
{
  return point_on_periphery(x) || point_inside(x);
}

polygon polygon::operator*(Eisenstein x) const
{
  return polygon(outline*x);
}

polygon polygon::operator+(Eisenstein x) const
{
  return polygon(outline+x);
}

ostream& operator<<(ostream& S, const polygon& P)
{
  S << make_pair(P.outline,P.reduced_outline);
  return S;
}

template<>
matrix3d coord3<double>::outer(const coord3<double>& y) const
{
  return matrix3d(x[0]*y[0],x[0]*y[1],x[0]*y[2],  x[1]*y[0],x[1]*y[1],x[1]*y[2],  x[2]*y[0],x[2]*y[1],x[2]*y[2]);
}

template<>
coord3<double> coord3<double>::operator*(const matrix3d& m) const
{
  return coord3<double>(x[0]*m(0,0)+x[1]*m(1,0)+x[2]*m(2,0),  x[0]*m(0,1)+x[1]*m(1,1)+x[2]*m(2,1),  x[0]*m(0,2)+x[1]*m(1,2)+x[2]*m(2,2));
}

// calculation of the angle beta at b(0,0,0)
template<>
double coord3<double>::angle(const coord3<double>& a, const coord3<double>& c)
{
  const double L2 = a.dot(a);
  const double R2 = c.dot(c);
  const double M2 = (c-a).dot(c-a);
  const double den = 2.0*sqrt(L2 * R2);
  double arg = (L2+R2-M2)/den;
  if(arg > 1)  arg = 1.0;
  if(arg < -1) arg = -1.0;
  return acos(arg);
}

// Nearest point of triangle t to p (metric projection onto the filled triangle):
// project p onto the triangle's plane; if the foot is inside, that is the answer,
// else clamp to the nearest edge. dist2 = |pos-p|^2, +inf for a degenerate face.
ClosestPoint closest_point_on_triangle(const coord3d& p, const tri_t& t,
                                       std::span<const coord3d> points)
{
  const Tri3D T(points, t);
  const double nn = T.n.dot(T.n);                          // |u x v|^2 == denominator
  if(nn < 1e-300)                                          // degenerate face: never nearest
    return {T.a, {1,0,0}, std::numeric_limits<double>::infinity()};

  const coord3d proj = p - T.n * (T.n.dot(p - T.a) / nn);  // foot of the perpendicular
  const coord3d w = proj - T.a;
  const double uu = T.u.dot(T.u), uv = T.u.dot(T.v), vv = T.v.dot(T.v),
               wu = w.dot(T.u), wv = w.dot(T.v);
  const double beta  = (vv*wu - uv*wv)/nn;
  const double gamma = (uu*wv - uv*wu)/nn;
  if(beta >= 0 && gamma >= 0 && beta + gamma <= 1)
    return {proj, {1-beta-gamma, beta, gamma}, (p-proj).dot(p-proj)};

  // Foot outside the triangle: clamp to the nearest edge.
  ClosestPoint best{T.a, {1,0,0}, std::numeric_limits<double>::infinity()};
  const coord3d corner[3] = {T.a, T.b, T.c};
  for(int e = 0; e < 3; e++){
    int i0 = e, i1 = (e+1)%3;
    coord3d p0 = corner[i0], edge = corner[i1] - p0;
    double s = std::clamp((p - p0).dot(edge) / std::max(edge.dot(edge), 1e-300), 0.0, 1.0);
    coord3d q = p0 + edge*s;
    double d2 = (p-q).dot(p-q);
    if(d2 < best.dist2){ best = {q, {0,0,0}, d2}; best.bary[i0] = 1-s; best.bary[i1] = s; }
  }
  return best;
}

// calculation of the derivative of angle beta at b(0,0,0) according to coordinates a and c with fixed b
template<>
void coord3<double>::dangle(const coord3<double>& a, const coord3<double>& c, coord3<double>& da, coord3<double>& dc)
{
  const double L2 = a.dot(a);
  const double R2 = c.dot(c);
  const double M2 = (c-a).dot(c-a);
  const double den = 2.0*sqrt(L2 * R2);
  const double arg = (L2+R2-M2)/den;

  const coord3d dM2__da = (a-c)*2.0;
  const coord3d dL2__da = a*2.0;
  const coord3d dden__da = dL2__da * R2/sqrt(L2*R2);
  const coord3d darg__da = dL2__da * 1.0/den - dM2__da * 1.0/den - dden__da * (L2+R2-M2)/(den*den);

  const coord3d dM2__dc = (c-a)*2.0;
  const coord3d dR2__dc = c*2.0;
  const coord3d dden__dc = dR2__dc * L2/sqrt(L2*R2);
  const coord3d darg__dc = dR2__dc * 1.0/den - dM2__dc * 1.0/den - dden__dc * (L2+R2-M2)/(den*den);

  da = -darg__da * 1.0/sqrt(1.0-arg*arg);
  dc = -darg__dc * 1.0/sqrt(1.0-arg*arg);
}


// calculation of the dihedral angle theta at a(0,0,0), b, c and d,  the result is an angle between -\pi and +\pi (in radians)
// rotate c-d around b-c to a-b mathematically positive
// in a polyhedron, the dihedral abcd, where a is in the centre, surrounded by b c and d:
//    clockwise, convex --> negative
//    counterclockwise, convex --> positive  ... this is the default case
//    clockwise, concave --> positive
//    counterclockwise, concave --> negative
template<>
double coord3<double>::dihedral(const coord3<double>& b, const coord3<double>& c, const coord3<double>& d)
{
  const coord3d ab = b; // a=0
  const coord3d bc = c-b;
  const coord3d cd = d-c;

  const coord3d abc = ab.cross(bc);
  const coord3d bcd = bc.cross(cd);

  const coord3d bc1 = bc/bc.norm();
  const coord3d abc1 = abc/abc.norm();
  const coord3d bcd1 = bcd/bcd.norm();
  const coord3d aux = abc1.cross(bc1);

  const double x = abc1.dot(bcd1);
  const double y = aux.dot(bcd1);

  return atan2(y,x);
}

// calculation of the derivative of dihedral angle theta at a(0,0,0), b, c and d  according to coordinates b, c and d with fixed a
template<>
void coord3<double>::ddihedral(const coord3<double>& b, const coord3<double>& c, const coord3<double>& d, coord3<double>& db, coord3<double>& dc, coord3<double>& dd)
{
  const coord3d ab = b; // a=0
  const coord3d bc = c-b;
  const coord3d cd = d-c;

  const double bc_length_inv = 1.0/bc.norm();
  const coord3d bc1 = bc*bc_length_inv;

  const coord3d abc = ab.cross(bc);
  const coord3d bcd = bc.cross(cd);

  const double abc_length_inv = 1.0/abc.norm();
  const double bcd_length_inv = 1.0/bcd.norm();
  const coord3d abc1 = abc * abc_length_inv;
  const coord3d bcd1 = bcd * bcd_length_inv;

  const coord3d aux = abc1.cross(bc1);

  const double x = abc1.dot(bcd1);
  const double y = aux.dot(bcd1);

//  const double dihedral_abcd = atan2(y,x);
//  cout << "D: "<< dihedral_abcd<<endl;

  const matrix3d dab__db = matrix3d::unit_matrix();
  const matrix3d dbc__db = -matrix3d::unit_matrix();
  const matrix3d dbc__dc = matrix3d::unit_matrix();
  const matrix3d dcd__dc = -matrix3d::unit_matrix();
  const matrix3d dcd__dd = matrix3d::unit_matrix();

  // bc_length_inv=1/dsqrt(bc_x**2 + bc_y**2 + bc_z**2)
  const coord3d dbc_length_inv__dbc = -bc * pow(bc_length_inv, 3);

  // bc_length_inv=1/dsqrt(bc_x**2 + bc_y**2 + bc_z**2)
  // vec = vec * mtx
  const coord3d dbc_length_inv__db = dbc_length_inv__dbc * dbc__db;
  const coord3d dbc_length_inv__dc = dbc_length_inv__dbc * dbc__dc;

  const matrix3d dbc1__dbc = matrix3d::unit_matrix() * bc_length_inv; 

  // bc1_x=bc_x*bc_length_inv
  // bc1_y=bc_y*bc_length_inv
  // bc1_z=bc_z*bc_length_inv
  // mtx = mtx * mtx + vec outer vec
  const matrix3d dbc1__db = dbc1__dbc*dbc__db + bc.outer(dbc_length_inv__db);
  const matrix3d dbc1__dc = dbc1__dbc*dbc__dc + bc.outer(dbc_length_inv__dc);

  // abc_x=ab_y*bc_z - ab_z*bc_y
  // abc_y=ab_z*bc_x - ab_x*bc_z
  // abc_z=ab_x*bc_y - ab_y*bc_x
  //FIXME is there a more elegant way of doing this?
  const matrix3d dabc__dab = matrix3d(0,bc[2],-bc[1], -bc[2],0,bc[0], bc[1],-bc[0],0); 
  const matrix3d dabc__dbc = matrix3d(0,-ab[2],ab[1], ab[2],0,-ab[0], -ab[1],ab[0],0); 

  // bcd_x=bc_y*cd_z - bc_z*cd_y
  // bcd_y=bc_z*cd_x - bc_x*cd_z
  // bcd_z=bc_x*cd_y - bc_y*cd_x
  //FIXME is there a more elegant way of doing this?
  const matrix3d dbcd__dbc = matrix3d(0,cd[2],-cd[1], -cd[2],0,cd[0], cd[1],-cd[0],0); 
  const matrix3d dbcd__dcd = matrix3d(0,-bc[2],bc[1], bc[2],0,-bc[0], -bc[1],bc[0],0); 

  // abc_x=-ab_y*bc_z + ab_z*bc_y
  // abc_y=-ab_z*bc_x + ab_x*bc_z
  // abc_z=-ab_x*bc_y + ab_y*bc_x
  // mtx = mtx * mtx + mtx * mtx
  const matrix3d dabc__db = dabc__dab*dab__db + dabc__dbc*dbc__db;
  const matrix3d dabc__dc =                     dabc__dbc*dbc__dc;

  // bcd_x=-bc_y*cd_z + bc_z*cd_y
  // bcd_y=-bc_z*cd_x + bc_x*cd_z
  // bcd_z=-bc_x*cd_y + bc_y*cd_x
  // mtx = mtx * mtx + mtx * mtx
  const matrix3d dbcd__db = dbcd__dbc*dbc__db;
  const matrix3d dbcd__dc = dbcd__dbc*dbc__dc + dbcd__dcd*dcd__dc;
  const matrix3d dbcd__dd =                     dbcd__dcd*dcd__dd;

  // abc_length_inv=1/dsqrt(abc_x**2 + abc_y**2 + abc_z**2)
  // bcd_length_inv=1/dsqrt(bcd_x**2 + bcd_y**2 + bcd_z**2)
  const coord3d dabc_length_inv__dabc = -abc*pow(abc_length_inv,3);
  const coord3d dbcd_length_inv__dbcd = -bcd*pow(bcd_length_inv,3);

  // abc_length_inv=1/dsqrt(abc_x**2 + abc_y**2 + abc_z**2)
  // vec = vec * mtx
  const coord3d dabc_length_inv__db = dabc_length_inv__dabc*dabc__db;
  const coord3d dabc_length_inv__dc = dabc_length_inv__dabc*dabc__dc;

  // bcd_length_inv=1/dsqrt(bcd_x**2 + bcd_y**2 + bcd_z**2)
  // vec = vec * mtx
  const coord3d dbcd_length_inv__db = dbcd_length_inv__dbcd * dbcd__db;
  const coord3d dbcd_length_inv__dc = dbcd_length_inv__dbcd * dbcd__dc;
  const coord3d dbcd_length_inv__dd = dbcd_length_inv__dbcd * dbcd__dd;

  // abc1_x=abc_x*abc_length_inv
  // abc1_y=abc_y*abc_length_inv
  // abc1_z=abc_z*abc_length_inv
  const matrix3d dabc1__dabc = matrix3d::unit_matrix() * abc_length_inv;

  // abc1_x=abc_x*abc_length_inv
  // abc1_y=abc_y*abc_length_inv
  // abc1_z=abc_z*abc_length_inv
  // mtx = mtx * mtx + vec outer vec
  const matrix3d dabc1__db = dabc1__dabc*dabc__db + abc.outer(dabc_length_inv__db);
  const matrix3d dabc1__dc = dabc1__dabc*dabc__dc + abc.outer(dabc_length_inv__dc);

  // bcd1_x=bcd_x*bcd_length_inv
  // bcd1_y=bcd_y*bcd_length_inv
  // bcd1_z=bcd_z*bcd_length_inv
  const matrix3d dbcd1__dbcd = matrix3d::unit_matrix() * bcd_length_inv;

  // bcd1_x=bcd_x*bcd_length_inv
  // bcd1_y=bcd_y*bcd_length_inv
  // bcd1_z=bcd_z*bcd_length_inv
  // mtx = mtx*mtx + vec outer vec
  const matrix3d dbcd1__db = dbcd1__dbcd * dbcd__db + bcd.outer(dbcd_length_inv__db);
  const matrix3d dbcd1__dc = dbcd1__dbcd * dbcd__dc + bcd.outer(dbcd_length_inv__dc);
  const matrix3d dbcd1__dd = dbcd1__dbcd * dbcd__dd + bcd.outer(dbcd_length_inv__dd);

  // aux_x=abc1_y*bc1_z-bc1_y*abc1_z
  // aux_y=abc1_z*bc1_x-bc1_z*abc1_x
  // aux_z=abc1_x*bc1_y-bc1_x*abc1_y
  //FIXME is there a more elegant way of doing this?
  const matrix3d daux__dabc1 = matrix3d(0,bc1[2],-bc1[1], -bc1[2],0,bc1[0], bc1[1],-bc1[0],0);
  const matrix3d daux__dbc1 = matrix3d(0,-abc1[2],abc1[1], abc1[2],0,-abc1[0], -abc1[1],abc1[0],0);

  // aux_x=abc1_y*bc1_z-bc1_y*abc1_z
  // aux_y=abc1_z*bc1_x-bc1_z*abc1_x
  // aux_z=abc1_x*bc1_y-bc1_x*abc1_y
  // mtx = mtx*mtx + mtx*mtx
  const matrix3d daux__db = daux__dabc1 * dabc1__db + daux__dbc1 * dbc1__db;
  const matrix3d daux__dc = daux__dabc1 * dabc1__dc + daux__dbc1 * dbc1__dc;

  // y=aux_x*bcd1_x + aux_y*bcd1_y + aux_z*bcd1_z
  // vec = vec * mtx
  const coord3d dy__db = bcd1 * daux__db + aux * dbcd1__db;
  const coord3d dy__dc = bcd1 * daux__dc + aux * dbcd1__dc;
  const coord3d dy__dd =                   aux * dbcd1__dd;

  // x=abc1_x*bcd1_x + abc1_y*bcd1_y + abc1_z*bcd1_z
  // vec = vec * mtx
  const coord3d dx__db = bcd1 * dabc1__db + abc1 * dbcd1__db;
  const coord3d dx__dc = bcd1 * dabc1__dc + abc1 * dbcd1__dc;
  const coord3d dx__dd =                    abc1 * dbcd1__dd;

  // df__dx=-y/(x**2 + y**2)
  // df__dy=x/(x**2 + y**2)
  const double df__dx = -y/(x*x + y*y);
  const double df__dy =  x/(x*x + y*y);

  // f=atan2(y,x)
  // vec = vec*sca + vec*sca
  db = dx__db*df__dx + dy__db*df__dy;
  dc = dx__dc*df__dx + dy__dc*df__dy;
  dd = dx__dd*df__dx + dy__dd*df__dy;
}


// assumes counter clockwise orientation and a convex singularity
template<>
double coord3<double>::ideal_dihedral(const int A, const int B, const int C, const double ur, const double us, const double ut)
{
//        t   B   s      
//          \   /
//        C   u   A
//            |
//            r
  const double eps = 1.0e-8;
  if (1.0/A + 1.0/B + 1.0/C > 0.5+eps){ // positive curvature // make sure 3 * 1/6 is recognised to be planar
    // neighbours are sorted CCW
    // assumption: Faces A, B and C are planar.
    // calculate squares of rs, st and rt via law of cosines
    const double ur_2 = ur*ur;
    const double us_2 = us*us;
    const double ut_2 = ut*ut;
    const double rs_2 = us_2 + ur_2 - 2.0*us*ur*cos(M_PI*(1.0-2.0/double(A)));
    const double st_2 = us_2 + ut_2 - 2.0*us*ut*cos(M_PI*(1.0-2.0/double(B)));
    const double tr_2 = ut_2 + ur_2 - 2.0*ut*ur*cos(M_PI*(1.0-2.0/double(C)));

    const double bx = 0;
    const double by = sqrt(ur_2);
    const double bz = 0;
    const double cx = 0.5 * sqrt((-pow(ur_2,2) - pow(us_2 - rs_2,2) + 2*ur_2*(us_2 + rs_2))/(ur_2));
    const double cy = (ur_2 + us_2 - rs_2)/(2.*sqrt(ur_2));
    const double cz = 0;
    const double dx = -(pow(ur_2,2) + (us_2 - rs_2)*(ut_2 - tr_2) - ur_2*(us_2 + ut_2 + rs_2 + tr_2 - 2*st_2))
                        / (2.*sqrt(ur_2)*sqrt(-pow(ur_2,2) - pow(us_2 - rs_2,2) + 2*ur_2*(us_2 + rs_2)));
    const double dy = (ur_2 + ut_2 - tr_2)/(2.*sqrt(ur_2));
    const double dz = sqrt((pow(us_2,2)*tr_2 + pow(ur_2,2)*st_2 + rs_2*(pow(ut_2,2) + ut_2*(rs_2 - tr_2 - st_2) + tr_2*st_2) 
                       - us_2*(ut_2*(rs_2 + tr_2 - st_2) + tr_2*(rs_2 - tr_2 + st_2))
                       - ur_2*((rs_2 + tr_2 - st_2)*st_2 + ut_2*(rs_2 - tr_2 + st_2) + us_2*(-rs_2 + tr_2 + st_2)))
                       / (pow(ur_2,2) + pow(us_2 - rs_2,2) - 2*ur_2*(us_2 + rs_2)));

   // cout << "bx,by,bz: " << bx<<" " <<by<<" "<<bz << endl;
   // cout << "cx,cy,cz: " << cx<<" " <<cy<<" "<<cz << endl;
   // cout << "dx,dy,dz: " << dx<<" " <<dy<<" "<<dz << endl;

    return coord3d::dihedral(coord3d(bx,by,bz), coord3d(cx,cy,cz), coord3d(dx,dy,dz));
  } else {
    //    cout << "planar or negative curvature\n";
    // in the case of negative curvature the adjacent faces cannot be planar
    // and the dihedral cannot be calculated (but 0 is a plausible value)
    return 0.0;
  }
}

