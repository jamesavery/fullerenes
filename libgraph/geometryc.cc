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

matrix3d coord3d::outer(const coord3d& y) const {
  return matrix3d(x[0]*y[0],x[0]*y[1],x[0]*y[2],  x[1]*y[0],x[1]*y[1],x[1]*y[2],  x[2]*y[0],x[2]*y[1],x[2]*y[2]);
}

// calculation of the angle beta at b(0,0,0)
double coord3d::angle(const coord3d& a, const coord3d& c)
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

// calculation of the derivative of angle beta at b(0,0,0) according to coordinates a and c with fixed b
void coord3d::dangle(const coord3d& a, const coord3d& c, coord3d& da, coord3d& dc)
{
  const double L2 = a.dot(a);
  const double R2 = c.dot(c);
  const double M2 = (c-a).dot(c-a);
  const double den = 2.0*sqrt(L2 * R2);
  const double arg = (L2+R2-M2)/den;

  const coord3d dM2__da = (a-c)*2;
  const coord3d dL2__da = a*2;
  const coord3d dden__da = dL2__da * R2/sqrt(L2*R2);
  const coord3d darg__da = dL2__da * 1.0/den - dM2__da * 1.0/den - dden__da * (L2+R2-M2)/(den*den);

  const coord3d dM2__dc = (c-a)*2.0;
  const coord3d dR2__dc = c*2.0;
  const coord3d dden__dc = dR2__dc * L2/sqrt(L2*R2);
  const coord3d darg__dc = dR2__dc * 1.0/den - dM2__dc * 1.0/den - dden__dc * (L2+R2-M2)/(den*den);

  da = -darg__da * 1.0/sqrt(1.0-arg*arg);
  dc = -darg__dc * 1.0/sqrt(1.0-arg*arg);
}


// calculation of the dihedral angle theta at a(0,0,0), b, c and d,  the result is an angel between -\pi and +\pi (in radians)
double coord3d::dihedral(const coord3d& b, const coord3d& c, const coord3d& d)
{
  const coord3d ab = b;
  const coord3d bc = c-b;
  const coord3d cd = d-c;
  cout << "ab,bc,cd" << ab<<bc<<cd<<endl;

  const coord3d bc1 = bc/bc.norm();
  cout << "bc1" << bc1 <<endl;

  const coord3d abc = ab.cross(bc);
  const coord3d bcd = bc.cross(cd);
  cout << "abc,bcd" << abc<<bcd <<endl;

  const coord3d abc1 = abc/abc.norm();
  const coord3d bcd1 = bcd/bcd.norm();
  cout << "abc1,bcd1" << abc1<<bcd1 <<endl;

  const coord3d aux = abc1.cross(bc1);
  cout << "aux" << aux <<endl;

  const double x = abc1.dot(bcd1);
  const double y = aux.dot(bcd1);
  const coord3d aux2 = bc1.dot(bcd1);
  cout << "aux2" << aux2 <<endl;
  cout << "x,y" << x<<y <<endl;

  return atan2(y,x);
}

// calculation of the derivative of dihedral angle theta at a(0,0,0), b, c and d  according to coordinates b, c and d with fixed a
void coord3d::ddihedral(const coord3d& b, const coord3d& c, const coord3d& d, coord3d& db, coord3d& dc, coord3d& dd)
{

// C at first the dihedral (copied from above)
// c vectors ab, bc and cd
//       ab_x=ax-bx
//       ab_y=ay-by
//       ab_z=az-bz
//       bc_x=bx-cx
//       bc_y=by-cy
//       bc_z=bz-cz
//       cd_x=cx-dx
//       cd_y=cy-dy
//       cd_z=cz-dz
// c vector bc normed to length 1
//       bc_length_inv=1/dsqrt(bc_x**2 + bc_y**2 + bc_z**2)
//       bc1_x=bc_x*bc_length_inv
//       bc1_y=bc_y*bc_length_inv
//       bc1_z=bc_z*bc_length_inv
// c normal vectors on abc and bcd
// c and the signs are this way because one of the two vectors points in the wrong direction
//       abc_x=-ab_y*bc_z + ab_z*bc_y
//       abc_y=-ab_z*bc_x + ab_x*bc_z
//       abc_z=-ab_x*bc_y + ab_y*bc_x
//       bcd_x=-bc_y*cd_z + bc_z*cd_y
//       bcd_y=-bc_z*cd_x + bc_x*cd_z
//       bcd_z=-bc_x*cd_y + bc_y*cd_x
// c their respective lengths
//       abc_length_inv=1/dsqrt(abc_x**2 + abc_y**2 + abc_z**2)
//       bcd_length_inv=1/dsqrt(bcd_x**2 + bcd_y**2 + bcd_z**2)
// c normal vectors (length 1) on abc and bcd
//       abc1_x=abc_x*abc_length_inv
//       abc1_y=abc_y*abc_length_inv
//       abc1_z=abc_z*abc_length_inv
//       bcd1_x=bcd_x*bcd_length_inv
//       bcd1_y=bcd_y*bcd_length_inv
//       bcd1_z=bcd_z*bcd_length_inv
// c abc \times bcd
//       aux_x=abc1_y*bc1_z-bc1_y*abc1_z
//       aux_y=abc1_z*bc1_x-bc1_z*abc1_x
//       aux_z=abc1_x*bc1_y-bc1_x*abc1_y
// c two auxiliary reals
// c     x=\vec abc1 \cdot \vec bcd_1
//       x=abc1_x*bcd1_x + abc1_y*bcd1_y + abc1_z*bcd1_z
// c     y=\vec aux  \cdot \vec bcd_1
//       y=aux_x*bcd1_x + aux_y*bcd1_y + aux_z*bcd1_z
// c the result
//       dihedral_abcd=datan2(y, x)

    const coord3d ab = b;
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

    const double dihedral_abcd = atan2(y,x);
    cout << "D: "<< dihedral_abcd<<endl;

// C THE DERIVATIVES
// c to be read from bottom to top
// 
// c derivatives of single vectors
// c all other combinations are zero
//       1=1
//       dab_x__dbx=-1
//       1=1
//       dab_y__dby=-1
//       1=1
//       dab_z__dbz=-1
// 
//       //dbc_x__dbx=1
//       dbc_x__dcx=-1
//       dbc_y__dby=1
//       dbc_y__dcy=-1
//       dbc_z__dbz=1
//       dbc_z__dcz=-1
// 
//       dcd_x__dcx=1
//       dcd_x__ddx=-1
//       dcd_y__dcy=1
//       dcd_y__ddy=-1
//       dcd_z__dcz=1
//       dcd_z__ddz=-1

//          dab__da = ab.dnorm();
//          dab__db = -ab.dnorm();
//          dab__dc = coord3d(0,0,0);
        const matrix3d dab__da = matrix3d(-1,0,0,0,-1,0,0,0,-1);
        const matrix3d dab__db = matrix3d(1,0,0,0,1,0,0,0,1);
        const matrix3d dab__dc = matrix3d();
        const matrix3d dab__dd = matrix3d();
// 
//          dbc__da = coord3d(0,0,0);
//          dbc__db = bc.dnorm();
//          dbc__dc = -bc.norm();
//          dbc__dd = coord3d(0,0,0);
        const matrix3d dbc__da = matrix3d();
        const matrix3d dbc__db = matrix3d(-1,0,0,0,-1,0,0,0,-1);
        const matrix3d dbc__dc = matrix3d(1,0,0,0,1,0,0,0,1);
        const matrix3d dbc__dd = matrix3d();
// 
//          dcd__da = coord3d(0,0,0);
//          dcd__db = coord3d(0,0,0);
//          dcd__dc = cd.dnorm();
//          dcd__dd = -cd.dnorm();
        const matrix3d dcd__da = matrix3d();
        const matrix3d dcd__db = matrix3d();
        const matrix3d dcd__dc = matrix3d(-1,0,0,0,-1,0,0,0,-1);
        const matrix3d dcd__dd = matrix3d(1,0,0,0,1,0,0,0,1);

// c bc_length_inv=1/dsqrt(bc_x**2 + bc_y**2 + bc_z**2)
//       bc_length_inv_cub=bc_length_inv**3
//       dbc_length_inv__dbc_x=-bc_x*bc_length_inv_cub
//       dbc_length_inv__dbc_y=-bc_y*bc_length_inv_cub
//       dbc_length_inv__dbc_z=-bc_z*bc_length_inv_cub

        const double bc_length_inv_cub = pow(bc_length_inv, 3);
        const coord3d dbc_length_inv__dbc = -bc * bc_length_inv_cub;

// c bc_length_inv=1/dsqrt(bc_x**2 + bc_y**2 + bc_z**2)
// c and the other terms are zero
//       dbc_length_inv__dbx = dbc_length_inv__dbc_x
//       dbc_length_inv__dby = dbc_length_inv__dbc_y
//       dbc_length_inv__dbz = dbc_length_inv__dbc_z
//       dbc_length_inv__dcx = -dbc_length_inv__dbc_x
//       dbc_length_inv__dcy = -dbc_length_inv__dbc_y
//       dbc_length_inv__dcz = -dbc_length_inv__dbc_z

       // vec = mtx * vec
       const coord3d dbc_length_inv__db = dbc__db * dbc_length_inv__dbc;
       const coord3d dbc_length_inv__dc = dbc__dc * dbc_length_inv__dbc;

// c bc1_x=bc_x*bc_length_inv
//        dbc1_x__dbx= bc_length_inv*dbc_x__dbx + bc_x*dbc_length_inv__dbx
//        dbc1_x__dby=                            bc_x*dbc_length_inv__dby
//        dbc1_x__dbz=                            bc_x*dbc_length_inv__dbz
//        dbc1_x__dcx= bc_length_inv*dbc_x__dcx + bc_x*dbc_length_inv__dcx
//        dbc1_x__dcy=                            bc_x*dbc_length_inv__dcy
//        dbc1_x__dcz=                            bc_x*dbc_length_inv__dcz
// c bc1_y=bc_y*bc_length_inv
//        dbc1_y__dbx=                            bc_y*dbc_length_inv__dbx
//        dbc1_y__dby= bc_length_inv*dbc_y__dby + bc_y*dbc_length_inv__dby
//        dbc1_y__dbz=                            bc_y*dbc_length_inv__dbz
//        dbc1_y__dcx=                            bc_y*dbc_length_inv__dcx
//        dbc1_y__dcy= bc_length_inv*dbc_y__dcy + bc_y*dbc_length_inv__dcy
//        dbc1_y__dcz=                            bc_y*dbc_length_inv__dcz
// c bc1_z=bc_z*bc_length_inv
//        dbc1_z__dbx=                            bc_z*dbc_length_inv__dbx
//        dbc1_z__dby=                            bc_z*dbc_length_inv__dby
//        dbc1_z__dbz= bc_length_inv*dbc_z__dbz + bc_z*dbc_length_inv__dbz
//        dbc1_z__dcx=                            bc_z*dbc_length_inv__dcx
//        dbc1_z__dcy=                            bc_z*dbc_length_inv__dcy
//        dbc1_z__dcz= bc_length_inv*dbc_z__dcz + bc_z*dbc_length_inv__dcz

      // mtx = mtx * sca + vec outer vec
      const matrix3d dbc1__da = matrix3d(); //FIXME remove
      const matrix3d dbc1__db = dbc__db*bc_length_inv + bc.outer(dbc_length_inv__db);
      const matrix3d dbc1__dc = dbc__dc*bc_length_inv + bc.outer(dbc_length_inv__dc);
      const matrix3d dbc1__dd = matrix3d(); //FIXME remove

// c abc_x=-ab_y*bc_z + ab_z*bc_y
//       dabc_x__dab_y=-bc_z
//       dabc_x__dbc_z=-ab_y
//       dabc_x__dab_z=bc_y
//       dabc_x__dbc_y=ab_z
// c abc_y=-ab_z*bc_x + ab_x*bc_z
//       dabc_y__dab_z=-bc_x
//       dabc_y__dbc_x=-ab_z
//       dabc_y__dab_x=bc_z
//       dabc_y__dbc_z=ab_x
// c abc_z=-ab_x*bc_y + ab_y*bc_x
//       dabc_z__dab_x=-bc_y
//       dabc_z__dbc_y=-ab_x
//       dabc_z__dab_y=bc_x
//       dabc_z__dbc_x=ab_y

    const matrix3d dabc__dab = matrix3d(0,-bc[2],bc[1],
                                        bc[2],0,-bc[0],
                                        -bc[1],bc[0],0); 
    const matrix3d dabc__dbc = matrix3d(0,ab[2],-ab[1],
                                        -ab[2],0,ab[0],
                                        ab[1],-ab[0],0); 
    const matrix3d dabc__dcd = matrix3d(); //FIXME remove

// c bcd_x=-bc_y*cd_z + bc_z*cd_y
//       dbcd_x__dbc_y=-cd_z
//       dbcd_x__dcd_z=-bc_y
//       dbcd_x__dbc_z=cd_y
//       dbcd_x__dcd_y=bc_z
// c bcd_y=-bc_z*cd_x + bc_x*cd_z
//       dbcd_y__dbc_z=-cd_x
//       dbcd_y__dcd_x=-bc_z
//       dbcd_y__dbc_x=cd_z
//       dbcd_y__dcd_z=bc_x
// c bcd_z=-bc_x*cd_y + bc_y*cd_x
//       dbcd_z__dbc_x=-cd_y
//       dbcd_z__dcd_y=-bc_x
//       dbcd_z__dbc_y=cd_x
//       dbcd_z__dcd_x=bc_y

    const matrix3d dbcd__dab = matrix3d(); //FIXME remove
    const matrix3d dbcd__dbc = matrix3d(0,-cd[2],cd[1],
                                        cd[2],0,-cd[0],
                                        -cd[1],cd[0],0); 
    const matrix3d dbcd__dcd = matrix3d(0,bc[2],-bc[1],
                                        -bc[2],0,bc[0],
                                        bc[1],-bc[0],0); 

// c abc_x=-ab_y*bc_z + ab_z*bc_y
// c      dabc_x__dax=0
//       dabc_x__day=dabc_x__dab_y*1
//       dabc_x__daz=dabc_x__dab_z*1
// c      dabc_x__dbx=0
//       dabc_x__dby=dabc_x__dab_y*dab_y__dby + dabc_x__dbc_y*dbc_y__dby
//       dabc_x__dbz=dabc_x__dbc_z*dbc_z__dbz + dabc_x__dab_z*dab_z__dbz
// c      dabc_x__dcx=0
//       dabc_x__dcy=dabc_x__dbc_y*dbc_y__dcy
//       dabc_x__dcz=dabc_x__dbc_z*dbc_z__dcz
// c abc_y=-ab_z*bc_x + ab_x*bc_z
//       dabc_y__dax=dabc_y__dab_x 
// c      dabc_y__day=0
//       dabc_y__daz=dabc_y__dab_z*1 
//       dabc_y__dbx=dabc_y__dbc_x + dabc_y__dab_x*dab_x__dbx
// c      dabc_y__dby=0
//       dabc_y__dbz=dabc_y__dab_z*dab_z__dbz + dabc_y__dbc_z*dbc_z__dbz
//       dabc_y__dcx=dabc_y__dbc_x*dbc_x__dcx
// c      dabc_y__dcy=0
//       dabc_y__dcz=dabc_y__dbc_z*dbc_z__dcz
// c abc_z=-ab_x*bc_y + ab_y*bc_x
//       dabc_z__dax=dabc_z__dab_x*1
//       dabc_z__day=dabc_z__dab_y*1 
//       dabc_z__daz=0
//       dabc_z__dbx=dabc_z__dbc_x + dabc_z__dab_x*dab_x__dbx
//       dabc_z__dby=dabc_z__dbc_y*dbc_y__dby + dabc_z__dab_y*dab_y__dby 
//       dabc_z__dbz=0 
//       dabc_z__dcx=dabc_z__dbc_x*dbc_x__dcx
//       dabc_z__dcy=dabc_z__dbc_y*dbc_y__dcy
//       dabc_z__dcz=0

         // mtx = mtx * mtx = mtx * mtx
         const matrix3d dabc__da = dabc__dbc*dbc__da + dabc__dcd*dcd__da;
         const matrix3d dabc__db = dabc__dbc*dbc__db + dabc__dcd*dcd__db;
         const matrix3d dabc__dc = dabc__dbc*dbc__dc + dabc__dcd*dcd__dc;
         const matrix3d dabc__dd = dabc__dbc*dbc__dd + dabc__dcd*dcd__dd;//FIXME remove

// c bcd_x=-bc_y*cd_z + bc_z*cd_y
// c      dbcd_x__dbx=0
//       dbcd_x__dby=dbcd_x__dbc_y*dbc_y__dby 
//       dbcd_x__dbz=dbcd_x__dbc_z*dbc_z__dbz 
// c      dbcd_x__dcx=0
//       dbcd_x__dcy=dbcd_x__dbc_y*dbc_y__dcy + dbcd_x__dcd_y*dcd_y__dcy
//       dbcd_x__dcz=dbcd_x__dcd_z*dcd_z__dcz + dbcd_x__dbc_z*dbc_z__dcz 
// c      dbcd_x__ddx=0
//       dbcd_x__ddy=dbcd_x__dcd_y*dcd_y__ddy
//       dbcd_x__ddz=dbcd_x__dcd_z*dcd_z__ddz
// c bcd_y=-bc_z*cd_x + bc_x*cd_z
//       dbcd_y__dbx=dbcd_y__dbc_x
// c      dbcd_y__dby=0
//       dbcd_y__dbz=dbcd_y__dbc_z*dbc_z__dbz
//       dbcd_y__dcx=dbcd_y__dcd_x*dcd_x__dcx + dbcd_y__dbc_x*dbc_x__dcx 
// c      dbcd_y__dcy=0
//       dbcd_y__dcz=dbcd_y__dbc_z*dbc_z__dcz + dbcd_y__dcd_z*dcd_z__dcz
//       dbcd_y__ddx=dbcd_y__dcd_x*dcd_x__ddx
// c      dbcd_y__ddy=0
//       dbcd_y__ddz=dbcd_y__dcd_z*dcd_z__ddz
// c bcd_z=-bc_x*cd_y + bc_y*cd_x
//       dbcd_z__dbx=dbcd_z__dbc_x 
//       dbcd_z__dby=dbcd_z__dbc_y*dbc_y__dby
// c      dbcd_z__dbz=0
//       dbcd_z__dcx=dbcd_z__dbc_x*dbc_x__dcx + dbcd_z__dcd_x*dcd_x__dcx
//       dbcd_z__dcy=dbcd_z__dcd_y*dcd_y__dcy + dbcd_z__dbc_y*dbc_y__dcy
// c      dbcd_z__dcz=0
//       dbcd_z__ddx=dbcd_z__dcd_x*dcd_x__ddx
//       dbcd_z__ddy=dbcd_z__dcd_y*dcd_y__ddy
// c      dbcd_z__ddz=0

         // mtx = mtx * mtx = mtx * mtx
         const matrix3d dbcd__da = dbcd__dbc*dbc__da + dbcd__dcd*dcd__da;//FIXME remove
         const matrix3d dbcd__db = dbcd__dbc*dbc__db + dbcd__dcd*dcd__db;
         const matrix3d dbcd__dc = dbcd__dbc*dbc__dc + dbcd__dcd*dcd__dc;
         const matrix3d dbcd__dd = dbcd__dbc*dbc__dd + dbcd__dcd*dcd__dd;

// c abc_length_inv=1/dsqrt(abc_x**2 + abc_y**2 + abc_z**2)
//       abc_length_inv_cub=abc_length_inv**3
//       dabc_length_inv__dabc_x=-abc_x*abc_length_inv_cub
//       dabc_length_inv__dabc_y=-abc_y*abc_length_inv_cub
//       dabc_length_inv__dabc_z=-abc_z*abc_length_inv_cub

        const double abc_length_inv_cub = pow(abc_length_inv,3);
        const coord3d dabc_length_inv__dabc = -abc*abc_length_inv_cub;

// c bcd_length_inv=1/dsqrt(bcd_x**2 + bcd_y**2 + bcd_z**2)
//       bcd_length_inv_cub=bcd_length_inv**3
//       dbcd_length_inv__dbcd_x=-bcd_x*bcd_length_inv_cub
//       dbcd_length_inv__dbcd_y=-bcd_y*bcd_length_inv_cub
//       dbcd_length_inv__dbcd_z=-bcd_z*bcd_length_inv_cub

        const double bcd_length_inv_cub = pow(bcd_length_inv,3);
        const coord3d dbcd_length_inv__dbcd = -bcd*bcd_length_inv_cub;

// c abc_length_inv=1/dsqrt(abc_x**2 + abc_y**2 + abc_z**2)
//       dabc_length_inv__dax=dabc_length_inv__dabc_y*dabc_y__dax
//      4 + dabc_length_inv__dabc_z*dabc_z__dax
//       dabc_length_inv__day=dabc_length_inv__dabc_x*dabc_x__day
//      4 + dabc_length_inv__dabc_z*dabc_z__day
//       dabc_length_inv__daz=dabc_length_inv__dabc_x*dabc_x__daz
//      3 + dabc_length_inv__dabc_y*dabc_y__daz
//       dabc_length_inv__dbx=dabc_length_inv__dabc_y*dabc_y__dbx
//      4 + dabc_length_inv__dabc_z*dabc_z__dbx
//       dabc_length_inv__dby=dabc_length_inv__dabc_x*dabc_x__dby
//      4 + dabc_length_inv__dabc_z*dabc_z__dby
//       dabc_length_inv__dbz=dabc_length_inv__dabc_x*dabc_x__dbz
//      3 + dabc_length_inv__dabc_y*dabc_y__dbz
//       dabc_length_inv__dcx=dabc_length_inv__dabc_y*dabc_y__dcx
//      4 + dabc_length_inv__dabc_z*dabc_z__dcx
//       dabc_length_inv__dcy=dabc_length_inv__dabc_x*dabc_x__dcy
//      4 + dabc_length_inv__dabc_z*dabc_z__dcy
//       dabc_length_inv__dcz=dabc_length_inv__dabc_x*dabc_x__dcz
//      3 + dabc_length_inv__dabc_y*dabc_y__dcz
// 
         // vec = mtx*vec + mtx*vec
         const coord3d dabc_length_inv__da = dabc__da*dabc_length_inv__dabc + dabc__da*dabc_length_inv__dabc;
         const coord3d dabc_length_inv__db = dabc__db*dabc_length_inv__dabc + dabc__db*dabc_length_inv__dabc;
         const coord3d dabc_length_inv__dc = dabc__dc*dabc_length_inv__dabc + dabc__dc*dabc_length_inv__dabc;
         const coord3d dabc_length_inv__dd = dabc__dd*dabc_length_inv__dabc + dabc__dd*dabc_length_inv__dabc; // FIXME remove

// c bcd_length_inv=1/dsqrt(bcd_x**2 + bcd_y**2 + bcd_z**2)
// c derivatives according to dax, day, daz
//       dbcd_length_inv__dbx=dbcd_length_inv__dbcd_y*dbcd_y__dbx
//      4 + dbcd_length_inv__dbcd_z*dbcd_z__dbx
//       dbcd_length_inv__dby=dbcd_length_inv__dbcd_x*dbcd_x__dby
//      4 + dbcd_length_inv__dbcd_z*dbcd_z__dby
//       dbcd_length_inv__dbz=dbcd_length_inv__dbcd_x*dbcd_x__dbz
//      3 + dbcd_length_inv__dbcd_y*dbcd_y__dbz
//       dbcd_length_inv__dcx=dbcd_length_inv__dbcd_y*dbcd_y__dcx
//      4 + dbcd_length_inv__dbcd_z*dbcd_z__dcx
//       dbcd_length_inv__dcy=dbcd_length_inv__dbcd_x*dbcd_x__dcy
//      4 + dbcd_length_inv__dbcd_z*dbcd_z__dcy
//       dbcd_length_inv__dcz=dbcd_length_inv__dbcd_x*dbcd_x__dcz
//      3 + dbcd_length_inv__dbcd_y*dbcd_y__dcz
//       dbcd_length_inv__ddx=dbcd_length_inv__dbcd_y*dbcd_y__ddx
//      4 + dbcd_length_inv__dbcd_z*dbcd_z__ddx
//       dbcd_length_inv__ddy=dbcd_length_inv__dbcd_x*dbcd_x__ddy
//      4 + dbcd_length_inv__dbcd_z*dbcd_z__ddy
//       dbcd_length_inv__ddz=dbcd_length_inv__dbcd_x*dbcd_x__ddz
//      3 + dbcd_length_inv__dbcd_y*dbcd_y__ddz

         // vec = mtx*vec + mtx*vec
         const coord3d dbcd_length_inv__da = dbcd__da*dbcd_length_inv__dbcd + dbcd__da*dbcd_length_inv__dbcd; // FIXME remove
         const coord3d dbcd_length_inv__db = dbcd__db*dbcd_length_inv__dbcd + dbcd__db*dbcd_length_inv__dbcd;
         const coord3d dbcd_length_inv__dc = dbcd__dc*dbcd_length_inv__dbcd + dbcd__dc*dbcd_length_inv__dbcd;
         const coord3d dbcd_length_inv__dd = dbcd__dd*dbcd_length_inv__dbcd + dbcd__dd*dbcd_length_inv__dbcd;

// c abc1_x=abc_x*abc_length_inv
// c abc1_y=abc_y*abc_length_inv
// c abc1_z=abc_z*abc_length_inv
//       dabc1_x__dabc_x=abc_length_inv
//       dabc1_y__dabc_y=abc_length_inv
//       dabc1_z__dabc_z=abc_length_inv
// 
//       dabc1_x__dabc_length_inv=abc_x
//       dabc1_y__dabc_length_inv=abc_y
//       dabc1_z__dabc_length_inv=abc_z

       const matrix3d dabc1__dabc=matrix3d(1,0,0,0,1,0,0,0,1) * abc_length_inv;
       const coord3d dabc1__dabc_length_inv=abc;

// c derivation of the components of the normals
// c abc1_x=abc_x*abc_length_inv
// c abc1_y=abc_y*abc_length_inv
// c abc1_z=abc_z*abc_length_inv
//       dabc1_x__dax=
//      2 dabc1_x__dabc_length_inv*dabc_length_inv__dax
//       dabc1_y__dax=dabc1_y__dabc_y*dabc_y__dax +
//      2 dabc1_y__dabc_length_inv*dabc_length_inv__dax
//       dabc1_z__dax=dabc1_z__dabc_z*dabc_z__dax +
//      2 dabc1_z__dabc_length_inv*dabc_length_inv__dax
//       dabc1_x__day=dabc1_x__dabc_x*dabc_x__day +
//      2 dabc1_x__dabc_length_inv*dabc_length_inv__day
//       dabc1_y__day=
//      2 dabc1_y__dabc_length_inv*dabc_length_inv__day
//       dabc1_z__day=dabc1_z__dabc_z*dabc_z__day +
//      2 dabc1_z__dabc_length_inv*dabc_length_inv__day
//       dabc1_x__daz=dabc1_x__dabc_x*dabc_x__daz +
//      2 dabc1_x__dabc_length_inv*dabc_length_inv__daz
//       dabc1_y__daz=dabc1_y__dabc_y*dabc_y__daz +
//      2 dabc1_y__dabc_length_inv*dabc_length_inv__daz
//       dabc1_z__daz=
//      2 dabc1_z__dabc_length_inv*dabc_length_inv__daz
// 
//       dabc1_x__dbx=
//      2 dabc1_x__dabc_length_inv*dabc_length_inv__dbx
//       dabc1_y__dbx=dabc1_y__dabc_y*dabc_y__dbx +
//      2 dabc1_y__dabc_length_inv*dabc_length_inv__dbx
//       dabc1_z__dbx=dabc1_z__dabc_z*dabc_z__dbx +
//      2 dabc1_z__dabc_length_inv*dabc_length_inv__dbx
//       dabc1_x__dby=dabc1_x__dabc_x*dabc_x__dby +
//      2 dabc1_x__dabc_length_inv*dabc_length_inv__dby
//       dabc1_y__dby=
//      2 dabc1_y__dabc_length_inv*dabc_length_inv__dby
//       dabc1_z__dby=dabc1_z__dabc_z*dabc_z__dby +
//      2 dabc1_z__dabc_length_inv*dabc_length_inv__dby
//       dabc1_x__dbz=dabc1_x__dabc_x*dabc_x__dbz +
//      2 dabc1_x__dabc_length_inv*dabc_length_inv__dbz
//       dabc1_y__dbz=dabc1_y__dabc_y*dabc_y__dbz +
//      2 dabc1_y__dabc_length_inv*dabc_length_inv__dbz
//       dabc1_z__dbz=
//      2 dabc1_z__dabc_length_inv*dabc_length_inv__dbz
// 
//       dabc1_x__dcx=
//      2 dabc1_x__dabc_length_inv*dabc_length_inv__dcx
//       dabc1_y__dcx=dabc1_y__dabc_y*dabc_y__dcx +
//      2 dabc1_y__dabc_length_inv*dabc_length_inv__dcx
//       dabc1_z__dcx=dabc1_z__dabc_z*dabc_z__dcx +
//      2 dabc1_z__dabc_length_inv*dabc_length_inv__dcx
//       dabc1_x__dcy=dabc1_x__dabc_x*dabc_x__dcy +
//      2 dabc1_x__dabc_length_inv*dabc_length_inv__dcy
//       dabc1_y__dcy=
//      2 dabc1_y__dabc_length_inv*dabc_length_inv__dcy
//       dabc1_z__dcy=dabc1_z__dabc_z*dabc_z__dcy +
//      2 dabc1_z__dabc_length_inv*dabc_length_inv__dcy
//       dabc1_x__dcz=dabc1_x__dabc_x*dabc_x__dcz +
//      2 dabc1_x__dabc_length_inv*dabc_length_inv__dcz
//       dabc1_y__dcz=dabc1_y__dabc_y*dabc_y__dcz +
//      2 dabc1_y__dabc_length_inv*dabc_length_inv__dcz
//       dabc1_z__dcz=
//      2 dabc1_z__dabc_length_inv*dabc_length_inv__dcz
// 
// c      dabc1_x__ddx=0
// c      dabc1_y__ddx=0
// c      dabc1_z__ddx=0
// c      dabc1_x__ddy=0
// c      dabc1_y__ddy=0
// c      dabc1_z__ddy=0
// c      dabc1_x__ddz=0
// c      dabc1_y__ddz=0
// c      dabc1_z__ddz=0

        // mtx = mtx * mtx + vec outer vec
        const matrix3d dabc1__da = dabc1__dabc*dabc__da + dabc1__dabc_length_inv.outer(dabc_length_inv__da); //FIXME cross is wrong
        const matrix3d dabc1__db = dabc1__dabc*dabc__db + dabc1__dabc_length_inv.outer(dabc_length_inv__db);
        const matrix3d dabc1__dc = dabc1__dabc*dabc__dc + dabc1__dabc_length_inv.outer(dabc_length_inv__dc);
        const matrix3d dabc1__dd = matrix3d();

// c bcd1_x=bcd_x*bcd_length_inv
// c bcd1_y=bcd_y*bcd_length_inv
// c bcd1_z=bcd_z*bcd_length_inv
//       dbcd1_x__dbcd_x=bcd_length_inv
//       dbcd1_y__dbcd_y=bcd_length_inv
//       dbcd1_z__dbcd_z=bcd_length_inv
// 
//       dbcd1_x__dbcd_length_inv=bcd_x
//       dbcd1_y__dbcd_length_inv=bcd_y
//       dbcd1_z__dbcd_length_inv=bcd_z

     const matrix3d dbcd1__dbcd = matrix3d(1,0,0,0,1,0,0,0,1) * bcd_length_inv;
     const coord3d dbcd1__dbcd_length_inv = bcd;

// c bcd1_x=bcd_x*bcd_length_inv
// c bcd1_y=bcd_y*bcd_length_inv
// c bcd1_z=bcd_z*bcd_length_inv
// 
// c      dbcd1_x__dax=0
// c      dbcd1_y__dax=0
// c      dbcd1_z__dax=0
// c      dbcd1_x__day=0
// c      dbcd1_y__day=0
// c      dbcd1_z__day=0
// c      dbcd1_x__daz=0
// c      dbcd1_y__daz=0
// c      dbcd1_z__daz=0
// 
//       dbcd1_x__dbx=
//      2 dbcd1_x__dbcd_length_inv*dbcd_length_inv__dbx
//       dbcd1_y__dbx=dbcd1_y__dbcd_y*dbcd_y__dbx +
//      2 dbcd1_y__dbcd_length_inv*dbcd_length_inv__dbx
//       dbcd1_z__dbx=dbcd1_z__dbcd_z*dbcd_z__dbx +
//      2 dbcd1_z__dbcd_length_inv*dbcd_length_inv__dbx
//       dbcd1_x__dby=dbcd1_x__dbcd_x*dbcd_x__dby +
//      2 dbcd1_x__dbcd_length_inv*dbcd_length_inv__dby
//       dbcd1_y__dby=
//      2 dbcd1_y__dbcd_length_inv*dbcd_length_inv__dby
//       dbcd1_z__dby=dbcd1_z__dbcd_z*dbcd_z__dby +
//      2 dbcd1_z__dbcd_length_inv*dbcd_length_inv__dby
//       dbcd1_x__dbz=dbcd1_x__dbcd_x*dbcd_x__dbz +
//      2 dbcd1_x__dbcd_length_inv*dbcd_length_inv__dbz
//       dbcd1_y__dbz=dbcd1_y__dbcd_y*dbcd_y__dbz +
//      2 dbcd1_y__dbcd_length_inv*dbcd_length_inv__dbz
//       dbcd1_z__dbz=
//      2 dbcd1_z__dbcd_length_inv*dbcd_length_inv__dbz
// 
//       dbcd1_x__dcx=
//      2 dbcd1_x__dbcd_length_inv*dbcd_length_inv__dcx
//       dbcd1_y__dcx=dbcd1_y__dbcd_y*dbcd_y__dcx +
//      2 dbcd1_y__dbcd_length_inv*dbcd_length_inv__dcx
//       dbcd1_z__dcx=dbcd1_z__dbcd_z*dbcd_z__dcx +
//      2 dbcd1_z__dbcd_length_inv*dbcd_length_inv__dcx
//       dbcd1_x__dcy=dbcd1_x__dbcd_x*dbcd_x__dcy +
//      2 dbcd1_x__dbcd_length_inv*dbcd_length_inv__dcy
//       dbcd1_y__dcy=
//      2 dbcd1_y__dbcd_length_inv*dbcd_length_inv__dcy
//       dbcd1_z__dcy=dbcd1_z__dbcd_z*dbcd_z__dcy +
//      2 dbcd1_z__dbcd_length_inv*dbcd_length_inv__dcy
//       dbcd1_x__dcz=dbcd1_x__dbcd_x*dbcd_x__dcz +
//      2 dbcd1_x__dbcd_length_inv*dbcd_length_inv__dcz
//       dbcd1_y__dcz=dbcd1_y__dbcd_y*dbcd_y__dcz +
//      2 dbcd1_y__dbcd_length_inv*dbcd_length_inv__dcz
//       dbcd1_z__dcz=
//      2 dbcd1_z__dbcd_length_inv*dbcd_length_inv__dcz
// 
//       dbcd1_x__ddx=
//      2 dbcd1_x__dbcd_length_inv*dbcd_length_inv__ddx
//       dbcd1_y__ddx=dbcd1_y__dbcd_y*dbcd_y__ddx +
//      2 dbcd1_y__dbcd_length_inv*dbcd_length_inv__ddx
//       dbcd1_z__ddx=dbcd1_z__dbcd_z*dbcd_z__ddx +
//      2 dbcd1_z__dbcd_length_inv*dbcd_length_inv__ddx
//       dbcd1_x__ddy=dbcd1_x__dbcd_x*dbcd_x__ddy +
//      2 dbcd1_x__dbcd_length_inv*dbcd_length_inv__ddy
//       dbcd1_y__ddy=
//      2 dbcd1_y__dbcd_length_inv*dbcd_length_inv__ddy
//       dbcd1_z__ddy=dbcd1_z__dbcd_z*dbcd_z__ddy +
//      2 dbcd1_z__dbcd_length_inv*dbcd_length_inv__ddy
//       dbcd1_x__ddz=dbcd1_x__dbcd_x*dbcd_x__ddz +
//      2 dbcd1_x__dbcd_length_inv*dbcd_length_inv__ddz
//       dbcd1_y__ddz=dbcd1_y__dbcd_y*dbcd_y__ddz +
//      2 dbcd1_y__dbcd_length_inv*dbcd_length_inv__ddz
//       dbcd1_z__ddz=
//      2 dbcd1_z__dbcd_length_inv*dbcd_length_inv__ddz

       // mtx = mtx*mtx + vec outer vec
       const matrix3d dbcd1__da = matrix3d();
       const matrix3d dbcd1__db = dbcd1__dbcd * dbcd__db + dbcd1__dbcd_length_inv.outer(dbcd_length_inv__db);//FIXME cross is wrong
       const matrix3d dbcd1__dc = dbcd1__dbcd * dbcd__dc + dbcd1__dbcd_length_inv.outer(dbcd_length_inv__dc);
       const matrix3d dbcd1__dd = dbcd1__dbcd * dbcd__dd + dbcd1__dbcd_length_inv.outer(dbcd_length_inv__dd);

// c aux_x=abc1_y*bc1_z-bc1_y*abc1_z
// c aux_y=abc1_z*bc1_x-bc1_z*abc1_x
// c aux_z=abc1_x*bc1_y-bc1_x*abc1_y
//       daux_x__dabc1_y=bc1_z
//       daux_x__dabc1_z=-bc1_y
//       daux_y__dabc1_z=bc1_x
//       daux_y__dabc1_x=-bc1_z
//       daux_z__dabc1_x=bc1_y
//       daux_z__dabc1_y=-bc1_x
// 
//       daux_x__dbc1_z=abc1_y
//       daux_x__dbc1_y=-abc1_z
//       daux_y__dbc1_x=abc1_z
//       daux_y__dbc1_z=-abc1_x
//       daux_z__dbc1_y=abc1_x
//       daux_z__dbc1_x=-abc1_y

         const matrix3d daux__dabc1 = matrix3d(0,bc1[2],-bc1[1],
                                               -bc1[2],0,bc1[0],
                                               bc1[1],-bc1[0],0);
         const matrix3d daux__dbc1 = matrix3d(0,-abc1[2],abc1[1],
                                              abc1[2],0,-abc1[0],
                                              -abc1[1],abc1[0],0);

// c aux_x=abc1_y*bc1_z-bc1_y*abc1_z
//       daux_x__dax=
//      2 daux_x__dabc1_y*dabc1_y__dax + daux_x__dabc1_z*dabc1_z__dax
//       daux_x__day=
//      2 daux_x__dabc1_y*dabc1_y__day + daux_x__dabc1_z*dabc1_z__day
//       daux_x__daz=
//      2 daux_x__dabc1_y*dabc1_y__daz + daux_x__dabc1_z*dabc1_z__daz
//       daux_x__dbx=
//      2 daux_x__dabc1_y*dabc1_y__dbx + daux_x__dbc1_z*dbc1_z__dbx
//      4 + daux_x__dbc1_y*dbc1_y__dbx + daux_x__dabc1_z*dabc1_z__dbx
//       daux_x__dby=
//      2 daux_x__dabc1_y*dabc1_y__dby + daux_x__dbc1_z*dbc1_z__dby
//      4 + daux_x__dbc1_y*dbc1_y__dby + daux_x__dabc1_z*dabc1_z__dby
//       daux_x__dbz=
//      2 daux_x__dabc1_y*dabc1_y__dbz + daux_x__dbc1_z*dbc1_z__dbz
//      4 + daux_x__dbc1_y*dbc1_y__dbz + daux_x__dabc1_z*dabc1_z__dbz
//       daux_x__dcx=
//      2 daux_x__dabc1_y*dabc1_y__dcx + daux_x__dbc1_z*dbc1_z__dcx
//      4 + daux_x__dbc1_y*dbc1_y__dcx + daux_x__dabc1_z*dabc1_z__dcx
//       daux_x__dcy=
//      2 daux_x__dabc1_y*dabc1_y__dcy + daux_x__dbc1_z*dbc1_z__dcy
//      4 + daux_x__dbc1_y*dbc1_y__dcy + daux_x__dabc1_z*dabc1_z__dcy
//       daux_x__dcz=
//      2 daux_x__dabc1_y*dabc1_y__dcz + daux_x__dbc1_z*dbc1_z__dcz
//      4 + daux_x__dbc1_y*dbc1_y__dcz + daux_x__dabc1_z*dabc1_z__dcz
// c      daux_x__ddx=0
// c      daux_x__ddy=0
// c      daux_x__ddz=0
// c aux_y=abc1_z*bc1_x-bc1_z*abc1_x
//       daux_y__dax=
//      2 daux_y__dabc1_z*dabc1_z__dax + daux_y__dabc1_x*dabc1_x__dax
//       daux_y__day=
//      2 daux_y__dabc1_z*dabc1_z__day + daux_y__dabc1_x*dabc1_x__day
//       daux_y__daz=
//      2 daux_y__dabc1_z*dabc1_z__daz + daux_y__dabc1_x*dabc1_x__daz
//       daux_y__dbx=
//      2 daux_y__dabc1_z*dabc1_z__dbx + daux_y__dbc1_x*dbc1_x__dbx
//      4 + daux_y__dbc1_z*dbc1_z__dbx + daux_y__dabc1_x*dabc1_x__dbx
//       daux_y__dby=
//      2 daux_y__dabc1_z*dabc1_z__dby + daux_y__dbc1_x*dbc1_x__dby
//      4 + daux_y__dbc1_z*dbc1_z__dby + daux_y__dabc1_x*dabc1_x__dby
//       daux_y__dbz=
//      2 daux_y__dabc1_z*dabc1_z__dbz + daux_y__dbc1_x*dbc1_x__dbz
//      4 + daux_y__dbc1_z*dbc1_z__dbz + daux_y__dabc1_x*dabc1_x__dbz
//       daux_y__dcx=
//      2 daux_y__dabc1_z*dabc1_z__dcx + daux_y__dbc1_x*dbc1_x__dcx
//      4 + daux_y__dbc1_z*dbc1_z__dcx + daux_y__dabc1_x*dabc1_x__dcx
//       daux_y__dcy=
//      2 daux_y__dabc1_z*dabc1_z__dcy + daux_y__dbc1_x*dbc1_x__dcy
//      4 + daux_y__dbc1_z*dbc1_z__dcy + daux_y__dabc1_x*dabc1_x__dcy
//       daux_y__dcz=
//      2 daux_y__dabc1_z*dabc1_z__dcz + daux_y__dbc1_x*dbc1_x__dcz
//      4 + daux_y__dbc1_z*dbc1_z__dcz + daux_y__dabc1_x*dabc1_x__dcz
// c      daux_y__ddx=0
// c      daux_y__ddy=0
// c      daux_y__ddz=0
// c aux_z=abc1_x*bc1_y-bc1_x*abc1_y
//       daux_z__dax=
//      2 daux_z__dabc1_x*dabc1_x__dax + daux_z__dabc1_y*dabc1_y__dax
//       daux_z__day=
//      2 daux_z__dabc1_x*dabc1_x__day + daux_z__dabc1_y*dabc1_y__day
//       daux_z__daz=
//      2 daux_z__dabc1_x*dabc1_x__daz + daux_z__dabc1_y*dabc1_y__daz
//       daux_z__dbx=
//      2 daux_z__dabc1_x*dabc1_x__dbx + daux_z__dbc1_y*dbc1_y__dbx
//      4 + daux_z__dbc1_x*dbc1_x__dbx + daux_z__dabc1_y*dabc1_y__dbx
//       daux_z__dby=
//      2 daux_z__dabc1_x*dabc1_x__dby + daux_z__dbc1_y*dbc1_y__dby
//      4 + daux_z__dbc1_x*dbc1_x__dby + daux_z__dabc1_y*dabc1_y__dby
//       daux_z__dbz=
//      2 daux_z__dabc1_x*dabc1_x__dbz + daux_z__dbc1_y*dbc1_y__dbz
//      4 + daux_z__dbc1_x*dbc1_x__dbz + daux_z__dabc1_y*dabc1_y__dbz
//       daux_z__dcx=
//      2 daux_z__dabc1_x*dabc1_x__dcx + daux_z__dbc1_y*dbc1_y__dcx
//      4 + daux_z__dbc1_x*dbc1_x__dcx + daux_z__dabc1_y*dabc1_y__dcx
//       daux_z__dcy=
//      2 daux_z__dabc1_x*dabc1_x__dcy + daux_z__dbc1_y*dbc1_y__dcy
//      4 + daux_z__dbc1_x*dbc1_x__dcy + daux_z__dabc1_y*dabc1_y__dcy
//       daux_z__dcz=
//      2 daux_z__dabc1_x*dabc1_x__dcz + daux_z__dbc1_y*dbc1_y__dcz
//      4 + daux_z__dbc1_x*dbc1_x__dcz + daux_z__dabc1_y*dabc1_y__dcz
// c      daux_z__ddx=0
// c      daux_z__ddy=0
// c      daux_z__ddz=0

        // mtx = mtx*mtx + mtx*mtx
         const matrix3d daux__da = daux__dabc1 * dabc1__da + daux__dbc1 * dbc1__da;
         const matrix3d daux__db = daux__dabc1 * dabc1__db + daux__dbc1 * dbc1__db;
         const matrix3d daux__dc = daux__dabc1 * dabc1__dc + daux__dbc1 * dbc1__dc;
         const matrix3d daux__dd = matrix3d(); //FIXME remove

// c y=aux_x*bcd1_x + aux_y*bcd1_y + aux_z*bcd1_z
// c x=abc1_x*bcd1_x + abc1_y*bcd1_y + abc1_z*bcd1_z
//       dy__daux_x=bcd1_x
//       dy__daux_y=bcd1_y
//       dy__daux_z=bcd1_z
// 
//       dy__dbcd1_x=aux_x
//       dy__dbcd1_y=aux_y
//       dy__dbcd1_z=aux_z
// 
//       dx__dabc1_x=bcd1_x
//       dx__dabc1_y=bcd1_y
//       dx__dabc1_z=bcd1_z
// 
//       dx__dbcd1_x=abc1_x
//       dx__dbcd1_y=abc1_y
//       dx__dbcd1_z=abc1_z

         const coord3d dy__daux = bcd1;
         const coord3d dy__dbcd1 = aux;
         const coord3d dx__dabc1 = bcd1;
         const coord3d dx__dbcd1 = abc1;

// c derivation of y
// c y=aux_x*bcd1_x + aux_y*bcd1_y + aux_z*bcd1_z
//       dy__dax=
//      2 dy__daux_x*daux_x__dax + 
//      3 dy__daux_y*daux_y__dax + 
//      4 dy__daux_z*daux_z__dax  
//       dy__day=
//      2 dy__daux_x*daux_x__day + 
//      3 dy__daux_y*daux_y__day + 
//      4 dy__daux_z*daux_z__day 
//       dy__daz=
//      2 dy__daux_x*daux_x__daz + 
//      3 dy__daux_y*daux_y__daz + 
//      4 dy__daux_z*daux_z__daz 
//       dy__dbx=
//      2 dy__daux_x*daux_x__dbx + dy__dbcd1_x*dbcd1_x__dbx +
//      3 dy__daux_y*daux_y__dbx + dy__dbcd1_y*dbcd1_y__dbx +
//      4 dy__daux_z*daux_z__dbx + dy__dbcd1_z*dbcd1_z__dbx
//       dy__dby=
//      2 dy__daux_x*daux_x__dby + dy__dbcd1_x*dbcd1_x__dby +
//      3 dy__daux_y*daux_y__dby + dy__dbcd1_y*dbcd1_y__dby +
//      4 dy__daux_z*daux_z__dby + dy__dbcd1_z*dbcd1_z__dby
//       dy__dbz=
//      2 dy__daux_x*daux_x__dbz + dy__dbcd1_x*dbcd1_x__dbz +
//      3 dy__daux_y*daux_y__dbz + dy__dbcd1_y*dbcd1_y__dbz +
//      4 dy__daux_z*daux_z__dbz + dy__dbcd1_z*dbcd1_z__dbz
//       dy__dcx=
//      2 dy__daux_x*daux_x__dcx + dy__dbcd1_x*dbcd1_x__dcx +
//      3 dy__daux_y*daux_y__dcx + dy__dbcd1_y*dbcd1_y__dcx +
//      4 dy__daux_z*daux_z__dcx + dy__dbcd1_z*dbcd1_z__dcx
//       dy__dcy=
//      2 dy__daux_x*daux_x__dcy + dy__dbcd1_x*dbcd1_x__dcy +
//      3 dy__daux_y*daux_y__dcy + dy__dbcd1_y*dbcd1_y__dcy +
//      4 dy__daux_z*daux_z__dcy + dy__dbcd1_z*dbcd1_z__dcy
//       dy__dcz=
//      2 dy__daux_x*daux_x__dcz + dy__dbcd1_x*dbcd1_x__dcz +
//      3 dy__daux_y*daux_y__dcz + dy__dbcd1_y*dbcd1_y__dcz +
//      4 dy__daux_z*daux_z__dcz + dy__dbcd1_z*dbcd1_z__dcz
//       dy__ddx=
//      2 dy__dbcd1_x*dbcd1_x__ddx +
//      3 dy__dbcd1_y*dbcd1_y__ddx +
//      4 dy__dbcd1_z*dbcd1_z__ddx
//       dy__ddy=
//      2 dy__dbcd1_x*dbcd1_x__ddy +
//      3 dy__dbcd1_y*dbcd1_y__ddy +
//      4 dy__dbcd1_z*dbcd1_z__ddy
//       dy__ddz=
//      2 dy__dbcd1_x*dbcd1_x__ddz +
//      3 dy__dbcd1_y*dbcd1_y__ddz +
//      4 dy__dbcd1_z*dbcd1_z__ddz


      //vec = mtx * vec
//      const coord3d dy__da = dy__daux.dot(daux__da);
      const coord3d dy__db = daux__db * dy__daux + dbcd1__db * dy__dbcd1;
      const coord3d dy__dc = daux__dc * dy__daux + dbcd1__dc * dy__dbcd1;
      const coord3d dy__dd =                       dbcd1__dd * dy__dbcd1;

// c derivation of x
// c x=abc1_x*bcd1_x + abc1_y*bcd1_y + abc1_z*bcd1_z
//       dx__dax=
//      2 dx__dabc1_x*dabc1_x__dax + 
//      3 dx__dabc1_y*dabc1_y__dax + 
//      4 dx__dabc1_z*dabc1_z__dax  
//       dx__day=
//      2 dx__dabc1_x*dabc1_x__day +
//      3 dx__dabc1_y*dabc1_y__day + 
//      4 dx__dabc1_z*dabc1_z__day 
//       dx__daz=
//      2 dx__dabc1_x*dabc1_x__daz + 
//      3 dx__dabc1_y*dabc1_y__daz + 
//      4 dx__dabc1_z*dabc1_z__daz 
//       dx__dbx=
//      2 dx__dabc1_x*dabc1_x__dbx + dx__dbcd1_x*dbcd1_x__dbx +
//      3 dx__dabc1_y*dabc1_y__dbx + dx__dbcd1_y*dbcd1_y__dbx +
//      4 dx__dabc1_z*dabc1_z__dbx + dx__dbcd1_z*dbcd1_z__dbx
//       dx__dby=
//      2 dx__dabc1_x*dabc1_x__dby + dx__dbcd1_x*dbcd1_x__dby +
//      3 dx__dabc1_y*dabc1_y__dby + dx__dbcd1_y*dbcd1_y__dby +
//      4 dx__dabc1_z*dabc1_z__dby + dx__dbcd1_z*dbcd1_z__dby
//       dx__dbz=
//      2 dx__dabc1_x*dabc1_x__dbz + dx__dbcd1_x*dbcd1_x__dbz +
//      3 dx__dabc1_y*dabc1_y__dbz + dx__dbcd1_y*dbcd1_y__dbz +
//      4 dx__dabc1_z*dabc1_z__dbz + dx__dbcd1_z*dbcd1_z__dbz
//       dx__dcx=
//      2 dx__dabc1_x*dabc1_x__dcx + dx__dbcd1_x*dbcd1_x__dcx +
//      3 dx__dabc1_y*dabc1_y__dcx + dx__dbcd1_y*dbcd1_y__dcx +
//      4 dx__dabc1_z*dabc1_z__dcx + dx__dbcd1_z*dbcd1_z__dcx
//       dx__dcy=
//      2 dx__dabc1_x*dabc1_x__dcy + dx__dbcd1_x*dbcd1_x__dcy +
//      3 dx__dabc1_y*dabc1_y__dcy + dx__dbcd1_y*dbcd1_y__dcy +
//      4 dx__dabc1_z*dabc1_z__dcy + dx__dbcd1_z*dbcd1_z__dcy
//       dx__dcz=
//      2 dx__dabc1_x*dabc1_x__dcz + dx__dbcd1_x*dbcd1_x__dcz +
//      3 dx__dabc1_y*dabc1_y__dcz + dx__dbcd1_y*dbcd1_y__dcz +
//      4 dx__dabc1_z*dabc1_z__dcz + dx__dbcd1_z*dbcd1_z__dcz
//       dx__ddx=
//      2 dx__dbcd1_x*dbcd1_x__ddx +
//      3 dx__dbcd1_y*dbcd1_y__ddx +
//      4 dx__dbcd1_z*dbcd1_z__ddx
//       dx__ddy=
//      2 dx__dbcd1_x*dbcd1_x__ddy +
//      3 dx__dbcd1_y*dbcd1_y__ddy +
//      4 dx__dbcd1_z*dbcd1_z__ddy
//       dx__ddz=
//      2 dx__dbcd1_x*dbcd1_x__ddz +
//      3 dx__dbcd1_y*dbcd1_y__ddz +
//      4 dx__dbcd1_z*dbcd1_z__ddz

      //vec = mtx * vec
//      const coord3d dx__da = dabc1__da * dx__dabc1;
      const coord3d dx__db = dabc1__db * dx__dabc1 + dbcd1__db * dx__dbcd1;
      const coord3d dx__dc = dabc1__dc * dx__dabc1 + dbcd1__dc * dx__dbcd1;
      const coord3d dx__dd =                         dbcd1__dd * dx__dbcd1;

// c derivation atan2(y,x) according to x and y
//       df__dx=-y/(x**2 + y**2)
//       df__dy=x/(x**2 + y**2)

      const double df__dx = -y/(x*x + y*y);
      const double df__dy =  x/(x*x + y*y);

// c derive f according to all 12 cartesion components of the four points
// c f=atan2(y,x)
//       df__dax=df__dx*dx__dax + df__dy*dy__dax
//       df__day=df__dx*dx__day + df__dy*dy__day
//       df__daz=df__dx*dx__daz + df__dy*dy__daz
//       df__dbx=df__dx*dx__dbx + df__dy*dy__dbx
//       df__dby=df__dx*dx__dby + df__dy*dy__dby
//       df__dbz=df__dx*dx__dbz + df__dy*dy__dbz
//       df__dcx=df__dx*dx__dcx + df__dy*dy__dcx
//       df__dcy=df__dx*dx__dcy + df__dy*dy__dcy
//       df__dcz=df__dx*dx__dcz + df__dy*dy__dcz
//       df__ddx=df__dx*dx__ddx + df__dy*dy__ddx
//       df__ddy=df__dx*dx__ddy + df__dy*dy__ddy
//       df__ddz=df__dx*dx__ddz + df__dy*dy__ddz

    // vec = vec*sca + vec*sca
    db=dx__db*df__dx + dy__db*df__dy;
    dc=dx__dc*df__dx + dy__dc*df__dy;
    dd=dx__dd*df__dx + dy__dd*df__dy;

}

double coord3d::ideal_dihedral(double lA, double lB, double lC){

  double theta_naught = 0.0;
  const double eps = 1.0e-10;
  if (1.0/lA + 1.0/lB + 1.0/lC > 0.5+eps){ // positive curvature // make sure 3 * 1/6 is recognised to be planar
    return 0.0; //FIXME: add actural functionality
  }
  else{ //planar or negative curvature
    return theta_naught;
  }

}

