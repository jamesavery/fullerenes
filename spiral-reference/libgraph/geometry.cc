#include <limits.h>
#include <iostream>
#include <cstdio>
#include <vector>
#include <algorithm>
#include <cassert>
#include <complex>
#include "geometry.hh"

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

string filename_extension(const string& filename)
{
  size_t i = filename.rfind(".");
  bool found = i != string::npos;
  if(found) 
    return filename.substr(i,filename.size());
  else 
    return "";
}



matrix3d coord3d::outer(const coord3d& y) const
{
  return matrix3d(x[0]*y[0],x[0]*y[1],x[0]*y[2],  x[1]*y[0],x[1]*y[1],x[1]*y[2],  x[2]*y[0],x[2]*y[1],x[2]*y[2]);
}

coord3d coord3d::operator*(const matrix3d& m) const
{
  return coord3d(x[0]*m(0,0)+x[1]*m(1,0)+x[2]*m(2,0),  x[0]*m(0,1)+x[1]*m(1,1)+x[2]*m(2,1),  x[0]*m(0,2)+x[1]*m(1,2)+x[2]*m(2,2));
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
double coord3d::dihedral(const coord3d& b, const coord3d& c, const coord3d& d)
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
void coord3d::ddihedral(const coord3d& b, const coord3d& c, const coord3d& d, coord3d& db, coord3d& dc, coord3d& dd)
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
double coord3d::ideal_dihedral(const int A, const int B, const int C, const double ur, const double us, const double ut)
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

//    cout << "bx,by,bz: " << bx<<" " <<by<<" "<<bz << endl;
//    cout << "cx,cy,cz: " << cx<<" " <<cy<<" "<<cz << endl;
//    cout << "dx,dy,dz: " << dx<<" " <<dy<<" "<<dz << endl;

    return coord3d::dihedral(coord3d(bx,by,bz), coord3d(cx,cy,cz), coord3d(dx,dy,dz));
  } else {
	// planar or negative curvature
    // in the case of negative curvature the adjacent faces cannot be planar
    // and the dihedral cannot be calculated (but 0 is a plausible value)
    return 0.0;
  }
}

