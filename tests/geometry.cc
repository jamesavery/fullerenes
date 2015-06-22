#include "libgraph/geometry.hh"

int main(int ac, char **av)
{

// test unit matrix and unary -
  cout << "unit matrix stuff" << endl;
  matrix3d m = matrix3d::unit_matrix();
  matrix3d mm = -m;
  cout << "unit matrix: " << m << endl;
  cout << "minus unit matrix: " << mm << endl;


  coord3d a, b, c, d, da, db, dc, dd;

// test angles
  cout << endl << "angle stuff" << endl;
  a = coord3d(1,0,0), b = coord3d(0,0,0), c = coord3d(0,1,0);
  coord3d::dangle(a-b,c-b,da,dc);
  cout << "angle: " << coord3d::angle(a-b,c-b) << " [pi/2]" << endl;
  cout << "da: " << da << " [0, -1, 0]" << endl;
  cout << "db: " << -(da+dc) << " [1, 1, 0]" << endl;
  cout << "dc: " << dc << " [-1, 0, 0]" << endl;

  a = coord3d(2,0,0), b = coord3d(0,0,0), c = coord3d(0,1,0);
  coord3d::dangle(a-b,c-b,da,dc);
  cout << "angle: " << coord3d::angle(a-b,c-b) << " [pi/2]" << endl;
  cout << "da: " << da << " [0, -0.5, 0]" << endl;
  cout << "db: " << -(da+dc) << " [1, 0.5, 0]" << endl;
  cout << "dc: " << dc << " [-1, 0, 0]" << endl;

  a = coord3d(-0.5, 2.16, -1.08), b = coord3d(0.93, -1.76, -1.63), c = coord3d(2.74, 0.7, 1.12);
  coord3d::dangle(a-b,c-b,da,dc);
  cout << "angle: " << coord3d::angle(a-b,c-b) << " [1.0526064221783733]" << endl;
  cout << "da: " << da       << " [ -0.1664835282, -0.0375385702, -0.1653095568]" << endl;
  cout << "db: " << -(da+dc) << " [  0.0102169758,  0.2153727597,  0.1090805960]" << endl;
  cout << "dc: " << dc       << " [  0.1562665524, -0.1778341895,  0.0562289607]" << endl;


// test dihedrals
  cout << endl << "dihedral stuff" << endl;
  a = coord3d(1.0,0.0,0.0), b = coord3d(0.0,0.0,0.0), c = coord3d(0.0,1.0,0.0), d = coord3d(0.0,1.0,1.0);
  cout << "a: " << a << " b: " << b << " c: " << c << " d: " << d << endl;
  coord3d::ddihedral(b-a,c-a,d-a,db,dc,dd);
  cout << "dihedral: " << coord3d::dihedral(b-a,c-a,d-a) << " [pi/2]" << endl;
  cout << "da: " << -(db+dc+dd) << " [0,0,-1]" << endl;
  cout << "db: " << db << " [0,0,1]" << endl;
  cout << "dc: " << dc << " [1,0,0]" <<endl;
  cout << "dd: " << dd << " [-1,0,0]" << endl;

  a = coord3d(1.0,1.0,0.0), b = coord3d(0.0,0.0,0.0), c = coord3d(-1.0,1.0,0.0), d = coord3d(-1.0,1.0,1.0);
  cout << "a: " << a << " b: " << b << " c: " << c << " d: " << d << endl;
  coord3d::ddihedral(b-a,c-a,d-a,db,dc,dd);
  cout << "dihedral: " << coord3d::dihedral(b-a,c-a,d-a) << " [pi/2]" << endl;
  cout << "da: " << -(db+dc+dd) << " [0,0,-0.707]" << endl;
  cout << "db: " << db << " [0,0,0.707]" << endl;
  cout << "dc: " << dc << " [0.707,0.707,0]" <<endl;
  cout << "dd: " << dd << " [-0.707,-0.707,0]" << endl;

  a = coord3d(-0.5, 2.16, -1.08), b = coord3d(0.93, -1.76, -1.63), c = coord3d(2.74, 0.7, 1.12), d = coord3d(-1.0,1.0,1.0);
  cout << "a: " << a << " b: " << b << " c: " << c << " d: " << d << endl;
  coord3d::ddihedral(b-a,c-a,d-a,db,dc,dd);
  cout << "dihedral: " << coord3d::dihedral(b-a,c-a,d-a) << " [-0.685]" << endl;
  cout << "da: " << -(db+dc+dd) << " [-0.172, -0.090, 0.193]" << endl;
  cout << "db: " << db << " [0.0757, -0.0347, -0.0188]" << endl;
  cout << "dc: " << dc << " [0.0725, -0.084, 0.0284]" <<endl;
  cout << "dd: " << dd << " [0.023, 0.209, -0.202]" << endl;


// ideal dihedrals
  cout << endl << "ideal dihedrals" << endl;
  cout << "ideal dihedral between three triangles, counter clockwise/convex --> positive: " << coord3d::ideal_dihedral(3,3,3) << " [1.23096]" << endl;
  cout << "ideal dihedral between three squares, counter clockwise/convex --> positive: "   << coord3d::ideal_dihedral(4,4,4) << " [0.95531]" << endl;
  cout << "ideal dihedral between three pentagons, counter clockwise/convex --> positive: " << coord3d::ideal_dihedral(5,5,5) << " []" << endl;
  cout << "ideal dihedral between three hexagons, counter clockwise/convex --> positive: "  << coord3d::ideal_dihedral(6,6,6) << " [0.0]" << endl;
  cout << "ideal dihedral between three heptagons, counter clockwise/convex --> positive: " << coord3d::ideal_dihedral(7,7,7) << " [0.0]" << endl;

  cout << "ideal dihedral between 5,6,6: " << coord3d::ideal_dihedral(5,6,6) << ", " << 180/M_PI * coord3d::ideal_dihedral(5,6,6) << " deg" << endl;
  cout << "ideal dihedral between 6,5,6: " << coord3d::ideal_dihedral(6,5,6) << ", " << 180/M_PI * coord3d::ideal_dihedral(6,5,6) << " deg" << endl;
  cout << "ideal dihedral between 6,6,5: " << coord3d::ideal_dihedral(6,6,5) << ", " << 180/M_PI * coord3d::ideal_dihedral(6,6,5) << " deg" << endl;

  cout << "ideal dihedral between 5,5,6: " << coord3d::ideal_dihedral(5,5,6) << ", " << 180/M_PI * coord3d::ideal_dihedral(5,5,6) << " deg" << endl;
  cout << "ideal dihedral between 5,6,5: " << coord3d::ideal_dihedral(5,6,5) << ", " << 180/M_PI * coord3d::ideal_dihedral(5,6,5) << " deg" << endl;
  cout << "ideal dihedral between 6,5,5: " << coord3d::ideal_dihedral(6,5,5) << ", " << 180/M_PI * coord3d::ideal_dihedral(6,5,5) << " deg" << endl;

  cout << "ideal dihedral between 5,5,5: " << coord3d::ideal_dihedral(5,5,5) << ", " << 180/M_PI * coord3d::ideal_dihedral(5,5,5) << " deg" << endl;

  cout << "ideal dihedral between 6,5,6 (corrected): " << coord3d::ideal_dihedral(6,5,6,1.401,1.458,1.458) << ", " << 180/M_PI * coord3d::ideal_dihedral(6,5,6,1.401,1.458,1.458) << " deg" << endl;
  cout << "ideal dihedral between 5,6,5 (corrected): " << coord3d::ideal_dihedral(5,6,5,1.479,1.458,1.458) << ", " << 180/M_PI * coord3d::ideal_dihedral(5,6,5,1.479,1.458,1.458) << " deg" << endl;

  return 0;
}

