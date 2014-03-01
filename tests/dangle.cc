#include "libgraph/fullerenegraph.hh"
#include "libgraph/polyhedron.hh"
#include "libgraph/geometry.hh"

int main(int ac, char **av)
{

  coord3d a, b, c, da, dc;

  a = coord3d(1,0,0), b = coord3d(0,0,0), c = coord3d(0,1,0);
  coord3d::dangle(a-b,c-b,da,dc);
  cout << "angle: " << coord3d::angle(a-b,c-b) << endl;
  cout << "da: " << da << endl;
  cout << "dc: " << dc << endl;

  a = coord3d(2,0,0), b = coord3d(0,0,0), c = coord3d(0,1,0);
  coord3d::dangle(a-b,c-b,da,dc);
  cout << "angle: " << coord3d::angle(a-b,c-b) << endl;
  cout << "da: " << da << endl;
  cout << "dc: " << dc << endl;

  a = coord3d(-0.5, 2.16, -1.08), b = coord3d(0.93, -1.76, -1.63), c = coord3d(2.74, 0.7, 1.12);
  coord3d::dangle(a-b,c-b,da,dc);
  cout << "angle: " << coord3d::angle(a-b,c-b) << "(ought to be 1.0526064221783733)" << endl;
  cout << "da: " << da       << "  (ought to be -0.1664835282, -0.0375385702, -0.1653095568)" << endl;
  cout << "db: " << -(da+dc) << "  (ought to be  0.0102169758,  0.2153727597,  0.1090805960)" << endl;
  cout << "dc: " << dc       << "  (ought to be  0.1562665524, -0.1778341895,  0.0562289607)" << endl;


  return 0;
}
