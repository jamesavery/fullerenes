#include "symmetry.hh"
using namespace std;

// Reference: Deza-2009, Theorem 2.1 (iii)
// For future reference: Theorem 2.2 lists symmetry groups for all triangulations
// of degree <= 6, arranged by signature (p3,p4,p5). 
PointGroup PointGroup::FullereneSymmetries[28] = {
  {C,1},      {C,2},       {C,REF_I},   {C,REF_S},
  {C,3},      {D,2},       {S,4},       {C,2,REF_V},
  {C,2,REF_H},{D,3},       {S,6},       {C,3,REF_V},
  {C,3,REF_H},{D,2,REF_H}, {D,2,REF_D}, {D,5},  
  {D,6},      {D,3,REF_H}, {D,3,REF_D}, {T},  
  {D,5,REF_H},{D,5,REF_D}, {D,6,REF_H}, {D,6,REF_D}, 
  {T,REF_D},  {T,REF_H},   {I},         {I,REF_H}
};


PointGroup::PointGroup(const Triangulation& T)
{
  // TODO: 
  // 1) Port the fortran symmetry detection routine, or
  // 2) Make the fortran routine callable from C++, or
  // 3) Implement the more efficient and general symmetry detection based on Brinkman's method.

  // Orders of fullerene point groups
  // 1: C1
  // 2: C2, Ci, Cs
  // 3: C3
  // 4: S4, D2, C2h, C2v
  // 6: S6, D3, C3h, C3v
  // 8: D2d, D2h
  // 10: D5
  // 12: T, D6, D3d, D3h
  // 20: D5d, D5h
  // 24: Th, Td, D6d, D6h
  // 60: I
  //120: Ih

}

string PointGroup::to_string() const {
  const char ts[7] = {'?','C','D','T','S','O','I'};
  const char rs[6] = {' ','v','h','d','i','s'};
  char result[4]   = {0,0,0,0};
  result[0] = ts[sym_type];
  result[1] = n>0? '0'+n : ' ';
  result[2] = rs[sym_reflection];
  return string(result);
}
