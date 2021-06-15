#pragma once

#include <array>
#include <fstream>
typedef double real_t;
typedef array<real_t,3> Coord3d;

inline Coord3d operator-(const Coord3d& a)                  { return {-a[0], -a[1], -a[2]};  }
inline Coord3d operator-(const Coord3d& a, const Coord3d& b){ return {a[0]-b[0], a[1]-b[1], a[2]-b[2]};  }
inline Coord3d operator+(const Coord3d& a, const Coord3d& b){ return {a[0]+b[0], a[1]+b[1], a[2]+b[2]};  }
inline Coord3d operator*(const Coord3d& a, const real_t s)  { return {a[0]*s, a[1]*s, a[2]*s};  }
inline Coord3d operator*(const real_t s, const Coord3d& a)  { return a*s; }
inline Coord3d operator*(const Coord3d& a, const Coord3d& b) { return {a[0]*b[0], a[1]*b[1], a[2]*b[2]};}
inline Coord3d operator/(const real_t s, const Coord3d& a)  { return a*(1/s); }
inline Coord3d operator/(const Coord3d& a, const real_t s)  { return a*(1/s); }
inline void operator+=(Coord3d& a, const Coord3d b) {a = a + b;}
inline void operator/=(Coord3d& a, const real_t s) {a = a / s;}

//5 FLOPs
inline real_t  dot(const Coord3d& a,  const Coord3d& b) { return a[0]*b[0] + a[1]*b[1] + a[2]*b[2]; }
inline real_t norm(const Coord3d& a)                    { return sqrt(dot(a,a)); }

inline Coord3d unit_vector(const Coord3d& a){
  real_t r = 1.0/sqrt(dot(a,a));
  return (a*r);
}

//10 FLOPs
inline Coord3d cross(const Coord3d& a, const Coord3d& b){ return {a[1]*b[2]-a[2]*b[1],
							   -a[0]*b[2]+a[2]*b[0],
							   a[0]*b[1]-a[1]*b[0]}; }
// $(a \otimes b) \cdot c$
inline Coord3d outer_dot(const Coord3d& a, const Coord3d& b, const Coord3d& c){
  return {a[0]*(b[0]*c[0] + b[1]*c[1] + b[2]*c[2]), 
	  a[1]*(b[0]*c[0] + b[1]*c[1] + b[2]*c[2]), 
	  a[2]*(b[0]*c[0] + b[1]*c[1] + b[2]*c[2])};
}

//7 FLOPs
inline tuple<real_t, Coord3d, Coord3d> split_norm(const Coord3d& a){
  real_t r = 1.0/sqrt(dot(a,a));
  return {r, a, a*r};
}

inline real_t bond_length(const Coord3d& ab){
    return 1.0/sqrt(dot(ab,ab));
}

inline void print_coord(const Coord3d& ab){

    printf("[%.16e, %.16e, %.16e]\n",ab[0],ab[1],ab[2]);
    //cout<< "[" << ab[0] << ", " << ab[1] << ", " << ab[2] << "]\n" ;
}

template <int N>
void write_to_file(const array<Coord3d,N>& a){
    FILE* pFile;
    pFile = fopen("test.bin","wb");
    fwrite(&a, sizeof(real_t), N*3, pFile);
    fclose(pFile);
}
