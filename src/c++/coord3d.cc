#pragma once

#include <array>
typedef double real_t;
typedef array<real_t,3> coord3d;

inline coord3d operator-(const coord3d& a)                  { return {-a[0], -a[1], -a[2]};  }
inline coord3d operator-(const coord3d& a, const coord3d& b){ return {a[0]-b[0], a[1]-b[1], a[2]-b[2]};  }
inline coord3d operator+(const coord3d& a, const coord3d& b){ return {a[0]+b[0], a[1]+b[1], a[2]+b[2]};  }
inline coord3d operator*(const coord3d& a, const real_t s)  { return {a[0]*s, a[1]*s, a[2]*s};  }
inline coord3d operator*(const real_t s, const coord3d& a)  { return a*s; }
inline coord3d operator/(const real_t s, const coord3d& a)  { return a*(1/s); }
inline coord3d operator/(const coord3d& a, const real_t s)  { return a*(1/s); }

inline real_t  dot(const coord3d& a,  const coord3d& b) { return a[0]*b[0] + a[1]*b[1] + a[2]*b[2]; }
inline real_t norm(const coord3d& a)                    { return sqrt(dot(a,a)); }

inline coord3d cross(const coord3d& a, const coord3d& b){ return {a[1]*b[2]-a[2]*b[1],
							   -a[0]*b[2]+a[2]*b[0],
							   a[0]*b[1]-a[1]*b[0]}; }
// $(a \otimes b) \cdot c$
inline coord3d outer_dot(const coord3d& a, const coord3d& b, const coord3d& c){
  return {a[0]*(b[0]*c[0] + b[1]*c[1] + b[2]*c[2]), 
	  a[1]*(b[0]*c[0] + b[1]*c[1] + b[2]*c[2]), 
	  a[2]*(b[0]*c[0] + b[1]*c[1] + b[2]*c[2])};
}

// 
inline pair<real_t, coord3d> split_norm(const coord3d& a){
  real_t r = sqrt(dot(a,a));
  return {r, {a[0], a[1], a[2]}};
}

inline real_t bond_length(const coord3d& ab){
    return sqrt(dot(ab,ab));
}

inline void print_coord(const coord3d& ab){
    cout<< "[" << ab[0] << ", " << ab[1] << ", " << ab[2] << "]\n" ;
}