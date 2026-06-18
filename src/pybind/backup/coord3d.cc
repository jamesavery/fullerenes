#pragma once

#include <array>
typedef double real_t;
typedef array<real_t,3> coord3d;

// inline is just a hint
#define INLINE __attribute__((always_inline)) inline

INLINE coord3d operator-(const coord3d& a)                  { return {-a[0], -a[1], -a[2]};  }
INLINE coord3d operator-(const coord3d& a, const coord3d& b){ return {a[0]-b[0], a[1]-b[1], a[2]-b[2]};  }
INLINE coord3d operator+(const coord3d& a, const coord3d& b){ return {a[0]+b[0], a[1]+b[1], a[2]+b[2]};  }
INLINE coord3d operator*(const coord3d& a, const real_t s)  { return {a[0]*s, a[1]*s, a[2]*s};  }
INLINE coord3d operator*(const real_t s, const coord3d& a)  { return a*s; }
INLINE coord3d operator/(const real_t s, const coord3d& a)  { return a*(1/s); }
INLINE coord3d operator/(const coord3d& a, const real_t s)  { return a*(1/s); }


INLINE real_t  dot(const coord3d& a,  const coord3d& b) { return a[0]*b[0] + a[1]*b[1] + a[2]*b[2]; }
INLINE real_t norm(const coord3d& a)                    { return sqrt(dot(a,a)); }

INLINE coord3d cross(const coord3d& a, const coord3d& b){ return {a[1]*b[2]-a[2]*b[1],
							   a[0]*b[2]-a[2]*b[0],
							   a[0]*b[1]-a[1]*b[0]}; }
// $(a \otimes b) \cdot c$
INLINE coord3d outer_dot(const coord3d& a, const coord3d& b, const coord3d& c){
  return {a[0]*(b[0]*c[0] + b[1]*c[1] + b[2]*c[2]), 
	  a[1]*(b[0]*c[0] + b[1]*c[1] + b[2]*c[2]), 
	  a[2]*(b[0]*c[0] + b[1]*c[1] + b[2]*c[2])};
}


// 
INLINE pair<real_t, coord3d> split_norm(const coord3d& a){
  real_t r = sqrt(dot(a,a));
  return {r, {a[0]/r, a[1]/r, a[2]/r}};
}
