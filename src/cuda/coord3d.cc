#pragma once
#include <array>
/** FLOAT OPERATIONS **/
typedef float real_t;
typedef array<real_t,3> coord3d;

inline coord3d coord3d_a_to_coord3d(const float4& b) { {b.x, b.y, b.z};}
inline coord3d operator-(const coord3d& a)                 { return {-a.x, -a.y, -a.z};  }
inline coord3d operator-(const coord3d& a, const coord3d& b){ return {a.x-b.x, a.y-b.y, a.z-b.z};  }
inline coord3d operator+(const coord3d& a, const coord3d& b){ return {a.x+b.x, a.y+b.y, a.z+b.z};  }
inline coord3d operator*(const coord3d& a, const real_t s)  { return {a.x*s, a.y*s, a.z*s};  }
inline coord3d operator*(const real_t s, const coord3d& a)  { return a*s; }
inline coord3d operator*(const coord3d& a, const coord3d& b) { return {a.x*b.x, a.y*b.y, a.z*b.z};}
inline coord3d operator/(const real_t s, const coord3d& a)  { return a*(1/s); }
inline coord3d operator/(const coord3d& a, const real_t s)  { return a*(1/s); }
inline void operator+=(coord3d& a, const coord3d& b) {a = a + b;}
inline void operator/=(coord3d& a, const real_t b) {a = a / b;}

inline void d_set(coord3d& a, const u_char j, real_t b){
  ((real_t*)&a)[j] = b; 
}

inline real_t d_get(const coord3d& a, const u_char j){
  return ((const real_t*)&a)[j]; 
}
//5 FLOPs
inline  real_t  dot(const coord3d& a,  const coord3d& b) { return a.x*b.x + a.y*b.y + a.z*b.z; }

//6 FLOPs
inline  real_t norm(const coord3d& a)                    { return sqrt(dot(a,a)); }

inline real_t sum(const coord3d& a) {return a.x + a.y + a.z;}

inline real_t max(const coord3d& a) {return max(max(a.x,a.y),a.z);}

//7 FLOPs
inline  coord3d unit_vector(const coord3d& a){
  real_t r = (real_t)1.0/sqrt(dot(a,a));
  return (a*r);
}
//10 FLOPs
inline  coord3d cross(const coord3d& a, const coord3d& b){ return {a.y*b.z-a.z*b.y,
							   -a.x*b.z+a.z*b.x,
							   a.x*b.y-a.y*b.x}; }
// $(a \otimes b) \cdot c$
inline  coord3d outer_dot(const coord3d& a, const coord3d& b, const coord3d& c){
  return {a.x*(b.x*c.x + b.y*c.y + b.z*c.z), 
  a.y*(b.x*c.x + b.y*c.y + b.z*c.z), 
  a.z*(b.x*c.x + b.y*c.y + b.z*c.z)};
}

//6 FLOPs
inline  real_t bond_length(const coord3d& ab){
    return (real_t)1.0/sqrtf(dot(ab,ab));
}

inline real_t non_resciprocal_bond_length(const coord3d& ab){
    return sqrt(dot(ab,ab));
}

inline void print_coord(const coord3d& ab){

    printf("[%.16e, %.16e, %.16e]\n",ab.x,ab.y,ab.z);
}