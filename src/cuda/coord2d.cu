#define INLINE __device__ __forceinline__
#include <inttypes.h>

INLINE device_coord2d operator-(const device_coord2d& a)                 { return {-a.x, -a.y};  }
INLINE device_coord2d operator-(const device_coord2d& a, const device_coord2d& b){ return {a.x-b.x, a.y-b.y};  }
INLINE device_coord2d operator+(const device_coord2d& a, const device_coord2d& b){ return {a.x+b.x, a.y+b.y};  }
INLINE device_coord2d operator*(const device_coord2d& a, const device_real_t s)  { return {a.x*s, a.y*s};  }
INLINE device_coord2d operator*(const device_real_t s, const device_coord2d& a)  { return a*s; }
INLINE device_coord2d operator*(const device_coord2d& a, const device_coord2d& b) { return {a.x*b.x, a.y*b.y};}
INLINE device_coord2d operator/(const device_real_t s, const device_coord2d& a)  { return a*(1/s); }
INLINE device_coord2d operator/(const device_coord2d& a, const device_real_t s)  { return a*(1/s); }
INLINE void operator+=(device_coord2d& a, const device_coord2d& b) {a = a + b;}
INLINE void operator/=(device_coord2d& a, const device_real_t b) {a = a / b;}

INLINE  device_real_t  dot(const device_coord2d& a,  const device_coord2d& b) { return a.x*b.x + a.y*b.y; }
INLINE  device_real_t norm(const device_coord2d& a)                    { return sqrt(dot(a,a)); }
INLINE device_real_t sum(const device_coord2d& a) {return a.x + a.y;}
INLINE device_real_t d_get(const device_coord2d& a, const uint8_t j){
  return ((const device_real_t*)&a)[j]; 
}

