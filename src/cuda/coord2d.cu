#define INLINE __device__ __forceinline__
#include <inttypes.h>

#define coord2d device_coord2d

INLINE coord2d operator-(const coord2d& a)                 { return {-a.x, -a.y};  }
INLINE coord2d operator-(const coord2d& a, const coord2d& b){ return {a.x-b.x, a.y-b.y};  }
INLINE coord2d operator+(const coord2d& a, const coord2d& b){ return {a.x+b.x, a.y+b.y};  }
INLINE coord2d operator*(const coord2d& a, const double s)  { return {a.x*s, a.y*s};  }
INLINE coord2d operator*(const double s, const coord2d& a)  { return a*s; }
INLINE coord2d operator*(const coord2d& a, const coord2d& b) { return {a.x*b.x, a.y*b.y};}
INLINE coord2d operator/(const double s, const coord2d& a)  { return a*(1/s); }
INLINE coord2d operator/(const coord2d& a, const double s)  { return a*(1/s); }
INLINE void operator+=(coord2d& a, const coord2d& b) {a = a + b;}
INLINE void operator/=(coord2d& a, const double b) {a = a / b;}

INLINE  double  dot(const coord2d& a,  const coord2d& b) { return a.x*b.x + a.y*b.y; }
INLINE  double norm(const coord2d& a)                    { return sqrt(dot(a,a)); }
INLINE double sum(const coord2d& a) {return a.x + a.y;}
