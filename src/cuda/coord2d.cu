#define INLINE __device__ __forceinline__
#include <inttypes.h>

INLINE float2 operator-(const float2& a)                 { return make_float2(-a.x, -a.y);  }
INLINE float2 operator-(const float2& a, const float2& b){ return make_float2(a.x-b.x, a.y-b.y);  }
INLINE float2 operator+(const float2& a, const float2& b){ return make_float2(a.x+b.x, a.y+b.y);  }
INLINE float2 operator*(const float2& a, const float s)  { return make_float2(a.x*s, a.y*s);  }
INLINE float2 operator*(const float s, const float2& a)  { return a*s; }
INLINE float2 operator*(const float2& a, const float2& b) { return make_float2(a.x*b.x, a.y*b.y);}
INLINE float2 operator/(const float s, const float2& a)  { return a*(1/s); }
INLINE float2 operator/(const float2& a, const float s)  { return a*(1/s); }
INLINE void operator+=(float2& a, const float2& b) {a = a + b;}
INLINE void operator/=(float2& a, const float b) {a = a / b;}

INLINE  float  dot(const float2& a,  const float2& b) { return a.x*b.x + a.y*b.y; }
INLINE  float norm(const float2& a)                    { return sqrt(dot(a,a)); }
INLINE float sum(const float2& a) {return a.x + a.y;}
INLINE float d_get(const float2& a, const uint8_t j){
  return ((const float*)&a)[j]; 
}

INLINE double2 operator-(const double2& a)                 { return make_double2(-a.x, -a.y);  }
INLINE double2 operator-(const double2& a, const double2& b){ return make_double2(a.x-b.x, a.y-b.y);  }
INLINE double2 operator+(const double2& a, const double2& b){ return make_double2(a.x+b.x, a.y+b.y);  }
INLINE double2 operator*(const double2& a, const double s)  { return make_double2(a.x*s, a.y*s);  }
INLINE double2 operator*(const double s, const double2& a)  { return a*s; }
INLINE double2 operator*(const double2& a, const double2& b) { return make_double2(a.x*b.x, a.y*b.y);}
INLINE double2 operator/(const double s, const double2& a)  { return a*(1/s); }
INLINE double2 operator/(const double2& a, const double s)  { return a*(1/s); }
INLINE void operator+=(double2& a, const double2& b) {a = a + b;}
INLINE void operator/=(double2& a, const double b) {a = a / b;}

INLINE  double  dot(const double2& a,  const double2& b) { return a.x*b.x + a.y*b.y; }
INLINE  double norm(const double2& a)                    { return sqrt(dot(a,a)); }
INLINE double sum(const double2& a) {return a.x + a.y;}