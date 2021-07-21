#pragma once

typedef float3 coord3d;
typedef float real_t;
typedef uint16_t node_t;
typedef float4 coord3d_a;

#define INLINE __device__ __host__ __forceinline__

INLINE float3 coord3d_a_to_coord3d(const float4& b) { make_float3(b.x, b.y, b.z);}
INLINE float3 operator-(const float3& a)                 { return make_float3(-a.x, -a.y, -a.z);  }
INLINE float3 operator-(const float3& a, const float3& b){ return make_float3(a.x-b.x, a.y-b.y, a.z-b.z);  }
INLINE float3 operator+(const float3& a, const float3& b){ return make_float3(a.x+b.x, a.y+b.y, a.z+b.z);  }
INLINE float3 operator*(const float3& a, const float s)  { return make_float3(a.x*s, a.y*s, a.z*s);  }
INLINE float3 operator*(const float s, const float3& a)  { return a*s; }
INLINE float3 operator*(const float3& a, const float3& b) { return make_float3(a.x*b.x, a.y*b.y, a.z*b.z);}
INLINE float3 operator/(const float s, const float3& a)  { return a*(1/s); }
INLINE float3 operator/(const float3& a, const float s)  { return a*(1/s); }
INLINE void operator+=(float3& a, const float3& b) {a = a + b;}
INLINE void operator/=(float3& a, const float b) {a = a / b;}

INLINE void set(float3& a, const uint8_t j, float b){
  ((float*)&a)[j] = b; 
}

INLINE float get(const float3& a, const uint8_t j){
  return ((const float*)&a)[j]; 
}
//5 FLOPs
INLINE  float  dot(const float3& a,  const float3& b) { return a.x*b.x + a.y*b.y + a.z*b.z; }

//6 FLOPs
INLINE  float norm(const float3& a)                    { return sqrt(dot(a,a)); }

//7 FLOPs
INLINE  float3 unit_vector(const float3& a){
  float r = rsqrt(dot(a,a));
  return (a*r);
}
//10 FLOPs
INLINE  float3 cross(const float3& a, const float3& b){ return make_float3(a.y*b.z-a.z*b.y,
							   -a.x*b.z+a.z*b.x,
							   a.x*b.y-a.y*b.x); }
// $(a \otimes b) \cdot c$
INLINE  float3 outer_dot(const float3& a, const float3& b, const float3& c){
  return make_float3(a.x*(b.x*c.x + b.y*c.y + b.z*c.z), 
  a.y*(b.x*c.x + b.y*c.y + b.z*c.z), 
  a.z*(b.x*c.x + b.y*c.y + b.z*c.z));
}

//6 FLOPs
INLINE  float bond_length(const float3& ab){
    return rsqrtf(dot(ab,ab));
}

INLINE float non_resciprocal_bond_length(const float3& ab){
    return sqrtf(dot(ab,ab));
}

__host__ __device__ void print_coord(const float3& ab){

    printf("[%.16e, %.16e, %.16e]\n",ab.x,ab.y,ab.z);
}



INLINE double3 operator-(const double3& a)                  { return make_double3(-a.x, -a.y, -a.z);  }
INLINE double3 operator-(const double3& a, const double3& b){ return make_double3(a.x-b.x, a.y-b.y, a.z-b.z);  }
INLINE double3 operator+(const double3& a, const double3& b){ return make_double3(a.x+b.x, a.y+b.y, a.z+b.z);  }
INLINE double3 operator*(const double3& a, const double s)  { return make_double3(a.x*s, a.y*s, a.z*s);  }
INLINE double3 operator*(const double s, const double3& a)  { return a*s; }
INLINE double3 operator*(const double3& a, const double3& b) { return make_double3(a.x*b.x, a.y*b.y, a.z*b.z);}
INLINE double3 operator/(const double s, const double3& a)  { return a*(1/s); }
INLINE double3 operator/(const double3& a, const double s)  { return a*(1/s); }
INLINE void operator+=(double3& a, const double3 b) {a = a + b;}
INLINE void operator/=(double3& a, const double b) {a = a / b;}

INLINE void set(double3& a, const uint8_t j, double b){
  ((double*)&a)[j] = b; 
}

INLINE double get(const double3& a, const uint8_t j){
  return ((const double*)&a)[j];
}

//5 FLOPs
INLINE  double  dot(const double3& a,  const double3& b) { return a.x*b.x + a.y*b.y + a.z*b.z; }

//6 FLOPs
INLINE  double norm(const double3& a)                    { return sqrt(dot(a,a)); }

//7 FLOPs
INLINE  double3 unit_vector(const double3& a){
  double r = rsqrt(dot(a,a));
  return (a*r);
}
//10 FLOPs
INLINE  double3 cross(const double3& a, const double3& b){ return make_double3(a.y*b.z-a.z*b.y,
							   -a.x*b.z+a.z*b.x,
							   a.x*b.y-a.y*b.x); }
// $(a \otimes b) \cdot c$
INLINE  double3 outer_dot(const double3& a, const double3& b, const double3& c){
  return make_double3(a.x*(b.x*c.x + b.y*c.y + b.z*c.z), 
  a.y*(b.x*c.x + b.y*c.y + b.z*c.z), 
  a.z*(b.x*c.x + b.y*c.y + b.z*c.z));
}

//6 FLOPs
INLINE  double bond_length(const double3& ab){
    return rsqrtf(dot(ab,ab));
}

INLINE double non_resciprocal_bond_length(const double3& ab){
    return sqrtf(dot(ab,ab));
}

__host__ __device__ void print_coord(const double3& ab){

    printf("[%.16e, %.16e, %.16e]\n",ab.x,ab.y,ab.z);
}
