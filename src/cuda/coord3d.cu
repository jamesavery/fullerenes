#pragma once

#ifndef DEVICE_COORD3D
#define DEVICE_COORD3D

#define INLINE __device__ __forceinline__

typedef float4 coord3d_a;



/** HALF OPERATIONS **/

__device__ __host__ __builtin_align__(8) struct half4
{
  __half2 x;
  __half2 y;
};



INLINE half4 coord3d_to_fp16(const float3& b) { return {{__float2half(b.x) , __float2half(b.y)}, { __float2half(b.z), __double2half(0.0)}};}
INLINE coord3d fp16_to_coord3d(const half4& b) { return {__half2float(b.x.x), __half2float(b.x.y), __half2float(b.y.x)};}
INLINE half4 operator-(const half4& a)                 { return {__hneg2(a.x), __hneg2(a.y)};  }
INLINE half4 operator-(const half4& a, const half4& b){ return { __hsub2(a.x,b.x),__hsub2(a.y,b.y)};  }
INLINE half4 operator+(const half4& a, const half4& b){ return { __hadd2(a.x,b.x), __hadd2(a.y,b.y)};  }
INLINE half4 operator*(const half4& a, const half s)  { return { __hmul2(a.x,make_half2(s,s)), __hmul2(a.y, make_half2(s,s)) };  }
INLINE half4 operator*(const half s, const half4& a)  { return a*s; }
INLINE half4 operator*(const half4& a, const half4& b) { return {__hmul2(a.x,b.x), __hmul2(a.y,b.y)};}
INLINE half4 operator/(const half s, const half4& a)  { return {__h2div(a.x,make_half2(s,s)), __h2div(a.y,make_half2(s,s))}; }
INLINE half4 operator/(const half4& a, const half s)  { return {__h2div(make_half2(s,s),a.x), __h2div(make_half2(s,s),a.y)}; }
INLINE void operator+=(half4& a, const half4& b) {a = a + b;}
INLINE void operator/=(half4& a, const half b) {a = a / b;}

INLINE void d_set(half4& a, const uint8_t j, half b){
  ((half*)&a)[j] = b; 
}

INLINE half d_get(const half4& a, const uint8_t j){
  return ((const half*)&a)[j]; 
}
//5 FLOPs
INLINE  half  dot(const half4& a,  const half4& b) { half2 temp = __hfma2(a.x,b.x,__hmul2(a.y,b.y)); return __hadd(temp.x,temp.y); }

//6 FLOPs
INLINE  half norm(const half4& a)                    { return hsqrt(dot(a,a)); }

//7 FLOPs
INLINE  half4 unit_vector(const half4& a){
  half r = hrsqrt(dot(a,a));
  return (a*r);
}
//10 FLOPs
INLINE  half4 cross(const half4& a, const half4& b){ return { {(a.x).y * (b.y).x-(a.y).x*(b.x).y,
							   -(a.x).x* (b.y).x+(a.y).x*(b.x).x},{
							   (a.x).x*(b.x).y-(a.x).y* (b.x).x, (a.y).y}}; }
// $(a \otimes b) \cdot c$

//6 FLOPs
INLINE  half bond_length(const half4& ab){
    return hrsqrt(dot(ab,ab));
}

INLINE half non_resciprocal_bond_length(const half4& ab){
    return hsqrt(dot(ab,ab));
}

INLINE void print_coord(const half4& ab){

    printf("[%.16e, %.16e, %.16e]\n",__half2float((ab.x).x),__half2float((ab.x).y),__half2float((ab.y).x));
}


/** 16BFLOAT OPERATIONS **/

#if defined(__CUDACC__) && (__CUDA_ARCH__ >= 800 || !defined(__CUDA_ARCH__))

typedef __nv_bfloat16 bhalf;
typedef __nv_bfloat162 bhalf2;
__device__ __host__ __builtin_align__(8) struct bhalf4
{
  bhalf2 x;
  bhalf2 y;
};



INLINE bhalf4 coord3d_to_bfp16(const float3& b) { return {{	__float2bfloat16(b.x) , 	__float2bfloat16(b.y)}, { 	__float2bfloat16(b.z), __double2bfloat16(0.0)}};}
INLINE coord3d bfp16_to_coord3d(const bhalf4& b) { return {		__bfloat162float(b.x.x), 		__bfloat162float(b.x.y), 		__bfloat162float(b.y.x)};}
INLINE bhalf4 operator-(const bhalf4& a)                 { return {__hneg2(a.x), __hneg2(a.y)};  }
INLINE bhalf4 operator-(const bhalf4& a, const bhalf4& b){ return { __hsub2(a.x,b.x),__hsub2(a.y,b.y)};  }
INLINE bhalf4 operator+(const bhalf4& a, const bhalf4& b){ return { __hadd2(a.x,b.x), __hadd2(a.y,b.y)};  }
INLINE bhalf4 operator*(const bhalf4& a, const bhalf s)  { return { __hmul2(a.x,{s,s}), __hmul2(a.y, {s,s}) };  }
INLINE bhalf4 operator*(const bhalf s, const bhalf4& a)  { return a*s; }
INLINE bhalf4 operator*(const bhalf4& a, const bhalf4& b) { return {__hmul2(a.x,b.x), __hmul2(a.y,b.y)};}
INLINE bhalf4 operator/(const bhalf s, const bhalf4& a)  { return {__h2div(a.x,{s,s}), __h2div(a.y,{s,s})}; }
INLINE bhalf4 operator/(const bhalf4& a, const bhalf s)  { return {__h2div({s,s},a.x), __h2div({s,s},a.y)}; }
INLINE void operator+=(bhalf4& a, const bhalf4& b) {a = a + b;}
INLINE void operator/=(bhalf4& a, const bhalf b) {a = a / b;}

INLINE void d_set(bhalf4& a, const uint8_t j, bhalf b){
  ((bhalf*)&a)[j] = b; 
}

INLINE bhalf d_get(const bhalf4& a, const uint8_t j){
  return ((const bhalf*)&a)[j]; 
}
//5 FLOPs
INLINE  bhalf  dot(const bhalf4& a,  const bhalf4& b) { bhalf2 temp = __hfma2(a.x,b.x,__hmul2(a.y,b.y)); return __hadd(temp.x,temp.y); }

//6 FLOPs
INLINE  bhalf norm(const bhalf4& a)                    { return hsqrt(dot(a,a)); }

//7 FLOPs
INLINE  bhalf4 unit_vector(const bhalf4& a){
  bhalf r = hrsqrt(dot(a,a));
  return (a*r);
}
//10 FLOPs
INLINE  bhalf4 cross(const bhalf4& a, const bhalf4& b){ return { {(a.x).y * (b.y).x-(a.y).x*(b.x).y,
							   -(a.x).x* (b.y).x+(a.y).x*(b.x).x},{
							   (a.x).x*(b.x).y-(a.x).y* (b.x).x, (a.y).y}}; }
// $(a \otimes b) \cdot c$

//6 FLOPs
INLINE  bhalf bond_length(const bhalf4& ab){
    return hrsqrt(dot(ab,ab));
}

INLINE bhalf non_resciprocal_bond_length(const bhalf4& ab){
    return hsqrt(dot(ab,ab));
}

INLINE void print_coord(const bhalf4& ab){

    printf("[%.16e, %.16e, %.16e]\n",__bfloat162float((ab.x).x),__bfloat162float((ab.x).y),__bfloat162float((ab.y).x));
}
#endif
/** FLOAT OPERATIONS **/

INLINE float3 coord3d_a_to_coord3d(const float4& b) { return make_float3(b.x, b.y, b.z); }
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

INLINE void d_set(float3& a, const u_char j, float b){
  ((float*)&a)[j] = b; 
}

INLINE float d_get(const float3& a, const u_char j){
  return ((const float*)&a)[j]; 
}
//5 FLOPs
INLINE  float  dot(const float3& a,  const float3& b) { return a.x*b.x + a.y*b.y + a.z*b.z; }

//6 FLOPs
INLINE  float norm(const float3& a)                    { return sqrt(dot(a,a)); }

INLINE float sum(const float3& a) {return a.x + a.y + a.z;}

INLINE float max(const float3& a) {return max(max(a.x,a.y),a.z);}

//7 FLOPs
INLINE  float3 unit_vector(const float3& a){
  float r = (float)1.0/sqrt(dot(a,a));
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
    return (float)1.0/sqrtf(dot(ab,ab));
}

INLINE float non_resciprocal_bond_length(const float3& ab){
    return sqrt(dot(ab,ab));
}

INLINE void print_coord(const float3& ab){

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

INLINE void d_set(double3& a, const u_char j, double b){
  ((double*)&a)[j] = b; 
}

INLINE double d_get(const double3& a, const u_char j){
  return ((const double*)&a)[j];
}

//5 FLOPs
INLINE  double  dot(const double3& a,  const double3& b) { return a.x*b.x + a.y*b.y + a.z*b.z; }

//6 FLOPs
INLINE  double norm(const double3& a)                    { return sqrt(dot(a,a)); }

INLINE double sum(const double3& a) {return a.x + a.y + a.z;}

INLINE double max(const double3& a) {return max(max(a.x,a.y),a.z);}

//7 FLOPs
INLINE  double3 unit_vector(const double3& a){
  double r = sqrt(dot(a,a));
  return (a/r);
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
    return 1.0/sqrt(dot(ab,ab));
}

INLINE double non_resciprocal_bond_length(const double3& ab){
    return sqrt(dot(ab,ab));
}

INLINE void print_coord(const double3& ab){

    printf("[%.16e, %.16e, %.16e]\n",ab.x,ab.y,ab.z);
}

INLINE device_node_t d_get(const device_node3& a, const uint8_t j){
  return ((const device_node_t*)&a)[j];
}
#endif