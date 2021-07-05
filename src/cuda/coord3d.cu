#pragma once
typedef float3 coord3d;


__device__ __forceinline__ float3 operator-(const float3& a)                  { return make_float3(-a.x, -a.y, -a.z);  }
__device__ __forceinline__ float3 operator-(const float3& a, const float3& b){ return make_float3(a.x-b.x, a.y-b.y, a.z-b.z);  }
__device__ __forceinline__ float3 operator+(const float3& a, const float3& b){ return make_float3(a.x+b.x, a.y+b.y, a.z+b.z);  }
__device__ __forceinline__ float3 operator*(const float3& a, const float s)  { return make_float3(a.x*s, a.y*s, a.z*s);  }
__device__ __forceinline__ float3 operator*(const float s, const float3& a)  { return a*s; }
__device__ __forceinline__ float3 operator*(const float3& a, const float3& b) { return make_float3(a.x*b.x, a.y*b.y, a.z*b.z);}
__device__ __forceinline__ float3 operator/(const float s, const float3& a)  { return a*(1/s); }
__device__ __forceinline__ float3 operator/(const float3& a, const float s)  { return a*(1/s); }
__device__ __forceinline__ void operator+=(float3& a, const float3 b) {a = a + b;}
__device__ __forceinline__ void operator/=(float3& a, const float b) {a = a / b;}

__device__ __forceinline__ void set(float3& a, const uint8_t j, float b){
  switch (j)
  {
  case 0:
    a.x = b;
    break;
  case 1:
    a.y = b;
    break;
  case 2:
    a.z = b;
    break;
  default:
    break;
  }
}

__device__ __forceinline__ float get(const float3& a, const uint8_t j){
  switch (j)
  {
  case 0:
    return a.x;
  case 1:
    return a.y;
  case 2:
    return a.z;
  default:
    break;
  }
}
//5 FLOPs
__device__ __forceinline__  float  dot(const float3& a,  const float3& b) { return a.x*b.x + a.y*b.y + a.z*b.z; }

//6 FLOPs
__device__ __forceinline__  float norm(const float3& a)                    { return sqrt(dot(a,a)); }

//7 FLOPs
__device__ __forceinline__  float3 unit_vector(const float3& a){
  float r = rsqrt(dot(a,a));
  return (a*r);
}
//10 FLOPs
__device__ __forceinline__  float3 cross(const float3& a, const float3& b){ return make_float3(a.y*b.z-a.z*b.y,
							   -a.x*b.z+a.z*b.x,
							   a.x*b.y-a.y*b.x); }
// $(a \otimes b) \cdot c$
__device__ __forceinline__  float3 outer_dot(const float3& a, const float3& b, const float3& c){
  return make_float3(a.x*(b.x*c.x + b.y*c.y + b.z*c.z), 
  a.y*(b.x*c.x + b.y*c.y + b.z*c.z), 
  a.z*(b.x*c.x + b.y*c.y + b.z*c.z));
}

//6 FLOPs
__device__ __forceinline__  float bond_length(const float3& ab){
    return rsqrtf(dot(ab,ab));
}

__host__ __device__ void print_coord(const float3& ab){

    printf("[%.16e, %.16e, %.16e]\n",ab.x,ab.y,ab.z);
    //cout<< "[" << ab.x << ", " << ab.y << ", " << ab.z << "]\n" ;
}
/*
template <int N>
void write_to_file(const array<double3,N>& a){
    FILE* pFile;
    pFile = fopen("test.bin","wb");
    fwrite(&a, sizeof(float), N*3, pFile);
    fclose(pFile);
}*/


__device__ __forceinline__ double3 operator-(const double3& a)                  { return make_double3(-a.x, -a.y, -a.z);  }
__device__ __forceinline__ double3 operator-(const double3& a, const double3& b){ return make_double3(a.x-b.x, a.y-b.y, a.z-b.z);  }
__device__ __forceinline__ double3 operator+(const double3& a, const double3& b){ return make_double3(a.x+b.x, a.y+b.y, a.z+b.z);  }
__device__ __forceinline__ double3 operator*(const double3& a, const double s)  { return make_double3(a.x*s, a.y*s, a.z*s);  }
__device__ __forceinline__ double3 operator*(const double s, const double3& a)  { return a*s; }
__device__ __forceinline__ double3 operator*(const double3& a, const double3& b) { return make_double3(a.x*b.x, a.y*b.y, a.z*b.z);}
__device__ __forceinline__ double3 operator/(const double s, const double3& a)  { return a*(1/s); }
__device__ __forceinline__ double3 operator/(const double3& a, const double s)  { return a*(1/s); }
__device__ __forceinline__ void operator+=(double3& a, const double3 b) {a = a + b;}
__device__ __forceinline__ void operator/=(double3& a, const double b) {a = a / b;}

__device__ __forceinline__ void set(double3& a, const uint8_t j, double b){
  switch (j)
  {
  case 0:
    a.x = b;
    break;
  case 1:
    a.y = b;
    break;
  case 2:
    a.z = b;
    break;
  default:
    break;
  }
}

__device__ __forceinline__ double get(const double3& a, const uint8_t j){
  switch (j)
  {
  case 0:
    return a.x;
  case 1:
    return a.y;
  case 2:
    return a.z;
  default:
    break;
  }
}

//5 FLOPs
__device__ __forceinline__  double  dot(const double3& a,  const double3& b) { return a.x*b.x + a.y*b.y + a.z*b.z; }

//6 FLOPs
__device__ __forceinline__  double norm(const double3& a)                    { return sqrt(dot(a,a)); }

//7 FLOPs
__device__ __forceinline__  double3 unit_vector(const double3& a){
  double r = rsqrt(dot(a,a));
  return (a*r);
}
//10 FLOPs
__device__ __forceinline__  double3 cross(const double3& a, const double3& b){ return make_double3(a.y*b.z-a.z*b.y,
							   -a.x*b.z+a.z*b.x,
							   a.x*b.y-a.y*b.x); }
// $(a \otimes b) \cdot c$
__device__ __forceinline__  double3 outer_dot(const double3& a, const double3& b, const double3& c){
  return make_double3(a.x*(b.x*c.x + b.y*c.y + b.z*c.z), 
  a.y*(b.x*c.x + b.y*c.y + b.z*c.z), 
  a.z*(b.x*c.x + b.y*c.y + b.z*c.z));
}

//6 FLOPs
__device__ __forceinline__  double bond_length(const double3& ab){
    return rsqrtf(dot(ab,ab));
}

__host__ __device__ void print_coord(const double3& ab){

    printf("[%.16e, %.16e, %.16e]\n",ab.x,ab.y,ab.z);
    //cout<< "[" << ab.x << ", " << ab.y << ", " << ab.z << "]\n" ;
}

/*
template <int N>
void write_to_file(const array<double3,N>& a){
    FILE* pFile;
    pFile = fopen("test.bin","wb");
    fwrite(&a, sizeof(real_t), N*3, pFile);
    fclose(pFile);
}*/