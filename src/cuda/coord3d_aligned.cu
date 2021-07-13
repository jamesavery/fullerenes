#pragma once

//16byte Aligned coord3d;
typedef float4 coord3d_a;


INLINE float4 operator-(const float4& a)                  { return make_float4(-a.x, -a.y, -a.z, 0);  }
INLINE float4 operator-(const float4& a, const float4& b){ return make_float4(a.x-b.x, a.y-b.y, a.z-b.z, 0);  }
INLINE float4 operator+(const float4& a, const float4& b){ return make_float4(a.x+b.x, a.y+b.y, a.z+b.z, 0);  }
INLINE float4 operator*(const float4& a, const float s)  { return make_float4(a.x*s, a.y*s, a.z*s, 0);  }
INLINE float4 operator*(const float s, const float4& a)  { return a*s; }
INLINE float4 operator*(const float4& a, const float4& b) { return make_float4(a.x*b.x, a.y*b.y, a.z*b.z, 0);}
INLINE float4 operator/(const float s, const float4& a)  { return a*(1/s); }
INLINE float4 operator/(const float4& a, const float s)  { return a*(1/s); }
INLINE void operator+=(float4& a, const float4 b) {a = a + b;}
INLINE void operator/=(float4& a, const float b) {a = a / b;}

INLINE void set(float4& a, const uint8_t j, float b){
  ((float*)&a)[j] = b;
  /* switch (j) */
  /* { */
  /* case 0: */
  /*   a.x = b; */
  /*   break; */
  /* case 1: */
  /*   a.y = b; */
  /*   break; */
  /* case 2: */
  /*   a.z = b; */
  /*   break; */
  /* default: */
  /*   break; */
  /* } */
}

INLINE float get(const float4& a, const uint8_t j){
  return ((const float*)&a)[j];
  /* switch (j) */
  /* { */
  /* case 0: */
  /*   return a.x; */
  /* case 1: */
  /*   return a.y; */
  /* case 2: */
  /*   return a.z; */
  /* default: */
  /*   break; */
  /* } */
}
//5 FLOPs
INLINE  float  dot(const float4& a,  const float4& b) { return a.x*b.x + a.y*b.y + a.z*b.z; }

//6 FLOPs
INLINE  float norm(const float4& a)                    { return sqrt(dot(a,a)); }

//7 FLOPs
INLINE  float4 unit_vector(const float4& a){
  float r = rsqrt(dot(a,a));
  return (a*r);
}
//10 FLOPs
INLINE  float4 cross(const float4& a, const float4& b){ return make_float4(a.y*b.z-a.z*b.y,
							   -a.x*b.z+a.z*b.x,
							   a.x*b.y-a.y*b.x, 0); }
// $(a \otimes b) \cdot c$
INLINE  float4 outer_dot(const float4& a, const float4& b, const float4& c){
  return make_float4(a.x*(b.x*c.x + b.y*c.y + b.z*c.z), 
  a.y*(b.x*c.x + b.y*c.y + b.z*c.z), 
  a.z*(b.x*c.x + b.y*c.y + b.z*c.z),0);
}

//6 FLOPs
INLINE  float bond_length(const float4& ab){
    return rsqrtf(dot(ab,ab));
}

__host__ __device__ void print_coord(const float4& ab){

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
