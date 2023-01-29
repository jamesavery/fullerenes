#ifndef DEVICE_COORD3D
#define DEVICE_COORD3D

INLINE device_coord3d operator-(const device_coord3d& a)                 { return {-a.x, -a.y, -a.z};  }
INLINE device_coord3d operator-(const device_coord3d& a, const device_coord3d& b){ return {a.x-b.x, a.y-b.y, a.z-b.z};  }
INLINE device_coord3d operator-(const device_real_t a, const device_coord3d& b){ return {a - b.x, a - b.y, a - b.z};  }
INLINE device_coord3d operator-(const device_coord3d& a, const device_real_t b){ return {a.x - b, a.y - b, a.z - b};  }
INLINE device_coord3d operator+(const device_coord3d& a, const device_real_t b){ return {a.x + b, a.y +b, a.z + b};}
INLINE device_coord3d operator+(const device_real_t a, const device_coord3d& b){ return b + a;}

INLINE device_coord3d operator+(const device_coord3d& a, const device_coord3d& b){ return {a.x+b.x, a.y+b.y, a.z+b.z};  }
INLINE device_coord3d operator*(const device_coord3d& a, const device_real_t s)  { return {a.x*s, a.y*s, a.z*s};  }
INLINE device_coord3d operator*(const device_real_t s, const device_coord3d& a)  { return a*s; }
INLINE device_coord3d operator*(const device_coord3d& a, const device_coord3d& b) { return {a.x*b.x, a.y*b.y, a.z*b.z};}
INLINE device_coord3d operator/(const device_real_t s, const device_coord3d& a)  { return a*(1/s); }
INLINE device_coord3d operator/(const device_coord3d& a, const device_real_t s)  { return a*(1/s); }
INLINE device_coord3d operator/(const device_coord3d& a, const device_coord3d& b)  { return {a.x/b.x, a.y/b.y, a.z/b.z}; }
INLINE device_coord3d operator<(const device_coord3d& a, const device_real_t& b) {return {(float)(a.x < b), (float)(a.y < b), (float)(a.z <b) };}
INLINE void operator+=(device_coord3d& a, const device_coord3d& b) {a = a + b;}
INLINE void operator-=(device_coord3d& a, const device_coord3d& b) {a = a - b;}
INLINE void operator/=(device_coord3d& a, const device_real_t b) {a = a / b;}
INLINE void operator*=(device_coord3d& a, const device_real_t b) {a = a * b;}

INLINE device_coord3d d_abs(const device_coord3d& a){ return {abs(a.x), abs(a.y), abs(a.z)};}
INLINE device_coord3d cos3(const device_coord3d& a){
  return {cos((double)a.x), cos((double)a.y), cos((double)a.z)};
}

INLINE void d_set(device_coord3d& a, const u_char j, device_real_t b){
  ((device_real_t*)&a)[j] = b; 
}

INLINE void d_set(uchar3& a, const u_char j, u_char b){
  ((u_char*)&a)[j] = b; 
}

INLINE void d_set(uchar4& a, const u_char j, u_char b){
  ((u_char*)&a)[j] = b; 
}

INLINE device_real_t d_get(const device_coord3d& a, const u_char j){
  return ((const device_real_t*)&a)[j]; 
}

INLINE constexpr uint8_t d_get(const uchar3& a, const u_char j){
  switch (j)
  {
  case 0:
    return a.x;
  case 1:
    return a.y;
  case 2:
    return a.z;
  }
}

INLINE constexpr uint8_t d_get(const uchar4& a, const u_char j){
  switch (j)
  {
  case 0:
    return a.x;
  case 1:
    return a.y;
  case 2:
    return a.z;
  }
}

//5 FLOPs
INLINE  device_real_t  dot(const device_coord3d& a,  const device_coord3d& b) { return a.x*b.x + a.y*b.y + a.z*b.z; }

//6 FLOPs
INLINE  device_real_t norm(const device_coord3d& a)                    { return sqrt(dot(a,a)); }

INLINE device_real_t sum(const device_coord3d& a) {return a.x + a.y + a.z;}

INLINE device_real_t max(const device_coord3d& a) {return d_max(d_max(a.x,a.y),a.z);}

INLINE device_real_t d_min(const device_real_t a, const device_real_t b){
  return a < b ? a : b;
}

INLINE device_real_t min(const device_coord3d& a) {return d_min(d_min(a.x,a.y),a.z);}

//7 FLOPs
INLINE  device_coord3d unit_vector(const device_coord3d& a){
  device_real_t r = (device_real_t)1.0/sqrt(dot(a,a));
  return (a*r);
}
//10 FLOPs
INLINE  device_coord3d cross(const device_coord3d& a, const device_coord3d& b){ return {a.y*b.z-a.z*b.y,
							   -a.x*b.z+a.z*b.x,
							   a.x*b.y-a.y*b.x}; }
// $(a \otimes b) \cdot c$
INLINE  device_coord3d outer_dot(const device_coord3d& a, const device_coord3d& b, const device_coord3d& c){
  return {a.x*(b.x*c.x + b.y*c.y + b.z*c.z), 
  a.y*(b.x*c.x + b.y*c.y + b.z*c.z), 
  a.z*(b.x*c.x + b.y*c.y + b.z*c.z)};
}

//6 FLOPs
INLINE  device_real_t bond_length(const device_coord3d& ab){
    return (device_real_t)1.0/sqrtf(dot(ab,ab));
}

INLINE device_real_t non_resciprocal_bond_length(const device_coord3d& ab){
    return sqrt(dot(ab,ab));
}

INLINE void print_coord(const device_coord3d& ab){

    printf("[%.16e, %.16e, %.16e]\n",ab.x,ab.y,ab.z);
}

INLINE device_node_t d_get(const device_node3& a, const uint8_t j){
  __builtin_assume(j < 3);
  return ((const device_node_t*)&a)[j];
}

struct UnitDyadic3{
  constexpr INLINE UnitDyadic3(){}
};

struct Mat3{
  device_real_t A[9] = {0,0,0,0,0,0,0,0,0};
  constexpr INLINE Mat3(device_real_t a, device_real_t b, device_real_t c, device_real_t d, device_real_t e, device_real_t f, device_real_t g, device_real_t h, device_real_t i){
    A[0] = a; A[1] = b; A[2] = c; A[3] = d; A[4] = e; A[5] = f; A[6] = g; A[7] = h; A[8] = i;
  }

  constexpr INLINE Mat3(){}

  INLINE Mat3(const Mat3& B){
     A[0] = B[0] ; A[1] = B[1] ; A[2] = B[2] ; A[3] = B[3] ; A[4] = B[4] ; A[5] = B[5] ; A[6] = B[6] ; A[7] = B[7] ; A[8] = B[8];
  }

  constexpr INLINE Mat3 zero_mat(){
    return Mat3(0,0,0,0,0,0,0,0,0);
  }

  
  INLINE device_real_t operator [](int i) const{
    return A[i];
  }

  INLINE device_real_t& operator [](int i){
    return A[i];
  }
};


INLINE Mat3 operator+ (const Mat3& A, const Mat3& B){
  return Mat3(A[0] + B[0], A[1] + B[1], A[2] + B[2], A[3] + B[3], A[4] + B[4], A[5] + B[5], A[6] + B[6], A[7] + B[7], A[8] + B[8]);
}

INLINE Mat3 operator- (const Mat3& A, const Mat3& B){
  return Mat3(A[0] - B[0], A[1] - B[1], A[2] - B[2], A[3] - B[3], A[4] - B[4], A[5] - B[5], A[6] - B[6], A[7] - B[7], A[8] - B[8]);
}

INLINE Mat3 operator* (const Mat3& A, const Mat3& B){
  return Mat3(A[0] * B[0], A[1] * B[1], A[2] * B[2], A[3] * B[3], A[4] * B[4], A[5] * B[5], A[6] * B[6], A[7] * B[7], A[8] * B[8]);
}

INLINE Mat3 operator/ (const Mat3& A, const Mat3& B){
  return Mat3(A[0] / B[0], A[1] / B[1], A[2] / B[2], A[3] / B[3], A[4] / B[4], A[5] / B[5], A[6] / B[6], A[7] / B[7], A[8] / B[8]);
}

INLINE Mat3 operator* (const UnitDyadic3& A , const Mat3& B){
  return Mat3(B);
}

INLINE Mat3 operator+ (const UnitDyadic3& A, const Mat3& B){
  return Mat3(1.f + B[0], B[1], B[2], B[3], 1.f + B[4], B[5], B[6], B[7], 1.f + B[8]);
}

INLINE Mat3 operator+ (const Mat3& A, const UnitDyadic3& B){
  return Mat3(1.f + A[0], A[1], A[2], A[3], 1.f + A[4], A[5], A[6], A[7], 1.f + A[8]);
}

INLINE Mat3 operator- (const UnitDyadic3& A, const Mat3& B){
  return Mat3(1.f - B[0], B[1], B[2], B[3], 1.f - B[4], B[5], B[6], B[7], 1.f - B[8]);
}

INLINE Mat3 operator- (const Mat3& A, const UnitDyadic3& B){
  return Mat3(A[0] - 1.f, A[1], A[2], A[3], A[4] - 1.f, A[5], A[6], A[7], A[8] - 1.f);
}

//Column wise cross product.
INLINE Mat3 cross(const Mat3& A, const device_coord3d& b){
  device_coord3d Aa = {A[0], A[3], A[6]},
                 Ab = {A[1], A[4], A[7]},
                 Ac = {A[2], A[5], A[8]};
  Aa = cross(Aa,b);
  Ab = cross(Ab,b);
  Ac = cross(Ac,b);
  return Mat3(Aa.x, Ab.x, Ac.x, Aa.y, Ab.y, Ac.y, Aa.z, Ab.z, Ac.z);
}

//Column wise cross product.
INLINE Mat3 cross(const device_coord3d& b, const Mat3& A){
  device_coord3d Aa = {A[0], A[3], A[6]},
                 Ab = {A[1], A[4], A[7]},
                 Ac = {A[2], A[5], A[8]};
  Aa = cross(b,Aa);
  Ab = cross(b,Ab);
  Ac = cross(b,Ac);
  return Mat3(Aa.x, Ab.x, Ac.x, Aa.y, Ab.y, Ac.y, Aa.z, Ab.z, Ac.z);
}

INLINE device_coord3d dot(const UnitDyadic3& A, const device_coord3d& b){
  return b;
}

INLINE device_coord3d dot(const device_coord3d& b, const UnitDyadic3& A){
  return b;
}

#endif
