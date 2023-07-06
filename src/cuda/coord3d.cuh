#ifndef DEVICE_COORD3D
#define DEVICE_COORD3D

INLINE device_coord3d operator-(const device_coord3d& a)                 { return {-a[0], -a[1], -a[2]};  }
INLINE device_coord3d operator-(const device_coord3d& a, const device_coord3d& b){ return {a[0]-b[0], a[1]-b[1], a[2]-b[2]};  }
INLINE device_coord3d operator-(const device_real_t a, const device_coord3d& b){ return {a - b[0], a - b[1], a - b[2]};  }
INLINE device_coord3d operator-(const device_coord3d& a, const device_real_t b){ return {a[0] - b, a[1] - b, a[2] - b};  }
INLINE device_coord3d operator+(const device_coord3d& a, const device_real_t b){ return {a[0] + b, a[1] +b, a[2] + b};}
INLINE device_coord3d operator+(const device_real_t a, const device_coord3d& b){ return b + a;}

INLINE device_coord3d operator+(const device_coord3d& a, const device_coord3d& b){ return {a[0]+b[0], a[1]+b[1], a[2]+b[2]};  }
INLINE device_coord3d operator*(const device_coord3d& a, const device_real_t s)  { return {a[0]*s, a[1]*s, a[2]*s};  }
INLINE device_coord3d operator*(const device_real_t s, const device_coord3d& a)  { return a*s; }
INLINE device_coord3d operator*(const device_coord3d& a, const device_coord3d& b) { return {a[0]*b[0], a[1]*b[1], a[2]*b[2]};}
INLINE device_coord3d operator/(const device_real_t s, const device_coord3d& a)  { return a*(device_real_t(1.)/s); }
INLINE device_coord3d operator/(const device_coord3d& a, const device_real_t s)  { return a*(device_real_t(1.)/s); }
INLINE device_coord3d operator/(const device_coord3d& a, const device_coord3d& b)  { return {a[0]/b[0], a[1]/b[1], a[2]/b[2]}; }
INLINE device_coord3d operator<(const device_coord3d& a, const device_real_t& b) {return {(device_real_t)(a[0] < b), (device_real_t)(a[1] < b), (device_real_t)(a[2] <b) };}
INLINE void operator+=(device_coord3d& a, const device_coord3d& b) {a = a + b;}
INLINE void operator-=(device_coord3d& a, const device_coord3d& b) {a = a - b;}
INLINE void operator/=(device_coord3d& a, const device_real_t b) {a = a / b;}
INLINE void operator*=(device_coord3d& a, const device_real_t b) {a = a * b;}
INLINE device_coord3d d_abs(const device_coord3d& a){ return {ABS(a[0]), ABS(a[1]), ABS(a[2])};}
INLINE void d_swap(device_real_t& a, device_real_t& b){
  device_real_t c = a;
  a = b;
  b = c;
}
INLINE device_coord3d d_sort(const device_coord3d& a){
  device_coord3d b = a;
  if(ABS(b[0]) > ABS(b[1])) d_swap(b[0], b[1]);
  if(ABS(b[1]) > ABS(b[2])) d_swap(b[1], b[2]);
  if(ABS(b[0]) > ABS(b[1])) d_swap(b[0], b[1]);
  return b;
}
INLINE device_coord3d rsqrt3(const device_coord3d& a){
  return {RSQRT(a[0]), RSQRT(a[1]), RSQRT(a[2])};
}
INLINE device_coord3d cos3(const device_coord3d& a){
  return {COS(a[0]), COS(a[1]), COS(a[2])};
}

INLINE void d_set(device_coord3d& a, const u_char j, device_real_t b){
  a[j] = b;
}

INLINE device_real_t d_get(const device_coord3d& a, const u_char j){
  return a[j];
}

//5 FLOPs
INLINE  device_real_t  dot(const device_coord3d& a,  const device_coord3d& b) { return a[0]*b[0] + a[1]*b[1] + a[2]*b[2]; }

//6 FLOPs
INLINE  device_real_t norm(const device_coord3d& a)                    { return SQRT(dot(a,a)); }

INLINE device_real_t sum(const device_coord3d& a) {return a[0] + a[1] + a[2];}

INLINE device_real_t max(const device_coord3d& a) {return d_max(d_max(a[0],a[1]),a[2]);}

INLINE device_real_t d_min(const device_real_t a, const device_real_t b){
  return a < b ? a : b;
}

INLINE device_real_t min(const device_coord3d& a) {return d_min(d_min(a[0],a[1]),a[2]);}

//7 FLOPs
INLINE  device_coord3d unit_vector(const device_coord3d& a){
  device_real_t r = RSQRT(dot(a,a));
  return (a*r);
}
//10 FLOPs
INLINE  device_coord3d cross(const device_coord3d& a, const device_coord3d& b){ return {a[1]*b[2]-a[2]*b[1],
							   -a[0]*b[2]+a[2]*b[0],
							   a[0]*b[1]-a[1]*b[0]}; }
// $(a \otimes b) \cdot c$
INLINE  device_coord3d outer_dot(const device_coord3d& a, const device_coord3d& b, const device_coord3d& c){
  return {a[0]*(b[0]*c[0] + b[1]*c[1] + b[2]*c[2]), 
  a[1]*(b[0]*c[0] + b[1]*c[1] + b[2]*c[2]), 
  a[2]*(b[0]*c[0] + b[1]*c[1] + b[2]*c[2])};
}

//6 FLOPs
INLINE  device_real_t bond_length(const device_coord3d& ab){
    return RSQRT(dot(ab,ab));
}

INLINE device_real_t non_resciprocal_bond_length(const device_coord3d& ab){
    return SQRT(dot(ab,ab));
}

INLINE void print_coord(const device_coord3d& ab){

    printf("[%.16e, %.16e, %.16e]\n",ab[0],ab[1],ab[2]);
}

INLINE device_node_t d_get(const device_node3& a, const uint8_t j){
  __builtin_assume(j < 3);
  return ((const device_node_t*)&a)[j];
}

#if FLOAT_TYPE != 3
INLINE void assign(device_coord3d& a, const std::array<double,3>& b){
  a[0] = b[0];
  a[1] = b[1];
  a[2] = b[2];
}
#endif

#if FLOAT_TYPE != 2
INLINE void assign(device_coord3d& a, const std::array<float,3>& b){
  a[0] = b[0];
  a[1] = b[1];
  a[2] = b[2];
}
#endif
INLINE void assign(std::array<float,3>& a, const device_coord3d& b){
  a[0] = b[0];
  a[1] = b[1];
  a[2] = b[2];
}

INLINE void assign(std::array<double,3>& a, const device_coord3d& b){
  a[0] = b[0];
  a[1] = b[1];
  a[2] = b[2];
}



INLINE void d_set(uchar3& a, const u_char j, u_char b){
  ((u_char*)&a)[j] = b; 
}

INLINE void d_set(uchar4& a, const u_char j, u_char b){
  ((u_char*)&a)[j] = b; 
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



#endif
