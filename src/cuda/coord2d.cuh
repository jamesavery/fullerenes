INLINE device_coord2d operator-(const device_coord2d& a)                 { return {-a[0], -a[1]};  }
INLINE device_coord2d operator-(const device_coord2d& a, const device_coord2d& b){ return {a[0]-b[0], a[1]-b[1]};  }
INLINE device_coord2d operator+(const device_coord2d& a, const device_coord2d& b){ return {a[0]+b[0], a[1]+b[1]};  }
INLINE device_coord2d operator*(const device_coord2d& a, const device_real_t s)  { return {a[0]*s, a[1]*s};  }
INLINE device_coord2d operator*(const device_real_t s, const device_coord2d& a)  { return a*s; }
INLINE device_coord2d operator*(const device_coord2d& a, const device_coord2d& b) { return {a[0]*b[0], a[1]*b[1]};}
INLINE device_coord2d operator/(const device_real_t s, const device_coord2d& a)  { return a*((device_real_t)1/s); }
INLINE device_coord2d operator/(const device_coord2d& a, const device_real_t s)  { return a*((device_real_t)1/s); }
INLINE void operator+=(device_coord2d& a, const device_coord2d& b) {a = a + b;}
INLINE void operator/=(device_coord2d& a, const device_real_t b) {a = a / b;}

INLINE  device_real_t  dot(const device_coord2d& a,  const device_coord2d& b) { return a[0]*b[0] + a[1]*b[1]; }
INLINE  device_real_t norm(const device_coord2d& a)                    { return SQRT(dot(a,a)); }
INLINE device_real_t sum(const device_coord2d& a) {return a[0] + a[1];}
INLINE device_real_t d_get(const device_coord2d& a, const uint8_t j){
  return ((const device_real_t*)&a)[j]; 
}
#if FLOAT_TYPE != 2
INLINE void assign(device_coord2d& a, const std::array<float,2>& b){
  a[0] = b[0];
  a[1] = b[1];
}
#endif
INLINE void assign(std::array<float,2>& a, const device_coord2d& b){
  a[0] = b[0];
  a[1] = b[1];
}
#if FLOAT_TYPE != 3
INLINE void assign(device_coord2d& a, const std::array<double,2>& b){
  a[0] = b[0];
  a[1] = b[1];
}
#endif
INLINE void assign(std::array<double,2>& a, const device_coord2d& b){
  a[0] = b[0];
  a[1] = b[1];
}