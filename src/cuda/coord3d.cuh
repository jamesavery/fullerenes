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

INLINE device_real_t d_get(const device_coord3d& a, const u_char j){
  return ((const device_real_t*)&a)[j]; 
}
//5 FLOPs
INLINE  device_real_t  dot(const device_coord3d& a,  const device_coord3d& b) { return a.x*b.x + a.y*b.y + a.z*b.z; }

//6 FLOPs
INLINE  device_real_t norm(const device_coord3d& a)                    { return sqrt(dot(a,a)); }

INLINE device_real_t sum(const device_coord3d& a) {return a.x + a.y + a.z;}

INLINE device_real_t max(const device_coord3d& a) {return d_max(d_max(a.x,a.y),a.z);}

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
  return ((const device_node_t*)&a)[j];
}
#endif
