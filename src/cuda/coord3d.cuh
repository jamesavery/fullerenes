#ifndef DEVICE_COORD3D
#define DEVICE_COORD3D

template <typename T> INLINE std::array<T,3> operator-(const std::array<T,3>& a)                 { return {-a[0], -a[1], -a[2]};  }
template <typename T> INLINE std::array<T,3> operator-(const std::array<T,3>& a, const std::array<T,3>& b){ return {a[0]-b[0], a[1]-b[1], a[2]-b[2]};  }
template <typename T> INLINE std::array<T,3> operator-(const T a, const std::array<T,3>& b){ return {a - b[0], a - b[1], a - b[2]};  }
template <typename T> INLINE std::array<T,3> operator-(const std::array<T,3>& a, const T b){ return {a[0] - b, a[1] - b, a[2] - b};  }
template <typename T> INLINE std::array<T,3> operator+(const std::array<T,3>& a, const T b){ return {a[0] + b, a[1] +b, a[2] + b};}
template <typename T> INLINE std::array<T,3> operator+(const T a, const std::array<T,3>& b){ return b + a;}

template <typename T> INLINE std::array<T,3> operator+(const std::array<T,3>& a, const std::array<T,3>& b){ return {a[0]+b[0], a[1]+b[1], a[2]+b[2]};  }
template <typename T> INLINE std::array<T,3> operator*(const std::array<T,3>& a, const T s)  { return {a[0]*s, a[1]*s, a[2]*s};  }
template <typename T> INLINE std::array<T,3> operator*(const T s, const std::array<T,3>& a)  { return a*s; }
template <typename T> INLINE std::array<T,3> operator*(const std::array<T,3>& a, const std::array<T,3>& b) { return {a[0]*b[0], a[1]*b[1], a[2]*b[2]};}
template <typename T> INLINE std::array<T,3> operator/(const T s, const std::array<T,3>& a)  { return a*(T(1.)/s); }
template <typename T> INLINE std::array<T,3> operator/(const std::array<T,3>& a, const T s)  { return a*(T(1.)/s); }
template <typename T> INLINE std::array<T,3> operator/(const std::array<T,3>& a, const std::array<T,3>& b)  { return {a[0]/b[0], a[1]/b[1], a[2]/b[2]}; }
template <typename T> INLINE std::array<T,3> operator<(const std::array<T,3>& a, const T& b) {return {(T)(a[0] < b), (T)(a[1] < b), (T)(a[2] <b) };}
template <typename T> INLINE void operator+=(std::array<T,3>& a, const std::array<T,3>& b) {a = a + b;}
template <typename T> INLINE void operator-=(std::array<T,3>& a, const std::array<T,3>& b) {a = a - b;}
template <typename T> INLINE void operator/=(std::array<T,3>& a, const T b) {a = a / b;}
template <typename T> INLINE void operator*=(std::array<T,3>& a, const T b) {a = a * b;}
template <typename T> INLINE std::array<T,3> d_abs(const std::array<T,3>& a){ return {ABS(a[0]), ABS(a[1]), ABS(a[2])};}
template <typename T> INLINE void d_swap(T& a, T& b){
  T c = a;
  a = b;
  b = c;
}
template <typename T> INLINE std::array<T,3> d_sort(const std::array<T,3>& a){
  std::array<T,3> b = a;
  if(ABS(b[0]) > ABS(b[1])) d_swap(b[0], b[1]);
  if(ABS(b[1]) > ABS(b[2])) d_swap(b[1], b[2]);
  if(ABS(b[0]) > ABS(b[1])) d_swap(b[0], b[1]);
  return b;
}
template <typename T> INLINE std::array<T,3> rsqrt3(const std::array<T,3>& a){
  return {RSQRT(a[0]), RSQRT(a[1]), RSQRT(a[2])};
}
template <typename T> INLINE std::array<T,3> cos3(const std::array<T,3>& a){
  return {COS(a[0]), COS(a[1]), COS(a[2])};
}

template <typename T, typename K> INLINE void d_set(std::array<T,3>& a, const u_char j, K b){
  a[j] = T(b);
}

template <typename T> INLINE T d_get(const std::array<T,3>& a, const u_char j){
  return a[j];
}

//5 FLOPs
template <typename T> INLINE  T  dot(const std::array<T,3>& a,  const std::array<T,3>& b) { return a[0]*b[0] + a[1]*b[1] + a[2]*b[2]; }

//6 FLOPs
template <typename T> INLINE  T norm(const std::array<T,3>& a)                    { return SQRT(dot(a,a)); }

template <typename T> INLINE T sum(const std::array<T,3>& a) {return a[0] + a[1] + a[2];}

template <typename T> INLINE T max(const std::array<T,3>& a) {return d_max(d_max(a[0],a[1]),a[2]);}

template <typename T> INLINE T d_min(const T a, const T b){
  return a < b ? a : b;
}

template <typename T> INLINE T min(const std::array<T,3>& a) {return d_min(d_min(a[0],a[1]),a[2]);}

//7 FLOPs
template <typename T> INLINE  std::array<T,3> unit_vector(const std::array<T,3>& a){
  T r = RSQRT(dot(a,a));
  return (a*r);
}
//10 FLOPs
template <typename T> INLINE  std::array<T,3> cross(const std::array<T,3>& a, const std::array<T,3>& b){ return {a[1]*b[2]-a[2]*b[1],
							   -a[0]*b[2]+a[2]*b[0],
							   a[0]*b[1]-a[1]*b[0]}; }
// $(a \otimes b) \cdot c$
template <typename T> INLINE  std::array<T,3> outer_dot(const std::array<T,3>& a, const std::array<T,3>& b, const std::array<T,3>& c){
  return {a[0]*(b[0]*c[0] + b[1]*c[1] + b[2]*c[2]), 
  a[1]*(b[0]*c[0] + b[1]*c[1] + b[2]*c[2]), 
  a[2]*(b[0]*c[0] + b[1]*c[1] + b[2]*c[2])};
}

//6 FLOPs
template <typename T> INLINE  T bond_length(const std::array<T,3>& ab){
    return RSQRT(dot(ab,ab));
}

template <typename T> INLINE T non_resciprocal_bond_length(const std::array<T,3>& ab){
    return SQRT(dot(ab,ab));
}

template <typename T> INLINE void print_coord(const std::array<T,3>& ab){

    printf("[%.16e, %.16e, %.16e]\n",ab[0],ab[1],ab[2]);
}

template <typename T1, typename T2> INLINE void assign(std::array<T1,3>& a, const std::array<T2,3>& b){
  a[0] = b[0];
  a[1] = b[1];
  a[2] = b[2];
}


template <typename T> INLINE void assign(std::array<float,3>& a, const std::array<T,3>& b){
  a[0] = b[0];
  a[1] = b[1];
  a[2] = b[2];
}

template <typename T> INLINE void assign(std::array<double,3>& a, const std::array<T,3>& b){
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
