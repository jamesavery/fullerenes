template <typename T> INLINE std::array<T,2> operator-(const std::array<T,2>& a)                 { return {-a[0], -a[1]};  }
template <typename T> INLINE std::array<T,2> operator-(const std::array<T,2>& a, const std::array<T,2>& b){ return {a[0]-b[0], a[1]-b[1]};  }
template <typename T> INLINE std::array<T,2> operator+(const std::array<T,2>& a, const std::array<T,2>& b){ return {a[0]+b[0], a[1]+b[1]};  }
template <typename T> INLINE std::array<T,2> operator*(const std::array<T,2>& a, const T s)  { return {a[0]*s, a[1]*s};  }
template <typename T> INLINE std::array<T,2> operator*(const T s, const std::array<T,2>& a)  { return a*s; }
template <typename T> INLINE std::array<T,2> operator*(const std::array<T,2>& a, const std::array<T,2>& b) { return {a[0]*b[0], a[1]*b[1]};}
template <typename T> INLINE std::array<T,2> operator/(const T s, const std::array<T,2>& a)  { return a*((T)1/s); }
template <typename T> INLINE std::array<T,2> operator/(const std::array<T,2>& a, const T s)  { return a*((T)1/s); }
template <typename T> INLINE void operator+=(std::array<T,2>& a, const std::array<T,2>& b) {a = a + b;}
template <typename T> INLINE void operator/=(std::array<T,2>& a, const T b) {a = a / b;}

template <typename T> INLINE T  dot(const std::array<T,2>& a,  const std::array<T,2>& b) { return a[0]*b[0] + a[1]*b[1]; }
template <typename T> INLINE T norm(const std::array<T,2>& a)                    { return SQRT(dot(a,a)); }
template <typename T> INLINE T sum(const std::array<T,2>& a) {return a[0] + a[1];}
template <typename T> INLINE T d_get(const std::array<T,2>& a, const uint8_t j){
  return ((const T*)&a)[j]; 
}



template <typename T1, typename T2> INLINE void assign(std::array<T1,2>& a, const std::array<T2,2>& b){
  a[0] = b[0];
  a[1] = b[1];
}