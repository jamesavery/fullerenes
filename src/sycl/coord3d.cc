#ifndef DEVICE_COORD3D
#define DEVICE_COORD3D
#include <array>
#include <limits>
#include <cmath>
#include <inttypes.h>


template <typename T> std::array<T,2> operator-(const std::array<T,2>& a)                 { return {-a[0], -a[1]};  }
template <typename T> std::array<T,2> operator-(const std::array<T,2>& a, const std::array<T,2>& b){ return {a[0]-b[0], a[1]-b[1]};  }
template <typename T> std::array<T,2> operator-(const T a, const std::array<T,2>& b){ return {a - b[0], a - b[1]};  }
template <typename T> std::array<T,2> operator-(const std::array<T,2>& a, const T b){ return {a[0] - b, a[1] - b};  }
template <typename T> std::array<T,2> operator+(const std::array<T,2>& a, const T b){ return {a[0] + b, a[1] +b};}
template <typename T> std::array<T,2> operator+(const T a, const std::array<T,2>& b){ return b + a;}
template <typename T> std::array<T,2> operator+(const std::array<T,2>& a, const std::array<T,2>& b){ return {a[0]+b[0], a[1]+b[1]};  }
template <typename T> std::array<T,2> operator*(const std::array<T,2>& a, const T s)  { return {a[0]*s, a[1]*s};  }
template <typename T> std::array<T,2> operator*(const T s, const std::array<T,2>& a)  { return a*s; }
template <typename T> std::array<T,2> operator*(const std::array<T,2>& a, const std::array<T,2>& b) { return {a[0]*b[0], a[1]*b[1]};}
template <typename T> std::array<T,2> operator/(const T s, const std::array<T,2>& a)  { return s/a[0], s/a[1]; }
template <typename T> std::array<T,2> operator/(const std::array<T,2>& a, const T s)  { return a*(T(1.)/s); }
template <typename T> void operator+=(std::array<T,2>& a, const std::array<T,2>& b) {a = a + b;}
template <typename T> void operator-=(std::array<T,2>& a, const std::array<T,2>& b) {a = a - b;}
template <typename T> void operator/=(std::array<T,2>& a, const T b) {a = a / b;}
template <typename T> void operator*=(std::array<T,2>& a, const T b) {a = a * b;}

template <typename T>  T dot(const std::array<T,2>& a,  const std::array<T,2>& b) { return a[0]*b[0] + a[1]*b[1]; }
template <typename T>  T norm(const std::array<T,2>& a)                    { return sycl::sqrt(dot(a,a)); }


template <typename T> std::array<T,3> operator-(const std::array<T,3>& a)                 { return {-a[0], -a[1], -a[2]};  }
template <typename T> std::array<T,3> operator-(const std::array<T,3>& a, const std::array<T,3>& b){ return {a[0]-b[0], a[1]-b[1], a[2]-b[2]};  }
template <typename T> std::array<T,3> operator-(const T a, const std::array<T,3>& b){ return {a - b[0], a - b[1], a - b[2]};  }
template <typename T> std::array<T,3> operator-(const std::array<T,3>& a, const T b){ return {a[0] - b, a[1] - b, a[2] - b};  }
template <typename T> std::array<T,3> operator+(const std::array<T,3>& a, const T b){ return {a[0] + b, a[1] +b, a[2] + b};}
template <typename T> std::array<T,3> operator+(const T a, const std::array<T,3>& b){ return b + a;}

template <typename T> std::array<T,3> operator+(const std::array<T,3>& a, const std::array<T,3>& b){ return {a[0]+b[0], a[1]+b[1], a[2]+b[2]};  }
template <typename T> std::array<T,3> operator*(const std::array<T,3>& a, const T s)  { return {a[0]*s, a[1]*s, a[2]*s};  }
template <typename T> std::array<T,3> operator*(const T s, const std::array<T,3>& a)  { return a*s; }
template <typename T> std::array<T,3> operator*(const std::array<T,3>& a, const std::array<T,3>& b) { return {a[0]*b[0], a[1]*b[1], a[2]*b[2]};}
template <typename T> std::array<T,3> operator/(const T s, const std::array<T,3>& a)  { return a*(T(1.)/s); }
template <typename T> std::array<T,3> operator/(const std::array<T,3>& a, const T s)  { return a*(T(1.)/s); }
template <typename T> std::array<T,3> operator/(const std::array<T,3>& a, const std::array<T,3>& b)  { return {a[0]/b[0], a[1]/b[1], a[2]/b[2]}; }
template <typename T> std::array<T,3> operator<(const std::array<T,3>& a, const T& b) {return {(T)(a[0] < b), (T)(a[1] < b), (T)(a[2] <b) };}
template <typename T> void operator+=(std::array<T,3>& a, const std::array<T,3>& b) {a = a + b;}
template <typename T> void operator-=(std::array<T,3>& a, const std::array<T,3>& b) {a = a - b;}
template <typename T> void operator/=(std::array<T,3>& a, const T b) {a = a / b;}
template <typename T> void operator*=(std::array<T,3>& a, const T b) {a = a * b;}
template <typename T> std::array<T,3> d_abs(const std::array<T,3>& a){ return {sycl::abs(a[0]), sycl::abs(a[1]), sycl::abs(a[2])};}
template <typename T> void d_swap(T& a, T& b){
  T c = a;
  a = b;
  b = c;
}
template <typename T> std::array<T,3> d_sort(const std::array<T,3>& a){
  std::array<T,3> b = a;
  if(sycl::abs(b[0]) > sycl::abs(b[1])) d_swap(b[0], b[1]);
  if(sycl::abs(b[1]) > sycl::abs(b[2])) d_swap(b[1], b[2]);
  if(sycl::abs(b[0]) > sycl::abs(b[1])) d_swap(b[0], b[1]);
  return b;
}
template <typename T> std::array<T,3> rsqrt3(const std::array<T,3>& a){
  return {sycl::rsqrt(a[0]), sycl::rsqrt(a[1]), sycl::rsqrt(a[2])};
}
template <typename T> std::array<T,3> cos3(const std::array<T,3>& a){
  return {sycl::cos(a[0]), sycl::cos(a[1]), sycl::cos(a[2])};
}

template <typename T, typename K> void d_set(std::array<T,3>& a, const u_char j, K b){
  a[j] = T(b);
}

template <typename T> T d_get(const std::array<T,3>& a, const u_char j){
  return a[j];
}

//5 FLOPs
template <typename T>  T  dot(const std::array<T,3>& a,  const std::array<T,3>& b) { return a[0]*b[0] + a[1]*b[1] + a[2]*b[2]; }

//6 FLOPs
template <typename T>  T norm(const std::array<T,3>& a)                    { return sycl::sqrt(dot(a,a)); }

template <typename T> T sum(const std::array<T,3>& a) {return a[0] + a[1] + a[2];}

template <typename T> T max(const std::array<T,3>& a) {return sycl::max(sycl::max(a[0],a[1]),a[2]);}

template <typename T> T d_min(const T a, const T b){
  return a < b ? a : b;
}

template <typename T> T min(const std::array<T,3>& a) {return sycl::min(sycl::min(a[0],a[1]),a[2]);}

//7 FLOPs
template <typename T>  std::array<T,3> unit_vector(const std::array<T,3>& a){
  T r = sycl::rsqrt(dot(a,a));
  return (a*r);
}
//10 FLOPs
template <typename T>  std::array<T,3> cross(const std::array<T,3>& a, const std::array<T,3>& b){ return {a[1]*b[2]-a[2]*b[1],
							   -a[0]*b[2]+a[2]*b[0],
							   a[0]*b[1]-a[1]*b[0]}; }
// $(a \otimes b) \cdot c$
template <typename T>  std::array<T,3> outer_dot(const std::array<T,3>& a, const std::array<T,3>& b, const std::array<T,3>& c){
  return {a[0]*(b[0]*c[0] + b[1]*c[1] + b[2]*c[2]), 
  a[1]*(b[0]*c[0] + b[1]*c[1] + b[2]*c[2]), 
  a[2]*(b[0]*c[0] + b[1]*c[1] + b[2]*c[2])};
}

//6 FLOPs
template <typename T>  T bond_length(const std::array<T,3>& ab){
    return sycl::rsqrt(dot(ab,ab));
}

template <typename T> T non_resciprocal_bond_length(const std::array<T,3>& ab){
    return sycl::sqrt(dot(ab,ab));
}

template <typename T> void print_coord(const std::array<T,3>& ab){

    printf("[%.16e, %.16e, %.16e]\n",ab[0],ab[1],ab[2]);
}

template <typename T1, typename T2> void assign(std::array<T1,3>& a, const std::array<T2,3>& b){
  a[0] = b[0];
  a[1] = b[1];
  a[2] = b[2];
}


template <typename T> void assign(std::array<float,3>& a, const std::array<T,3>& b){
  a[0] = b[0];
  a[1] = b[1];
  a[2] = b[2];
}

template <typename T> void assign(std::array<double,3>& a, const std::array<T,3>& b){
  a[0] = b[0];
  a[1] = b[1];
  a[2] = b[2];
}

template <typename T> T safe_rsqrt(T a){
  static_assert(std::is_floating_point<T>::value, "T must be floating point");
  return sycl::rsqrt(a + std::numeric_limits<T>::epsilon()*T(1e2));
}

#endif
