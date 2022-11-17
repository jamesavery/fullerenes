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

struct symMat3
{
  device_real_t a, b, c, d, e, f;
  
  //[[a , b,  c]
  // [b,  d,  e]
  // [c,  e,  f]]
  INLINE symMat3(){}
  INLINE symMat3(device_real_t a, device_real_t b, device_real_t c, device_real_t d, device_real_t e, device_real_t f) : a(a), b(b), c(c), d(d), e(e), f(f){}
  
  //Approx 107 FLOPS
  INLINE device_coord3d eigenvalues() const{
    DEVICE_TYPEDEFS
     // Coefficients of characteristic polynomial, calculated with Mathematica
    real_t 
      A = -1.f,
      B = a+d+f,
      C = b*b + c*c - a*d + e*e - a*f - d*f,
      D = -c*c*d + (real_t)2.f*b*c*e - a*e*e - b*b*f + a*d*f;

    if(abs(D) < 1e-12){
      auto temp = B*B - real_t(4.f)*A*C;
      real_t Disc = temp > (real_t)0. ? sqrt(B*B - real_t(4.f)*A*C) : 0;

      return {0.f, (-B-Disc)/( real_t(2.f)*A),(-B+Disc)/( real_t(2.f)*A)};
    }

    // Depress characteristic polynomial - see http://en.wikipedia.org/wiki/Cubic_equation#Reduction_to_a_depressed_cubic
    real_t
      p  = ( (real_t)3.f*A*C - B*B)/( (real_t)3.f*A*A),
      q  = ( (real_t)2.f*B*B*B - (real_t)9.f*A*B*C + (real_t)27.f*A*A*D)/( (real_t)27.f*A*A*A),
      xc = B/( (real_t)3.f*A);

    // François Viète's solution to cubic polynomials with three real roots. 
    device_coord3d t;
    if(abs(p) < 1e-12) {
      t = {(real_t)0., (real_t)0., (real_t)0.};
      return t - xc;}

    //For numerical stability we must ensure that acos doesn't receive an arugment which is outside [-1,1]
    auto frac = ( (real_t)3.f*q)/( (real_t)2.f*p)*sqrt((real_t)-3.f/p);
    frac = d_max((real_t)-1.,d_min((real_t)1., frac));

    real_t K = (real_t)2.f*sqrt(-p/ (real_t)3.f), 
                  theta0 = ((real_t)1.f/ (real_t)3.f)*acos(frac);
    for(int k=0;k<3;k++) d_set(t,k,K*cos(theta0-k* (real_t)2.f* (real_t)M_PI/ (real_t)3.f) );
    // lambda = t - B/(3A)
    return t - xc;
    
  }
  //Best case 25 FLOPS
  INLINE device_coord3d eigenvector(const device_real_t lambda) const{
    // using the first two eqs
    // [ a_12 * a_23 - a_13 * (a_22 - r) ]
    // [ a_12 * a_13 - a_23 * (a_11 - r) ]
    // [ (a_11 - r) * (a_22 - r) - a_12^2 ]
    device_real_t normx;
    device_coord3d x = {b*e - c*(d-lambda),
                 b*c - e*(a-lambda),
                 (a-lambda)*(d-lambda) - b*b };
    normx = norm(x);
    if (normx / (a + d + f) > 1.e-12){ // not zero-ish
      return x/normx;
    }
  
    // using the first+last eqs
    // [ a_12 * (a_33 - r) - a_13 * a_23 ]
    // [ a_13^2 - (a_11 - r) * (a_33 - r) ]
    // [ a_23 * (a_11 - r) - a_12 * a_13 ]
    x = { b*(f-lambda) - c*e,
                 c*c - (a-lambda)*(f-lambda),
                 e*(a-lambda) - b*c };
    normx = norm(x);
    if (normx / (a + d + f) > 1.e-12){ // not zero-ish
      return x/normx;
    }

    // using the last two eqs
    // [ a_23^2 - (a_22 - r) * (a_33 - r) ]
    // [ a_12 * (a_33 - r) - a_13 * a_23 ]
    // [ a_13 * (a_22 - r) - a_12 * a_23 ]
    x ={ e*e - (d-lambda)*(f-lambda),
                 b*(f-lambda) - c*e,
                 c*(d-lambda) - b*e };
    normx = norm(x);
    if (normx / (a + d + f) > 1.e-12){ // not zero-ish
      return x/normx;
    } 
    //assert(false); // Something went wrong possibly two degenerate evals.
    return device_coord3d();
  }
};

#endif
