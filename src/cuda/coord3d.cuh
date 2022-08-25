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
  return ((const device_node_t*)&a)[j];
}

struct symMat3
{
  device_real_t a, b, c, d, e, f;
  
  //[[a , b,  c]
  // [b,  d,  e]
  // [c,  e,  f]]

  INLINE symMat3(device_real_t a, device_real_t b, device_real_t c, device_real_t d, device_real_t e, device_real_t f) : a(a), b(b), c(c), d(d), e(e), f(f){}
  INLINE device_coord3d eigenvalues() const{
     // Coefficients of characteristic polynomial, calculated with Mathematica
    device_real_t 
      A = -1.f,
      B = a+d+f,
      C = b*b + c*c - a*d + e*e - a*f - d*f,
      D = -c*c*d + (device_real_t)2.f*b*c*e - a*e*e - b*b*f + a*d*f;

    if(abs(D) < 1e-12){
      device_real_t Disc = sqrt(B*B - device_real_t(4.f)*A*C);
      return {0.f, (-B-Disc)/( device_real_t(2.f)*A),(-B+Disc)/( device_real_t(2.f)*A)};
    }

    // Depress characteristic polynomial - see http://en.wikipedia.org/wiki/Cubic_equation#Reduction_to_a_depressed_cubic
    device_real_t
      p  = ( (device_real_t)3.f*A*C - B*B)/( (device_real_t)3.f*A*A),
      q  = ( (device_real_t)2.f*B*B*B - (device_real_t)9.f*A*B*C + (device_real_t)27.f*A*A*D)/( (device_real_t)27.f*A*A*A),
      xc = B/( (device_real_t)3.f*A);

    // François Viète's solution to cubic polynomials with three real roots. 
    device_coord3d t;
    device_real_t K = (device_real_t)2.f*sqrt(-p/ (device_real_t)3.f), 
                  theta0 = ((device_real_t)1.f/ (device_real_t)3.f)*acos(( (device_real_t)3.f*q)/( (device_real_t)2.f*p)*sqrt((device_real_t)-3.f/p));
    for(int k=0;k<3;k++) d_set(t,k,K*cos(theta0-k* (device_real_t)2.f* (device_real_t)M_PI/ (device_real_t)3.f) );

    // lambda = t - B/(3A)
    return t - xc;
    
  }

  INLINE device_coord3d eigenvector(const device_real_t lambda){
    // using the first two eqs
    // [ a_12 * a_23 - a_13 * (a_22 - r) ]
    // [ a_12 * a_13 - a_23 * (a_11 - r) ]
    // [ (a_11 - r) * (a_22 - r) - a_12^2 ]
    device_coord3d x = {b*e - c*(d-lambda),
                 b*c - e*(a-lambda),
                 (a-lambda)*(d-lambda) - b*b };
    if (norm(x) / (a + d + f) > 1.e-12){ // not zero-ish
      return x/norm(x);
    }
  
    // using the first+last eqs
    // [ a_12 * (a_33 - r) - a_13 * a_23 ]
    // [ a_13^2 - (a_11 - r) * (a_33 - r) ]
    // [ a_23 * (a_11 - r) - a_12 * a_13 ]
    x = { b*(f-lambda) - c*e,
                 c*c - (a-lambda)*(f-lambda),
                 e*(a-lambda) - b*c };
    if (norm(x) / (a + d + f) > 1.e-12){ // not zero-ish
      return x/norm(x);
    }

    // using the last two eqs
    // [ a_23^2 - (a_22 - r) * (a_33 - r) ]
    // [ a_12 * (a_33 - r) - a_13 * a_23 ]
    // [ a_13 * (a_22 - r) - a_12 * a_23 ]
    x ={ e*e - (d-lambda)*(f-lambda),
                 b*(f-lambda) - c*e,
                 c*(d-lambda) - b*e };
    if (norm(x) / (a + d + f) > 1.e-12){ // not zero-ish
      return x/norm(x);
    } 
    
    assert(false); // Something went wrong possibly two degenerate evals.
    return device_coord3d();
  }
};

#endif
