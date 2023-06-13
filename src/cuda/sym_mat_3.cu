#include <limits>
#include <algorithm>

// TODO: Hvor skal disse bo?
INLINE device_real_t __device__ relative_error(device_real_t x, device_real_t xref)
{
  return (x-xref)/(xref + (xref==0));
}

INLINE bool __device__ fp_close(device_real_t x, device_real_t y, device_real_t tolerance=100*std::numeric_limits<device_real_t>::epsilon())
{
  return std::abs(relative_error(x,y)) < tolerance;
}

template <typename t>
std::array<t,3> __device__ eig_sort(t l1, t l2, t l3)
{
  t a = std::min(l1,l2), b = std::max(l1,l2);
  if(l3<a) return {{l3,a,b}};
  if(l3<b) return {{a,l3,b}};
  else     return {{a,b,l3}};
}

struct symMat3
{
  static constexpr device_real_t epsilon = std::numeric_limits<device_real_t>::epsilon();
  static constexpr double      epsilon64 = std::numeric_limits<double>::epsilon();  
  
  device_real_t a, b, c, d, e, f;
  
  //[[a , b,  c]
  // [b,  d,  e]
  // [c,  e,  f]]
  INLINE symMat3(): a(0), b(0), c(0), d(0), e(0), f(0){}
  INLINE symMat3(device_real_t a, device_real_t b, device_real_t c, device_real_t d, device_real_t e, device_real_t f) : a(a), b(b), c(c), d(d), e(e), f(f){}
  
  INLINE device_coord3d eigenvalues() const{
    DEVICE_TYPEDEFS

      // Combined home-derived eigenvalue w/ Wikipedia method
    long double 
      A = -1.L,
      B = a+d+f,
      C = b*b + c*c - a*d + e*e - a*f - d*f,
      D = -c*c*d + 2*b*c*e - a*e*e - b*b*f + a*d*f;

    if(fabs(D) < 100*epsilon){			
      long double Disc = SQRT(B*B-4*A*C); // TODO: Kahan's formula for b^2-4ac.
      device_real_t lam1 = -0.5L*(-B-Disc), lam2 = -0.5L*(-B+Disc);
      return eig_sort<device_real_t>(0,lam1,lam2);
    }
      
    device_real_t p1 = b*b + c*c + e*e;
    if(fabs(p1) < 100*epsilon){
      return eig_sort<device_real_t>(a,d,f);
     }

    device_real_t q = (a+d+f)/3.L; // q=Tr(M)/3
    device_real_t p2 = (a-q)*(a-q) + (d-q)*(d-q) + (f-q)*(f-q) + 2.L*p1;
    device_real_t p = SQRT(p2/6.L);

    device_real_t detBxp3 = -c*c*(d-q) + 2*b*c*e - (a-q)*e*e - b*b*(f-q) + (a-q)*(d-q)*(f-q);
    device_real_t r = detBxp3/(2*p*p*p);

    device_real_t phi = r <= -1? M_PI/3.L : (r >= 1? 0.L : ACOS(r)/3.L);
    device_real_t lam1 = q + 2*p*COS(phi);
    device_real_t lam3 = q + 2*p*COS(phi +  (2.L/3.L)*M_PI);
    device_real_t lam2 = 3*q - lam1 - lam3;

    return {{lam1,lam2,lam3}};
  }

  INLINE std::array<double,3> eigenvalues_fp64() const {
    DEVICE_TYPEDEFS

      // Combined home-derived eigenvalue w/ Wikipedia method
    long double 
      A = -1.L,
      B = a+d+f,
      C = b*b + c*c - a*d + e*e - a*f - d*f,
      D = -c*c*d + 2*b*c*e - a*e*e - b*b*f + a*d*f;

    if(fabs(D) < 100*epsilon64){			
      long double Disc = sqrt(B*B-4*A*C); // TODO: Kahan's formula for b^2-4ac.
      double lam1 = -0.5L*(-B-Disc), lam2 = -0.5L*(-B+Disc);
      return eig_sort<double>(0,lam1,lam2);
    }
      
    double p1 = b*b + c*c + e*e;
    if(fabs(p1) < 100*epsilon64){
      return eig_sort<double>(a,d,f);
     }

    double q = (a+d+f)/3.L; // q=Tr(M)/3
    double p2 = (a-q)*(a-q) + (d-q)*(d-q) + (f-q)*(f-q) + 2.L*p1;
    double p = sqrt(p2/6.L);

    double detBxp3 = -c*c*(d-q) + 2*b*c*e - (a-q)*e*e - b*b*(f-q) + (a-q)*(d-q)*(f-q);
    double r = detBxp3/(2*p*p*p);

    double phi = r <= -1? M_PI/3.L : (r >= 1? 0.L : acos(r)/3.L);
    double lam1 = q + 2*p*cos(phi);
    double lam3 = q + 2*p*cos(phi +  (2.L/3.L)*M_PI);
    double lam2 = 3*q - lam1 - lam3;

    return {{lam1,lam2,lam3}};      
  }

  template <typename t>
  INLINE device_real_t xTAx(const std::array<t,3> &x) const {
    return x[0]*(a*x[0] + b*x[1] + c*x[2])
         + x[1]*(b*x[0] + d*x[1] + e*x[2])
         + x[2]*(c*x[0] + e*x[1] + f*x[2]);
  }
  
  //Best case 25 FLOPS
  INLINE device_coord3d eigenvector(const device_real_t lambda) const{
    DEVICE_TYPEDEFS
    device_real_t  aa = a-lambda, dd = d-lambda, ff = f-lambda;
    device_coord3d xs[3];
    // using the first two eqs
    // [ b * e - c * (d - r) ]
    // [ b * c - e * (a - r) ]
    // [ (a - r)*(d - r)-b^2 ]
    xs[0] = {b*e - c*dd, b*c - e*aa, aa*dd - b*b };
    // using the first+last eqs
    // [ b * (f - r) - c * e ]
    // [ c^2-(a - r) * (f - r) ]
    // [ e * (a - r) - b * c ]
    xs[1] = { b*ff - c*e,  c*c - aa*ff, e*aa - b*c };

    // using the last two eqs
    // [ e^2-(d - r) * (f - r) ]
    // [ b * (f - r) - c * e ]
    // [ c * (d - r) - b * e ]
    xs[2] ={ e*e - dd*ff, b*ff - c*e, c*dd - b*e };

    device_real_t ns[3] = {dot(xs[0],xs[0]), dot(xs[1],xs[1]), dot(xs[2],xs[2])};
    int imax=0;			
    for(int i=0;i<3;i++) if(ns[i]>ns[imax]) imax = i; // ns[imax] = max(ns)

    return xs[imax]/sqrt(ns[imax]);
    
    // {
    //   device_real_t ns[3] = {dot(xs[0],xs[0]), dot(xs[1],xs[1]), dot(xs[2],xs[2])};
    //   int imax=0;			
    //   for(int i=0;i<3;i++) if(ns[i]>ns[imax]) imax = i; // ns[imax] = max(ns)
      
    //   device_real_t lambda1 = xTAx(xs[imax])/ns[imax];
    //   aa = a-lambda1, dd = d-lambda1, ff = f-lambda1;
    //   xs[0] = { b*e - c*dd, b*c - e*aa, aa*dd - b*b };
    //   xs[1] = { b*ff - c*e,  c*c - aa*ff, e*aa - b*c };
    //   xs[2] = { e*e - dd*ff, b*ff - c*e, c*dd - b*e };
    // }

    // {
    //   device_real_t ns[3] = {dot(xs[0],xs[0]), dot(xs[1],xs[1]), dot(xs[2],xs[2])};
    //   int imax=0;			
    //   for(int i=0;i<3;i++) if(ns[i]>ns[imax]) imax = i; // ns[imax] = max(ns)
      
    //   return xs[imax]/SQRT(ns[imax]);
    // }
    
    // x1 /= n1; x2 /= n2; x3 /= n3;
    // device_real_t l1 = xTAx(x1), l2 = xTAx(x2), l3 = xTAx(x3);
    // device_coord3d
    //    r1 = mul(x1) - l1*x1,
    //    r2 = mul(x2) - l2*x2,
    //    r3 = mul(x3) - l3*x3;

    // device_real_t N1 = norm(r1), N2 = norm(r2), N3 = norm(r3);
    
    // if(N1<N2) return (N1<N3? x1 : x3);
    //    else      return (N2<N3? x2 : x3);
  }

  INLINE std::array<device_coord3d,3> eigenvectors(const device_coord3d lambdas) const {
    std::array<device_coord3d,3> v;
    bool close01 = fp_close(lambdas[0],lambdas[1]), close12 = fp_close(lambdas[1],lambdas[2]);
    if(close01 && close12) return {{{1,0,0},{0,1,0},{0,0,1}}};

    bool largest_gap = fabs(lambdas[1]-lambdas[0]) < fabs(lambdas[2]-lambdas[1]);
    if(largest_gap==0){
      v[0] = eigenvector(lambdas[0]);
      v[1] = eigenvector(lambdas[1]);
      v[2] = cross(v[0],v[1]);
    } else {
      v[1] = eigenvector(lambdas[1]);      
      v[2] = eigenvector(lambdas[2]);
      v[0] = cross(v[1],v[2]);      
    }

    return v;
  }  
  
  INLINE device_coord3d mul(const device_coord3d &x) const{
    return {a*x[0] + b*x[1] + c*x[2],
            b*x[0] + d*x[1] + e*x[2],
            c*x[0] + e*x[1] + f*x[2]};
  }

  INLINE device_coord3d eigenvector_brute_force(const device_real_t lambda) const{
    DEVICE_TYPEDEFS
    // using the first two eqs
    // [ a_12 * a_23 - a_13 * (a_22 - r) ]
    // [ a_12 * a_13 - a_23 * (a_11 - r) ]
    // [ (a_11 - r) * (a_22 - r) - a_12^2 ]
    coord3d x1 = {b*e - c*(d-lambda),
                 b*c - e*(a-lambda),
                 (a-lambda)*(d-lambda) - b*b };
    
    coord3d x2 = { b*(f-lambda) - c*e,
                 c*c - (a-lambda)*(f-lambda),
                 e*(a-lambda) - b*c };

    coord3d x3 ={ e*e - (d-lambda)*(f-lambda),
                 b*(f-lambda) - c*e,
                 c*(d-lambda) - b*e };

    device_real_t err1 = norm(mul(x1)/x1 - lambda);
    device_real_t err2 = norm(mul(x2)/x2 - lambda);
    device_real_t err3 = norm(mul(x3)/x3 - lambda);
    if( threadIdx.x + blockIdx.x == 0){
      coord3d lambdax1 = mul(x1)/x1; 
      coord3d lambdax2 = mul(x2)/x2;
      coord3d lambdax3 = mul(x3)/x3;
      printf("Errors: %f, %f, %f\n", lambdax1[0], lambdax1[1], lambdax1[2]);
      printf("Errors: %f, %f, %f\n", lambdax2[0], lambdax2[1], lambdax2[2]);
      printf("Errors: %f, %f, %f\n", lambdax3[0], lambdax3[1], lambdax3[2]);

    }
    if (err1 <= err2 && err1 <= err3) return x1/norm(x1);
    if (err2 <= err1 && err2 <= err3) return x2/norm(x2);
    if (err3 <= err1 && err3 <= err2) return x3/norm(x3);
    //assert(false); // Something went wrong possibly two degenerate evals.
    return device_coord3d();
  }
};
