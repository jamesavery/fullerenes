#include <limits>
#include <algorithm>
struct symMat3
{
  static constexpr device_real_t epsilon = std::numeric_limits<device_real_t>::epsilon();
  
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
      double Disc = SQRT(B*B-4*A*C); // TODO: Kahan's formula for b^2-4ac.
      return {{0,real_t(-0.5L*(-B-Disc)), real_t(-0.5L*(-B+Disc))}};
    }
      
    device_real_t p1 = b*b + c*c + e*e;
    if(fabs(p1) < 100*epsilon){
      return {{a,d,f}};
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

    if(fabs(D) < 100*epsilon){			
      double Disc = SQRT(B*B-4*A*C); // TODO: Kahan's formula for b^2-4ac.
      return {{0,double(-0.5L*(-B-Disc)), double(-0.5L*(-B+Disc))}};
    }
      
    double p1 = b*b + c*c + e*e;
    if(fabs(p1) < 100*epsilon){
      return {{a,d,f}};
     }
      
    double p1 = b*b + c*c + e*e;
    if(fabs(p1) < 100*epsilon){
      return {{a,d,f}};
     }

    double q = (a+d+f)/3L; // q=Tr(M)/3
    double p2 = (a-q)*(a-q) + (d-q)*(d-q) + (f-q)*(f-q) + 2*p1;
    double p = SQRT(p2/6);

    double detBxp3 = -c*c*(d-q) + 2*b*c*e - (a-q)*e*e - b*b*(f-q) + (a-q)*(d-q)*(f-q);
    double r = detBxp3/(2*p*p*p);

    double phi = r <= -1? M_PI/3L : (r >= 1? 0L : ACOS(r)/3L);
    double lam1 = q + 2*p*COS(phi);
    double lam3 = q + 2*p*COS(phi +  (2L/3L)*M_PI);
    double lam2 = 3*q - lam1 - lam3;

    return {{lam1,lam2,lam3}};      
  }
  //Best case 25 FLOPS
  INLINE device_coord3d eigenvector(const device_real_t lambda) const{
    DEVICE_TYPEDEFS
    // using the first two eqs
    // [ a_12 * a_23 - a_13 * (a_22 - r) ]
    // [ a_12 * a_13 - a_23 * (a_11 - r) ]
    // [ (a_11 - r) * (a_22 - r) - a_12^2 ]
    device_real_t mat_norm = SQRT(a*a + b*b + c*c + d*d + e*e + f*f);
    device_real_t normx;
    coord3d x = {b*e - c*(d-lambda),
                 b*c - e*(a-lambda),
                 (a-lambda)*(d-lambda) - b*b };
    normx = norm(x);
    if (normx / mat_norm > (device_real_t)1e4*std::numeric_limits<device_real_t>::epsilon()){ // not zero-ish
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
    if (normx / mat_norm > (device_real_t)1e4*std::numeric_limits<device_real_t>::epsilon()){ // not zero-ish
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
    if (normx / mat_norm > (device_real_t)1e4*std::numeric_limits<device_real_t>::epsilon()){ // not zero-ish
      return x/normx;
    } 
    //assert(false); // Something went wrong possibly two degenerate evals.
    return device_coord3d();
  }

  INLINE device_coord3d mul(const device_coord3d x) const{
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
    device_real_t normx;
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
