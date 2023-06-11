#include <limits>
#include <algorithm>
struct symMat3
{
  device_real_t a, b, c, d, e, f;
  
  //[[a , b,  c]
  // [b,  d,  e]
  // [c,  e,  f]]
  INLINE symMat3(): a(0), b(0), c(0), d(0), e(0), f(0){}
  INLINE symMat3(device_real_t a, device_real_t b, device_real_t c, device_real_t d, device_real_t e, device_real_t f) : a(a), b(b), c(c), d(d), e(e), f(f){}
  
  //Approx 107 FLOPS
  INLINE device_coord3d eigenvalues() const{
    DEVICE_TYPEDEFS
     // Coefficients of characteristic polynomial, calculated with Mathematica
    real_t 
      A = -1.f,
      B = a+d+f,
      C = b*b + c*c - a*d + e*e - a*f - d*f,
      D = -c*c*d + real_t(2.)*b*c*e - a*e*e - b*b*f + a*d*f;

    if(ABS(D) < real_t(1e3*std::numeric_limits<real_t>::epsilon())){
      auto temp = B*B - real_t(4.)*A*C;
      real_t Disc = temp > real_t(0.) ? SQRT(B*B - real_t(4.)*A*C) : real_t(0.);

      return {0., (-B-Disc)/( real_t(2.)*A),(-B+Disc)/( real_t(2.)*A)};
    }

    // Depress characteristic polynomial - see http://en.wikipedia.org/wiki/Cubic_equation#Reduction_to_a_depressed_cubic
    real_t
      p  = ( real_t(3.)*A*C - B*B)/( real_t(3.)*A*A),
      q  = ( real_t(2.)*B*B*B - real_t(9.)*A*B*C + real_t(27.)*A*A*D)/( real_t(27.)*A*A*A),
      xc = B/( real_t(3.)*A);

    // François Viète's solution to cubic polynomials with three real roots. 
    device_coord3d t;
    if(ABS(p) < real_t(1e3*std::numeric_limits<real_t>::epsilon())) {
      t = {real_t(0.), real_t(0.), real_t(0.)};
      return t - xc;}

    //For numerical stability we must ensure that acos doesn't receive an arugment which is outside [-1,1]
    auto frac = ( real_t(3.)*q)/( real_t(2.)*p)*SQRT(real_t(-3.)/p);
    frac = d_max(real_t(-1.),d_min(real_t(1.), frac));

    real_t K = real_t(2.)*SQRT(-p/ real_t(3.)), 
                  theta0 = (real_t(1.)/ real_t(3.))*real_t(ACOS(frac));
    for(int k=0;k<3;k++) d_set(t,k,K*COS(theta0- k* real_t(2.)* real_t(M_PI)/ real_t(3.)) );
    // lambda = t - B/(3A)
    return t - xc;
    
  }

  INLINE std::array<double,3> eigenvalues_fp64() const {
    DEVICE_TYPEDEFS
     // Coefficients of characteristic polynomial, calculated with Mathematica
    double 
      A = -1.f,
      B = a+d+f,
      C = b*b + c*c - a*d + e*e - a*f - d*f,
      D = -c*c*d + double(2.)*b*c*e - a*e*e - b*b*f + a*d*f;

    if(ABS(D) < double(1e3*std::numeric_limits<double>::epsilon())){
      double temp = B*B - double(4.)*A*C;
      double Disc = temp > double(0.) ? SQRT(B*B - double(4.)*A*C) : double(0.);

      return {0., (-B-Disc)/( double(2.)*A),(-B+Disc)/( double(2.)*A)};
    }

    // Depress characteristic polynomial - see http://en.wikipedia.org/wiki/Cubic_equation#Reduction_to_a_depressed_cubic
    double
      p  = ( double(3.)*A*C - B*B)/( double(3.)*A*A),
      q  = ( double(2.)*B*B*B - double(9.)*A*B*C + double(27.)*A*A*D)/( double(27.)*A*A*A),
      xc = B/( double(3.)*A);

    // François Viète's solution to cubic polynomials with three real roots. 
    std::array<double,3> t;
    if(ABS(p) < double(1e3*std::numeric_limits<double>::epsilon())) {
      t = {double(0.), double(0.), double(0.)};
      return {t[0]- xc, t[1] - xc, t[2] - xc};}

    //For numerical stability we must ensure that acos doesn't receive an arugment which is outside [-1,1]
    double frac = ( double(3.)*q)/( double(2.)*p)*sqrt(double(-3.)/p);
    frac = std::clamp(frac, double(-1), double(1.0));

    double K = double(2.)*sqrt(-p/ double(3.)), 
                  theta0 = (double(1.)/ double(3.))*double(acos(frac));
    for(int k=0;k<3;k++) t[k] = K*cos(theta0- k* double(2.)* double(M_PI)/ double(3.));
    // lambda = t - B/(3A)
    return {t[0]- xc, t[1] - xc, t[2] - xc};
  }
  //Best case 25 FLOPS
  INLINE device_coord3d eigenvector(const device_real_t lambda) const{
    DEVICE_TYPEDEFS
    // using the first two eqs
    // [ a_12 * a_23 - a_13 * (a_22 - r) ]
    // [ a_12 * a_13 - a_23 * (a_11 - r) ]
    // [ (a_11 - r) * (a_22 - r) - a_12^2 ]
    real_t mat_norm = SQRT(a*a + b*b + c*c + d*d + e*e + f*f);
    real_t normx;
    coord3d x = {b*e - c*(d-lambda),
                 b*c - e*(a-lambda),
                 (a-lambda)*(d-lambda) - b*b };
    normx = norm(x);
    if (normx / mat_norm > (real_t)1e4*std::numeric_limits<real_t>::epsilon()){ // not zero-ish
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
    if (normx / mat_norm > (real_t)1e4*std::numeric_limits<real_t>::epsilon()){ // not zero-ish
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
    if (normx / mat_norm > (real_t)1e4*std::numeric_limits<real_t>::epsilon()){ // not zero-ish
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
    real_t normx;
    coord3d x1 = {b*e - c*(d-lambda),
                 b*c - e*(a-lambda),
                 (a-lambda)*(d-lambda) - b*b };
    
    coord3d x2 = { b*(f-lambda) - c*e,
                 c*c - (a-lambda)*(f-lambda),
                 e*(a-lambda) - b*c };

    coord3d x3 ={ e*e - (d-lambda)*(f-lambda),
                 b*(f-lambda) - c*e,
                 c*(d-lambda) - b*e };

    real_t err1 = norm(mul(x1)/x1 - lambda);
    real_t err2 = norm(mul(x2)/x2 - lambda);
    real_t err3 = norm(mul(x3)/x3 - lambda);
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