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

    if(ABS(D) < (real_t)1e-12){
      auto temp = B*B - real_t(4.f)*A*C;
      real_t Disc = temp > (real_t)0. ? SQRT(B*B - real_t(4.f)*A*C) : (real_t)0;

      return {0.f, (-B-Disc)/( real_t(2.f)*A),(-B+Disc)/( real_t(2.f)*A)};
    }

    // Depress characteristic polynomial - see http://en.wikipedia.org/wiki/Cubic_equation#Reduction_to_a_depressed_cubic
    real_t
      p  = ( (real_t)3.f*A*C - B*B)/( (real_t)3.f*A*A),
      q  = ( (real_t)2.f*B*B*B - (real_t)9.f*A*B*C + (real_t)27.f*A*A*D)/( (real_t)27.f*A*A*A),
      xc = B/( (real_t)3.f*A);

    // François Viète's solution to cubic polynomials with three real roots. 
    device_coord3d t;
    if(ABS(p) < (real_t)1e-12) {
      t = {(real_t)0., (real_t)0., (real_t)0.};
      return t - xc;}

    //For numerical stability we must ensure that acos doesn't receive an arugment which is outside [-1,1]
    auto frac = ( (real_t)3.f*q)/( (real_t)2.f*p)*SQRT((real_t)-3.f/p);
    frac = d_max((real_t)-1.,d_min((real_t)1., frac));

    real_t K = (real_t)2.f*SQRT(-p/ (real_t)3.f), 
                  theta0 = ((real_t)1.f/ (real_t)3.f)*(real_t)acos((float)frac);
    for(int k=0;k<3;k++) d_set(t,k,K*COS(theta0-(real_t)k* (real_t)2.f* (real_t)M_PI/ (real_t)3.f) );
    // lambda = t - B/(3A)
    return t - xc;
    
  }
  //Best case 25 FLOPS
  INLINE device_coord3d eigenvector(const device_real_t lambda) const{
    DEVICE_TYPEDEFS
    // using the first two eqs
    // [ a_12 * a_23 - a_13 * (a_22 - r) ]
    // [ a_12 * a_13 - a_23 * (a_11 - r) ]
    // [ (a_11 - r) * (a_22 - r) - a_12^2 ]
    real_t normx;
    coord3d x = {b*e - c*(d-lambda),
                 b*c - e*(a-lambda),
                 (a-lambda)*(d-lambda) - b*b };
    normx = norm(x);
    if (normx / (a + d + f) > (real_t)1.e-6){ // not zero-ish
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
    if (normx / (a + d + f) > (real_t)1.e-6){ // not zero-ish
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
    if (normx / (a + d + f) > (real_t)1.e-6){ // not zero-ish
      return x/normx;
    } 
    //assert(false); // Something went wrong possibly two degenerate evals.
    return device_coord3d();
  }
};