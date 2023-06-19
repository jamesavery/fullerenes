#include <limits>
#include <algorithm>

// TODO: Hvor skal disse bo?
template <typename scalar> INLINE scalar __device__
relative_error(scalar x, scalar xref)
{
  return (x-xref)/(xref + (xref==0));
}

template <typename scalar> INLINE bool __device__
fp_close(scalar x, scalar y, scalar tolerance=100*std::numeric_limits<scalar>::epsilon())
{
  return std::abs(relative_error(x,y)) < tolerance;
}

template <typename scalar> INLINE void __device__
swap(scalar& a, scalar&b){ scalar c = a; a = b; b = c; }

template <typename scalar> INLINE array<scalar,3> __device__
cross(const array<scalar,3>& a, const array<scalar,3>& b){
  return {a[1]*b[2]-a[2]*b[1],
         -a[0]*b[2]+a[2]*b[0],
	  a[0]*b[1]-a[1]*b[0]}; }

template <typename scalar> INLINE scalar  __device__ 
dot(const array<scalar,3>& a, const array<scalar,3>& b){ return a[0]*b[0]+a[1]*b[1]+a[2]*b[2]; }

template <typename scalar> INLINE array<scalar,3> __device__
dot(const array<array<scalar,3>,3>& A, const array<scalar,3>& x){
  return
    {A[0][0]*x[0] + A[0][1]*x[1] + A[0][2]*x[2],
     A[1][0]*x[0] + A[1][1]*x[1] + A[1][2]*x[2],
     A[2][0]*x[0] + A[2][1]*x[1] + A[2][2]*x[2]};
}

template <typename scalar> INLINE array<scalar,3> __device__
operator-(const array<scalar,3>& a, const array<scalar,3>& b){ return {a[0]-b[0],a[1]-b[1],a[2]-b[2]}; }

template <typename scalar> INLINE array<scalar,3> __device__
operator+(const array<scalar,3>& a, const array<scalar,3>& b){ return {a[0]+b[0],a[1]+b[1],a[2]+b[2]}; }

template <typename scalar> INLINE array<scalar,3> __device__
operator*(const scalar& a, const array<scalar,3>& b){ return {a*b[0],a*b[1],a*b[2]}; }

template <typename scalar> INLINE array<scalar,3> __device__
operator*(const array<scalar,3>& b, const scalar& a){ return {a*b[0],a*b[1],a*b[2]}; }


template <typename scalar>
array<scalar,3> __device__ eig_sort(scalar l1, scalar l2, scalar l3)
{
  scalar a = std::min(l1,l2), b = std::max(l1,l2);
  if(l3<a) return {{l3,a,b}};
  if(l3<b) return {{a,l3,b}};
  else     return {{a,b,l3}};
}

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
    DEVICE_TYPEDEFS;

      // Combined home-derived eigenvalue w/ Wikipedia method
    long double 
      A = -1.L,
      B = a+d+f,
      C = b*b + c*c - a*d + e*e - a*f - d*f,
      D = -c*c*d + 2*b*c*e - a*e*e - b*b*f + a*d*f;

    if(fabs(D) < 100*epsilon){			
      long double Disc = SQRT(B*B-4*A*C); // TODO: Kahan's formula for b^2-4ac.
      device_real_t lam1 = -0.5L*(-B-Disc), lam2 = -0.5L*(-B+Disc);
      return eig_sort<real_t>(0,lam1,lam2);
    }
      
    device_real_t p1 = b*b + c*c + e*e;
    if(fabs(p1) < 100*epsilon)
      return eig_sort(a,d,f);

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

  device_real_t __device__
  xTAx(const device_coord3d &x) const {
    return x[0]*(a*x[0] + b*x[1] + c*x[2])
         + x[1]*(b*x[0] + d*x[1] + e*x[2])
         + x[2]*(c*x[0] + e*x[1] + f*x[2]);
  }
  
  
// Compute eigenvector associated with lambda, directly from 3x3 equations.
// NB: lambda is input/output and updated to Rayleigh quotient.  
device_coord3d __device__
eigenvector3x3(device_real_t &lambda) const{
    DEVICE_TYPEDEFS;
    
    array<real_t,3> xs[3];    
    int imax=0;
    real_t ns[3];
    
    real_t aa = a-lambda, dd = d-lambda, ff = f-lambda;
    // using the first two eqs yields determinant
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

    // Choose largest determinant (avoid using two near-linearly dependent equations of the three)
    imax = 0;      
    for(int i=0;i<3;i++) ns[i] = dot(xs[i],xs[i]);
    for(int i=0;i<3;i++) if(ns[i]>ns[imax]) imax = i; // ns[imax] = max(ns)      

    real_t n = 1.L/SQRT(ns[imax]);
    for(int i=0;i<3;i++) xs[imax][i] *= n;
			     
    lambda = xTAx(xs[imax]);
    
    return xs[imax];
  }

  // Eigenvector asssociated with lambda1 and orthogonal to v0. Calculated by deflating to 2x2 against v0.
  // NB: lambda1 is input/output and updated to with improved accuracy.
  device_coord3d __device__
  eigenvector2x2(const device_coord3d &v0, device_real_t &lambda1) const {
    DEVICE_TYPEDEFS;
    
    // Vælg vilkårlige u og v som er orthogonale til v0:
    real_t n = (v0[0]+v0[1]+v0[2]);

    device_coord3d u{1-n*v0[0],1-n*v0[1],1-n*v0[2]};  u = u*real_t(1/sqrt(dot(u,u)));
    device_coord3d v = cross(v0,u);

    // Symmetrisk 2x2 matrix i u,v-basis
    device_coord3d Au = dot(mat(),u), Av = dot(mat(),v); 

    // std::cout << "v0 = " << v0 << "\n";
    // std::cout << "dot(u,u)  = " << dot(u,u)  << ", dot(v,v)  = " << dot(v,v) << "\n";
    // std::cout << "dot(u,v0) = " << dot(u,v0) << ", dot(v,v0) = " << dot(v,v0) << ", dot(u,v) = " << dot(u,v) << "\n";
    // std::cout << "dot(Au,v0) = " << dot(Au,v0) << ", dot(Av,v0) = " << dot(Av,v0) <<"\n";    

    real_t aa = dot(u,Au)-lambda1, bb = dot(u,Av), cc = dot(v,Av)-lambda1;

    // Characteristic polynmoial is chi(l) = det(A-l*I) = (a-l)*(c-l) - b^2 = l^2 - (a+c)l +ac-b^2
    double A = 1, B = -(aa+cc), C = aa*cc-bb*bb;
    double Disc = sqrt(B*B-4*A*C); // TODO: Kahan's formula for b^2-4ac instead of extra precision
    double lam1 = 0.5L*(-B-Disc), lam2 = 0.5L*(-B+Disc);

    real_t lam = std::abs(lam1) < std::abs(lam2)? lam1 : lam2;

    // std::cout << "a,b,c = "<<device_coord3d{{aa,bb,cc}} << "; {A,B,C} = " << array<long double,3>{{A,B,C}} << "\n";
    // std::cout << "{lambda1,lam1,lam1} = "<< array<long double,3>{{lambda1,lam1,lam2}} << "\n";
    // //    real_t x[2] = {bb,lam-aa};
    real_t x[2] = {lam-cc,bb};
    // real_t y[2] = {aa*x[0]+bb*x[1], bb*x[0]+cc*x[1]};
    // real_t lam2x2 = (x[0]*y[0] + x[1]*y[2])/(x[0]*x[0]+x[1]*x[1]);
    // std::cout << "lam2x2 = " << lam2x2 << "\n";
    // std::cout << "r2x2 =   " << (y[0] - lam2x2*x[0]) << ", " << (y[1]-lam2x2*x[1]) << "\n";
    // std::cout << "q2x2 =   " << (y[0]/x[0]) << ", " << (y[1]/x[1]) << "\n";
    
    device_coord3d v1 = x[0]*u + x[1]*v;
    v1 = v1*real_t(1/sqrt(dot(v1,v1)));

    lambda1 += lam;		// Update lambda1 with correction
    return v1;
  }


  pair<array<device_coord3d,3>,device_coord3d> __device__
  eigensystem() const {  return eigensystem(eigenvalues());  }
  
  pair<array<device_coord3d,3>,device_coord3d> __device__
  eigensystem(const device_coord3d &lambdas_) const {
    DEVICE_TYPEDEFS;
    
    array<coord3d,3> v;
    // Compute eigenvalues using closed-form expressions
    coord3d lambdas = lambdas_;

    // If all the eigenvalues are close withing numerical accuracy, eigenvectors form identity matrix
    bool close01 = fp_close(lambdas[0],lambdas[1]), close12 = fp_close(lambdas[1],lambdas[2]);
    if(close01 && close12) return {{{{1,0,0},{0,1,0},{0,0,1}}},lambdas};

    // If there are at least two distinct eigenvalues, we start by directly computing the eigenvector
    // for the most isolated eigenvalue (eigenvector3x3). 
    // Then we compute 
    bool largest_gap = fabs(lambdas[1]-lambdas[0]) < fabs(lambdas[2]-lambdas[1]);
    if(largest_gap==0){			       // Smallest eigenvalue is most isolated.
      v[0] = eigenvector3x3(lambdas[0]);       // 1. Direct 3x3 eigenvector computation. lambdas[0] is updated with Rayleigh quotient.
      v[0] = eigenvector3x3(lambdas[0]);       //    Iterate once with Rayleigh coefficient to increase accuracy.
      v[1] = eigenvector2x2(v[0],lambdas[1]);  // 2. Deflate to 2x2 and solve to get second eigenvector (robust against degeneracy).

      v[2] = cross(v[0],v[1]);	               // 3. Final eigenvector is simply cross product, as it must be orthogonal to first two.
      lambdas[2] = xTAx(v[2]);		       //    Compute Rayleigh coefficient to increase eigenvalue accuracy.
    } else {
      v[2] = eigenvector3x3(lambdas[2]);       // Same, but largest eigenvalue is most isolated.
      v[2] = eigenvector3x3(lambdas[2]);    
      v[1] = eigenvector2x2(v[2],lambdas[1]);

      v[0] = cross(v[1],v[2]);
      lambdas[0] = xTAx(v[0]);
    }

    return {v,lambdas};
  }
  
  array<device_coord3d,3> __device__
  mat() const { return {{{{a,b,c}},{{b,d,e}},{{c,e,f}}}}; }  
};