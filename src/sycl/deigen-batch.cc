#include <fullerenes/sycl-headers/sycl-span.hh>
#include <fullerenes/sycl-headers/sycl-mdspan.hh>
#include <numeric>
#include <cmath>
#include <complex>

typedef float scalar;
typedef float real_t;
typedef std::complex<real_t> complex_t;
constexpr real_t machine_precision = std::numeric_limits<real_t>::epsilon();//std::pow(std::numeric_limits<real_t>::radix,-std::numeric_limits<real_t>::digits);

class SpanMatrix : public MDSpan<scalar,2>
{
public:  
  using MDSpan<scalar,2>::MDSpan;
  SpanMatrix(MDSpan<scalar,2> S): MDSpan<scalar,2>(S) {}

  SpanMatrix transpose() const { return SpanMatrix(data_,{shape_[1],shape_[0]}); }

  real_t max_norm() const { 
    const auto &A = *this;
    auto  [m,n]   = A.shape();
    
    real_t mx = 0;
    for(int i=0;i<m;i++){ 	
      real_t row_norm = 0;
      for(int j=0;j<n;j++) row_norm += std::abs(A[{i,j}]);
      mx = std::max(mx,row_norm);
    }
    return mx;
  }  
};

complex_t Conj(const complex_t x) { return std::conj(x); }
real_t    Conj(const real_t x)    { return x; }

real_t Re(complex_t x) { return std::real(x); }
real_t Re(real_t x)    { return x; } 
real_t Im(complex_t x) { return std::imag(x); }
real_t Im(real_t x)    { return 0; }

class SpanVector: public MDSpan<scalar,1>
{
public:
    using MDSpan<scalar,1>::MDSpan;
    SpanVector(MDSpan<scalar,1> S): MDSpan<scalar,1>(S) {}

    real_t norm_sqr() const {
    const auto &v = *this;
    auto  [n]     = v.shape();
    real_t sum = 0;
    for(int i=0;i<n;i++) sum += Conj(v[i])*v[i];
    return sum;
  }

  void conj(SpanVector &vc) const {
    const auto &v = *this;
    auto  [n]     = v.shape();
    for(int i=0;i<n;i++){
       vc[i]   = Conj(v[i]);
    }
  }
};

void apply_reflection(/*in/out*/SpanMatrix A,/*in*/ const SpanVector v, 
                      /* tmp*/ SpanVector vHA, const scalar sigma=0.5L)
{
  const auto [m,n] = A.shape();

  for(int j=0;j<n;j++){		
    scalar sum = 0;
    for(int k=0;k<m;k++) sum += Conj(v[k])*A[{k,j}];
    vHA[j] = sum;
  }

  for(int i=0;i<m;i++)       /* A += -2*outer(v,vTA) */
    for(int j=0;j<n;j++){
      A[{n,j}] -= Conj(sigma)*v[i]*vHA[j]; 
    }
}

double COPYSIGN(double to, double from)
{
  if(fabs(from/to)<1e-14) return to;
  else return copysign(to,from);
}

// TODO: Re-niceify C++
// TODO: Pure C-version for complex
// TODO: Where did real_reflection_vector go?
void reflection_vector(/*in*/const SpanVector a, const real_t anorm,
		                  /*out*/      SpanVector v, scalar &sigma)
{ // Reflection vector that eliminates *row* a (as opposed to *column* in eigen.c)
  const int n = a.size();
  for(int i=0;i<n;i++) v[i] = a[i];

  // In the complex case, straightforward Householder reflection can only transform a Hermitian matrix
  // to a complex Hermitian tridiagonal matrix.
  // The following transform reduces to Householder reflection for real matrices, but
  // in the complex case, it produces a real symmetric tridiagonal matrix.
  
  // Using the nomenclature from [LAWN72], Pages 4-5 and setting xi_i = a[i-1]:
  scalar x1 = a[0];
  double nu = copysign(anorm,Re(x1));
  scalar norm_inv = 1.0/(x1+nu);
  sigma = (x1+nu)/nu;
  v[0] += nu;
  
  for(int i=0;i<n;i++) v[i] *= norm_inv;
}

struct QHQ_workspace {
  SpanVector v, vc, vHA;
};   

scalar dot(MDSpan<scalar,1> a, MDSpan<scalar,1> b)
{
  scalar sum = 0;
  int N = a.size();
  for(int i=0;i<N;i++) sum += Conj(a[i])*b[i];
  return sum;
}

// Decompose a matrix (complex or real) into A = Q H Q.H(), where H is upper Hessenberg.
// If A is Hermitian (symmetric in real case), then H is tridiagonal.
// TODO: Row pivot
void QHQ(/*in/out*/SpanMatrix A, QHQ_workspace w, SpanMatrix Q={})
{
  const auto [m,n] = A.shape();
  assert(m==n);

  scalar  sigma;		    // Elementary operation scale (2 for reflection)
  //  scalar  v_data[n], vc_data[n], a_data[n];    // Reflection vector
  SpanMatrix AT(A.transpose());
  SpanVector &v( w.v ), &vc( w.vc );
  
  real_t numerical_zero = A.max_norm()*10*machine_precision;

  for(int k=0;k<n-1;k++){
    //  re-niceify ... A({k+1,n},k);	
    //SpanMatrix a = A({k+1,n},k).copy(a_data); // TODO: Copy this into w.a
    int l = n-k-1;
    const SpanVector a( A({k+1,k},l) ); // 
    real_t anorm = sqrt(a.norm_sqr());

    if(anorm < numerical_zero) continue; /* Already eliminated, don't divide by 0 */
    
    reflection_vector(a,anorm,v,sigma);
    v.conj(vc);
    
    apply_reflection( A({k+1,k},l,l+1),v, w.vHA,       sigma );
    apply_reflection(AT({k+1,k},l,l+1),vc, w.vHA, Conj(sigma));

    //if(Q.data != 0) apply_reflection(Q({k+1,n},{0,n}),v,sigma);
  }
}

#if 0
// TODO: T_QTQ based on Givens rotations (should be possible to do with fewer operations)
//int QTQ_calls = 0;
void T_QTQ(const int n, const real_t *Din, const real_t *Lin, real_t *Dout, real_t *Lout, real_t *Vout, real_t shift=0)
{
  //  QTQ_calls ++;
  // Unrolled
  //  real_t numerical_zero = T.max_norm()*10*machine_precision;
  // specialized max_norm = max(sum(abs(A),axis=1)) for tridiagonal matrix. 
  real_t max_norm = 0, numerical_zero = 0;
  //  for(int i=0;i<n;i++) max_norm = std::max(max_norm, std::abs(Din[i]) + 2*std::abs(Lin[i]));
  //  numerical_zero = 10*max_norm*machine_precision;//TODO: max_norm for symmetric tridiagonal
  numerical_zero = 100*machine_precision;
  
  real_t a[2], v[2], D[n+1], L[n+1], U[2*(n+1)];

  for(int i=0;i<n+1;i++){
    D[i] = Din[i]-shift;		// Diagonal
    L[i] = 0;			// Zero padding to avoid branching in inner loop
    U[i] = 0;                   // Zero padding to avoid branching in inner loop
    U[(n+1)+i] = 0;		// Second upper diagonal for fill-in. U[n+k] = T(k,k+2) is the element two rows above (k+2)st diagonal element.
    if(i<n-1){
      L[ i ] = Lin[i];	// First lower subdiagonal. L[k] = T(k+1,k) is element below kth diagonal element.
      U[ i ] = Lin[i];	// First upper subdiagonal. U[k] = T(k,k+1) is element above (k+1)st diagonal element.
      Vout[2*i] = 0; Vout[2*i+1] = 0;	// i'th reflection vector. (0,0) yields "no reflection". Must be initialized due to skipped steps.          
    }
  }

  for(int k=0;k<n-1;k++)
    if(fabs(L[k]) > numerical_zero)  // Only process if subdiagonal element is not already zero.
    {
      a[0] = D[k]; a[1] = L[k];       // a = T[k:k+2,k] is the vector of nonzeros in kth subdiagonal column.
      
      real_t anorm = sqrt(a[0]*a[0] + a[1]*a[1]); 

      // // Udrullet
      // //    reflection_vector(a,anorm,v);
      v[0] = D[k]; v[1] = L[k];
      real_t alpha = -copysign(anorm,a[0]); // Koster ingenting
      v[0] -= alpha;

      real_t vnorm = sqrt(v[0]*v[0]+v[1]*v[1]);
      real_t norm_inv = 1/vnorm;               /* Normalize */
      v[0] *= norm_inv;  v[1] *= norm_inv;

      Vout[2*k] = v[0]; Vout[2*k+1] = v[1];
      
      // // Udrullet 
      // //    apply_reflection(T({k,k+2},{k,k+3}),v);
      // //      if(k+1<n){			// k=n-1 case handled by padding with zeros
      real_t   vTA[3] = {D[ k ]*v[0] + L[ k ]*v[1],  // T(k,k  )*v[0] + T(k+1,k  )*v[1]
      			 U[ k ]*v[0] + D[k+1]*v[1],  // T(k,k+1)*v[0] + T(k+1,k+1)*v[1]
      			 U[(n+1)+k]*v[0] + U[k+1]*v[1]}; // T(k,k+2)*v[0] + T(k+1,k+2)*v[1]

      D[ k ]     -= 2*v[0]*vTA[0];
      L[ k ]     -= 2*v[1]*vTA[0];
      U[ k ]     -= 2*v[0]*vTA[1];
      D[k+1]     -= 2*v[1]*vTA[1];
      U[(n+1)+k] -= 2*v[0]*vTA[2];
      U[k+1]     -= 2*v[1]*vTA[2];
    }

  // Transform from the right = transform columns of the transpose.
  {
    int k = 0;
    const real_t *v = &Vout[0];
    real_t   vTA[2] = {D[ k ]*v[0] + U[  k  ]*v[1],  // T(k,k  )*v[0] + T(k,  k+1)*v[1]
  		          0        + D[ k+1 ]*v[1]}; // T(k+1,k)*v[0] + T(k+1,k+1)*v[1]. Lower subdiagonal is zero at this stage.
    
    D[k]       -= 2*v[0]*vTA[0]; // T(k,k)     -= 2*v[0]*vTA[0]
    U[k]       -= 2*v[1]*vTA[0]; // T(k,k+1)   -= 2*v[1]*vTA[0]
    L[k]       -= 2*v[0]*vTA[1]; // T(k+1,k)   -= 2*v[0]*vTA[1]
    D[k+1]     -= 2*v[1]*vTA[1]; // T(k+1,k+1) -= 2*v[1]*vTA[1]        
  }    

  for(int k=1;k<n-1;k++){
    const real_t *v = &Vout[2*k];

    real_t   vTA[3] = {U[k-1]*v[0] + U[(n+1)+k-1]*v[1], // T(k-1,k)*v[0] + T(k-1,k+1)*v[1]  
  		       D[ k ]*v[0] + U[  k  ]*v[1],     // T(k,k  )*v[0] + T(k,  k+1)*v[1]
  		       L[ k ]*v[0] + D[ k+1 ]*v[1]};    // T(k+1,k)*v[0] + T(k+1,k+1)*v[1]. Lower subdiagonal is zero at this stage

    U[k-1]     -= 2*v[0]*vTA[0];     // T(k-1,k)   -= 2*v[0]*vTA[0]
    U[(n+1)+(k-1)] -= 2*v[1]*vTA[0]; // T(k-1,k+1) -= 2*v[1]*vTA[0]
    D[k]       -= 2*v[0]*vTA[1];     // T(k,  k)     -= 2*v[0]*vTA[1]
    U[k]       -= 2*v[1]*vTA[1];     // T(k,  k+1)   -= 2*v[1]*vTA[1]
    L[k]       -= 2*v[0]*vTA[2];     // T(k+1,k)   -= 2*v[0]*vTA[2]
    D[k+1]     -= 2*v[1]*vTA[2];     // T(k+1,k+1) -= 2*v[1]*vTA[2]        
  } 
  
  // Copy working diagonals to output
  for(int i=0;i<n;i++){
    Dout[i] = D[i] + shift;	  // Diagonal
    if(i<n-1){
      Lout[i] = U[i];  // First lower subdiagonal. L[k] = T(k+1,k) is element below kth diagonal element.
    }
  }
}

template <typename scalar> 
void apply_all_reflections(const real_t *V, const int n, SpanMatrix& Q)
{
  if(Q.data != 0){       // Do we want eigenvectors?
    int m = Q.shape[0];
    
    for(int k=0;k<n;k++){
      const real_t &v0 = V[2*k], &v1 = V[2*k+1];      
      // Udrullet:
      //       apply_reflection(Q({k,k+2},{0,m}), v);
      for(int l=0;l<m;l++){
	scalar &q0 = Q.data[k*m+l], &q1 = Q.data[(k+1)*m+l];
	scalar vTA = q0*v0 + q1*v1;
	q0 -= 2*v0*vTA;
	q1 -= 2*v1*vTA;
      }      
    }  
  }  
}

array<real_t,2> eigvalsh2x2(const array<real_t,4> &A){
  auto [a,b,c,d] = A;
  real_t D = sqrt(4*b*c+(a-d)*(a-d));
  return {(a+d-D)/2, (a+d+D)/2};
}

void SymmetrizeT(real_matrix &A)
{
  auto [m,n] = A.shape;
  assert(m==n);
  
  for(int i=0;i<n-1;i++){
    real_t mean = 0.5*(A(i+1,i)+A(i,i+1));
    A(i+1,i) = mean;
    A(i,i+1) = mean;
  }
}

int nth_time = 0;

// TODO: Til tridiagonale matricer er Givens-rotation nemmere/hurtigere (kun een sqrt)
// TODO: Assumes all different eigenvalues. Does this break with multiples?
// TODO: Stop after max_steps for fixed k. Return max Gershgorin radius as convergence -- or max Rayleigh quotient residual?
// TODO: Implement implicit QR iteration using Francis' Q theorem/bulge chasing
template <typename scalar>
std::pair<real_t,int> eigensystem_hermitian(const SpanMatrix& A,
					       matrix<real_t>& lambdas, SpanMatrix Qt={0},
					       const real_t tolerance=1e4*machine_precision,
					       const int max_iterations=40)
{
  auto [m,n] = A.shape;
  //  assert(A.is_hermitian());  
  assert(lambdas.shape[0]*lambdas.shape[1] == n);
  
  scalar T_data[n*n];
  real_t Qhat_data[n*n], tmp_data[n*n];
  SpanMatrix T = A.copy(T_data);
  SpanMatrix Q = Qt.T();
  matrix<real_t> Qhat, tmp(tmp_data,{n,n});

  real_t max_error    = tolerance;
  int n_iterations = 0;

  if(Q.data !=0){ // Do we want to compute eigenvectors?
    Q    = identity(Q.data,{n,n}); // Yes, so initialize Q and Qhat to identity    
    Qhat = identity(Qhat_data,{n,n}); // Yes, so initialize Q and Qhat to identity
  }

  nth_time++;
  // 1. Initial O(n^3) decomposition A = Q T Q.T to tridiagonal form
  QHQ(T,Q);	
  if(!T.is_hermitian(sqrt(machine_precision))){
    //    printf("T is not Hermitian within absolute error " G ":\n%s\n+i*%s\n\n",sqrt(machine_precision ),string(Re(T,tmp)).c_str(),string(Im(T,tmp)).c_str());
    abort();
  }
  
  real_t D[n], L[n], V[2*(n-1)];
  for(int i=0;i<n;i++){
    D[i] = Re(T(i,i));
    L[i] = (i+1<n)? Re(T(i,i+1)) : 0;
  }

  // 2. After tridiagonal decomposition, we can do an eigenvalue
  //    QR-iteration step in O(n), and an eigenvector QR-iteration
  //    step in O(n^2).
  for(int k=n-1;k>=0;k--){
    // We start by targeting the (n,n)-eigenvalue, and gradually
    // deflate, working on smaller and smaller submatrices.
    real_t d = D[k];		// d = T(k,k)
    real_t shift = d;

    // The Gershgorin disk radius is defined by just the row-sums of
    // absolute off-diagonal elements, since T is symmetric. As T is
    // tridiagonal, only T(k,k-1),T(k,k), and T(k,k+1) are nonzero.
    // Thus, the k'th Gershgorin radius is just |T(k,k-1)| +
    // |T(k,k+1)| = |T(k,k-1)| + |T(k+1,k)| = |L[k-1]|+|L[k]|.
    int i=0;
    real_t GR = (k>0?fabs(L[k-1]):0)+fabs(L[k]);
    while(GR > tolerance){
      i++;   
      T_QTQ(k+1, D,L, D,L, V, shift);	
      apply_all_reflections(V,k,Qhat);
      
      GR = (k>0?fabs(L[k-1]):0)+(k+1<n?fabs(L[k]):0);      

      // Best guess to eigenvalue in position n-1,n-1.
      if(k>0){
	auto [l0,l1]  = eigvalsh2x2({D[k-1],L[k-1], /* T[(k-1):k, (k-1):k] */
				     L[k-1],D[k]  });

	shift    = fabs(l0-d) < fabs(l1-d)? l0 : l1; // Pick closest eigenvalue
      } else
	shift    = D[k];

      if(i>max_iterations){
	printf("%dth run: Cannot converge eigenvalue %d to tolerance " G
	       " using machine precision " G " (d=" G ", shift=" G ", G=" G ")\n"
	       "D[k] = " G ", L[k-1] = " G ", L[k] = " G "\n",
	       nth_time,k,tolerance,
	       machine_precision,d,shift,GR,
	       D[k], (k>0)?L[k-1]:0, (k+1<n)?L[k]:0);
	
	max_error = std::max(max_error,GR);
	break;
      }
      n_iterations++;
    }
  }
  for(int k=0;k<n;k++) lambdas[k] = D[k];

  if(Q.data != 0){
    //    std::cout << "Q = \n" << Q.str() << "\n\nQhat = \n" << Qhat.str() << "\n\n";
    Q.transpose() *= Qhat.transpose();
  }
  
  return {max_error,n_iterations};
}



#endif