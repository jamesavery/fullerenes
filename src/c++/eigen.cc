// TODO: T_QTQ based on Givens rotations (should be possible to do with fewer operations)
//size_t QTQ_calls = 0;
void T_QTQ(const int n, const real_t *Din, const real_t *Lin, real_t *Dout, real_t *Lout, real_t *Vout, real_t shift=0)
{
  //  QTQ_calls ++;
  // Unrolled
  //  real_t numerical_zero = T.max_norm()*10*std::numeric_limits<real_t>::epsilon();
  // specialized max_norm = max(sum(abs(A),axis=1)) for tridiagonal matrix. 
  real_t max_norm = 0, numerical_zero = 0;
  for(int i=0;i<n;i++) max_norm = std::max(max_norm, std::abs(Din[i]) + 2*std::abs(Lin[i]));
  numerical_zero = 10*max_norm*std::numeric_limits<real_t>::epsilon();
  
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
void apply_all_reflections(const real_t *V, const int n, matrix<scalar>& Q)
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
std::pair<real_t,size_t> eigensystem_hermitian(const matrix<scalar>& A,
					       matrix<real_t>& lambdas, matrix<scalar> Qt={0},
					       const real_t tolerance=1e4*std::numeric_limits<real_t>::epsilon(),
					       const int max_iterations=40)
{
  auto [m,n] = A.shape;
  //  assert(A.is_hermitian());  
  assert(lambdas.shape[0]*lambdas.shape[1] == n);
  
  scalar T_data[n*n];
  real_t Qhat_data[n*n], tmp_data[n*n];
  matrix<scalar> T = A.copy(T_data);
  matrix<scalar> Q = Qt.T();
  matrix<real_t> Qhat, tmp(tmp_data,{n,n});

  real_t max_error    = tolerance;
  size_t n_iterations = 0;

  if(Q.data !=0){ // Do we want to compute eigenvectors?
    Q    = identity(Q.data,{n,n}); // Yes, so initialize Q and Qhat to identity    
    Qhat = identity(Qhat_data,{n,n}); // Yes, so initialize Q and Qhat to identity
  }

  nth_time++;
  // 1. Initial O(n^3) decomposition A = Q T Q.T to tridiagonal form
  QHQ(T,Q);	
  if(!T.is_hermitian(sqrt(std::numeric_limits<real_t>::epsilon()))){
    //    printf("T is not Hermitian within absolute error " G ":\n%s\n+i*%s\n\n",sqrt(std::numeric_limits<real_t>::epsilon() ),string(Re(T,tmp)).c_str(),string(Im(T,tmp)).c_str());
    abort();
  }

  //@Jonas: Herfra arbejder vi med en tridiagonal reel matrix. 
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
    int not_done = 1;    
    while(not_done > 0){	// GPU NB: Kan erstattes med fornuftig konstant antal iterationer, f.eks. 4-5 stykker.
      i++;   
      T_QTQ(k+1, D,L, D,L, V, shift);  // 
      apply_all_reflections(V,k,Qhat);
      
      GR = (k>0?fabs(L[k-1]):0)+(k+1<n?fabs(L[k]):0);      

      // Best guess to eigenvalue in position n-1,n-1.
      if(k>0){
	auto [l0,l1]  = eigvalsh2x2({D[k-1],L[k-1],   /* Diagonalize T[(k-1):k, (k-1):k] 2x2 submatrix */
				     L[k-1],D[k]  });

	shift    = fabs(l0-d) < fabs(l1-d)? l0 : l1; // Pick closest eigenvalue
      } else
	shift    = D[k];
      
      if(GR <= tolerance) not_done--; // Do one (or optionally more) steps after reaching tolerance, to get all off-diagonals below.
                                      // GPU NB: Se GPU NB ovenfor.
      if(i>max_iterations){
	printf("%dth run: Cannot converge eigenvalue %d to tolerance " G " using machine precision %g (d=%g, shift=%g, G=%g)\n"
	       "D[k] = %g, L[k-1] = %g, L[k] = %g\n", nth_time,k,tolerance, std::numeric_limits<real_t>::epsilon(),d,shift,GR, D[k], (k>0)?L[k-1]:0, (k+1<n)?L[k]:0);
	
	max_error = std::max(max_error,GR);
	break;
      }
      n_iterations++;
    }
  }
  for(int k=0;k<n;k++) lambdas[k] = D[k]; // Extract eigenvalues into result.

  if(Q.data != 0)		// Are we doing eigenvectors? Then combine A->T and T->D transforms to produce eigenvectors.
    Q.transpose() *= Qhat.transpose();
  
  return {max_error,n_iterations};
}



