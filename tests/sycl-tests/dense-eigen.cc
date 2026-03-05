#include <fullerenes/sycl-isomer-batch.hh>
#include <numeric>
#include <array>
#include "sycl/sycl.hpp"
#include <oneapi/dpl/random>
using namespace sycl;

template <typename T>
void apply_all_reflections(const sycl::group<1> &cta, const sycl::local_accessor<T,1>& V, const int n, const int m, multi_ptr<T, sycl::access::address_space::global_space, sycl::access::decorated::legacy> Q)
{   
    static_assert(std::is_floating_point<T>::value, "T must be floating point");
    auto tid = cta.get_local_linear_id();
    auto bdim = cta.get_local_linear_range();
    for(int k=0;k<n;k++){
        const T &v0 = V[2*k], &v1 = V[2*k+1];      
        // Udrullet:
        //       apply_reflection(Q({k,k+2},{0,m}), v);
        for(int l=tid;l<m; l+=bdim){
            T &q0 = Q[k*m+l], &q1 = Q[(k+1)*m+l];
            T vTA = q0*v0 + q1*v1;
            q0 -= 2*v0*vTA;
            q1 -= 2*v1*vTA;
        }      
    }  
}

template <typename T>
void T_QTQ(sycl::group<1>& cta, const int n, const sycl::local_accessor<T,1>& D, const sycl::local_accessor<T,1>& L, const sycl::local_accessor<T,1>& U, const sycl::local_accessor<T,1>& Vout, T shift=0){
    int tix = cta.get_local_linear_id();
    int bdim = cta.get_local_linear_range();
    FLOAT_TYPEDEFS(T);
    //  QTQ_calls ++;
    // Unrolled
    //  T numerical_zero = T.max_norm()*10*std::numeric_limits<T>::epsilon();
    // specialized max_norm = max(sum(abs(A),axis=1)) for tridiagonal matrix. 
    T local_max = T(0.);
    for (int i = tix; i < n; i += bdim){
        local_max = std::max(local_max, abs(D[i]) + 2*abs(L[i]));
    }
    T max_norm = reduce_over_group(cta, local_max, sycl::maximum<T>());
    T numerical_zero = 10*std::numeric_limits<T>::epsilon();
    T d_n, l_n, l_nm1;
    d_n = D[n]; l_n = L[n]; l_nm1 = L[n-1];
    sycl::group_barrier(cta);
    //T a[2], v[2], D[n+1], L[n+1], U[2*(n+1)];
    T a[2], v[2];//, D[n+1], L[n+1], U[2*(n+1)];
    for(int k = tix; k < n + 1; k += bdim){
        D[k] -= shift;
        U[n+1 + k] = T(0.);
        if(k < n-1){
            U[k] = L[k];
            Vout[2*k] = T(0.); Vout[2*k+1] = T(0.);
        } else {
            L[k] = T(0.);
            U[k] = T(0.);
        }
    }

    sycl::group_barrier(cta);
    if(tix == 0)
        for(int k=0;k<n-1;k++){
            if (sycl::abs(L[k]) > numerical_zero){
            a[0] = D[k]; a[1] = L[k];       // a = T[k:k+2,k] is the vector of nonzeros in kth subdiagonal column.
            
            T anorm = sqrt(a[0]*a[0] + a[1]*a[1]); 

            // // Udrullet
            // //    reflection_vector(a,anorm,v);
            v[0] = D[k]; v[1] = L[k];
            T alpha = -copysign(anorm,a[0]); // Koster ingenting
            v[0] -= alpha;

            T vnorm = sqrt(v[0]*v[0]+v[1]*v[1]);
            T norm_inv = T(1.)/vnorm;               //Normalize
            v[0] *= norm_inv;  v[1] *= norm_inv;

            Vout[2*k] = v[0]; Vout[2*k+1] = v[1];
            
            // // Udrullet 
            // //    apply_reflection(T({k,k+2},{k,k+3}),v);
            // //      if(k+1<n){			// k=n-1 case handled by padding with zeros
            coord3d vTA = { D[ k ]*v[0] + L[ k ]*v[1],  // T(k,k  )*v[0] + T(k+1,k  )*v[1]
                            U[ k ]*v[0] + D[k+1]*v[1],  // T(k,k+1)*v[0] + T(k+1,k+1)*v[1]
                            U[(n+1)+k]*v[0] + U[k+1]*v[1]}; // T(k,k+2)*v[0] + T(k+1,k+2)*v[1]

        
            D[k]     -= T(2.)*v[0]*vTA[0];
            L[k]     -= T(2.)*v[1]*vTA[0];
            U[k]     -= T(2.)*v[0]*vTA[1];
            D[k+1]     -= T(2.)*v[1]*vTA[1];
            U[(n+1)+k] -= T(2.)*v[0]*vTA[2];
            U[k+1]     -= T(2.)*v[1]*vTA[2];
            }
        }

    if(tix == 0)
    { // Transform from the right = transform columns of the transpose.
        int k = 0;
        const T *v = &Vout[0];
        T   vTA[2] = {D[ k ]*v[0] + U[  k  ]*v[1],  // T(k,k  )*v[0] + T(k,  k+1)*v[1]
                    0        + D[ k+1 ]*v[1]}; // T(k+1,k)*v[0] + T(k+1,k+1)*v[1]. Lower subdiagonal is zero at this stage.
        
        D[k]       -= T(2.)*v[0]*vTA[0]; // T(k,k)     -= 2*v[0]*vTA[0]
        U[k]       -= T(2.)*v[1]*vTA[0]; // T(k,k+1)   -= 2*v[1]*vTA[0]
        L[k]       -= T(2.)*v[0]*vTA[1]; // T(k+1,k)   -= 2*v[0]*vTA[1]
        D[k+1]     -= T(2.)*v[1]*vTA[1]; // T(k+1,k+1) -= 2*v[1]*vTA[1]        
    }
    sycl::group_barrier(cta);


    if(tix == 0){
        for(int k=1;k<n-1;k++){
            const T *v = &Vout[2*k];
            coord3d vTA = {U[k-1]*v[0] + U[(n+1)+k-1]*v[1], // T(k-1,k)*v[0] + T(k-1,k+1)*v[1]  
                            D[ k ]*v[0] + U[  k  ]*v[1],     // T(k,k  )*v[0] + T(k,  k+1)*v[1]
                            L[ k ]*v[0] + D[ k+1 ]*v[1]};    // T(k+1,k)*v[0] + T(k+1,k+1)*v[1]. Lower subdiagonal is zero at this stage

            U[k-1]     -= T(2.)*v[0]*vTA[0];     // T(k-1,k)   -= 2*v[0]*vTA[0]
            U[(n+1)+(k-1)] -= T(2.)*v[1]*vTA[0]; // T(k-1,k+1) -= 2*v[1]*vTA[0]
            U[k]       -= T(2.)*v[1]*vTA[1];     // T(k,  k+1)   -= 2*v[1]*vTA[1]
            D[k]       -= T(2.)*v[0]*vTA[1];     // T(k,  k)     -= 2*v[0]*vTA[1]
            L[k]       -= T(2.)*v[0]*vTA[2];     // T(k+1,k)   -= 2*v[0]*vTA[2]
            D[k+1]     -= T(2.)*v[1]*vTA[2];     // T(k+1,k+1) -= 2*v[1]*vTA[2]        
        }
    }


    sycl::group_barrier(cta);
    for (int k = tix; k<n; k+=bdim){  // Copy working diagonals to output
        D[k] += shift;
        if(k < n-1){
            L[k] = U[k];
        }
    }
    sycl::group_barrier(cta);
    if (tix==0){
        D[n] = d_n;
        L[n-1] = l_nm1;
        L[n] = l_n;
    }
    sycl::group_barrier(cta);

}

template <typename T>
std::array<T,2> eigvalsh2x2(const std::array<T,4> &A){
    auto [a,b,c,d] = A;
    T D = sqrt(4*b*c+(a-d)*(a-d));
    return {(a+d-D)/2, (a+d+D)/2};
}

template <typename T>
T max_norm(const sycl::local_accessor<T,1>& A, const int m, n){
{
  T mx = 0;
  int j = cta.get_local_id(0);
  for(int i=0;i<m;i++){ 	
    T row_norm = sycl::reduce_over_group(cta, sycl::abs(A[i*n+j]), sycl::plus<T>());
    mx = sycl::max(mx,row_norm);
  }
  return mx;  
}

template <typename T>
T reflection_vector(const T a,const T anorm)
{
    int tid = cta.get_local_id(0);
    T alpha = -copysign(anorm,a);
    T v = a + (tid==0)*alpha; // TODO: Check fortegn 
    T vnorm = sycl::sqrt(sycl::reduce_over_group(cta, v*v, sycl::plus<T>()));
    return v / vnorm;
}

/* TODO: Kan jeg håndtere a[0] =~= 0 mere robust? Kan jeg inkludere pivots? */
/* Note: Restricted to real matrix A, redo for complex. */
template <typename T>
void QHQ(sycl::group<1>& cta, /*in/out*/T *A, int n, T *Q)
{
  T numerical_zero = max_norm(A,n,n)*10*std::numeric_limits<T>::epsilon();
  
  int tid = cta.get_local_id(0);
  for(int k=0;k<n-1;k++){ // Sequential loop
    int l = n-(k+1);		         /* length of kth subdiagonal column */
    T a     = tid<l? A[k*n+(k+1)+j_tid] : 0;      /* a = A[k,(k+1):n], kth postdiagonal row. */
    T anorm = sycl::sqrt(sycl::reduce_over_group(cta, a*a, sycl::plus<T>())); /* Norm of a */

    if(anorm < numerical_zero) continue; /* Already eliminated, don't divide by 0 */

    v =  reflection_vector(a,anorm);    /* Vector definining elimination operations */
    
    reflect_region(A,n,/* reg. start:*/k+1,k, /* reg. size:*/l,l+1, v,2, ROWS);/* TODO: fix */
    reflect_region(A,n,/* reg. start:*/k+1,k, /* reg. size:*/l,l+1, v,2, COLS);

    if(Q != 0) reflect_region(Q,n, k+1,0, l,n, v, sigma, ROWS); 
  }
}

//A[i,j] = A[i,j] - sum(2*v[i]*A[j,k]*v[k],k)
template <typename T>
void reflect_region(sycl::local_accessor<T,1>& A, int cols)
{
  T vHA_j;
  int stride[2] = {N*(!cols) + 1*cols, N*cols + 1*(!cols)};

  int i_tid = cta.get_local_id(0);
  for(int j=0;j<n;j++){
    size_t IJ = (i_tid+i0)*stride[0] + (j+j0)*stride[1];
    T vAij = (i_tid<m)? v[i]*A[IJ] : 0;
    T val = sycl::reduce_over_group(cta, vAij, sycl::plus<T>());
    if(i_tid==j) vHA_j = val;
  }

  for(size_t i=0;i<m;i++){       /* A += -2*outer(v,vTA) */
    int j_tid = cta.get_local_id(0);
    size_t IJ = (i+i0)*stride[0] + (j_tid+j0)*stride[1];
    if(tid<n) A[IJ] -= 2*v[i]*vHA_j; 
  }  
}

