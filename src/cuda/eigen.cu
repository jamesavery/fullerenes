#include <cuda.h>
#include <cooperative_groups.h>
#include <cooperative_groups/reduce.h>
#include <cooperative_groups/scan.h>
#include "cuda_runtime.h"
#include <cuda_runtime_api.h>
#include "fullerenes/config.h"
#include "fullerenes/gpu/cu_array.hh"
#include "fullerenes/isomer_batch.hh"
#include "fullerenes/gpu/kernels.hh"

#if(CUSOLVER)
# include <cusolverDn.h>
# include <cusparse.h>
#endif

#define N_STREAMS 16


// TODO: T_QTQ based on Givens rotations (should be possible to do with fewer operations)
//size_t QTQ_calls = 0;
void T_QTQ(const int n, device_real_t *D, device_real_t *L, device_real_t *Vout, device_real_t shift=0)
{
  //  QTQ_calls ++;
  // Unrolled
  //  device_real_t numerical_zero = T.max_norm()*10*std::numeric_limits<device_real_t>::epsilon();
  // specialized max_norm = max(sum(abs(A),axis=1)) for tridiagonal matrix. 
  device_real_t max_norm = 0, numerical_zero = 10*max_norm*std::numeric_limits<device_real_t>::epsilon();
  for(int i=0;i<n;i++) max_norm = std::max(max_norm, std::abs(D[i]) + 2*std::abs(L[i]));
  
  device_real_t a[2], v[2], U[2*(n+1)];//, D[n+1], L[n+1];
  device_real_t d_n, l_n, l_nm1;
    d_n = D[n]; l_n = L[n]; l_nm1 = L[n-1];
  for(int i=0;i<n+1;i++){
    D[i] -= shift;		// Diagonal
    //L[i] = 0;			// Zero padding to avoid branching in inner loop
    //U[i] = 0;                   // Zero padding to avoid branching in inner loop
    U[(n+1)+i] = 0;		// Second upper diagonal for fill-in. U[n+k] = T(k,k+2) is the element two rows above (k+2)st diagonal element.
    if(i<n-1){
      L[ i ] = L[i];	// First lower subdiagonal. L[k] = T(k+1,k) is element below kth diagonal element.
      U[ i ] = L[i];	// First upper subdiagonal. U[k] = T(k,k+1) is element above (k+1)st diagonal element.
      Vout[2*i] = 0; Vout[2*i+1] = 0;	// i'th reflection vector. (0,0) yields "no reflection". Must be initialized due to skipped steps.          
    } else {
        L[ i ] = 0;		// Zero padding to avoid branching in inner loop
        U[ i ] = 0;		// Zero padding to avoid branching in inner loop
    }
  }
   
  for(int k=0;k<n-1;k++)
    if(fabs(L[k]) > numerical_zero)  // Only process if subdiagonal element is not already zero.
    {
      a[0] = D[k]; a[1] = L[k];       // a = T[k:k+2,k] is the vector of nonzeros in kth subdiagonal column.
      
      device_real_t anorm = sqrt(a[0]*a[0] + a[1]*a[1]); 

      // // Udrullet
      // //    reflection_vector(a,anorm,v);
      v[0] = D[k]; v[1] = L[k];
      device_real_t alpha = -copysign(anorm,a[0]); // Koster ingenting
      v[0] -= alpha;

      device_real_t vnorm = sqrt(v[0]*v[0]+v[1]*v[1]);
      device_real_t norm_inv = 1/vnorm;               /* Normalize */
      v[0] *= norm_inv;  v[1] *= norm_inv;

      Vout[2*k] = v[0]; Vout[2*k+1] = v[1];
      
      // // Udrullet 
      // //    apply_reflection(T({k,k+2},{k,k+3}),v);
      // //      if(k+1<n){			// k=n-1 case handled by padding with zeros
      device_real_t   vTA[3] = {D[ k ]*v[0] + L[ k ]*v[1],  // T(k,k  )*v[0] + T(k+1,k  )*v[1]
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
    const device_real_t *v = &Vout[0];
    device_real_t   vTA[2] = {D[ k ]*v[0] + U[  k  ]*v[1],  // T(k,k  )*v[0] + T(k,  k+1)*v[1]
  		          0        + D[ k+1 ]*v[1]}; // T(k+1,k)*v[0] + T(k+1,k+1)*v[1]. Lower subdiagonal is zero at this stage.
    
    D[k]       -= 2*v[0]*vTA[0]; // T(k,k)     -= 2*v[0]*vTA[0]
    U[k]       -= 2*v[1]*vTA[0]; // T(k,k+1)   -= 2*v[1]*vTA[0]
    L[k]       -= 2*v[0]*vTA[1]; // T(k+1,k)   -= 2*v[0]*vTA[1]
    D[k+1]     -= 2*v[1]*vTA[1]; // T(k+1,k+1) -= 2*v[1]*vTA[1]        
  }    

  for(int k=1;k<n-1;k++){
    const device_real_t *v = &Vout[2*k];

    device_real_t   vTA[3] = {U[k-1]*v[0] + U[(n+1)+k-1]*v[1], // T(k-1,k)*v[0] + T(k-1,k+1)*v[1]  
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
    D[i] = D[i] + shift;	  // Diagonal
    if(i<n-1){
      L[i] = U[i];  // First lower subdiagonal. L[k] = T(k+1,k) is element below kth diagonal element.
    }
  }
  D[n] = d_n;
  L[n-1] = l_nm1;
  L[n] = l_n;
}

void apply_all_reflections(const device_real_t *V, const int n, const int m, vector<device_real_t>& Q)
{
    for(int k=0;k<n;k++){
        const device_real_t &v0 = V[2*k], &v1 = V[2*k+1];      
        // Udrullet:
        //       apply_reflection(Q({k,k+2},{0,m}), v);
        for(int l=0;l<m;l++){
            device_real_t &q0 = Q[k*m+l], &q1 = Q[(k+1)*m+l];
            device_real_t vTA = q0*v0 + q1*v1;
            q0 -= 2*v0*vTA;
            q1 -= 2*v1*vTA;
        }      
    }  
}

array<device_real_t,2> eigvalsh2x2(const array<device_real_t,4> &A){
  auto [a,b,c,d] = A;
  device_real_t D = sqrt(4*b*c+(a-d)*(a-d));
  return {(a+d-D)/2, (a+d+D)/2};
}


int nth_time = 0;

// TODO: Til tridiagonale matricer er Givens-rotation nemmere/hurtigere (kun een sqrt)
// TODO: Assumes all different eigenvalues. Does this break with multiples?
// TODO: Stop after max_steps for fixed k. Return max Gershgorin radius as convergence -- or max Rayleigh quotient residual?
// TODO: Implement implicit QR iteration using Francis' Q theorem/bulge chasing
std::pair<device_real_t,size_t> eigensystem_hermitian(const int n, const 
                            vector<device_real_t>& D_, 
                            const vector<device_real_t>& L_, 
                            vector<device_real_t>& Q, 
					        vector<device_real_t>& lambdas,
					        const device_real_t tolerance=1e4*std::numeric_limits<device_real_t>::epsilon(),
					        const int max_iterations=5)
{
  device_real_t max_error = 0;
  int n_iterations = 0;

  //@Jonas: Herfra arbejder vi med en tridiagonal reel matrix. 
  device_real_t D[n + 1], L[n + 1], V[2*(n-1)];
  for(int i=0;i<n;i++){
    D[i] = D_[i];
    L[i] = (i+1<n)? L_[i] : 0;
  }

    for (int i = 0; i < n; i++) {
        Q[i*n+i] = device_real_t(1);
    }

  // 2. After tridiagonal decomposition, we can do an eigenvalue
  //    QR-iteration step in O(n), and an eigenvector QR-iteration
  //    step in O(n^2).
  for(int k=n-1;k>=0;k--){
    // We start by targeting the (n,n)-eigenvalue, and gradually
    // deflate, working on smaller and smaller submatrices.
    device_real_t d = D[k];		// d = T(k,k)
    device_real_t shift = d;

    // The Gershgorin disk radius is defined by just the row-sums of
    // absolute off-diagonal elements, since T is symmetric. As T is
    // tridiagonal, only T(k,k-1),T(k,k), and T(k,k+1) are nonzero.
    // Thus, the k'th Gershgorin radius is just |T(k,k-1)| +
    // |T(k,k+1)| = |T(k,k-1)| + |T(k+1,k)| = |L[k-1]|+|L[k]|.
    int i=0;
    device_real_t GR = (k>0?fabs(L[k-1]):0)+fabs(L[k]);
    int not_done = 1;    
    while(not_done > 0){	// GPU NB: Kan erstattes med fornuftig konstant antal iterationer, f.eks. 4-5 stykker.
      i++;   
      T_QTQ(k+1, D,L, V, shift);  // 
      apply_all_reflections(V,k,n,Q);
      
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
	//printf("%dth run: Cannot converge eigenvalue %d to tolerance using machine precision %g (d=%g, shift=%g, G=%g)\n D[k] = %g, L[k-1] = %g, L[k] = %g\n", nth_time,k,tolerance, std::numeric_limits<device_real_t>::epsilon(),d,shift,GR, D[k], (k>0)?L[k-1]:0, (k+1<n)?L[k]:0);
	
	max_error = std::max(max_error,GR);
	break;
      }
      n_iterations++;
    }
  }
  for(int k=0;k<n;k++) lambdas[k] = D[k]; // Extract eigenvalues into result.
  
  return {max_error,n_iterations};
}





namespace gpu_kernels{
    namespace isomerspace_eigen{
        template void spectrum_ends<GPU, float, uint16_t>(const IsomerBatch<GPU>& B, const CuArray<float>& hessians, const CuArray<uint16_t>& cols, CuArray<float>& lambda_mins, CuArray<float>& lambda_maxs, int m_lanczos, const LaunchCtx& ctx, const LaunchPolicy policy);
        template void spectrum_ends<GPU, double, uint16_t>(const IsomerBatch<GPU>& B, const CuArray<double>& hessians, const CuArray<uint16_t>& cols, CuArray<double>& lambda_mins, CuArray<double>& lambda_maxs, int m_lanczos, const LaunchCtx& ctx, const LaunchPolicy policy);
        template void spectrum_ends<GPU, float, uint16_t>(const IsomerBatch<GPU>& B, const CuArray<float>& hessians, const CuArray<uint16_t>& cols, CuArray<float>& lambda_mins, CuArray<float>& lambda_maxs, CuArray<float>& eigvect_mins, CuArray<float>& eigvect_maxs, int m_lanczos, const LaunchCtx& ctx, const LaunchPolicy policy);
        template void spectrum_ends<GPU, double, uint16_t>(const IsomerBatch<GPU>& B, const CuArray<double>& hessians, const CuArray<uint16_t>& cols, CuArray<double>& lambda_mins, CuArray<double>& lambda_maxs, CuArray<double>& eigvect_mins, CuArray<double>& eigvect_maxs, int m_lanczos, const LaunchCtx& ctx, const LaunchPolicy policy);
        template void eigensolve<GPU, float, uint16_t>(const IsomerBatch<GPU>& B, CuArray<float>& Q, const CuArray<float>& hessians, const CuArray<uint16_t>& cols, CuArray<float>& eigenvalues, const LaunchCtx& ctx, const LaunchPolicy policy);
        template void eigensolve<GPU, double, uint16_t>(const IsomerBatch<GPU>& B, CuArray<double>& Q, const CuArray<double>& hessians, const CuArray<uint16_t>& cols, CuArray<double>& eigenvalues, const LaunchCtx& ctx, const LaunchPolicy policy);
        template void eigensolve_special<GPU, float, uint16_t>(const IsomerBatch<GPU>& B, CuArray<float>& Q, const CuArray<float>& hessians, const CuArray<uint16_t>& cols, CuArray<float>& eigenvalues, const LaunchCtx& ctx, const LaunchPolicy policy);

        #include "device_includes.cu"
        enum class EigensolveMode {NO_VECTORS, VECTORS, ENDS, FULL_SPECTRUM, FULL_SPECTRUM_MOLECULE}; 
#if(CUSOLVER)      
        void eigensolve_cusolver(const IsomerBatch& B, const CuArray<T>& hessians, const CuArray<K>& cols, CuArray<T>& eigenvalues, const LaunchCtx& ctx, const LaunchPolicy policy){
            static std::vector<cusolverDnHandle_t> cusolverHs(N_STREAMS, NULL);
            static int m = B.n_atoms*3;
            static int lda = m;
            static bool initialized = false;
            static std::vector<CuArray<device_real_t>> As(N_STREAMS);
            static std::vector<LaunchCtx> ctxs(N_STREAMS);
            static int nisomers = B.isomer_capacity;
            static std::vector<CuArray<device_real_t>> d_works(N_STREAMS);
            static std::vector<std::vector<device_real_t>> h_works(N_STREAMS);
            static std::vector<int*> d_infos(N_STREAMS, nullptr);
            static std::vector<int> infos(N_STREAMS, 1);
            static std::vector<int> lworks(N_STREAMS, 0);
            static std::vector<size_t> workspaceInBytesOnDevice(N_STREAMS, 0);
            static std::vector<size_t> workspaceInBytesOnHost(N_STREAMS, 0);
            int ncols = 10*3;
            int nn = m*ncols;
            cusolverEigMode_t jobz = CUSOLVER_EIG_MODE_NOVECTOR; // compute eigenvalues and eigenvectors.
            cublasFillMode_t uplo = CUBLAS_FILL_MODE_LOWER;
            auto T0 = std::chrono::steady_clock::now();
            if (!initialized){
                for (size_t I = 0; I < N_STREAMS ; I++){
                    As[I] = CuArray<device_real_t>(m*m, 0.);
                    ctxs[I] = LaunchCtx(0);
                    cusolverDnCreate(&cusolverHs[I]);
                    cusolverDnSetStream(cusolverHs[I], ctxs[I].stream);   
                    cudaMalloc(reinterpret_cast<void **>(&(d_infos[I])), sizeof(int));
                  
                    #if FLOAT_TYPE == 3
                    //    cusolverDnDsyevd_bufferSize(cusolverHs[I], jobz, uplo, m, As[I].data, lda, &(eigenvalues.data[I*m]), &(lworks[I]));
                        cusolverDnXsyevd_bufferSize(cusolverHs[I], NULL, jobz, uplo, m, CUDA_R_64F, (const void*)As[I].data, lda, CUDA_R_64F, (const void*)&(eigenvalues.data[I*m]), CUDA_R_64F, &workspaceInBytesOnDevice[I], &workspaceInBytesOnHost[I]);  
                    #elif FLOAT_TYPE == 2
                        cusolverDnXsyevd_bufferSize(cusolverHs[I], NULL, jobz, uplo, m, CUDA_R_32F, (const void*)As[I].data, lda, CUDA_R_32F, (const void*)&(eigenvalues.data[I*m]), CUDA_R_32F, &workspaceInBytesOnDevice[I], &workspaceInBytesOnHost[I]);  
                    //  cusolverDnSsyevd_bufferSize(cusolverHs[I], jobz, uplo, m, As[I].data, lda, &(eigenvalues.data[I*m]), &(lworks[I]));
                    #endif
                    d_works[I] = CuArray<device_real_t>(workspaceInBytesOnDevice[I]/sizeof(device_real_t));
                    h_works[I].resize(workspaceInBytesOnHost[I]/sizeof(device_real_t));
                    //cudaMalloc(reinterpret_cast<void **>(&(d_works[I])), sizeof(device_real_t) * workspaceInBytesOnDevice[I]);
                    //cudaMallocHost(reinterpret_cast<void **>(&(h_works[I])), sizeof(device_real_t) * workspaceInBytesOnHost[I]);
                    printLastCudaError("eigensolve");
                }    
                initialized = true;
            }
            auto T1 = std::chrono::steady_clock::now();
            std::cout << "Initializing the eigensolver took " << std::chrono::duration_cast<std::chrono::microseconds>(T1 - T0).count() / (device_real_t)nisomers << " us / isomer" << std::endl;
            int counter = 0;
            //Loading the sparse hessians into dense matrices; Might be a bottleneck.
            auto start = std::chrono::steady_clock::now();
            auto fill_As = [&] (){
                int I = 0;
                auto start_counter = counter;
                for (int II = start_counter; II < std::min(start_counter + N_STREAMS, nisomers)  ; II++){
                counter++;
                for (size_t i = 0; i < m; i++){ //Number of rows in the hessian
                    for (size_t j = 0; j < ncols; j++){ //Number of columns in the hessian, it is 10*3 because we have a term for the node itself, it's 3 neighbours and 6 outer neighbours, each with 3 dx, dy, dz terms
                        int col = cols.data[int(II*nn + i*ncols + j)];
                        As[I][int(i*m + col)] = hessians.data[int(II*nn + i*ncols + j)];
                    }
                }
                ++I;
            }
            };
            
            auto end = std::chrono::steady_clock::now();
            std::cout << "Loading the sparse hessians into dense matrices took " << std::chrono::duration_cast<std::chrono::microseconds>(end - start).count() / (device_real_t)nisomers << " us / isomer" << std::endl;   

            if (policy == LaunchPolicy::SYNC){
                ctx.wait();
            }

            start = std::chrono::steady_clock::now();
            for (size_t I = 0; I < nisomers ; I++){
                int idx = I % N_STREAMS;
                if (idx == 0) fill_As();
                //std::cout  << "Inum: " << counter << std::endl;
                #if FLOAT_TYPE == 3
                    cusolverDnXsyevd( cusolverHs[idx], NULL, jobz, uplo, m, CUDA_R_64F, As[idx].data, lda, CUDA_R_64F, &eigenvalues.data[I*m], CUDA_R_64F, (void *) d_works[idx].data, workspaceInBytesOnDevice[idx], (void *) h_works[idx].data(), workspaceInBytesOnHost[idx], d_infos[idx]);
                //    cusolverDnDsyevd(cusolverHs[idx], jobz, uplo, m, As[idx].data,    lda, &eigenvalues.data[I*m], d_works[idx],   lworks[idx], d_infos[idx]);
                #elif FLOAT_TYPE == 2
                    cusolverDnXsyevd( cusolverHs[idx], NULL, jobz, uplo, m, CUDA_R_32F, As[idx].data, lda, CUDA_R_32F, &eigenvalues.data[I*m], CUDA_R_32F, (void *) d_works[idx].data, workspaceInBytesOnDevice[idx], (void *) h_works[idx].data(), workspaceInBytesOnHost[idx], d_infos[idx]);
                //    cusolverDnSsyevd(cusolverHs[idx], jobz, uplo, m, As[idx].data,    lda, &eigenvalues.data[I*m], d_works[idx],   lworks[idx], d_infos[idx]);
                #endif
                //cudaMemcpyAsync(&infos[idx], d_infos[idx], sizeof(int), cudaMemcpyDeviceToHost, ctx.stream);
            }

            if (policy == LaunchPolicy::SYNC){
                for (size_t I = 0; I < N_STREAMS ; I++){
                    ctxs[I].wait();
                }
            }
            end = std::chrono::steady_clock::now();
            std::cout << "eigensolve time: " << std::chrono::duration_cast<std::chrono::microseconds>(end - start).count() / (device_real_t)nisomers <<  " us / isomer" <<std::endl;
        
        }
#endif
      
        /* template <typename T>
        struct array_wrapper{
            T* data;
            int stride;
            __device__ array_wrapper(T* data, int stride) : data(data), stride(stride) {}
            T& __device__ operator[](int i) const{
                return data[i*stride];
            }
        };

        typedef array_wrapper<device_real_t> real_wrap; */
        template <typename T>
        void __device__ apply_all_reflections(const T *V, const int n, const int m, T* Q)
        {   
            static_assert(std::is_floating_point<T>::value, "T must be floating point");

            for(int k=0;k<n;k++){
                const T &v0 = V[2*k], &v1 = V[2*k+1];      
                // Udrullet:
                //       apply_reflection(Q({k,k+2},{0,m}), v);
                for(int l=threadIdx.x;l<m; l+=blockDim.x){
                    T &q0 = Q[k*m+l], &q1 = Q[(k+1)*m+l];
                    T vTA = q0*v0 + q1*v1;
                    q0 -= 2*v0*vTA;
                    q1 -= 2*v1*vTA;
                }      
            }  
        }
        //Customized diagonalization routine for symmetric tridiagonal matrices
        template <typename T>
        void __device__ T_QTQ(const int n, T* D, T* L, T* U, T* Vout, T shift=0)
        {
        int tix = threadIdx.x;
        FLOAT_TYPEDEFS(T);
        SMEM(T);
        //  QTQ_calls ++;
        // Unrolled
        //  real_t numerical_zero = T.max_norm()*10*std::numeric_limits<real_t>::epsilon();
        // specialized max_norm = max(sum(abs(A),axis=1)) for tridiagonal matrix. 
        real_t local_max = real_t(0.);
        for (int i = tix; i < n; i += blockDim.x){
            local_max = std::max(local_max, ABS(D[i]) + 2*ABS(L[i]));
        }
        real_t max_norm = reduction_max(smem, local_max);
        real_t numerical_zero = 10*std::numeric_limits<real_t>::epsilon();
        device_real_t d_n, l_n, l_nm1;
        d_n = D[n]; l_n = L[n]; l_nm1 = L[n-1];
        BLOCK_SYNC
        //real_t a[2], v[2], D[n+1], L[n+1], U[2*(n+1)];
        real_t a[2], v[2];//, D[n+1], L[n+1], U[2*(n+1)];
        for(int k = tix; k < n + 1; k += blockDim.x){
            D[k] -= shift;
            U[n+1 + k] = real_t(0.);
            if(k < n-1){
                U[k] = L[k];
                Vout[2*k] = real_t(0.); Vout[2*k+1] = real_t(0.);
            } else {
                L[k] = real_t(0.);
                U[k] = real_t(0.);
            }
        }
        
        BLOCK_SYNC
        if(tix == 0)
            for(int k=0;k<n-1;k++){
                if (ABS(L[k]) > numerical_zero){
                a[0] = D[k]; a[1] = L[k];       // a = T[k:k+2,k] is the vector of nonzeros in kth subdiagonal column.
                
                real_t anorm = SQRT(a[0]*a[0] + a[1]*a[1]); 

                // // Udrullet
                // //    reflection_vector(a,anorm,v);
                v[0] = D[k]; v[1] = L[k];
                real_t alpha = -copysign(anorm,a[0]); // Koster ingenting
                v[0] -= alpha;

                real_t vnorm = SQRT(v[0]*v[0]+v[1]*v[1]);
                real_t norm_inv = real_t(1.)/vnorm;               //Normalize
                v[0] *= norm_inv;  v[1] *= norm_inv;

                Vout[2*k] = v[0]; Vout[2*k+1] = v[1];
                
                // // Udrullet 
                // //    apply_reflection(T({k,k+2},{k,k+3}),v);
                // //      if(k+1<n){			// k=n-1 case handled by padding with zeros
                coord3d vTA = { D[ k ]*v[0] + L[ k ]*v[1],  // T(k,k  )*v[0] + T(k+1,k  )*v[1]
                                U[ k ]*v[0] + D[k+1]*v[1],  // T(k,k+1)*v[0] + T(k+1,k+1)*v[1]
                                U[(n+1)+k]*v[0] + U[k+1]*v[1]}; // T(k,k+2)*v[0] + T(k+1,k+2)*v[1]

            
                D[k]     -= real_t(2.)*v[0]*vTA[0];
                L[k]     -= real_t(2.)*v[1]*vTA[0];
                U[k]     -= real_t(2.)*v[0]*vTA[1];
                D[k+1]     -= real_t(2.)*v[1]*vTA[1];
                U[(n+1)+k] -= real_t(2.)*v[0]*vTA[2];
                U[k+1]     -= real_t(2.)*v[1]*vTA[2];
                }
            }
        
        if(tix == 0)
        { // Transform from the right = transform columns of the transpose.
            int k = 0;
            const real_t *v = &Vout[0];
            real_t   vTA[2] = {D[ k ]*v[0] + U[  k  ]*v[1],  // T(k,k  )*v[0] + T(k,  k+1)*v[1]
                        0        + D[ k+1 ]*v[1]}; // T(k+1,k)*v[0] + T(k+1,k+1)*v[1]. Lower subdiagonal is zero at this stage.
            
            D[k]       -= real_t(2.)*v[0]*vTA[0]; // T(k,k)     -= 2*v[0]*vTA[0]
            U[k]       -= real_t(2.)*v[1]*vTA[0]; // T(k,k+1)   -= 2*v[1]*vTA[0]
            L[k]       -= real_t(2.)*v[0]*vTA[1]; // T(k+1,k)   -= 2*v[0]*vTA[1]
            D[k+1]     -= real_t(2.)*v[1]*vTA[1]; // T(k+1,k+1) -= 2*v[1]*vTA[1]        
        }
        BLOCK_SYNC

        
        if(tix == 0){
            for(int k=1;k<n-1;k++){
                const real_t *v = &Vout[2*k];
                coord3d vTA = {U[k-1]*v[0] + U[(n+1)+k-1]*v[1], // T(k-1,k)*v[0] + T(k-1,k+1)*v[1]  
                                D[ k ]*v[0] + U[  k  ]*v[1],     // T(k,k  )*v[0] + T(k,  k+1)*v[1]
                                L[ k ]*v[0] + D[ k+1 ]*v[1]};    // T(k+1,k)*v[0] + T(k+1,k+1)*v[1]. Lower subdiagonal is zero at this stage

                U[k-1]     -= real_t(2.)*v[0]*vTA[0];     // T(k-1,k)   -= 2*v[0]*vTA[0]
                U[(n+1)+(k-1)] -= real_t(2.)*v[1]*vTA[0]; // T(k-1,k+1) -= 2*v[1]*vTA[0]
                U[k]       -= real_t(2.)*v[1]*vTA[1];     // T(k,  k+1)   -= 2*v[1]*vTA[1]
                D[k]       -= real_t(2.)*v[0]*vTA[1];     // T(k,  k)     -= 2*v[0]*vTA[1]
                L[k]       -= real_t(2.)*v[0]*vTA[2];     // T(k+1,k)   -= 2*v[0]*vTA[2]
                D[k+1]     -= real_t(2.)*v[1]*vTA[2];     // T(k+1,k+1) -= 2*v[1]*vTA[2]        
            }
        }
       

        BLOCK_SYNC
        for (int k = tix; k<n; k+=blockDim.x){  // Copy working diagonals to output
            D[k] += shift;
            if(k < n-1){
                L[k] = U[k];
            }
        }
        BLOCK_SYNC
        if (tix==0){
         D[n] = d_n;
         L[n-1] = l_nm1;
         L[n] = l_n;
        }
        BLOCK_SYNC
        
        }

        template <typename T>
        array<T,2> INLINE eigvalsh2x2(const array<T,4> &A){
            auto [a,b,c,d] = A;
            T D = SQRT(4*b*c+(a-d)*(a-d));
            return {(a+d-D)/2, (a+d+D)/2};
        }

        //Takes a set of tridiagonal matrices and solves them

        template<EigensolveMode mode, Device DEV, typename T, typename K>
        void __global__ eigensolve_(const IsomerBatch<DEV> B, CuArray<T> D_, CuArray<T> L_, CuArray<T> U_, CuArray<T> Q_, int n){
            FLOAT_TYPEDEFS(T);
            SMEM(T);
            T *D = smem + blockDim.x*2, *L = D + (n+1), *U = L + (n+1), *V = U + (n+1)*2;
            //Expected layout is that each thread reads the (threadIdx.x + blockIdx.x*blockDim.x)^th column of D and L, in that way reads should be coalesced.
            for (int I = blockIdx.x; I < B.isomer_capacity; I += gridDim.x) if(B.statuses[I] == IsomerStatus::CONVERGED){
                for(int i = threadIdx.x; i < n; i += blockDim.x){
                    D[i] = D_.data[n*I + i];
                    L[i] = L_.data[n*I + i];
                    U[i] = L_.data[n*I + i];
                }
                for (int i = threadIdx.x; i < n; i += blockDim.x){
                    Q_.data[n*n*I + i*(n+1)] = real_t(1.); //Set Q to the identity matrix
                }
            
            BLOCK_SYNC
                

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
                real_t GR = (k>0?ABS(L[k-1]):0)+ABS(L[k]);
                int not_done = 1;    
                while(not_done > 0){	// GPU NB: Kan erstattes med fornuftig konstant antal iterationer, f.eks. 4-5 stykker.
                i++;   
                T_QTQ(k+1, D,L, U, V, shift);  // 
                if(mode == EigensolveMode::VECTORS || mode == EigensolveMode::FULL_SPECTRUM_MOLECULE){
                    apply_all_reflections(V,k,n,Q_.data + n*n*I);
                }
                
                GR = (k>0?ABS(L[k-1]):0)+(k+1<n?ABS(L[k]):0);      

                // Best guess to eigenvalue in position n-1,n-1.
                if(k>0){
                    std::array<T,4> args = {D[k-1],L[k-1],L[k-1],D[k]};
                    auto [l0,l1]  = eigvalsh2x2(args);

                shift    = ABS(l0-d) < ABS(l1-d)? l0 : l1; // Pick closest eigenvalue
                } else
                shift    = D[k];
                
                if(GR <= std::numeric_limits<real_t>::epsilon()*real_t(10.)) not_done--; // Do one (or optionally more) steps after reaching tolerance, to get all off-diagonals below.
                                                // GPU NB: Se GPU NB ovenfor.
                if(i>5){
                //printf("%dth run: Cannot converge eigenvalue %d to tolerance " G " using machine precision %g (d=%g, shift=%g, G=%g)\n" "D[k] = %g, L[k-1] = %g, L[k] = %g\n", nth_time,k,tolerance, std::numeric_limits<real_t>::epsilon(),d,shift,GR, D[k], (k>0)?L[k-1]:0, (k+1<n)?L[k]:0);
                
                auto max_error = std::max(std::numeric_limits<real_t>::epsilon()*real_t(10.),GR);
                break;
                }

                }
            }
            BLOCK_SYNC
            //Copy back to global memory
            for (int i = threadIdx.x; i < n; i += blockDim.x){
                if( mode == EigensolveMode::FULL_SPECTRUM_MOLECULE){
                    if(i < 6) {D_.data[(n+6)*I + i] = 0;}
                    D_.data[(n+6)*I + 6 + i] = D[i];
                } else {
                    D_.data[n*I + i] = D[i];
                }
            }
        }
        }
        
        template<EigensolveMode mode, Device DEV, typename T, typename K>
        void __global__ eigensolve_min_max_(const IsomerBatch<DEV> B, CuArray<T> D_, CuArray<T> L_, CuArray<T> U_, CuArray<T> Q_, CuArray<T> EigMin_, CuArray<T> EigMax_, CuArray<int> MinIdx_, CuArray<int> MaxIdx_, int n){
            FLOAT_TYPEDEFS(T);
            SMEM(T);
            T *D = smem + blockDim.x*2, *L = D + (n+1), *U = L + (n+1), *V = U + (n+1)*2;
            //Expected layout is that each thread reads the (threadIdx.x + blockIdx.x*blockDim.x)^th column of D and L, in that way reads should be coalesced.
            for (int I = blockIdx.x; I < B.isomer_capacity; I += gridDim.x) if(B.statuses[I] == IsomerStatus::CONVERGED){
                for(int i = threadIdx.x; i < n; i += blockDim.x){
                    D[i] = D_.data[n*I + i];
                    L[i] = L_.data[n*I + i];
                    U[i] = L_.data[n*I + i];
                }
                if (mode == EigensolveMode::VECTORS)
                    for (int i =  threadIdx.x; i < n; i += blockDim.x){
                        Q_.data[n*n*I + i*(n+1)] = real_t(1.); //Set Q to the identity matrix
                    }
            
            BLOCK_SYNC
                

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
                real_t GR = (k>0?ABS(L[k-1]):0)+ABS(L[k]);
                int not_done = 1;    
                while(not_done > 0){	// GPU NB: Kan erstattes med fornuftig konstant antal iterationer, f.eks. 4-5 stykker.
                i++;   
                T_QTQ(k+1, D,L, U, V, shift);  // 
                if(mode == EigensolveMode::VECTORS) apply_all_reflections(V,k,n,Q_.data + n*n*I);
                
                GR = (k>0?ABS(L[k-1]):0)+(k+1<n?ABS(L[k]):0);      

                // Best guess to eigenvalue in position n-1,n-1.
                if(k>0){
                std::array<T,4> args = {D[k-1],L[k-1],
                            L[k-1],D[k]};
                auto [l0,l1]  = eigvalsh2x2(args   /* Diagonalize T[(k-1):k, (k-1):k] 2x2 submatrix */);

                shift    = ABS(l0-d) < ABS(l1-d)? l0 : l1; // Pick closest eigenvalue
                } else
                shift    = D[k];
                
                if(GR <= std::numeric_limits<real_t>::epsilon()*real_t(10.)) not_done--; // Do one (or optionally more) steps after reaching tolerance, to get all off-diagonals below.
                                                // GPU NB: Se GPU NB ovenfor.
                if(i>5){
                //printf("%dth run: Cannot converge eigenvalue %d to tolerance " G " using machine precision %g (d=%g, shift=%g, G=%g)\n" "D[k] = %g, L[k-1] = %g, L[k] = %g\n", nth_time,k,tolerance, std::numeric_limits<real_t>::epsilon(),d,shift,GR, D[k], (k>0)?L[k-1]:0, (k+1<n)?L[k]:0);
                
                auto max_error = std::max(std::numeric_limits<real_t>::epsilon()*real_t(10.),GR);
                break;
                }

                }
            }
            BLOCK_SYNC
            //Copy back to global memory
            for (int i = threadIdx.x; i < n; i += blockDim.x){
                D_.data[n*I + i] = D[i];}
            smem[threadIdx.x] = real_t(0.);
            real_t local_max = real_t(0.);
            real_t local_min = numeric_limits<real_t>::max();
            int local_min_idx = 0;
            int local_max_idx = 0;
            for (int i = threadIdx.x; i < n; i += blockDim.x){
                local_max = ISNAN(D[i]) ? NAN : std::max(local_max, ABS(D[i]));
                local_min = ISNAN(D[i]) ? NAN : std::min(local_min, ABS(D[i]));
                local_min_idx = ISNAN(D[i]) ? NAN : (local_min == ABS(D[i]) ? i : local_min_idx);
                local_max_idx = ISNAN(D[i]) ? NAN : (local_max == ABS(D[i]) ? i : local_max_idx);
            }
            real_t max_eig = reduction_max(smem, local_max);
            smem[threadIdx.x] = numeric_limits<real_t>::max();
            real_t min_eig = reduction_min(smem, local_min > 1e-1 ? local_min : numeric_limits<real_t>::max());
            if(threadIdx.x == 0){
                EigMax_.data[I] = max_eig;
                EigMin_.data[I] = min_eig;
            }
            BLOCK_SYNC
            //Argmax and argmin
            if (min_eig == D[local_min_idx]){
                //If by some miracle multiple eigenvalues are exactly equal, we just pick one of them at random using atomicExch_block
                atomicExch_block(MinIdx_.data + I, local_min_idx);
            }
            if (max_eig == D[local_max_idx]){
                //If by some miracle multiple eigenvalues are exactly equal, we just pick one of them at random using atomicExch_block
                atomicExch_block(MaxIdx_.data + I, local_max_idx);
            }
        }
        }

        template <EigensolveMode mode, Device DEV, typename T, typename K>
        void __global__ lanczos_(const IsomerBatch<DEV> B, CuArray<T> V_, CuArray<T> U, CuArray<T> D, const CuArray<T> H, const CuArray<K> cols, int m){
            TEMPLATE_TYPEDEFS(T,K);
            SMEM(T);
            int N = B.n_atoms * 3; //Number of rows in the hessian
            int atom_idx = threadIdx.x/3; //Atom index
            constexpr int M = 10*3;          //Number of columns in the hessian
            real_t* betas = smem + N;
            real_t* alphas = betas + m;
            real_t A[M]; //Hessian matrix, threadIdx.x'th row
            node_t C[M]; //Column indices of the threadIdx.x'th row 3-fold degenerate
            real_t* V;
            float* X_ptr = B.X + N*blockIdx.x; //WARNING float 
            real_t Z[6]; //Eigenvectors spanning the kernel of the hessian (Rotations, Translations)
            if (mode == EigensolveMode::ENDS || mode == EigensolveMode::FULL_SPECTRUM_MOLECULE){
                Z[0] = real_t(threadIdx.x%3 == 0)/SQRT(B.n_atoms); Z[1] = real_t(threadIdx.x%3 == 1)/SQRT(B.n_atoms); Z[2] = real_t(threadIdx.x%3 == 2)/SQRT(B.n_atoms); // Translation eigenvectors
                
            }
            BLOCK_SYNC
      
            clear_cache(smem, N);


            auto mat_vect = [&](const real_t x){
                real_t result = real_t(0);
                smem[threadIdx.x] = x;
                BLOCK_SYNC
                #pragma unroll
                for (int j = 0; j < M; j++){
                    int col = C[j];
                    result += A[j] * smem[col];
                }
                BLOCK_SYNC
                return result;
            };
            
            //Modified Gram-Schmidt, Also orthogonalizes against the deflation space
            auto MGS = [&](int index){
                BLOCK_SYNC
                real_t result = V[index*N];
                smem[threadIdx.x] = 0;
                if (mode == EigensolveMode::ENDS || mode == EigensolveMode::FULL_SPECTRUM_MOLECULE){
                    #pragma unroll
                    for (int j = 0; j < 6; j++){
                        auto proj = reduction(smem, result * Z[j]) * Z[j];
                        result -= proj; //Remove the component along Z[j] from result
                    }
                }

                #pragma unroll
                for (int j = 0; j < index; j++){
                    auto proj = reduction(smem, result * V[j*N]) * V[j*N];
                    result -= proj; //Remove the component along V[j*N] from result
                }
                result /= sqrt(reduction(smem, result * result));
                return result;
            };
           
            curandState state;            
            curand_init(42 + threadIdx.x, 0, 0, &state);

            for (int I = blockIdx.x; I < B.isomer_capacity; I += gridDim.x) if(B.statuses[I] == IsomerStatus::CONVERGED){
                V = V_.data + I * m * N + threadIdx.x;
                if(mode == EigensolveMode::ENDS || mode == EigensolveMode::FULL_SPECTRUM_MOLECULE){
                    X_ptr = B.X + N*I;
                    if (threadIdx.x%3 == 0) {
                        Z[3] = real_t(0.);
                        Z[4] = -X_ptr[atom_idx*3 + 2];
                        Z[5] = -X_ptr[atom_idx*3 + 1];
                    } else if (threadIdx.x%3 == 1) {
                        Z[3] = -X_ptr[atom_idx*3 + 2];
                        Z[4] = real_t(0.);
                        Z[5] = X_ptr[atom_idx*3 + 0];
                    } else {
                        Z[3] = X_ptr[atom_idx*3 + 1];
                        Z[4] = X_ptr[atom_idx*3 + 0];
                        Z[5] = real_t(0.);
                    }
                        clear_cache(smem, N);
                        Z[3] /= SQRT(reduction(smem, Z[3]*Z[3]));
                        clear_cache(smem, N);
                        Z[4] /= SQRT(reduction(smem, Z[4]*Z[4]));
                        clear_cache(smem, N);
                        Z[5] /= SQRT(reduction(smem, Z[5]*Z[5]));

                }
                //Load the hessian and cols into local memory
                memcpy(A, &H.data[I*N*M + threadIdx.x*M], M*sizeof(real_t));
                for (int j = 0; j < M; j++){ 
                    A[j] = H.data[I*N*M + threadIdx.x*M + j];
                    C[j] = cols.data[I*N*M + threadIdx.x*M + j];
                }
                BLOCK_SYNC
                /* for(int j = 0; j< 6; j++){
                clear_cache(smem, N);
                real_t test_ = mat_vect(Z[j]);
                clear_cache (smem, N);
                real_t rayleigh_ = reduction(smem, test_ * Z[j]);
                real_t resid_ = test_ - rayleigh_ * Z[j];

                print_single("\nRayleigh: \n");
                print_single(rayleigh_);
                print_single("\nResidual: \n");
                print_single(resid_);
                print_single("\n");
                }
 */
                //Lanczos algorithm 
                if(threadIdx.x == 0) betas[0] = real_t(0);
                real_t beta = real_t(0);
                real_t alpha = real_t(0);
                V[0*N] = curand_uniform(&state);
                smem[threadIdx.x] = real_t(0); //Clear the shared memory
                V[0*N] /= SQRT(reduction(smem, V[0*N] * V[0*N]));
                V[0*N] = MGS(0);
                for (int i = 0; i < m; i++){
                    if (i % 2 == 0 && i > 1){
                        V[(i-1)*N] = MGS(i-1);
                        V[i*N] = MGS(i);
                    }
                    real_t v = mat_vect(V[i*N]);
                    smem[threadIdx.x] = real_t(0); //Clear the shared memory
                    alpha = reduction(smem, v * V[i*N]);
                    if (threadIdx.x == i) alphas[i] = alpha;
                    if (i == 0){
                        v -= alpha * V[i*N];
                    } else {
                        v -= betas[i-1] * V[(i-1)*N] + alpha * V[i*N];
                    }
                    smem[threadIdx.x] = real_t(0); //Clear the shared memory
                    beta = SQRT(reduction(smem, v * v));
                    if (threadIdx.x == i) betas[i] = beta;
                    if (i < m-1) V[(i+1)*N] = v / beta;
                    //if (i < N-1) V[(i+1)*N] = beta;
                }
                if (threadIdx.x < m){
                    D.data[I*m + threadIdx.x] = alphas[threadIdx.x];
                    U.data[I*m + threadIdx.x] = betas[threadIdx.x];
                }
            }   
        }
        //Assumes that N = B.n_atoms * 3
        template <EigensolveMode mode, Device DEV,typename T, typename K>
        void __global__ compute_eigenvectors_(const IsomerBatch<DEV> B, CuArray<T> Q, CuArray<T> V, CuArray<T> E, int m){
            int atom_idx = threadIdx.x/3; //Atom index (Integer division so the result is rounded down)
            TEMPLATE_TYPEDEFS(T,K);
            SMEM(T);
            int n = B.n_atoms * 3;
            int offset = 0;
            if (mode == EigensolveMode::ENDS || mode == EigensolveMode::FULL_SPECTRUM_MOLECULE){
                offset = 6;
            }
            for (int I = blockIdx.x; I < B.isomer_capacity; I += gridDim.x){
                real_t* v = V.data + I * m * n;
                real_t* e = E.data + I * n * n;
                real_t* q = Q.data + I * m * m;
                if (mode == EigensolveMode::FULL_SPECTRUM_MOLECULE) {
                    for(int i = 0; i < 3; i++){
                        e[i*n + threadIdx.x] = real_t(threadIdx.x%3 == i)/SQRT(B.n_atoms); 
                    }
                    float* X_ptr = B.X + n*I; //WARNING float hardcoded here to match the type of B.X (IsomerBatch)
                    if (threadIdx.x%3 == 0) {
                        e[3*n + threadIdx.x] = real_t(0.);
                        e[4*n + threadIdx.x] = -X_ptr[atom_idx*3 + 2];
                        e[5*n + threadIdx.x] = -X_ptr[atom_idx*3 + 1];
                    } else if (threadIdx.x%3 == 1) {
                        e[3*n + threadIdx.x] = -X_ptr[atom_idx*3 + 2];
                        e[4*n + threadIdx.x] = real_t(0.);
                        e[5*n + threadIdx.x] = X_ptr[atom_idx*3 + 0];
                    } else {
                        e[3*n + threadIdx.x] = X_ptr[atom_idx*3 + 1];
                        e[4*n + threadIdx.x] = X_ptr[atom_idx*3 + 0];
                        e[5*n + threadIdx.x] = real_t(0.);
                    }
                        clear_cache(smem,  n + warpSize);
                        e[3*n + threadIdx.x] /= SQRT(reduction(smem, e[3*n + threadIdx.x]*e[3*n + threadIdx.x]));
                        clear_cache(smem,  n + warpSize);
                        e[4*n + threadIdx.x] /= SQRT(reduction(smem, e[4*n + threadIdx.x]*e[4*n + threadIdx.x]));
                        clear_cache(smem,  n + warpSize);
                        e[5*n + threadIdx.x] /= SQRT(reduction(smem, e[5*n + threadIdx.x]*e[5*n + threadIdx.x]));
                }
                
                if(threadIdx.x < n)
                for (int k = offset; k < n; k++){
                    e = E.data + I * n * n + k * n;
                    e[threadIdx.x] = real_t(0.);
                    q = Q.data + I * m * m + (k-offset) * m;
                    for (int i = 0; i < m; i++){
                        e[threadIdx.x] += v[i*n + threadIdx.x] * q[i];
                    }
                }
            }
        }
        
        template <Device DEV,typename T, typename K>
        void __global__ compute_eigenvectors_ends_(const IsomerBatch<DEV> B, CuArray<T> Q, CuArray<T> V, CuArray<T> Emin, CuArray<T> Emax, CuArray<int> MinIdx, CuArray<int> MaxIdx, int m){
            TEMPLATE_TYPEDEFS(T,K);
            int n = B.n_atoms * 3;
            for (int I = blockIdx.x; I < B.isomer_capacity; I += gridDim.x){
                int minidx = MinIdx.data[I];
                int maxidx = MaxIdx.data[I];
                real_t* emin = Emin.data + I * n;
                real_t* emax = Emax.data + I * n;
                real_t* v = V.data + I * m * n;
                real_t* qmin = Q.data + I * m * m + minidx * m;
                real_t* qmax = Q.data + I * m * m + maxidx * m;
                emin[threadIdx.x] = real_t(0.);
                emax[threadIdx.x] = real_t(0.);
                for (int i = 0; i < m; i++){
                    emin[threadIdx.x] += v[i*n + threadIdx.x] * qmin[i];
                    emax[threadIdx.x] += v[i*n + threadIdx.x] * qmax[i];
                }
            }
        }

        template <Device DEV,typename T, typename K>
        void eigensolve(const IsomerBatch<DEV>& B, CuArray<T>& Q, const CuArray<T>& hessians, const CuArray<K>& cols, CuArray<T>& eigenvalues, const LaunchCtx& ctx, const LaunchPolicy policy){
            if (policy == LaunchPolicy::SYNC) ctx.wait();
            FLOAT_TYPEDEFS(T);
            cudaSetDevice(B.get_device_id());
            auto dev = B.get_device_id();
            static int Nd = LaunchCtx::get_device_count();
            static std::vector<bool> init(Nd, false);
            if(policy == LaunchPolicy::SYNC) {ctx.wait();}
            static std::vector<CuArray<T>> Us(Nd); //Upper diagonals
            static std::vector<CuArray<T>> Ls(Nd); //Lower diagonals
            static std::vector<CuArray<T>> Vs(Nd); //Store Lanczos vectors
            static std::vector<CuArray<T>> Qs(Nd); //Store transformation matrices Q DEV,T,K Q^DEV,T,K = H
            static int m_natoms = B.n_atoms;
            static int m_isomer_capacity = B.isomer_capacity;
            static int n_deflation = 0;
            if (!init[dev] || m_natoms != B.n_atoms || m_isomer_capacity != B.isomer_capacity){
                init[dev] = true;
                Us[dev].resize(B.n_atoms*B.isomer_capacity*3); Us[dev].fill(0.); Us[dev].to_device(dev);
                Ls[dev].resize(B.n_atoms*B.isomer_capacity*3); Ls[dev].fill(0.); Ls[dev].to_device(dev);
                Vs[dev].resize(B.n_atoms*B.isomer_capacity*3*3*B.n_atoms); Vs[dev].fill(0.); Vs[dev].to_device(dev);
                Qs[dev].resize(B.n_atoms*B.isomer_capacity*3*3*B.n_atoms); Qs[dev].fill(0.); Qs[dev].to_device(dev);
                m_natoms = B.n_atoms;
                m_isomer_capacity = B.isomer_capacity;
            }

            size_t smem = sizeof(coord3d)*B.n_atoms*3 + sizeof(T)*Block_Size_Pow_2;
            size_t smem_qr = sizeof(T)*(B.n_atoms*3+1)*6 + sizeof(T)*(B.n_atoms*3)*2;
            size_t smem_eig = 0;

            static LaunchDims dims((void*)lanczos_<EigensolveMode::FULL_SPECTRUM,DEV,T,K>, B.n_atoms*3, smem, B.isomer_capacity);
            static LaunchDims qr_dims((void*)eigensolve_<EigensolveMode::VECTORS,DEV,T,K>, 64, smem_qr, B.isomer_capacity);
            static LaunchDims eig_dims((void*)compute_eigenvectors_<EigensolveMode::FULL_SPECTRUM,DEV,T,K>, B.n_atoms*3, smem_eig, B.isomer_capacity);
            dims.update_dims((void*)lanczos_<EigensolveMode::FULL_SPECTRUM,DEV,T,K>, B.n_atoms*3, smem, B.isomer_capacity);
            qr_dims.update_dims((void*)eigensolve_<EigensolveMode::VECTORS,DEV,T,K>, 64, smem_qr, B.isomer_capacity);
            eig_dims.update_dims((void*)compute_eigenvectors_<EigensolveMode::FULL_SPECTRUM,DEV,T,K>, B.n_atoms*3, smem_eig, B.isomer_capacity);
            
            //The hessian has 6 degrees of freedom, so in order to find the smallest eigenvalue we must find the 6 eigenvectors corresponding to these lambda=0
            //for (int i = 0; i <  6; i++){ 
            int Nlanczos = B.n_atoms*3;
            void* kargs[]{(void*)&B, (void*)&Vs[dev], (void*)&Ls[dev], (void*)&eigenvalues, (void*)&hessians, (void*)&cols, (void*)&Nlanczos};
            void* kargs_qr[]{(void*)&B, (void*)&eigenvalues, (void*)&Ls[dev], (void*)&Us[dev], (void*)&Qs[dev], (void*)&Nlanczos};
            void* kargs_vector[]{(void*)&B, (void*)&Qs[dev], (void*)&Vs[dev], (void*)&Q, (void*)&Nlanczos};
            safeCudaKernelCall((void*)lanczos_<EigensolveMode::FULL_SPECTRUM,DEV,T,K>, dims.get_grid(), dims.get_block(), kargs, smem, ctx.stream);
            //ctx.wait();
            //vector<T> Qhost(B.n_atoms*3*3*B.n_atoms);
            //vector<T> Diags = vector<T>(eigenvalues.data, eigenvalues.data + B.n_atoms*3);
            //vector<T> OffDiags = vector<T>(Ls[dev].data, Ls[dev].data + B.n_atoms*3);
            //vector<T> LambdaHost(B.n_atoms*3);
            //eigensystem_hermitian(Nlanczos, Diags, OffDiags, Qhost, LambdaHost);
            //ofstream out("Qhost.float32", ios::binary); out.write((char*)Qhost.data(), Qhost.size()*sizeof(T)); out.close();
            //ofstream out2("D.float32", ios::binary); out2.write((char*)eigenvalues.data, eigenvalues.size()*sizeof(T)); out2.close();
            //ofstream out3("L.float32", ios::binary); out3.write((char*)Ls[dev].data, Ls[dev].size()*sizeof(T)); out3.close();

            safeCudaKernelCall((void*)eigensolve_<EigensolveMode::VECTORS,DEV,T,K>, qr_dims.get_grid(), qr_dims.get_block(), kargs_qr, smem_qr, ctx.stream);
            safeCudaKernelCall((void*)compute_eigenvectors_<EigensolveMode::FULL_SPECTRUM,DEV,T,K>, eig_dims.get_grid(), eig_dims.get_block(), kargs_vector, smem_eig, ctx.stream);

            if (policy == LaunchPolicy::SYNC) ctx.wait();
            printLastCudaError("Full Spectrum Eigensolver Failed : ");
        }

        template <Device DEV,typename T, typename K>
        void eigensolve(const IsomerBatch<DEV>& B, const CuArray<T>& hessians, const CuArray<K>& cols, CuArray<T>& eigenvalues, const LaunchCtx& ctx, const LaunchPolicy policy){
            if (policy == LaunchPolicy::SYNC) ctx.wait();
            cudaSetDevice(B.get_device_id());
            FLOAT_TYPEDEFS(T);
            auto dev = B.get_device_id();
            static int Nd = LaunchCtx::get_device_count();
            static std::vector<bool> init(Nd, false);
            if(policy == LaunchPolicy::SYNC) {ctx.wait();}
            static std::vector<CuArray<T>> Us(Nd); //Upper diagonals
            static std::vector<CuArray<T>> Ls(Nd); //Lower diagonals
            static std::vector<CuArray<T>> Vs(Nd); //Store Lanczos vectors
            static std::vector<CuArray<T>> Qs(Nd); //Store transformation matrices Q DEV,T,K Q^DEV,T,K = H
            static int m_natoms = B.n_atoms;
            static int m_isomer_capacity = B.isomer_capacity;
            static int n_deflation = 0;
            if (!init[dev] || m_natoms != B.n_atoms || m_isomer_capacity != B.isomer_capacity){
                init[dev] = true;
                Us[dev].resize(B.n_atoms*B.isomer_capacity*3); Us[dev].fill(0.); Us[dev].to_device(dev);
                Ls[dev].resize(B.n_atoms*B.isomer_capacity*3); Ls[dev].fill(0.); Ls[dev].to_device(dev);
                Vs[dev].resize(B.n_atoms*B.isomer_capacity*3*3*B.n_atoms); Vs[dev].fill(0.); Vs[dev].to_device(dev);
                Qs[dev].resize(B.n_atoms*B.isomer_capacity*3*3*B.n_atoms); Qs[dev].fill(0.); Qs[dev].to_device(dev);
                m_natoms = B.n_atoms;
                m_isomer_capacity = B.isomer_capacity;
            }

            size_t smem = sizeof(coord3d)*B.n_atoms*3 + sizeof(T)*Block_Size_Pow_2;
            size_t smem_qr = sizeof(T)*(B.n_atoms*3+1)*6 + sizeof(T)*(B.n_atoms*3)*2;
            size_t smem_eig = 0;
            static LaunchDims dims((void*)lanczos_<EigensolveMode::FULL_SPECTRUM, DEV,T,K>, B.n_atoms*3, smem, B.isomer_capacity);
            static LaunchDims qr_dims((void*)eigensolve_<EigensolveMode::NO_VECTORS, DEV,T,K>, 64, smem_qr, B.isomer_capacity);

            dims.update_dims((void*)lanczos_<EigensolveMode::FULL_SPECTRUM, DEV,T,K>, B.n_atoms*3, smem, B.isomer_capacity);
            qr_dims.update_dims((void*)eigensolve_<EigensolveMode::NO_VECTORS, DEV,T,K>, 64, smem_qr, B.isomer_capacity);

            
            //The hessian has 6 degrees of freedom, so in order to find the smallest eigenvalue we must find the 6 eigenvectors corresponding to these lambda=0
            //for (int i = 0; i <  6; i++){ 
            int Nlanczos = B.n_atoms*3;
            void* kargs[]{(void*)&B, (void*)&Vs[dev], (void*)&Ls[dev], (void*)&eigenvalues, (void*)&hessians, (void*)&cols, (void*)&Nlanczos};
            void* kargs_qr[]{(void*)&B, (void*)&eigenvalues, (void*)&Ls[dev], (void*)&Us[dev], (void*)&Qs[dev], (void*)&Nlanczos};
            safeCudaKernelCall((void*)lanczos_<EigensolveMode::FULL_SPECTRUM, DEV,T,K>, dims.get_grid(), dims.get_block(), kargs, smem, ctx.stream);
            //ctx.wait();
            //vector<T> Qhost(B.n_atoms*3*3*B.n_atoms);
            //vector<T> Diags = vector<T>(eigenvalues.data, eigenvalues.data + B.n_atoms*3);
            //vector<T> OffDiags = vector<T>(Ls[dev].data, Ls[dev].data + B.n_atoms*3);
            //vector<T> LambdaHost(B.n_atoms*3);
            //eigensystem_hermitian(Nlanczos, Diags, OffDiags, Qhost, LambdaHost);
            //ofstream out("Qhost.float32", ios::binary); out.write((char*)Qhost.data(), Qhost.size()*sizeof(T)); out.close();
            //ofstream out2("D.float32", ios::binary); out2.write((char*)eigenvalues.data, eigenvalues.size()*sizeof(T)); out2.close();
            //ofstream out3("L.float32", ios::binary); out3.write((char*)Ls[dev].data, Ls[dev].size()*sizeof(T)); out3.close();
            safeCudaKernelCall((void*)eigensolve_<EigensolveMode::NO_VECTORS, DEV,T,K>, qr_dims.get_grid(), qr_dims.get_block(), kargs_qr, smem_qr, ctx.stream);
            if (policy == LaunchPolicy::SYNC) ctx.wait();
            printLastCudaError("Full Spectrum Eigensolver Failed : ");
        }

        template <Device DEV,typename T, typename K>
        void spectrum_ends(const IsomerBatch<DEV>& B, const CuArray<T>& hessians, const CuArray<K>& cols, CuArray<T>& lambda_mins, CuArray<T>& lambda_maxs, int m_lanczos, const LaunchCtx& ctx, const LaunchPolicy policy){
            if (policy == LaunchPolicy::SYNC) ctx.wait();
            FLOAT_TYPEDEFS(T);
            cudaSetDevice(B.get_device_id());
            auto dev = B.get_device_id();
            static int Nd = LaunchCtx::get_device_count();
            static std::vector<bool> init(Nd, false);
            if(policy == LaunchPolicy::SYNC) {ctx.wait();}
            static std::vector<CuArray<T>> Us(Nd); //Upper diagonals
            static std::vector<CuArray<T>> Ls(Nd); //Lower diagonals
            static std::vector<CuArray<T>> Ds(Nd); //Diagonals
            static std::vector<CuArray<T>> Vs(Nd); //Lanczos vectors
            static std::vector<CuArray<T>> Qs(Nd); //Q matrix for QR decomposition
            static std::vector<CuArray<int>> EigMinIdxs(Nd); //Indices of the smallest eigenvalues
            static std::vector<CuArray<int>> EigMaxIdxs(Nd); //Indices of the largest eigenvalues
            static int m_natoms = B.n_atoms;
            static int m_isomer_capacity = B.isomer_capacity;
            static int m = m_lanczos;

            if (!init[dev] || m != m_lanczos || m_natoms != B.n_atoms || m_isomer_capacity != B.isomer_capacity){
                init[dev] = true;
                m_natoms = B.n_atoms;
                m_isomer_capacity = B.isomer_capacity;
                m = m_lanczos;
                Us[dev].resize(B.isomer_capacity*m); Us[dev].fill(0.); Us[dev].to_device(dev);
                Ls[dev].resize(B.isomer_capacity*m); Ls[dev].fill(0.); Ls[dev].to_device(dev);
                Ds[dev].resize(B.isomer_capacity*m); Ds[dev].fill(0.); Ds[dev].to_device(dev);
                Vs[dev].resize(B.isomer_capacity*B.n_atoms*3*m); Vs[dev].fill(0.); Vs[dev].to_device(dev);
                Qs[dev].resize(B.isomer_capacity*m*m); Qs[dev].fill(0.); Qs[dev].to_device(dev);
                EigMinIdxs[dev].resize(B.isomer_capacity); EigMinIdxs[dev].fill(0); EigMinIdxs[dev].to_device(dev);
                EigMaxIdxs[dev].resize(B.isomer_capacity); EigMaxIdxs[dev].fill(0); EigMaxIdxs[dev].to_device(dev);
            }

            
            size_t smem = sizeof(T)*B.n_atoms*3*2 + m*2*sizeof(T);
            size_t smem_qr = sizeof(T)*(m+1)*6 + sizeof(T)*(64)*2;
            //size_t smem_transform = sizeof(T)*B.n_atoms*3*2;

            static LaunchDims dims((void*)lanczos_<EigensolveMode::ENDS,DEV,T,K>, B.n_atoms*3, smem, B.isomer_capacity);
            static LaunchDims qr_dims((void*)eigensolve_min_max_<EigensolveMode::NO_VECTORS,DEV,T,K>, 64, smem_qr, B.isomer_capacity);
            //static LaunchDims transform_dims((void*)compute_eigenvectors_, B.n_atoms*3, smem_transform, B.isomer_capacity);
            dims.update_dims((void*)lanczos_<EigensolveMode::ENDS,DEV,T,K>, B.n_atoms*3, smem, B.isomer_capacity);
            qr_dims.update_dims((void*)eigensolve_min_max_<EigensolveMode::NO_VECTORS,DEV,T,K>, 64, smem_qr, B.isomer_capacity);

            //transform_dims.update_dims((void*)compute_eigenvectors_, B.n_atoms*3, smem_transform, B.isomer_capacity);
            
            //The hessian has 6 degrees of freedom, so in order to find the smallest eigenvalue we must find the 6 eigenvectors corresponding to these lambda=0
            int n_deflate = 0;
            void* kargs[]{(void*)&B, (void*)&Vs[dev], (void*)&Ls[dev], (void*)&Ds[dev], (void*)&hessians, (void*)&cols, (void*)&m};
            void* kargs_qr[]{(void*)&B, (void*)&Ds[dev], (void*)&Ls[dev], (void*)&Us[dev], (void*)&Qs[dev], (void*)&lambda_mins, (void*)&lambda_maxs, (void*)&EigMinIdxs[dev], (void*)&EigMaxIdxs[dev], (void*)&m};
            //void* kargs_transform[]{(void*)&B, (void*)&Ds[dev], (void*)&Qs[dev], (void*)&Vs[dev], (void*)&Es[dev], (void*)&m, (void*)&n_deflate};
            safeCudaKernelCall((void*)lanczos_<EigensolveMode::ENDS,DEV,T,K>, dims.get_grid(), dims.get_block(), kargs, smem, ctx.stream);
            safeCudaKernelCall((void*)eigensolve_min_max_<EigensolveMode::NO_VECTORS,DEV,T,K>, qr_dims.get_grid(), qr_dims.get_block(), kargs_qr, smem_qr, ctx.stream);

            if (policy == LaunchPolicy::SYNC) ctx.wait();
            printLastCudaError("Spectrum Ends Failed: ");
        }

        template <Device DEV,typename T, typename K>
        void spectrum_ends(const IsomerBatch<DEV>& B, const CuArray<T>& hessians, const CuArray<K>& cols, CuArray<T>& lambda_mins, CuArray<T>& lambda_maxs, CuArray<T>& eigvect_mins, CuArray<T>& eigvect_maxs, int m_lanczos, const LaunchCtx& ctx, const LaunchPolicy policy){
            if (policy == LaunchPolicy::SYNC) ctx.wait();
            cudaSetDevice(B.get_device_id());
            auto dev = B.get_device_id();
            static int Nd = LaunchCtx::get_device_count();
            static std::vector<bool> init(Nd, false);
            if(policy == LaunchPolicy::SYNC) {ctx.wait();}
            static std::vector<CuArray<T>> Us(Nd); //Upper diagonals
            static std::vector<CuArray<T>> Ls(Nd); //Lower diagonals
            static std::vector<CuArray<T>> Ds(Nd); //Diagonals
            static std::vector<CuArray<T>> Vs(Nd); //Lanczos vectors
            static std::vector<CuArray<T>> Qs(Nd); //Q matrix for QR decomposition
            static std::vector<CuArray<int>> EigMinIdxs(Nd); //Indices of the smallest eigenvalues
            static std::vector<CuArray<int>> EigMaxIdxs(Nd); //Indices of the largest eigenvalues
            static int m_natoms = B.n_atoms;
            static int m_isomer_capacity = B.isomer_capacity;
            static int m = m_lanczos;

            if (!init[dev] || m != m_lanczos || m_natoms != B.n_atoms || m_isomer_capacity != B.isomer_capacity){
                init[dev] = true;
                m_natoms = B.n_atoms;
                m_isomer_capacity = B.isomer_capacity;
                m = m_lanczos;
                Us[dev].resize(B.isomer_capacity*m); Us[dev].fill(0.); Us[dev].to_device(dev);
                Ls[dev].resize(B.isomer_capacity*m); Ls[dev].fill(0.); Ls[dev].to_device(dev);
                Ds[dev].resize(B.isomer_capacity*m); Ds[dev].fill(0.); Ds[dev].to_device(dev);
                Vs[dev].resize(B.isomer_capacity*B.n_atoms*3*m); Vs[dev].fill(0.); Vs[dev].to_device(dev);
                Qs[dev].resize(B.isomer_capacity*m*m); Qs[dev].fill(0.); Qs[dev].to_device(dev);
                EigMinIdxs[dev].resize(B.isomer_capacity); EigMinIdxs[dev].fill(0); EigMinIdxs[dev].to_device(dev);
                EigMaxIdxs[dev].resize(B.isomer_capacity); EigMaxIdxs[dev].fill(0); EigMaxIdxs[dev].to_device(dev);
            }

            
            size_t smem = sizeof(T)*B.n_atoms*3*2 + m*2*sizeof(T);
            size_t smem_qr = sizeof(T)*(m+1)*6 + sizeof(T)*(64)*2;
            size_t smem_eigs = 0;


            static LaunchDims dims((void*)lanczos_<EigensolveMode::ENDS,DEV,T,K>, B.n_atoms*3, smem, B.isomer_capacity);
            static LaunchDims qr_dims((void*)eigensolve_min_max_<EigensolveMode::VECTORS,DEV,T,K>, 64, smem_qr, B.isomer_capacity);
            static LaunchDims eigs_dims((void*)compute_eigenvectors_ends_<DEV,T,K>, B.n_atoms*3, smem_eigs, B.isomer_capacity);

            dims.update_dims((void*)lanczos_<EigensolveMode::ENDS,DEV,T,K>, B.n_atoms*3, smem, B.isomer_capacity);
            qr_dims.update_dims((void*)eigensolve_min_max_<EigensolveMode::VECTORS,DEV,T,K>, 64, smem_qr, B.isomer_capacity);
            eigs_dims.update_dims((void*)compute_eigenvectors_ends_<DEV,T,K>, B.n_atoms*3, smem_eigs, B.isomer_capacity);
            
            //The hessian has 6 degrees of freedom, so in order to find the smallest eigenvalue we must find the 6 eigenvectors corresponding to these lambda=0
            int n_deflate = 0;
            void* kargs[]{(void*)&B, (void*)&Vs[dev], (void*)&Ls[dev], (void*)&Ds[dev], (void*)&hessians, (void*)&cols, (void*)&m};
            void* kargs_qr[]{(void*)&B, (void*)&Ds[dev], (void*)&Ls[dev], (void*)&Us[dev], (void*)&Qs[dev], (void*)&lambda_mins, (void*)&lambda_maxs, (void*)&EigMinIdxs[dev], (void*)&EigMaxIdxs[dev], (void*)&m};
            void* kargs_eigs[]{(void*)&B, (void*)&Qs[dev], (void*)&Vs[dev], (void*)&eigvect_mins, (void*)&eigvect_maxs, (void*)&EigMinIdxs[dev], (void*)&EigMaxIdxs[dev], (void*)&m};
            //void* kargs_transform[]{(void*)&B, (void*)&Ds[dev], (void*)&Qs[dev], (void*)&Vs[dev], (void*)&Es[dev], (void*)&m, (void*)&n_deflate};

            safeCudaKernelCall((void*)lanczos_<EigensolveMode::ENDS,DEV,T,K>, dims.get_grid(), dims.get_block(), kargs, smem, ctx.stream);
            safeCudaKernelCall((void*)eigensolve_min_max_<EigensolveMode::VECTORS,DEV,T,K>, qr_dims.get_grid(), qr_dims.get_block(), kargs_qr, smem_qr, ctx.stream);
            safeCudaKernelCall((void*)compute_eigenvectors_ends_<DEV,T,K>, eigs_dims.get_grid(), eigs_dims.get_block(), kargs_eigs, smem_eigs, ctx.stream);

            if (policy == LaunchPolicy::SYNC) ctx.wait();
            printLastCudaError("Spectrum Ends Failed: ");
        }

        template <Device U, typename T, typename K>        //Compute the full spectrum of the hessian matrix in the special case where the hessian has 6 zero eigenvalues pertaining to the 3 translational and 3 rotational degrees of freedom.
        void eigensolve_special(const IsomerBatch<U>& B, CuArray<T>& Q, const CuArray<T>& hessians, const CuArray<K>& cols, CuArray<T>& eigenvalues, const LaunchCtx& ctx, const LaunchPolicy policy){
            if (policy == LaunchPolicy::SYNC) ctx.wait();
            FLOAT_TYPEDEFS(T);
            cudaSetDevice(B.get_device_id());
            auto dev = B.get_device_id();
            static int Nd = LaunchCtx::get_device_count();
            static std::vector<bool> init(Nd, false);
            if(policy == LaunchPolicy::SYNC) {ctx.wait();}
            static std::vector<CuArray<T>> Us(Nd); //Upper diagonals
            static std::vector<CuArray<T>> Ls(Nd); //Lower diagonals
            static std::vector<CuArray<T>> Vs(Nd); //Store Lanczos vectors
            static std::vector<CuArray<T>> Qs(Nd); //Store transformation matrices Q DEV,T,K Q^DEV,T,K = H
            int m_natoms = B.n_atoms;
            int m_isomer_capacity = B.isomer_capacity;
            int n_deflation = 6;
            int n_eigs = B.n_atoms *3;
            if (!init[dev] || m_natoms != B.n_atoms || m_isomer_capacity != B.isomer_capacity){
                init[dev] = true;
                Us[dev].resize((n_eigs)*B.isomer_capacity); Us[dev].fill(0.); Us[dev].to_device(dev);
                Ls[dev].resize((n_eigs)*B.isomer_capacity); Ls[dev].fill(0.); Ls[dev].to_device(dev);
                Vs[dev].resize(B.isomer_capacity*(n_eigs) * n_eigs); Vs[dev].fill(0.); Vs[dev].to_device(dev);
                Qs[dev].resize(B.isomer_capacity*(n_eigs) * n_eigs); Qs[dev].fill(0.); Qs[dev].to_device(dev);
            }

            size_t smem = sizeof(coord3d)*n_eigs + sizeof(T)*Block_Size_Pow_2;
            size_t smem_qr = sizeof(T)*(n_eigs+1)*6 + sizeof(T)*(n_eigs)*2;
            size_t smem_eig = sizeof(T)*n_eigs*2; //We just need enough memory to compute reductions of eigenvectors such that we can normalize them.

            static LaunchDims dims((void*)lanczos_<EigensolveMode::FULL_SPECTRUM_MOLECULE,U,T,K>, n_eigs, smem, B.isomer_capacity);
            static LaunchDims qr_dims((void*)eigensolve_<EigensolveMode::FULL_SPECTRUM_MOLECULE,U,T,K>, 64, smem_qr, B.isomer_capacity);
            static LaunchDims eig_dims((void*)compute_eigenvectors_<EigensolveMode::FULL_SPECTRUM_MOLECULE,U,T,K>, n_eigs, smem_eig, B.isomer_capacity);
            dims.update_dims((void*)lanczos_<EigensolveMode::FULL_SPECTRUM_MOLECULE,U,T,K>, n_eigs, smem, B.isomer_capacity);
            qr_dims.update_dims((void*)eigensolve_<EigensolveMode::FULL_SPECTRUM_MOLECULE,U,T,K>, 64, smem_qr, B.isomer_capacity);
            eig_dims.update_dims((void*)compute_eigenvectors_<EigensolveMode::FULL_SPECTRUM_MOLECULE,U,T,K>, n_eigs, smem_eig, B.isomer_capacity);
            
            //The hessian has 6 degrees of freedom, so in order to find the smallest eigenvalue we must find the 6 eigenvectors corresponding to these lambda=0
            //for (int i = 0; i <  6; i++){ 
            int Nlanczos = n_eigs - n_deflation;
            void* kargs[]{(void*)&B, (void*)&Vs[dev], (void*)&Ls[dev], (void*)&eigenvalues, (void*)&hessians, (void*)&cols, (void*)&Nlanczos};
            void* kargs_qr[]{(void*)&B, (void*)&eigenvalues, (void*)&Ls[dev], (void*)&Us[dev], (void*)&Qs[dev], (void*)&Nlanczos};
            void* kargs_vector[]{(void*)&B, (void*)&Qs[dev], (void*)&Vs[dev], (void*)&Q, (void*)&Nlanczos};
            safeCudaKernelCall((void*)lanczos_<EigensolveMode::FULL_SPECTRUM_MOLECULE,U,T,K>, dims.get_grid(), dims.get_block(), kargs, smem, ctx.stream);
            //ctx.wait();
            //vector<T> Qhost(B.n_atoms*3*3*B.n_atoms);
            //vector<T> Diags = vector<T>(eigenvalues.data, eigenvalues.data + B.n_atoms*3);
            //vector<T> OffDiags = vector<T>(Ls[dev].data, Ls[dev].data + B.n_atoms*3);
            //vector<T> LambdaHost(B.n_atoms*3);
            //eigensystem_hermitian(Nlanczos, Diags, OffDiags, Qhost, LambdaHost);
            //ofstream out("Qhost.float32", ios::binary); out.write((char*)Qhost.data(), Qhost.size()*sizeof(T)); out.close();
            //ofstream out2("D.float32", ios::binary); out2.write((char*)eigenvalues.data, eigenvalues.size()*sizeof(T)); out2.close();
            //ofstream out3("L.float32", ios::binary); out3.write((char*)Ls[dev].data, Ls[dev].size()*sizeof(T)); out3.close();

            safeCudaKernelCall((void*)eigensolve_<EigensolveMode::FULL_SPECTRUM_MOLECULE,U,T,K>, qr_dims.get_grid(), qr_dims.get_block(), kargs_qr, smem_qr, ctx.stream);
            safeCudaKernelCall((void*)compute_eigenvectors_<EigensolveMode::FULL_SPECTRUM_MOLECULE,U,T,K>, eig_dims.get_grid(), eig_dims.get_block(), kargs_vector, smem_eig, ctx.stream);

            if (policy == LaunchPolicy::SYNC) ctx.wait();
            printLastCudaError("Full Spectrum Eigensolver Failed : ");

        }


    }
}
