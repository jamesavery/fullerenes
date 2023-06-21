#include <cuda.h>
#include <cooperative_groups.h>
#include <cooperative_groups/reduce.h>
#include <cooperative_groups/scan.h>
#include "cuda_runtime.h"
#include <cuda_runtime_api.h>
#include "fullerenes/gpu/cuda_definitions.h"
#include "fullerenes/gpu/cu_array.hh"
#include "fullerenes/gpu/isomer_batch.hh"
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
        #include "device_includes.cu"

#if(CUSOLVER)      
        void eigensolve_cusolver(const IsomerBatch& B, const CuArray<device_real_t>& hessians, const CuArray<device_node_t>& cols, CuArray<device_real_t>& eigenvalues, const LaunchCtx& ctx, const LaunchPolicy policy){
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
            std::cout << "Initializing the eigensolver took " << std::chrono::duration_cast<std::chrono::microseconds>(T1 - T0).count() / (float)nisomers << " us / isomer" << std::endl;
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
            std::cout << "Loading the sparse hessians into dense matrices took " << std::chrono::duration_cast<std::chrono::microseconds>(end - start).count() / (float)nisomers << " us / isomer" << std::endl;   

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
            std::cout << "eigensolve time: " << std::chrono::duration_cast<std::chrono::microseconds>(end - start).count() / (float)nisomers <<  " us / isomer" <<std::endl;
        
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
        void __device__ apply_all_reflections(const device_real_t *V, const int n, const int m, device_real_t* Q)
        {
            if (threadIdx.x == 0)
            for(int k=0;k<n;k++){
                const device_real_t &v0 = V[2*k], &v1 = V[2*k+1];      
                // Udrullet:
                //       apply_reflection(Q({k,k+2},{0,m}), v);
                for(int l=0;l<m; l++){
                    device_real_t &q0 = Q[k*m+l], &q1 = Q[(k+1)*m+l];
                    device_real_t vTA = q0*v0 + q1*v1;
                    q0 -= 2*v0*vTA;
                    q1 -= 2*v1*vTA;
                }      
            }  
        }
        //Customized diagonalization routine for symmetric tridiagonal matrices
        void __device__ T_QTQ(const int n, device_real_t* D, device_real_t* L, device_real_t* U, device_real_t* Vout, device_real_t shift=0)
        {
        int tix = threadIdx.x;
        DEVICE_TYPEDEFS;
        extern __shared__ real_t shared[];
        //  QTQ_calls ++;
        // Unrolled
        //  real_t numerical_zero = T.max_norm()*10*std::numeric_limits<real_t>::epsilon();
        // specialized max_norm = max(sum(abs(A),axis=1)) for tridiagonal matrix. 
        real_t local_max = real_t(0.);
        for (int i = tix; i < n; i += blockDim.x){
            local_max = std::max(local_max, ABS(D[i]) + 2*ABS(L[i]));
        }
        real_t max_norm = reduction_max(shared, local_max);
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

        array<device_real_t,2> INLINE eigvalsh2x2(const array<device_real_t,4> &A){
            auto [a,b,c,d] = A;
            device_real_t D = SQRT(4*b*c+(a-d)*(a-d));
            return {(a+d-D)/2, (a+d+D)/2};
        }

        //Takes a set of tridiagonal matrices and solves them
        void __global__ eigensolve_(const IsomerBatch B, CuArray<device_real_t> D_, CuArray<device_real_t> L_, CuArray<device_real_t> U_, CuArray<device_real_t> Q_, int n){
            DEVICE_TYPEDEFS;
            extern __shared__ device_real_t smem[];
            device_real_t *D = smem + blockDim.x*2, *L = D + (n+1), *U = L + (n+1), *V = U + (n+1)*2;
            //Expected layout is that each thread reads the (threadIdx.x + blockIdx.x*blockDim.x)^th column of D and L, in that way reads should be coalesced.
            for (int I = blockIdx.x; I < B.isomer_capacity; I += gridDim.x) if(B.statuses[I] != IsomerStatus::EMPTY){
                for(int i = threadIdx.x; i < n; i += blockDim.x){
                    D[i] = D_.data[n*I + i];
                    L[i] = L_.data[n*I + i];
                    U[i] = L_.data[n*I + i];
                }
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
                apply_all_reflections(V,k,n,Q_.data + n*n*I);
                
                GR = (k>0?ABS(L[k-1]):0)+(k+1<n?ABS(L[k]):0);      

                // Best guess to eigenvalue in position n-1,n-1.
                if(k>0){
                auto [l0,l1]  = eigvalsh2x2({D[k-1],L[k-1],   /* Diagonalize T[(k-1):k, (k-1):k] 2x2 submatrix */
                                L[k-1],D[k]  });

                shift    = ABS(l0-d) < ABS(l1-d)? l0 : l1; // Pick closest eigenvalue
                } else
                shift    = D[k];
                
                if(GR <= std::numeric_limits<real_t>::epsilon()*real_t(10.)) not_done--; // Do one (or optionally more) steps after reaching tolerance, to get all off-diagonals below.
                                                // GPU NB: Se GPU NB ovenfor.
                if(i>20){
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
        }
        }

        void __global__ eigensolve_largest_(const IsomerBatch B, CuArray<device_real_t> D_, CuArray<device_real_t> L_, CuArray<device_real_t> U_, CuArray<device_real_t> Q_, CuArray<device_real_t> EigMax_, int n){
            DEVICE_TYPEDEFS;
            extern __shared__ device_real_t smem[];
            device_real_t *D = smem + blockDim.x*2, *L = D + (n+1), *U = L + (n+1), *V = U + (n+1)*2;
            //Expected layout is that each thread reads the (threadIdx.x + blockIdx.x*blockDim.x)^th column of D and L, in that way reads should be coalesced.
            for (int I = blockIdx.x; I < B.isomer_capacity; I += gridDim.x) if(B.statuses[I] != IsomerStatus::EMPTY){
                for(int i = threadIdx.x; i < n; i += blockDim.x){
                    D[i] = D_.data[n*I + i];
                    L[i] = L_.data[n*I + i];
                    U[i] = L_.data[n*I + i];
                }
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
                //apply_all_reflections(V,k,n,Q_.data + n*n*I);
                
                GR = (k>0?ABS(L[k-1]):0)+(k+1<n?ABS(L[k]):0);      

                // Best guess to eigenvalue in position n-1,n-1.
                if(k>0){
                auto [l0,l1]  = eigvalsh2x2({D[k-1],L[k-1],   /* Diagonalize T[(k-1):k, (k-1):k] 2x2 submatrix */
                                L[k-1],D[k]  });

                shift    = ABS(l0-d) < ABS(l1-d)? l0 : l1; // Pick closest eigenvalue
                } else
                shift    = D[k];
                
                if(GR <= std::numeric_limits<real_t>::epsilon()*real_t(10.)) not_done--; // Do one (or optionally more) steps after reaching tolerance, to get all off-diagonals below.
                                                // GPU NB: Se GPU NB ovenfor.
                if(i>20){
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
            for (int i = threadIdx.x; i < n; i += blockDim.x){
                local_max = std::max(local_max, D[i]);
            }
            real_t max_eig = reduction_max(smem, local_max);
            if(threadIdx.x == 0){
                EigMax_.data[I] = max_eig;
            }
        }
        }


        void __global__ lanczos_(const IsomerBatch B, CuArray<device_real_t> E_, CuArray<device_real_t> V_, CuArray<device_real_t> U, CuArray<device_real_t> D, const CuArray<device_real_t> H, const CuArray<device_node_t> cols, int m, int n_deflation){
            DEVICE_TYPEDEFS;
            extern __shared__ real_t smem[];
            int N = B.n_atoms * 3; //Number of rows in the hessian
            constexpr int M = 10*3;          //Number of columns in the hessian
            real_t* betas = smem + N;
            real_t* alphas = betas + m;
            real_t A[M]; //Hessian matrix, threadIdx.x'th row
            node_t C[M]; //Column indices of the threadIdx.x'th row 3-fold degenerate
            real_t* V;
            real_t* E;
            
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
                #pragma unroll
                for (int j = 0; j < index; j++){
                    auto proj = reduction(smem, result * V[j*N]) * V[j*N];
                    result -= proj; //Remove the component along V[j*N] from result
                }
                #pragma unroll
                for (int j = 0; j < n_deflation; j++){
                    auto proj = reduction(smem, result * E[j*N]) * E[j*N];
                    result -= proj; //Remove the component along E[j*N] from result 
                }
                result /= sqrt(reduction(smem, result * result));
                return result;
            };
           
            curandState state;            
            curand_init(42 + threadIdx.x, 0, 0, &state);

            for (int I = blockIdx.x; I < B.isomer_capacity; I += gridDim.x) if(B.statuses[I] != IsomerStatus::EMPTY){
                V = V_.data + I * m * N + threadIdx.x;
                E = E_.data + I * 6 * N + threadIdx.x; //6 eigenvectors per isomer (Related to 0 eigenvalues, 6 degrees of freedom)
                //Load the hessian and cols into local memory
                memcpy(A, &H.data[I*N*M + threadIdx.x*M], M*sizeof(real_t));
                for (int j = 0; j < M; j++){ 
                    A[j] = H.data[I*N*M + threadIdx.x*M + j];
                    C[j] = cols.data[I*N*M + threadIdx.x*M + j];
                }

                //Lanczos algorithm 
                if(threadIdx.x == 0) betas[0] = real_t(0);
                real_t beta = real_t(0);
                real_t alpha = real_t(0);
                V[0*N] = curand_uniform(&state);
                V[0*N] /= SQRT(reduction(smem, V[0*N] * V[0*N]));
                for (int i = 0; i < m; i++){
                    if (i % 2 == 0 && i > 1){
                        V[(i-1)*N] = MGS(i-1);
                        V[i*N] = MGS(i);
                        //if(threadIdx.x + blockIdx.x == 0) printf("i = %d, N = %d, V[i*N] = %f\n", i, N, V[i*N]);
                    }
                    real_t v = mat_vect(V[i*N]);
                    alpha = reduction(smem, v * V[i*N]);
                    if (threadIdx.x == i) alphas[i] = alpha;
                    if (i == 0){
                        v -= alpha * V[i*N];
                    } else {
                        v -= betas[i-1] * V[(i-1)*N] + alpha * V[i*N];
                    }
                    beta = SQRT(reduction(smem, v * v));
                    if (threadIdx.x == i) betas[i] = beta;
                    if (i < N-1) V[(i+1)*N] = v / beta;
                    //if (i < N-1) V[(i+1)*N] = beta;
                }
                if (threadIdx.x < m){
                    D.data[I*m + threadIdx.x] = ISNAN(alphas[threadIdx.x])  ? real_t(0) : alphas[threadIdx.x];
                    U.data[I*m + threadIdx.x] = ISNAN(betas[threadIdx.x]) ? real_t(0) : betas[threadIdx.x];
                }
            }   
        }
        //Assumes that N = B.n_atoms * 3
        void __global__ compute_eigenvector_(const IsomerBatch B, CuArray<device_real_t> Q, CuArray<device_real_t> V, CuArray<device_real_t> E, int m, int i_th_vector){
            DEVICE_TYPEDEFS;
            int n = B.n_atoms * 3;
            for (int I = blockIdx.x; I < B.isomer_capacity; I += gridDim.x){
                real_t* v = V.data + I * m * n;
                real_t* e = E.data + (I * 6 + i_th_vector) * n;
                real_t* q = Q.data + I * m * m; //Always access the 1st eigenvector in Q

                e[threadIdx.x] = real_t(0.);
                for (int i = 0; i < m; i++){
                    e[threadIdx.x] += v[i*n + threadIdx.x] * q[i];
                }

            }
        }

        void eigensolve(const IsomerBatch& B, CuArray<device_real_t>& Q, const CuArray<device_real_t>& hessians, const CuArray<device_node_t>& cols, CuArray<device_real_t>& eigenvalues, const LaunchCtx& ctx, const LaunchPolicy policy){
            if (policy == LaunchPolicy::SYNC) ctx.wait();
            cudaSetDevice(B.get_device_id());
            auto dev = B.get_device_id();
            static int Nd = LaunchCtx::get_device_count();
            static std::vector<bool> init(Nd, false);
            if(policy == LaunchPolicy::SYNC) {ctx.wait();}
            static std::vector<CuArray<device_real_t>> Us(Nd); //Upper diagonals
            static std::vector<CuArray<device_real_t>> Ls(Nd); //Lower diagonals
            static std::vector<CuArray<device_real_t>> Ds(Nd); //Diagonals
            static std::vector<CuArray<device_real_t>> Vs(Nd); //Diagonals
            static int m_natoms = B.n_atoms;
            static int m_isomer_capacity = B.isomer_capacity;

            if (!init[dev] || m_natoms != B.n_atoms || m_isomer_capacity != B.isomer_capacity){
                init[dev] = true;
                Us[dev].resize(B.n_atoms*B.isomer_capacity*3); Us[dev].to_device(dev);
                Ls[dev].resize(B.n_atoms*B.isomer_capacity*3); Ls[dev].to_device(dev);
                Ds[dev].resize(B.n_atoms*B.isomer_capacity*3); Ds[dev].to_device(dev);
                Vs[dev].resize(B.n_atoms*B.isomer_capacity*3*3*B.n_atoms); Vs[dev].to_device(dev);
                m_natoms = B.n_atoms;
                m_isomer_capacity = B.isomer_capacity;
            }

            size_t smem = sizeof(device_coord3d)*B.n_atoms*3 + sizeof(device_real_t)*Block_Size_Pow_2;
            size_t smem_qr = sizeof(device_real_t)*(B.n_atoms*3+1)*6 + sizeof(device_real_t)*(B.n_atoms*3)*2;
            static LaunchDims dims((void*)lanczos_, B.n_atoms*3, smem, B.isomer_capacity);
            static LaunchDims qr_dims((void*)eigensolve_, 64, smem_qr, B.isomer_capacity);
            dims.update_dims((void*)lanczos_, B.n_atoms*3, smem, B.isomer_capacity);
            qr_dims.update_dims((void*)eigensolve_, 64, smem_qr, B.isomer_capacity);
            
            //The hessian has 6 degrees of freedom, so in order to find the smallest eigenvalue we must find the 6 eigenvectors corresponding to these lambda=0
            //for (int i = 0; i <  6; i++){ 
            int Nlanczos = B.n_atoms*3;
            void* kargs[]{(void*)&B, (void*)&Vs[dev], (void*)&Ls[dev], (void*)&Ds[dev], (void*)&hessians, (void*)&cols, (void*)&Nlanczos};
            void* kargs_qr[]{(void*)&B, (void*)&Ds[dev], (void*)&Ls[dev], (void*)&Us[dev], (void*)&Vs[dev], (void*)&Nlanczos};
            safeCudaKernelCall((void*)lanczos_, dims.get_grid(), dims.get_block(), kargs, smem, ctx.stream);
            safeCudaKernelCall((void*)eigensolve_, qr_dims.get_grid(), qr_dims.get_block(), kargs_qr, smem_qr, ctx.stream);
            if (policy == LaunchPolicy::SYNC) ctx.wait();
            printLastCudaError("Full Spectrum Eigensolver Failed : ");
        }

        void spectrum_ends(const IsomerBatch& B, const CuArray<device_real_t>& hessians, const CuArray<device_node_t>& cols, CuArray<device_real_t>& lambda_mins, CuArray<device_real_t>& lambda_maxs, const LaunchCtx& ctx, const LaunchPolicy policy){
            if (policy == LaunchPolicy::SYNC) ctx.wait();
            cudaSetDevice(B.get_device_id());
            auto dev = B.get_device_id();
            static int Nd = LaunchCtx::get_device_count();
            static std::vector<bool> init(Nd, false);
            if(policy == LaunchPolicy::SYNC) {ctx.wait();}
            static std::vector<CuArray<device_real_t>> Us(Nd); //Upper diagonals
            static std::vector<CuArray<device_real_t>> Ls(Nd); //Lower diagonals
            static std::vector<CuArray<device_real_t>> Ds(Nd); //Diagonals
            static std::vector<CuArray<device_real_t>> Vs(Nd); //Lanczos vectors
            static std::vector<CuArray<device_real_t>> Qs(Nd); //Q matrix for QR decomposition
            static std::vector<CuArray<device_real_t>> Es(Nd); //6 Eigenvectors to deflate 'H' with
            static int m_natoms = B.n_atoms;
            static int m_isomer_capacity = B.isomer_capacity;
            static int m = B.n_atoms*3;

            if (!init[dev] || m_natoms != B.n_atoms || m_isomer_capacity != B.isomer_capacity){
                init[dev] = true;
                Us[dev].resize(B.isomer_capacity*m); Us[dev].fill(0.); Us[dev].to_device(dev);
                Ls[dev].resize(B.isomer_capacity*m); Ls[dev].fill(0.); Ls[dev].to_device(dev);
                Ds[dev].resize(B.isomer_capacity*m); Ds[dev].fill(0.); Ds[dev].to_device(dev);
                Vs[dev].resize(B.isomer_capacity*B.n_atoms*3*m); Vs[dev].fill(0.); Vs[dev].to_device(dev);
                Qs[dev].resize(B.isomer_capacity*m*m); Qs[dev].fill(0.); Qs[dev].to_device(dev);
                Es[dev].resize(B.isomer_capacity*B.n_atoms*3*6); Es[dev].fill(0.); Es[dev].to_device(dev); 
                m_natoms = B.n_atoms;
                m_isomer_capacity = B.isomer_capacity;
            }

            
            size_t smem = sizeof(device_real_t)*B.n_atoms*3*2 + m*2*sizeof(device_real_t);
            size_t smem_qr = sizeof(device_real_t)*(m+1)*6 + sizeof(device_real_t)*(64)*2;
            size_t smem_transform = 0;
            static LaunchDims dims((void*)lanczos_, B.n_atoms*3, smem, B.isomer_capacity);
            static LaunchDims qr_dims((void*)eigensolve_, 64, smem_qr, B.isomer_capacity);
            static LaunchDims transform_dims((void*)compute_eigenvector_, B.n_atoms*3, smem_transform, B.isomer_capacity);
            dims.update_dims((void*)lanczos_, B.n_atoms*3, smem, B.isomer_capacity);
            qr_dims.update_dims((void*)eigensolve_, 64, smem_qr, B.isomer_capacity);
            transform_dims.update_dims((void*)compute_eigenvector_, B.n_atoms*3, smem_transform, B.isomer_capacity);
            
            //The hessian has 6 degrees of freedom, so in order to find the smallest eigenvalue we must find the 6 eigenvectors corresponding to these lambda=0
            for (int i = 0; i <  1; i++){ 
                void* kargs[]{(void*)&B, (void*)&Es[dev], (void*)&Vs[dev], (void*)&Ls[dev], (void*)&Ds[dev], (void*)&hessians, (void*)&cols, (void*)&m, (void*)&i};
                void* kargs_qr[]{(void*)&B, (void*)&Ds[dev], (void*)&Ls[dev], (void*)&Us[dev], (void*)&Qs[dev], (void*)&m};
                void* kargs_transform[]{(void*)&B, (void*)&Qs[dev], (void*)&Vs[dev], (void*)&Es[dev], (void*)&m, (void*)&i};
                safeCudaKernelCall((void*)lanczos_, dims.get_grid(), dims.get_block(), kargs, smem, ctx.stream);
                ctx.wait();
                std::vector<device_real_t> Dhost = vector<device_real_t>(Ds[dev].data, Ds[dev].data + Ds[dev].size());
                std::vector<device_real_t> Lhost = vector<device_real_t>(Ls[dev].data, Ls[dev].data + Ls[dev].size());
                safeCudaKernelCall((void*)eigensolve_, qr_dims.get_grid(), qr_dims.get_block(), kargs_qr, smem_qr, ctx.stream);
                safeCudaKernelCall((void*)compute_eigenvector_, transform_dims.get_grid(), transform_dims.get_block(), kargs_transform, smem_transform, ctx.stream);
                

                ctx.wait();
                
                ofstream diag("diag.float32", std::ios::binary); diag.write((char*)Ds[dev].data, Ds[dev].size()*sizeof(device_real_t)); 
                ofstream lower_diag("lower_diag.float32", std::ios::binary); lower_diag.write((char*)Ls[dev].data, Ls[dev].size()*sizeof(device_real_t));
                ofstream Qvects ("Qvects.float32", std::ios::binary); Qvects.write((char*)Qs[dev].data, Qs[dev].size()*sizeof(device_real_t));
                ofstream Evects ("Evects.float32", std::ios::binary); Evects.write((char*)Es[dev].data, Es[dev].size()*sizeof(device_real_t));
                std::vector<device_real_t> Qhost(m*m);
                std::vector<device_real_t> Lambdas(m);
                eigensystem_hermitian(m, Dhost, Lhost, Qhost, Lambdas);
                std::sort(Ds[dev].data, Ds[dev].data + Ds[dev].size());
                std::sort(Lambdas.begin(), Lambdas.end());
                std::vector<device_real_t> error(m); 
                std::transform(Lambdas.begin(), Lambdas.end(), Ds[dev].data, error.begin(), [](device_real_t lambda, device_real_t d) { return (d - lambda) / std::abs(lambda); });
                std::cout << "Error: " << vector(error.data(), error.data() + m) << std::endl;
                std::cout << "Qvects: " << vector(Qs[dev].data, Qs[dev].data + m) << std::endl;
                std::cout << "Qhost: " << vector(Qhost.data(), Qhost.data() + m) << std::endl;
            }
            //std::cout << "QR Time: " << std::chrono::duration_cast<std::chrono::microseconds>(T1-T0).count()/(float)B.isomer_capacity << std::endl;
            
            

            //std::sort(lambdas.begin(), lambdas.end());
            //ofstream eigs("eigs.float32", std::ios::binary); 
            //eigs.write((char*)Ds[dev].data, Ds[dev].size()*sizeof(device_real_t));
            //std::cout << lambdas << std::endl;
            //}

            //safeCudaKernelCall((void*)lanczos_, dims.get_grid(), dims.get_block(), kargs, smem, ctx.stream);  
            //safeCudaKernelCall((void*)eigensolve_, std::min((size_t)ceil(m_isomer_capacity/qr_dims.get_block().x), (size_t)qr_dims.get_grid().x), qr_dims.get_block(), kargs, smem_qr, ctx.stream);

            
            

            if (policy == LaunchPolicy::SYNC) ctx.wait();
            printLastCudaError("Spectrum Ends Failed: ");
        }

        void lambda_max(const IsomerBatch& B, const CuArray<device_real_t>& hessians, const CuArray<device_node_t>& cols, CuArray<device_real_t>& lambda_maxs, int lanczos_steps, const LaunchCtx& ctx, const LaunchPolicy policy){
            if (policy == LaunchPolicy::SYNC) ctx.wait();
            cudaSetDevice(B.get_device_id());
            auto dev = B.get_device_id();
            static int Nd = LaunchCtx::get_device_count();
            static std::vector<bool> init(Nd, false);
            if(policy == LaunchPolicy::SYNC) {ctx.wait();}
            static std::vector<CuArray<device_real_t>> Us(Nd); //Upper diagonals
            static std::vector<CuArray<device_real_t>> Ls(Nd); //Lower diagonals
            static std::vector<CuArray<device_real_t>> Ds(Nd); //Diagonals
            static std::vector<CuArray<device_real_t>> Vs(Nd); //Lanczos vectors
            static std::vector<CuArray<device_real_t>> Qs(Nd); //Q matrix for QR decomposition
            static std::vector<CuArray<device_real_t>> Es(Nd); //6 Eigenvectors to deflate 'H' with
            static int m_natoms = B.n_atoms;
            static int m_isomer_capacity = B.isomer_capacity;
            static int m = lanczos_steps;

            if (!init[dev] || m_natoms != B.n_atoms || m_isomer_capacity != B.isomer_capacity || m != lanczos_steps){
                init[dev] = true;
                m = lanczos_steps;
                m_natoms = B.n_atoms;
                m_isomer_capacity = B.isomer_capacity;
                Us[dev].resize(B.isomer_capacity*m); Us[dev].fill(0.); Us[dev].to_device(dev);
                Ls[dev].resize(B.isomer_capacity*m); Ls[dev].fill(0.); Ls[dev].to_device(dev);
                Ds[dev].resize(B.isomer_capacity*m); Ds[dev].fill(0.); Ds[dev].to_device(dev);
                Vs[dev].resize(B.isomer_capacity*B.n_atoms*3*m); Vs[dev].fill(0.); Vs[dev].to_device(dev);
                Qs[dev].resize(B.isomer_capacity*m*m); Qs[dev].fill(0.); Qs[dev].to_device(dev);
                Es[dev].resize(B.isomer_capacity*B.n_atoms*3*6); Es[dev].fill(0.); Es[dev].to_device(dev); 
            }
            size_t smem = sizeof(device_real_t)*B.n_atoms*3*2 + m*2*sizeof(device_real_t);
            size_t smem_qr = sizeof(device_real_t)*(m+1)*6 + sizeof(device_real_t)*(64)*2;
            static LaunchDims dims((void*)lanczos_, B.n_atoms*3, smem, B.isomer_capacity);
            static LaunchDims qr_dims((void*)eigensolve_largest_, 64, smem_qr, B.isomer_capacity);    
            dims.update_dims((void*)lanczos_, B.n_atoms*3, smem, B.isomer_capacity);
            qr_dims.update_dims((void*)eigensolve_largest_, 64, smem_qr, B.isomer_capacity);
            int n_deflate = 0;
            void* kargs[]{(void*)&B, (void*)&Es[dev], (void*)&Vs[dev], (void*)&Ls[dev], (void*)&Ds[dev], (void*)&hessians, (void*)&cols, (void*)&m, (void*)&n_deflate};
            void* kargs_qr[]{(void*)&B, (void*)&Ds[dev], (void*)&Ls[dev], (void*)&Us[dev], (void*)&Qs[dev], (void*)&lambda_maxs, (void*)&m};
            safeCudaKernelCall((void*)lanczos_, dims.get_grid(), dims.get_block(), kargs, smem, ctx.stream);
            safeCudaKernelCall((void*)eigensolve_largest_, qr_dims.get_grid(), qr_dims.get_block(), kargs_qr, smem_qr, ctx.stream);
            if (policy == LaunchPolicy::SYNC) ctx.wait();
            //std::cout << "Lambda Max: " << vector(lambda_maxs.data, lambda_maxs.data + B.isomer_capacity) << std::endl;
            /* for(int i = 0; i < B.isomer_capacity; i++){
                std::cout << "Diagonals[" << i << vector(Ds[dev].data + i*m, Ds[dev].data + (i+1)*m) << std::endl;
            } */

        }


    }
}