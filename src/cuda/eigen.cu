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
#include <cusolverDn.h>
#include <cusparse.h>

#define N_STREAMS 16

namespace gpu_kernels{
    namespace isomerspace_eigen{
        #include "device_includes.cu"

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
        void __global__ eigensolve_(const IsomerBatch B, CuArray<device_real_t> Qglobal, const CuArray<device_real_t> H, const CuArray<device_node_t> cols, CuArray<device_real_t> eigenvalues){
            DEVICE_TYPEDEFS
            extern __shared__ real_t smem[];
            int N = B.n_atoms * 3; //Number of rows in the hessian
            constexpr int M = 10*3;          //Number of columns in the hessian
            real_t* betas = smem + N;
            real_t* alphas = smem + N*2;
            real_t A[M]; //Hessian matrix, threadIdx.x'th row
            node_t C[M]; //Column indices of the threadIdx.x'th row 3-fold degenerate
            real_t* Q;
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
            //Modified Gram-Schmidt
            auto MGS = [&](int index){
                BLOCK_SYNC
                real_t result = Q[index*N];
                smem[threadIdx.x] = 0;
                #pragma unroll
                for (int j = 0; j < index; j++){
                    auto proj = reduction(smem, result * Q[j*N]) * Q[j*N];
                    result -= proj; //Remove the component along Q[j*N] from result
                }
                result /= sqrt(reduction(smem, result * result));
                return result;
            };
           
            curandState state;            
            curand_init(42 + threadIdx.x, 0, 0, &state);

            for (int I = blockIdx.x; I < B.isomer_capacity; I += gridDim.x){
                Q = Qglobal.data + I * N * N + threadIdx.x;
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
                Q[0*N] = curand_uniform(&state);
                Q[0*N] /= SQRT(reduction(smem, Q[0*N] * Q[0*N]));
                for (int i = 0; i < N; i++){
                    if (i % 2 == 0 && i > 1){
                        Q[(i-1)*N] = MGS(i-1);
                        Q[i*N] = MGS(i);
                        //if(threadIdx.x + blockIdx.x == 0) printf("i = %d, N = %d, Q[i*N] = %f\n", i, N, Q[i*N]);
                    }
                    real_t v = mat_vect(Q[i*N]);
                    alpha = reduction(smem, v * Q[i*N]);
                    if (threadIdx.x == i) alphas[i] = alpha;
                    if (i == 0){
                        v -= alpha * Q[i*N];
                    } else {
                        v -= betas[i-1] * Q[(i-1)*N] + alpha * Q[i*N];
                    }
                    beta = SQRT(reduction(smem, v * v));
                    if (threadIdx.x == i) betas[i] = beta;
                    if (i < N-1) Q[(i+1)*N] = v / beta;
                    //if (i < N-1) Q[(i+1)*N] = beta;
                }
                eigenvalues.data[I*N*2 + threadIdx.x] = ISNAN(alphas[threadIdx.x])  ? real_t(0) : alphas[threadIdx.x];
                eigenvalues.data[I*N*2 + N + threadIdx.x] = ISNAN(betas[threadIdx.x]) ? real_t(0) : betas[threadIdx.x];
            }   
        }

        void eigensolve(const IsomerBatch& B, CuArray<device_real_t>& Q, const CuArray<device_real_t>& hessians, const CuArray<device_node_t>& cols, CuArray<device_real_t>& eigenvalues, const LaunchCtx& ctx, const LaunchPolicy policy){
            if (policy == LaunchPolicy::SYNC) ctx.wait();
            cudaSetDevice(B.get_device_id());
            auto dev = B.get_device_id();

            if(policy == LaunchPolicy::SYNC) {ctx.wait();}

            size_t smem = sizeof(device_coord3d)*B.n_atoms*3 + sizeof(device_real_t)*Block_Size_Pow_2;
            static LaunchDims dims((void*)eigensolve_, B.n_atoms*3, smem, B.isomer_capacity);
            dims.update_dims((void*)eigensolve_, B.n_atoms*3, smem, B.isomer_capacity);
            void* kargs[]{(void*)&B, (void*)&Q, (void*)&hessians, (void*)&cols, (void*)&eigenvalues};
            safeCudaKernelCall((void*)eigensolve_, dims.get_grid(), dims.get_block(), kargs, smem, ctx.stream);  
            
            
            printLastCudaError("Cubic Layout: ");

            if (policy == LaunchPolicy::SYNC) ctx.wait();
        }


    }
}