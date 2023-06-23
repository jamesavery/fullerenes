#include <cuda.h>
#include <cooperative_groups.h>
#include <cooperative_groups/reduce.h>
#include <cooperative_groups/scan.h>
#include "cuda_runtime.h"
#include <cuda_runtime_api.h>
//#define FLOAT_TYPE 2
#include "fullerenes/gpu/cuda_definitions.h"
#include "fullerenes/gpu/cu_array.hh"
#include "fullerenes/gpu/isomer_batch.hh"
#include "fullerenes/gpu/kernels.hh"

namespace gpu_kernels{
    namespace isomerspace_properties{
        #include "device_includes.cu"
        symMat3 __device__ inertia_matrix(const device_coord3d* X){
            DEVICE_TYPEDEFS;
            extern __shared__ real_t smem[];
            clear_cache(smem, blockDim.x);
            int tid = threadIdx.x;
            symMat3 I;
            real_t diag = reduction(smem, dot(X[tid], X[tid]));
            I.a = diag;
            I.d = diag;
            I.f = diag;
            I.a -= reduction(smem, X[tid][0]*X[tid][0]);
            I.b -= reduction(smem, X[tid][0]*X[tid][1]);
            I.c -= reduction(smem, X[tid][0]*X[tid][2]);
            I.d -= reduction(smem, X[tid][1]*X[tid][1]);
            I.e -= reduction(smem, X[tid][1]*X[tid][2]);
            I.f -= reduction(smem, X[tid][2]*X[tid][2]);
            return I;
        }

     
      std::array<device_coord3d,3> __device__ principal_axes(const device_coord3d* X){
            DEVICE_TYPEDEFS;
            auto I = inertia_matrix(X);
	    auto [V,lambdas] = I.eigensystem();
	    return V;
        }

        //Returns the best ellipsoid for the given coordinates, lambda0 = a, lambda1 = b, lambda2 = c.
        device_coord3d __device__ best_ellipsoid (const device_coord3d* X){
            DEVICE_TYPEDEFS;
            auto I = inertia_matrix(X);
            return rsqrt3(d_sort(d_abs(I.eigenvalues()))); 
        }

        void __global__ transform_coordinates_(IsomerBatch B){
            DEVICE_TYPEDEFS;
            extern __shared__ real_t shared_memory[];
            const int tid = threadIdx.x;
            auto limit = ((B.isomer_capacity + gridDim.x - 1) / gridDim.x ) * gridDim.x;  //Fast ceiling integer division.
            for (int isomer_idx = blockIdx.x; isomer_idx < limit; isomer_idx += gridDim.x){
                if (isomer_idx < B.isomer_capacity) if(B.statuses[isomer_idx] != IsomerStatus::EMPTY)
                {
                clear_cache(shared_memory, blockDim.x);
                BLOCK_SYNC
                size_t offset = isomer_idx * blockDim.x;
                coord3d* X              = reinterpret_cast<coord3d*>(shared_memory) + B.n_atoms;
                assign(X[threadIdx.x],reinterpret_cast<std::array<float,3>*>(B.X+offset*3)[threadIdx.x]);
                BLOCK_SYNC
                coord3d centroid = {reduction(shared_memory, X[tid][0]), reduction(shared_memory, X[tid][1]), reduction(shared_memory, X[tid][2])};
                X[tid] -= centroid/real_t(B.n_atoms);
                BLOCK_SYNC
		            mat3 P{principal_axes(X)};
                if (ISNAN(P[0][0]) || ISNAN(P[0][1]) || ISNAN(P[0][2]) || ISNAN(P[1][0]) || ISNAN(P[1][1]) || ISNAN(P[1][2]) || ISNAN(P[2][0]) || ISNAN(P[2][1]) || ISNAN(P[2][2])) {
                    //assert(false);
                    return;
                } 
                BLOCK_SYNC
                X[tid] = dot(P, X[tid]);
                BLOCK_SYNC
                assign(reinterpret_cast<std::array<float,3>*>(B.X)[offset + threadIdx.x], X[threadIdx.x]);
                //sequential_print(P[0][0],0);
                }
            }
            
        }

        void __global__ moments_of_inertia_(IsomerBatch B, CuArray<device_real_t> lambdas){
            DEVICE_TYPEDEFS;
            extern __shared__ real_t shared_memory[];
            const int tid = threadIdx.x;
            auto limit = ((B.isomer_capacity + gridDim.x - 1) / gridDim.x ) * gridDim.x;  //Fast ceiling integer division.
            for (int isomer_idx = blockIdx.x; isomer_idx < limit; isomer_idx += gridDim.x){
	      //TODO: simplify
                if (isomer_idx < B.isomer_capacity) if(B.statuses[isomer_idx] == IsomerStatus::CONVERGED)
                {
		BLOCK_SYNC;
                size_t offset = isomer_idx * blockDim.x;
                coord3d* X           = reinterpret_cast<coord3d*>(shared_memory) + B.n_atoms;
                assign(X[threadIdx.x],reinterpret_cast<std::array<float,3>*>(B.X+offset*3)[threadIdx.x]);
                BLOCK_SYNC;
		auto I = inertia_matrix(X);
		device_coord3d lams = I.eigenvalues(); 
                BLOCK_SYNC;
                if (tid == 0){ lambdas.data[3*isomer_idx] = lams[0];lambdas.data[3*isomer_idx+1] = lams[0];lambdas.data[3*isomer_idx+2] = lams[2]; }
                }
            }
        }
        
      
      
        void __global__ eccentricities_(IsomerBatch B, CuArray<device_real_t> ecce){
            DEVICE_TYPEDEFS;
            extern __shared__ real_t shared_memory[];
            const int tid = threadIdx.x;
            auto limit = ((B.isomer_capacity + gridDim.x - 1) / gridDim.x ) * gridDim.x;  //Fast ceiling integer division.
            for (int isomer_idx = blockIdx.x; isomer_idx < limit; isomer_idx += gridDim.x){
                if (isomer_idx < B.isomer_capacity) if(B.statuses[isomer_idx] == IsomerStatus::CONVERGED)
                {
                BLOCK_SYNC
                size_t offset = isomer_idx * blockDim.x;
                coord3d* X            = reinterpret_cast<coord3d*>(shared_memory) + B.n_atoms;
                assign(X[threadIdx.x],reinterpret_cast<std::array<float,3>*>(B.X+offset*3)[threadIdx.x]);
                BLOCK_SYNC
                auto result = best_ellipsoid(X);
                BLOCK_SYNC
                if (tid == 0) ecce.data[isomer_idx] = result[0] / result[2];
                }
            }
        }
        
        void __global__ volume_divergences_(IsomerBatch B, CuArray<device_real_t> vd){
            DEVICE_TYPEDEFS;
            typedef device_node3 tri_t;
            extern __shared__ real_t smems[];
            const int tid = threadIdx.x;
            auto limit = ((B.isomer_capacity + gridDim.x - 1) / gridDim.x ) * gridDim.x;  //Fast ceiling integer division.
            for (int isomer_idx = blockIdx.x; isomer_idx < limit; isomer_idx += gridDim.x){
                if (isomer_idx < B.isomer_capacity) if(B.statuses[isomer_idx] == IsomerStatus::CONVERGED)
                {
                BLOCK_SYNC
                size_t offset = isomer_idx * blockDim.x;
                NodeNeighbours node_graph    = NodeNeighbours(B, isomer_idx, smems);
                coord3d* X              = reinterpret_cast<coord3d*>(smems) + B.n_atoms;
                assign(X[threadIdx.x],reinterpret_cast<std::array<float,3>*>(B.X+offset*3)[threadIdx.x]);
                BLOCK_SYNC
                real_t V = 0.;
                if (tid < B.n_faces) {
                    coord3d face_center = (coord3d){0.,0.,0.};
                    for (int i = 0; i < node_graph.face_size; i++) face_center += X[node_graph.face_nodes[i]];
                    face_center /= node_graph.face_size; //The center of the threadIdx.x-th face.

                    for (int i = 0; i < node_graph.face_size; i++){
                        coord3d a = X[node_graph.face_nodes[i]];
                        coord3d b = X[node_graph.face_nodes[(i+1)%node_graph.face_size]];
                        coord3d c = face_center;
                        coord3d u = b - a;
                        coord3d v = c - a;
                        coord3d n = cross(u,v);
                        V += dot(a,n) / real_t(2.0);
                    }
                }
                clear_cache(smems, blockDim.x);
                auto result = reduction(smems, V)/real_t(3.0);
                if (tid == 0) vd.data[isomer_idx] = result;
                }
            }
        }

        void __global__ surface_areas_(IsomerBatch B, CuArray<device_real_t> sa){
            DEVICE_TYPEDEFS;
            typedef device_node3 tri_t;
            extern __shared__ real_t smems[];
            const int tid = threadIdx.x;
            auto limit = ((B.isomer_capacity + gridDim.x - 1) / gridDim.x ) * gridDim.x;  //Fast ceiling integer division.
            for (int isomer_idx = blockIdx.x; isomer_idx < limit; isomer_idx += gridDim.x){
                if (isomer_idx < B.isomer_capacity) if(B.statuses[isomer_idx] != IsomerStatus::EMPTY)
                {
                BLOCK_SYNC
                size_t offset = isomer_idx * blockDim.x;
                NodeNeighbours node_graph    = NodeNeighbours(B, isomer_idx, smems);
                coord3d* X              = reinterpret_cast<coord3d*>(smems) + B.n_atoms;
                assign(X[threadIdx.x],reinterpret_cast<std::array<float,3>*>(B.X+offset*3)[threadIdx.x]);
                BLOCK_SYNC
                real_t A = real_t(0.);
                if (tid < B.n_faces) {
                    coord3d face_center = (coord3d){0.,0.,0.};
                    for (int i = 0; i < node_graph.face_size; i++) face_center += X[node_graph.face_nodes[i]];
                    face_center /= node_graph.face_size; //The center of the threadIdx.x-th face.

                    for (int i = 0; i < node_graph.face_size; i++){
                        coord3d a = X[node_graph.face_nodes[i]];
                        coord3d b = X[node_graph.face_nodes[(i+1)%node_graph.face_size]];
                        coord3d c = face_center;
                        coord3d u = b - a;
                        coord3d v = c - a;
                        coord3d n = cross(u,v);
                        A += norm(n);
                    }
                }
                clear_cache(smems, blockDim.x);
                auto result = reduction(smems, A)/real_t(2.0);
                if (tid == 0) sa.data[isomer_idx] = result;
                }
            }
        }

        void __global__ debug_function_(IsomerBatch B, CuArray<device_real_t> eigenvalues, CuArray<device_real_t> eigenvectors, CuArray<device_real_t> inertia_matrices, CuArray<device_real_t> orthogonality){
            DEVICE_TYPEDEFS;
            extern __shared__ real_t shared_memory[];
            const int tid = threadIdx.x;
            auto limit = ((B.isomer_capacity + gridDim.x - 1) / gridDim.x ) * gridDim.x;  //Fast ceiling integer division.
            for (int isomer_idx = blockIdx.x; isomer_idx < limit; isomer_idx += gridDim.x){
                if (isomer_idx < B.isomer_capacity)
                if (B.statuses[isomer_idx] == IsomerStatus::CONVERGED)
                {
                BLOCK_SYNC
                size_t offset = isomer_idx * blockDim.x;
                Constants constants          = Constants(B, isomer_idx);
                NodeNeighbours node_graph    = NodeNeighbours(B, isomer_idx, shared_memory);
                coord3d* X              = reinterpret_cast<coord3d*>(shared_memory) + B.n_atoms;
                assign(X[threadIdx.x],reinterpret_cast<std::array<float,3>*>(B.X+offset*3)[threadIdx.x]);
                BLOCK_SYNC
                auto I = inertia_matrix(X);
                BLOCK_SYNC
                if (tid == 0){
		  auto [P,eigs] = I.eigensystem();		  

		  reinterpret_cast<coord3d*>(eigenvalues.data) [isomer_idx] = eigs;
		  reinterpret_cast<coord3d*>(eigenvectors.data)[isomer_idx*3+0] = P[0];
		  reinterpret_cast<coord3d*>(eigenvectors.data)[isomer_idx*3+1] = P[1];
		  reinterpret_cast<coord3d*>(eigenvectors.data)[isomer_idx*3+2] = P[2];
		
		  real_t orthog = real_t(1.0);

		  if(ISNAN(eigs[0]) || ISNAN(eigs[1]) || ISNAN(eigs[2]) || ISNAN(I.a) || ISNAN(I.b) || ISNAN(I.c) || ISNAN(I.d) || ISNAN(I.e) || ISNAN(I.f)){
                    orthog = real_t(2.0);
		  } else if ((ABS(eigs[0]- eigs[1])/ABS(eigs[0]) < 1e-5) && (ABS(eigs[1]- eigs[2])/ABS(eigs[0]) < 1e-5)){
                    orthog = real_t(3.0);
		  } else if (ABS(eigs[0] - eigs[1])/ABS(eigs[0]) < 1e-5) {
                    orthog = ABS(dot(P[0], P[2]));
		  } else {
                    orthog = ABS(dot(P[0], P[1]));
		  }

		  orthogonality.data[isomer_idx] = orthog;
		  reinterpret_cast<coord3d*>(inertia_matrices.data)[isomer_idx*3] = {I.a, I.b, I.c};
		  reinterpret_cast<coord3d*>(inertia_matrices.data)[isomer_idx*3+1] = {I.b, I.d, I.e};
		  reinterpret_cast<coord3d*>(inertia_matrices.data)[isomer_idx*3+2] = {I.c, I.e, I.f};
                }
                }
            }
        }



        cudaError_t transform_coordinates(IsomerBatch& B, const LaunchCtx& ctx, const LaunchPolicy policy){
            cudaSetDevice(B.get_device_id());
                
            //If launch ploicy is synchronous then wait.
            if(policy == LaunchPolicy::SYNC) ctx.wait();

            size_t smem = sizeof(device_coord3d)* (3*B.n_atoms + 4) + sizeof(device_real_t)*Block_Size_Pow_2;
            static LaunchDims dims((void*)transform_coordinates_, B.n_atoms, smem, B.isomer_capacity);
            dims.update_dims((void*)transform_coordinates_, B.n_atoms, smem, B.isomer_capacity);
            void* kargs[]{(void*)&B};
            auto error = safeCudaKernelCall((void*)transform_coordinates_, dims.get_grid(), dims.get_block(), kargs, smem, ctx.stream);
    
            if(policy == LaunchPolicy::SYNC) ctx.wait();
            printLastCudaError("Transformation of Coordinates Failed: ");
            return error;
        }


        cudaError_t eccentricities(const IsomerBatch& B, CuArray<device_real_t>& eccentricities, const LaunchCtx& ctx, const LaunchPolicy policy){
            cudaSetDevice(B.get_device_id());

            //If launch ploicy is synchronous then wait.
            if(policy == LaunchPolicy::SYNC) ctx.wait();

            size_t smem = sizeof(device_coord3d)* (3*B.n_atoms + 4) + sizeof(device_real_t)*Block_Size_Pow_2;
            static LaunchDims dims((void*)eccentricities_, B.n_atoms, smem, B.isomer_capacity);
            dims.update_dims((void*)eccentricities_, B.n_atoms, smem, B.isomer_capacity);
            void* kargs[]{(void*)&B, (void*)&eccentricities};
            auto error = safeCudaKernelCall((void*)eccentricities_, dims.get_grid(), dims.get_block(), kargs, smem, ctx.stream);

            if(policy == LaunchPolicy::SYNC) ctx.wait();
            printLastCudaError("Calculation of Eccentricities Failed: ");
            return error;
        }

      cudaError_t moments_of_inertia(const IsomerBatch& B, CuArray<device_real_t>& lambdas, const LaunchCtx& ctx, const LaunchPolicy policy){
            cudaSetDevice(B.get_device_id());

            //If launch ploicy is synchronous then wait.
            if(policy == LaunchPolicy::SYNC) ctx.wait();

            size_t smem = sizeof(device_coord3d)* (3*B.n_atoms + 4) + sizeof(device_real_t)*Block_Size_Pow_2;
            static LaunchDims dims((void*)moments_of_inertia_, B.n_atoms, smem, B.isomer_capacity);
            dims.update_dims((void*)moments_of_inertia_, B.n_atoms, smem, B.isomer_capacity);
            void* kargs[]{(void*)&B, (void*)&lambdas};
            auto error = safeCudaKernelCall((void*)moments_of_inertia_, dims.get_grid(), dims.get_block(), kargs, smem, ctx.stream);

            if(policy == LaunchPolicy::SYNC) ctx.wait();
            printLastCudaError("Calculation of moments of inertia Failed: ");
            return error;
        }

        cudaError_t surface_areas(const IsomerBatch& B, CuArray<device_real_t>& surface_areas, const LaunchCtx& ctx, const LaunchPolicy policy){
            cudaSetDevice(B.get_device_id());

            //If launch ploicy is synchronous then wait.
            if(policy == LaunchPolicy::SYNC) ctx.wait();

            size_t smem = sizeof(device_coord3d)* (3*B.n_atoms + 4) + sizeof(device_real_t)*Block_Size_Pow_2;
            static LaunchDims dims((void*)surface_areas_, B.n_atoms, smem, B.isomer_capacity);
            dims.update_dims((void*)surface_areas_, B.n_atoms, smem, B.isomer_capacity);
            void* kargs[]{(void*)&B, (void*)&surface_areas};
            auto error = safeCudaKernelCall((void*)surface_areas_, dims.get_grid(), dims.get_block(), kargs, smem, ctx.stream);

            if(policy == LaunchPolicy::SYNC) ctx.wait();
            printLastCudaError("Calculation of Volume Divergences Failed: ");
            return error;
        }

        cudaError_t volume_divergences(const IsomerBatch& B, CuArray<device_real_t>& volumes, const LaunchCtx& ctx, const LaunchPolicy policy){
            cudaSetDevice(B.get_device_id());

            //If launch ploicy is synchronous then wait.
            if(policy == LaunchPolicy::SYNC) ctx.wait();

            size_t smem = sizeof(device_coord3d)* (3*B.n_atoms + 4) + sizeof(device_real_t)*Block_Size_Pow_2;
            static LaunchDims dims((void*)volume_divergences_, B.n_atoms, smem, B.isomer_capacity);
            dims.update_dims((void*)volume_divergences_, B.n_atoms, smem, B.isomer_capacity);
            void* kargs[]{(void*)&B, (void*)&volumes};
            auto error = safeCudaKernelCall((void*)volume_divergences_, dims.get_grid(), dims.get_block(), kargs, smem, ctx.stream);

            if(policy == LaunchPolicy::SYNC) ctx.wait();
            printLastCudaError("Calculation of Volume Divergences Failed: ");
            return error;
        }

        cudaError_t debug_function(const IsomerBatch& B, CuArray<device_real_t>& eigenvalues, CuArray<device_real_t>& eigenvectors, CuArray<device_real_t>& inertia_matrices, CuArray<device_real_t>& orthogonality, const LaunchCtx& ctx, const LaunchPolicy policy){
            cudaSetDevice(B.get_device_id());

            //If launch ploicy is synchronous then wait.
            if(policy == LaunchPolicy::SYNC) ctx.wait();

            size_t smem = sizeof(device_coord3d)* (3*B.n_atoms + 4) + sizeof(device_real_t)*Block_Size_Pow_2;
            static LaunchDims dims((void*)debug_function_, B.n_atoms, smem, B.isomer_capacity);
            dims.update_dims((void*)debug_function_, B.n_atoms, smem, B.isomer_capacity);
            void* kargs[]{(void*)&B, (void*)&eigenvalues, (void*)&eigenvectors, (void*)&inertia_matrices, (void*)&orthogonality};
            auto error = safeCudaKernelCall((void*)debug_function_, dims.get_grid(), dims.get_block(), kargs, smem, ctx.stream);

            if(policy == LaunchPolicy::SYNC) ctx.wait();
            printLastCudaError("Debug Function Failed: ");
            return error;
        }





    }
}
