#ifndef ISOMERSPACE_KERNELS_H
#define ISOMERSPACE_KERNELS_H
#include "fullerenes/gpu/cuda_definitions.h"
#include "fullerenes/gpu/launch_ctx.hh"
#include "fullerenes/gpu/isomer_batch.hh"
#include "launch_ctx.hh"
#include "cu_array.hh"
#include <chrono>

enum ForcefieldType {WIRZ, PEDERSEN, FLATNESS_ENABLED, FLAT_BOND, BOND, ANGLE, DIH, ANGLE_M, ANGLE_P, DIH_A, DIH_M, DIH_P};

namespace gpu_kernels{
    namespace isomerspace_forcefield{
        int optimal_batch_size(const int N, const int device_id = 0);

        //Returns total time spent in kernel in milliseconds. First call is not timed as that would require synchronous execution.
        std::chrono::microseconds time_spent(); 

        //Resets the time spent in kernel to 0
        void reset_time();

        //Uses forcefield optimization to relax the positions of atoms in all isomers of a batch.
        template <ForcefieldType T, BufferType U>
        cudaError_t optimise(IsomerBatch<U>& B, const size_t iterations, const size_t max_iterations, const LaunchCtx& ctx = LaunchCtx(), const LaunchPolicy policy = LaunchPolicy::SYNC);

        template <BufferType U>
        cudaError_t test_fun(IsomerBatch<U>& B, CuArray<float>& output);
    
        //ForcefieldType agnostic functions since the forcefield type doesnt change the way bonds are computed
        //Retrieves all interatomic bond lengths, angles and dihedrals
        template <BufferType U> cudaError_t get_bonds           (const IsomerBatch<U>& B, CuArray<float>& bonds);          //N x M x 3
        template <BufferType U> cudaError_t get_angles          (const IsomerBatch<U>& B, CuArray<float>& angles);         //N x M x 3
        template <BufferType U> cudaError_t get_dihedrals       (const IsomerBatch<U>& B, CuArray<float>& dihedrals);      //N x M x 3


        //MAE = Mean Absolute Error
        //RMSE = Root Mean Squared Error
        //RRMSE = Relative Root Mean Squared Error

        //Stat functions:
        template <ForcefieldType T, BufferType U> cudaError_t get_dihedral_mean   (const IsomerBatch<U>& B, CuArray<float>& dihedral_mean);  //M x 1
        template <ForcefieldType T, BufferType U> cudaError_t get_dihedral_mae    (const IsomerBatch<U>& B, CuArray<float>& dihedral_mae);   //M x 1 
        template <ForcefieldType T, BufferType U> cudaError_t get_dihedral_rmse   (const IsomerBatch<U>& B, CuArray<float>& dihedral_rmse);   //M x 1
        template <ForcefieldType T, BufferType U> cudaError_t get_dihedral_rrmse  (const IsomerBatch<U>& B, CuArray<float>& dihedral_rrmse);  //M x 1
        template <ForcefieldType T, BufferType U> cudaError_t get_dihedral_max    (const IsomerBatch<U>& B, CuArray<float>& dihedral_max);   //M x 1
        template <ForcefieldType T, BufferType U> cudaError_t get_energies        (const IsomerBatch<U>& B, CuArray<float>& energies);       //M x 1

        template <ForcefieldType T, BufferType U> cudaError_t get_energy          (const IsomerBatch<U>& B, CuArray<float>& energy);         //M x 1
        template <ForcefieldType T, BufferType U> cudaError_t get_bond_mean       (const IsomerBatch<U>& B, CuArray<float>& bond_mean);      //M x 1
        template <ForcefieldType T, BufferType U> cudaError_t get_bond_mae        (const IsomerBatch<U>& B, CuArray<float>& bond_mae);       //M x 1
        template <ForcefieldType T, BufferType U> cudaError_t get_bond_rmse       (const IsomerBatch<U>& B, CuArray<float>& bond_rmse);       //M x 1
        template <ForcefieldType T, BufferType U> cudaError_t get_bond_rrmse      (const IsomerBatch<U>& B, CuArray<float>& bond_rrmse);      //M x 1
        template <ForcefieldType T, BufferType U> cudaError_t get_bond_max        (const IsomerBatch<U>& B, CuArray<float>& bond_max);       //M x 1
        template <ForcefieldType T, BufferType U> cudaError_t get_angle_mean      (const IsomerBatch<U>& B, CuArray<float>& angle_mean);     //M x 1
        template <ForcefieldType T, BufferType U> cudaError_t get_angle_mae       (const IsomerBatch<U>& B, CuArray<float>& angle_mae);      //M x 1
        template <ForcefieldType T, BufferType U> cudaError_t get_angle_rmse      (const IsomerBatch<U>& B, CuArray<float>& angle_rms);      //M x 1
        template <ForcefieldType T, BufferType U> cudaError_t get_angle_rrmse     (const IsomerBatch<U>& B, CuArray<float>& angle_rrmse);     //M x 1
        template <ForcefieldType T, BufferType U> cudaError_t get_angle_max       (const IsomerBatch<U>& B, CuArray<float>& angle_max);      //M x 1
        
        template <ForcefieldType T, BufferType U> cudaError_t get_flat_mean       (const IsomerBatch<U>& B, CuArray<float>& flat_mean);      //M x 1
        template <ForcefieldType T, BufferType U> cudaError_t get_flat_max        (const IsomerBatch<U>& B, CuArray<float>& flat_max);      //M x 1
        template <ForcefieldType T, BufferType U> cudaError_t get_flat_rmse       (const IsomerBatch<U>& B, CuArray<float>& flat_rms);      //M x 1
        
        
        template <ForcefieldType T, BufferType U> cudaError_t get_gradient_mean   (const IsomerBatch<U>& B, CuArray<float>& gradient_mean);  //M x 1
        template <ForcefieldType T, BufferType U> cudaError_t get_gradient_norm   (const IsomerBatch<U>& B, CuArray<float>& gradient_norm);  //M x 1
        template <ForcefieldType T, BufferType U> cudaError_t get_gradient_rms    (const IsomerBatch<U>& B, CuArray<float>& gradient_rms);   //M x 1
        template <ForcefieldType T, BufferType U> cudaError_t get_gradient_max    (const IsomerBatch<U>& B, CuArray<float>& gradient_max);   //M x 1

    }

    namespace isomerspace_X0{
        //Returns total time spent in kernel in milliseconds. First call is not timed as that would require synchronous execution.
        std::chrono::microseconds time_spent(); 

        void reset_time(); //Resets time spent in kernel to 0.

        //Uses spherical projection and 'wrapping a sphere' technique to generate starting coordinates in 3D space for a batch of fullerenes.
        template <BufferType U> 
        cudaError_t zero_order_geometry(IsomerBatch<U>& B, const float scalerad, const LaunchCtx& ctx = LaunchCtx(), const LaunchPolicy policy = LaunchPolicy::SYNC);
    }

    namespace isomerspace_tutte{
        //Returns total time spent in kernel in milliseconds. First call is not timed as that would require synchronous execution.
        std::chrono::microseconds time_spent(); 

        void reset_time(); //Resets time spent in kernel to 0.
        //Computes tutte embedding for an entire batch of fullernes.
        template <BufferType U>
        cudaError_t tutte_layout(IsomerBatch<U>& B, const size_t max_iterations = 10000000, const LaunchCtx& ctx = LaunchCtx(), const LaunchPolicy policy = LaunchPolicy::SYNC);
    }

    namespace isomerspace_dual{
        //Returns the optimal batch size for a given number of atoms, for the dualisation kernel.   
        int optimal_batch_size(const int N, const int device_id = 0); 

        int optimal_batch_size_2(const int N, const int device_id = 0);
        
        //Returns total time spent in kernel in milliseconds. First call is not timed as that would require synchronous execution.
        std::chrono::microseconds time_spent(); 

        void reset_time(); //Resets time spent in kernel to 0.
        //Computes the cubic neighbour list from the dual neighbour list
        template <BufferType U>
        cudaError_t dualise(IsomerBatch<U>& B, const LaunchCtx& ctx = LaunchCtx(), const LaunchPolicy policy = LaunchPolicy::SYNC);

        template <BufferType U>
        cudaError_t dualise_2(IsomerBatch<U>& B, const LaunchCtx& ctx = LaunchCtx(), const LaunchPolicy policy = LaunchPolicy::SYNC);

        template <BufferType U>
        void dualise_3(IsomerBatch<U>& B);

        template <BufferType U>
        void dualise_4(IsomerBatch<U>& B);
    }

    namespace isomerspace_eigen{
        template <BufferType U>
        void eigensolve_cusolver(const IsomerBatch<U>& B, const CuArray<device_real_t>& hessians, const CuArray<device_node_t>& cols, CuArray<device_real_t>& eigenvalues, const LaunchCtx& ctx = LaunchCtx(), const LaunchPolicy policy = LaunchPolicy::SYNC);

        template <BufferType U>
        void eigensolve(const IsomerBatch<U>& B, CuArray<device_real_t>& Q, const CuArray<device_real_t>& hessians, const CuArray<device_node_t>& cols, CuArray<device_real_t>& eigenvalues, const LaunchCtx& ctx = LaunchCtx(), const LaunchPolicy policy = LaunchPolicy::SYNC);

        template <BufferType U>
        void spectrum_ends(const IsomerBatch<U>& B, const CuArray<device_real_t>& hessians, const CuArray<device_node_t>& cols, CuArray<device_real_t>& lambda_mins, CuArray<device_real_t>& lambda_maxs, int m_lanczos=40, const LaunchCtx& ctx = LaunchCtx(), const LaunchPolicy policy = LaunchPolicy::SYNC);
        template <BufferType U>
        void spectrum_ends(const IsomerBatch<U>& B, const CuArray<device_real_t>& hessians, const CuArray<device_node_t>& cols, CuArray<device_real_t>& lambda_mins, CuArray<device_real_t>& lambda_maxs, CuArray<device_real_t>& eigvect_mins, CuArray<device_real_t>& eigvect_maxs, int m_lanczos=40, const LaunchCtx& ctx = LaunchCtx(), const LaunchPolicy policy = LaunchPolicy::SYNC);
    }

    namespace isomerspace_hessian{
        template <ForcefieldType T, BufferType U>
        cudaError_t compute_hessians(const IsomerBatch<U>& B, CuArray<device_real_t>& hessians, CuArray<device_node_t>& cols, const LaunchCtx& ctx = LaunchCtx(), const LaunchPolicy policy = LaunchPolicy::SYNC);
        template <ForcefieldType T, BufferType U>
        cudaError_t compute_hessians_fd(const IsomerBatch<U>& B, CuArray<device_real_t>& hessians, CuArray<device_node_t>& cols, const float rel_delta = 1e-5, const LaunchCtx& ctx = LaunchCtx(), const LaunchPolicy policy = LaunchPolicy::SYNC);
    }

    namespace isomerspace_properties{

        //Computes the volume of each molecule in the batch. See sequential code in polyhedron.cc
        template <BufferType U>
        cudaError_t volume_divergences(const IsomerBatch<U>& B, CuArray<device_real_t>& volumes, const LaunchCtx& ctx = LaunchCtx(), const LaunchPolicy policy = LaunchPolicy::SYNC);

        //Computes the surface area of each molecule in the batch. See sequential code in polyhedron.cc
        template <BufferType U>
        cudaError_t surface_areas(const IsomerBatch<U>& B, CuArray<device_real_t>& surface_areas, const LaunchCtx& ctx = LaunchCtx(), const LaunchPolicy policy = LaunchPolicy::SYNC);

      //Computes the moments of inertia (eigenvalues to the moment of inertia matrix)
        template <BufferType U>
        cudaError_t moments_of_inertia(const IsomerBatch<U>& B, CuArray<device_real_t>& lambdas, const LaunchCtx& ctx = LaunchCtx(), const LaunchPolicy policy = LaunchPolicy::SYNC);      
      //Computes the best elipsoid fit for each molecule in the batch and returns the eccentricities of the ellipsoids.
        template <BufferType U>
        cudaError_t eccentricities(const IsomerBatch<U>& B, CuArray<device_real_t>& eccentricities, const LaunchCtx& ctx = LaunchCtx(), const LaunchPolicy policy = LaunchPolicy::SYNC);

        //Moves coordinates to the centre of mass of the molecule and aligns the molecule to the principal axes of inertia. See sequential code in polyhedron.cc
        template <BufferType U>
        cudaError_t transform_coordinates(IsomerBatch<U>& B, const LaunchCtx& ctx = LaunchCtx(), const LaunchPolicy policy = LaunchPolicy::SYNC);

        //cudaError_t debug_function(const IsomerBatch<U>& B, CuArray<double>& eigenvalues, CuArray<double>& eigenvectors, CuArray<double>& inertia_matrices, CuArray<double>& orthogonality, const LaunchCtx& ctx = LaunchCtx(), const LaunchPolicy policy = LaunchPolicy::SYNC);
        template <BufferType U>
        cudaError_t debug_function(const IsomerBatch<U>& B, CuArray<device_real_t>& eigenvalues, CuArray<device_real_t>& eigenvectors, CuArray<device_real_t>& inertia_matrices, CuArray<device_real_t>& orthogonality, const LaunchCtx& ctx = LaunchCtx(), const LaunchPolicy policy = LaunchPolicy::SYNC);
    }
}
#endif
