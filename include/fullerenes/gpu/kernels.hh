#ifndef ISOMERSPACE_KERNELS_H
#define ISOMERSPACE_KERNELS_H
#include "fullerenes/gpu/cuda_definitions.h"
#include "fullerenes/gpu/launch_ctx.hh"
#include "fullerenes/isomer_batch.hh"
#include "launch_ctx.hh"
#include "cu_array.hh"
#include <chrono>

enum ForcefieldType {WIRZ, PEDERSEN, FLATNESS_ENABLED, FLAT_BOND, BOND, ANGLE, DIH, ANGLE_M, ANGLE_P, DIH_A, DIH_M, DIH_P};

namespace gpu_kernels{
    namespace isomerspace_forcefield{
        template <ForcefieldType FFT = PEDERSEN, Device U = GPU, typename T = float, typename K = uint16_t>
        int optimal_batch_size(const int N, const int device_id = 0);

        //Returns total time spent in kernel in milliseconds. First call is not timed as that would require synchronous execution.
        std::chrono::microseconds time_spent(); 

        //Resets the time spent in kernel to 0
        void reset_time();

        //Uses forcefield optimization to relax the positions of atoms in all isomers of a batch.
        
        template <ForcefieldType FFT = PEDERSEN, Device U = GPU, typename T = float, typename K = uint16_t>
        cudaError_t optimise(IsomerBatch<U>& B, const size_t iterations, const size_t max_iterations, const LaunchCtx& ctx = LaunchCtx(), const LaunchPolicy policy = LaunchPolicy::SYNC);

        template <Device U>
        cudaError_t test_fun(IsomerBatch<U>& B, CuArray<float>& output);
    
        //ForcefieldType agnostic functions since the forcefield type doesnt change the way bonds are computed
        //Retrieves all interatomic bond lengths, angles and dihedrals
        template <Device U, typename T = float, typename K = uint16_t > cudaError_t get_bonds           (const IsomerBatch<U>& B, CuArray<T>& bonds);          //N x M x 3
        template <Device U, typename T = float, typename K = uint16_t > cudaError_t get_angles          (const IsomerBatch<U>& B, CuArray<T>& angles);         //N x M x 3
        template <Device U, typename T = float, typename K = uint16_t > cudaError_t get_dihedrals       (const IsomerBatch<U>& B, CuArray<T>& dihedrals);      //N x M x 3


        //MAE = Mean Absolute Error
        //RMSE = Root Mean Squared Error
        //RRMSE = Relative Root Mean Squared Error

        //Stat functions:
        template <ForcefieldType FFT = PEDERSEN, Device U = GPU, typename T = float, typename K = uint16_t> cudaError_t get_dihedral_mean   (const IsomerBatch<U>& B, CuArray<T>& dihedral_mean);  //M x 1
        template <ForcefieldType FFT = PEDERSEN, Device U = GPU, typename T = float, typename K = uint16_t> cudaError_t get_dihedral_mae    (const IsomerBatch<U>& B, CuArray<T>& dihedral_mae);   //M x 1 
        template <ForcefieldType FFT = PEDERSEN, Device U = GPU, typename T = float, typename K = uint16_t> cudaError_t get_dihedral_rmse   (const IsomerBatch<U>& B, CuArray<T>& dihedral_rmse);   //M x 1
        template <ForcefieldType FFT = PEDERSEN, Device U = GPU, typename T = float, typename K = uint16_t> cudaError_t get_dihedral_rrmse  (const IsomerBatch<U>& B, CuArray<T>& dihedral_rrmse);  //M x 1
        template <ForcefieldType FFT = PEDERSEN, Device U = GPU, typename T = float, typename K = uint16_t> cudaError_t get_dihedral_max    (const IsomerBatch<U>& B, CuArray<T>& dihedral_max);   //M x 1
        template <ForcefieldType FFT = PEDERSEN, Device U = GPU, typename T = float, typename K = uint16_t> cudaError_t get_energies        (const IsomerBatch<U>& B, CuArray<T>& energies);       //M x 1

        template <ForcefieldType FFT = PEDERSEN, Device U = GPU, typename T = float, typename K = uint16_t> cudaError_t get_energy          (const IsomerBatch<U>& B, CuArray<T>& energy);         //M x 1
        template <ForcefieldType FFT = PEDERSEN, Device U = GPU, typename T = float, typename K = uint16_t> cudaError_t get_bond_mean       (const IsomerBatch<U>& B, CuArray<T>& bond_mean);      //M x 1
        template <ForcefieldType FFT = PEDERSEN, Device U = GPU, typename T = float, typename K = uint16_t> cudaError_t get_bond_mae        (const IsomerBatch<U>& B, CuArray<T>& bond_mae);       //M x 1
        template <ForcefieldType FFT = PEDERSEN, Device U = GPU, typename T = float, typename K = uint16_t> cudaError_t get_bond_rmse       (const IsomerBatch<U>& B, CuArray<T>& bond_rmse);       //M x 1
        template <ForcefieldType FFT = PEDERSEN, Device U = GPU, typename T = float, typename K = uint16_t> cudaError_t get_bond_rrmse      (const IsomerBatch<U>& B, CuArray<T>& bond_rrmse);      //M x 1
        template <ForcefieldType FFT = PEDERSEN, Device U = GPU, typename T = float, typename K = uint16_t> cudaError_t get_bond_max        (const IsomerBatch<U>& B, CuArray<T>& bond_max);       //M x 1
        template <ForcefieldType FFT = PEDERSEN, Device U = GPU, typename T = float, typename K = uint16_t> cudaError_t get_angle_mean      (const IsomerBatch<U>& B, CuArray<T>& angle_mean);     //M x 1
        template <ForcefieldType FFT = PEDERSEN, Device U = GPU, typename T = float, typename K = uint16_t> cudaError_t get_angle_mae       (const IsomerBatch<U>& B, CuArray<T>& angle_mae);      //M x 1
        template <ForcefieldType FFT = PEDERSEN, Device U = GPU, typename T = float, typename K = uint16_t> cudaError_t get_angle_rmse      (const IsomerBatch<U>& B, CuArray<T>& angle_rms);      //M x 1
        template <ForcefieldType FFT = PEDERSEN, Device U = GPU, typename T = float, typename K = uint16_t> cudaError_t get_angle_rrmse     (const IsomerBatch<U>& B, CuArray<T>& angle_rrmse);     //M x 1
        template <ForcefieldType FFT = PEDERSEN, Device U = GPU, typename T = float, typename K = uint16_t> cudaError_t get_angle_max       (const IsomerBatch<U>& B, CuArray<T>& angle_max);      //M x 1
        
        template <ForcefieldType FFT = PEDERSEN, Device U = GPU, typename T = float, typename K = uint16_t> cudaError_t get_flat_mean       (const IsomerBatch<U>& B, CuArray<T>& flat_mean);      //M x 1
        template <ForcefieldType FFT = PEDERSEN, Device U = GPU, typename T = float, typename K = uint16_t> cudaError_t get_flat_max        (const IsomerBatch<U>& B, CuArray<T>& flat_max);      //M x 1
        template <ForcefieldType FFT = PEDERSEN, Device U = GPU, typename T = float, typename K = uint16_t> cudaError_t get_flat_rmse       (const IsomerBatch<U>& B, CuArray<T>& flat_rms);      //M x 1
        
        
        template <ForcefieldType FFT = PEDERSEN, Device U = GPU, typename T = float, typename K = uint16_t> cudaError_t get_gradient_mean   (const IsomerBatch<U>& B, CuArray<T>& gradient_mean);  //M x 1
        template <ForcefieldType FFT = PEDERSEN, Device U = GPU, typename T = float, typename K = uint16_t> cudaError_t get_gradient_norm   (const IsomerBatch<U>& B, CuArray<T>& gradient_norm);  //M x 1
        template <ForcefieldType FFT = PEDERSEN, Device U = GPU, typename T = float, typename K = uint16_t> cudaError_t get_gradient_rms    (const IsomerBatch<U>& B, CuArray<T>& gradient_rms);   //M x 1
        template <ForcefieldType FFT = PEDERSEN, Device U = GPU, typename T = float, typename K = uint16_t> cudaError_t get_gradient_max    (const IsomerBatch<U>& B, CuArray<T>& gradient_max);   //M x 1

    }

    namespace isomerspace_X0{
        //Returns total time spent in kernel in milliseconds. First call is not timed as that would require synchronous execution.
        std::chrono::microseconds time_spent(); 

        void reset_time(); //Resets time spent in kernel to 0.

        //Uses spherical projection and 'wrapping a sphere' technique to generate starting coordinates in 3D space for a batch of fullerenes.
        template <Device U = GPU, typename T = float, typename K = uint16_t> 
        cudaError_t zero_order_geometry(IsomerBatch<U>& B, const float scalerad, const LaunchCtx& ctx = LaunchCtx(), const LaunchPolicy policy = LaunchPolicy::SYNC);
    }

    namespace isomerspace_tutte{
        //Returns total time spent in kernel in milliseconds. First call is not timed as that would require synchronous execution.
        std::chrono::microseconds time_spent(); 

        void reset_time(); //Resets time spent in kernel to 0.
        //Computes tutte embedding for an entire batch of fullernes.
        template <Device U = GPU, typename T = float, typename K = uint16_t> 
        cudaError_t tutte_layout(IsomerBatch<U>& B, const size_t max_iterations = 10000000, const LaunchCtx& ctx = LaunchCtx(), const LaunchPolicy policy = LaunchPolicy::SYNC);
    }

    namespace isomerspace_dual{
        //Returns the optimal batch size for a given number of atoms, for the dualisation kernel.   
        template <Device U = GPU, typename K = uint16_t>
        int optimal_batch_size(const int N, const int device_id = 0); 
        template <Device U = GPU, typename K = uint16_t>
        int optimal_batch_size_2(const int N, const int device_id = 0);
        
        //Returns total time spent in kernel in milliseconds. First call is not timed as that would require synchronous execution.
        std::chrono::microseconds time_spent(); 

        void reset_time(); //Resets time spent in kernel to 0.
        //Computes the cubic neighbour list from the dual neighbour list
        template <Device U = GPU, typename K = uint16_t>
        cudaError_t dualise(IsomerBatch<U>& B, const LaunchCtx& ctx = LaunchCtx(), const LaunchPolicy policy = LaunchPolicy::SYNC);

        template <Device U = GPU, typename K = uint16_t>
        cudaError_t dualise_2(IsomerBatch<U>& B, const LaunchCtx& ctx = LaunchCtx(), const LaunchPolicy policy = LaunchPolicy::SYNC);

        template <Device U = GPU, typename K = uint16_t>
        void dualise_3(IsomerBatch<U>& B);

        template <Device U = GPU, typename K = uint16_t>
        void dualise_4(IsomerBatch<U>& B);
    }

    namespace isomerspace_eigen{
        template <Device U = GPU, typename T = float, typename K = uint16_t>
        void eigensolve_cusolver(const IsomerBatch<U>& B, const CuArray<T>& hessians, const CuArray<K>& cols, CuArray<T>& eigenvalues, const LaunchCtx& ctx = LaunchCtx(), const LaunchPolicy policy = LaunchPolicy::SYNC);

        template <Device U = GPU, typename T = float, typename K = uint16_t>        //Compute the full spectrum and eigenvectors of the hessian matrix for a batch of fullerenes.
        void eigensolve(const IsomerBatch<U>& B, CuArray<T>& Q, const CuArray<T>& hessians, const CuArray<K>& cols, CuArray<T>& eigenvalues, const LaunchCtx& ctx = LaunchCtx(), const LaunchPolicy policy = LaunchPolicy::SYNC);
        
        template <Device U = GPU, typename T = float, typename K = uint16_t>        //Compute the full spectrum of the hessian matrix for a batch of fullerenes.
        void eigensolve(const IsomerBatch<U>& B, const CuArray<T>& hessians, const CuArray<K>& cols, CuArray<T>& eigenvalues, const LaunchCtx& ctx = LaunchCtx(), const LaunchPolicy policy = LaunchPolicy::SYNC);  
      
        template <Device U = GPU, typename T = float, typename K = uint16_t>         //Computes the smallest and largest eigenvalues of the hessian matrix for a batch of fullerenes.
        void spectrum_ends(const IsomerBatch<U>& B, const CuArray<T>& hessians, const CuArray<K>& cols, CuArray<T>& lambda_mins, CuArray<T>& lambda_maxs, int m_lanczos=40, const LaunchCtx& ctx = LaunchCtx(), const LaunchPolicy policy = LaunchPolicy::SYNC);
      
        template <Device U = GPU, typename T = float, typename K = uint16_t>        //Computes the smallest and largest eigenvalues and corresponding eigenvectors of the hessian matrix for a batch of fullerenes.
        void spectrum_ends(const IsomerBatch<U>& B, const CuArray<T>& hessians, const CuArray<K>& cols, CuArray<T>& lambda_mins, CuArray<T>& lambda_maxs, CuArray<T>& eigvect_mins, CuArray<T>& eigvect_maxs, int m_lanczos=40, const LaunchCtx& ctx = LaunchCtx(), const LaunchPolicy policy = LaunchPolicy::SYNC);
    }

    namespace isomerspace_hessian{
        template <ForcefieldType FFT = PEDERSEN, Device U = GPU, typename T = float, typename K = uint16_t>
        cudaError_t compute_hessians(const IsomerBatch<U>& B, CuArray<T>& hessians, CuArray<K>& cols, const LaunchCtx& ctx = LaunchCtx(), const LaunchPolicy policy = LaunchPolicy::SYNC);
        template <ForcefieldType FFT = PEDERSEN, Device U = GPU, typename T = float, typename K = uint16_t>
        cudaError_t compute_hessians_fd(const IsomerBatch<U>& B, CuArray<T>& hessians, CuArray<K>& cols, const T rel_delta = 1e-5, const LaunchCtx& ctx = LaunchCtx(), const LaunchPolicy policy = LaunchPolicy::SYNC);
    }

    namespace isomerspace_properties{

        //Computes the volume of each molecule in the batch. See sequential code in polyhedron.cc
        template <Device U = GPU, typename T = float, typename K = uint16_t> 
        cudaError_t volume_divergences(const IsomerBatch<U>& B, CuArray<T>& volumes, const LaunchCtx& ctx = LaunchCtx(), const LaunchPolicy policy = LaunchPolicy::SYNC);

        //Computes the surface area of each molecule in the batch. See sequential code in polyhedron.cc
        template <Device U = GPU, typename T = float, typename K = uint16_t> 
        cudaError_t surface_areas(const IsomerBatch<U>& B, CuArray<T>& surface_areas, const LaunchCtx& ctx = LaunchCtx(), const LaunchPolicy policy = LaunchPolicy::SYNC);

      //Computes the moments of inertia (eigenvalues to the moment of inertia matrix)
        template <Device U = GPU, typename T = float, typename K = uint16_t> 
        cudaError_t moments_of_inertia(const IsomerBatch<U>& B, CuArray<T>& lambdas, const LaunchCtx& ctx = LaunchCtx(), const LaunchPolicy policy = LaunchPolicy::SYNC);      
      //Computes the best elipsoid fit for each molecule in the batch and returns the eccentricities of the ellipsoids.
        template <Device U = GPU, typename T = float, typename K = uint16_t> 
        cudaError_t eccentricities(const IsomerBatch<U>& B, CuArray<T>& eccentricities, const LaunchCtx& ctx = LaunchCtx(), const LaunchPolicy policy = LaunchPolicy::SYNC);

        //Moves coordinates to the centre of mass of the molecule and aligns the molecule to the principal axes of inertia. See sequential code in polyhedron.cc
        template <Device U = GPU, typename T = float, typename K = uint16_t> 
        cudaError_t transform_coordinates(IsomerBatch<U>& B, const LaunchCtx& ctx = LaunchCtx(), const LaunchPolicy policy = LaunchPolicy::SYNC);

        //cudaError_t debug_function(const IsomerBatch<U>& B, CuArray<double>& eigenvalues, CuArray<double>& eigenvectors, CuArray<double>& inertia_matrices, CuArray<double>& orthogonality, const LaunchCtx& ctx = LaunchCtx(), const LaunchPolicy policy = LaunchPolicy::SYNC);
        template <Device U = GPU, typename T = float, typename K = uint16_t> 
        cudaError_t debug_function(const IsomerBatch<U>& B, CuArray<T>& eigenvalues, CuArray<T>& eigenvectors, CuArray<T>& inertia_matrices, CuArray<T>& orthogonality, const LaunchCtx& ctx = LaunchCtx(), const LaunchPolicy policy = LaunchPolicy::SYNC);
    }
}
#endif
