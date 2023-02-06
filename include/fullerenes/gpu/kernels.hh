#ifndef ISOMERSPACE_KERNELS_H
#define ISOMERSPACE_KERNELS_H
#include "fullerenes/gpu/cuda_definitions.h"
#include "fullerenes/gpu/launch_ctx.hh"
#include "fullerenes/gpu/isomer_batch.hh"
#include "launch_ctx.hh"
#include "cu_array.hh"
#include <chrono>

enum ForcefieldType {WIRZ, PEDERSEN, FLATNESS_ENABLED, FLAT_BOND};

namespace gpu_kernels{
    namespace isomerspace_forcefield{
        int optimal_batch_size(const int N, const int device_id = 0);

        //Returns total time spent in kernel in milliseconds. First call is not timed as that would require synchronous execution.
        std::chrono::microseconds time_spent(); 

        //Resets the time spent in kernel to 0
        void reset_time();

        //Uses forcefield optimization to relax the positions of atoms in all isomers of a batch.
        template <ForcefieldType T>
        cudaError_t optimise(IsomerBatch& B, const size_t iterations, const size_t max_iterations, const LaunchCtx& ctx = LaunchCtx(), const LaunchPolicy policy = LaunchPolicy::SYNC);
        cudaError_t test_fun(IsomerBatch& B, CuArray<device_real_t>& output);
        //ForcefieldType agnostic functions since the forcefield type doesnt change the way bonds are computed
        //Retrieves all interatomic bond lengths, angles and dihedrals
        cudaError_t get_bonds           (const IsomerBatch& B, CuArray<device_real_t>& bonds);          //N x M x 3
        cudaError_t get_angles          (const IsomerBatch& B, CuArray<device_real_t>& angles);         //N x M x 3
        cudaError_t get_dihedrals       (const IsomerBatch& B, CuArray<device_real_t>& dihedrals);      //N x M x 3


        //MAE = Mean Absolute Error
        //RMSE = Root Mean Squared Error
        //RRMSE = Relative Root Mean Squared Error

        //Stat functions:
        template <ForcefieldType T> cudaError_t get_dihedral_mean   (const IsomerBatch& B, CuArray<device_real_t>& dihedral_mean);  //M x 1
        template <ForcefieldType T> cudaError_t get_dihedral_mae    (const IsomerBatch& B, CuArray<device_real_t>& dihedral_mae);   //M x 1 
        template <ForcefieldType T> cudaError_t get_dihedral_rmse   (const IsomerBatch& B, CuArray<device_real_t>& dihedral_rmse);   //M x 1
        template <ForcefieldType T> cudaError_t get_dihedral_rrmse  (const IsomerBatch& B, CuArray<device_real_t>& dihedral_rrmse);  //M x 1
        template <ForcefieldType T> cudaError_t get_dihedral_max    (const IsomerBatch& B, CuArray<device_real_t>& dihedral_max);   //M x 1
        template <ForcefieldType T> cudaError_t get_energies        (const IsomerBatch& B, CuArray<device_real_t>& energies);       //M x 1

        template <ForcefieldType T> cudaError_t get_energy          (const IsomerBatch& B, CuArray<device_real_t>& energy);         //M x 1
        template <ForcefieldType T> cudaError_t get_bond_mean       (const IsomerBatch& B, CuArray<device_real_t>& bond_mean);      //M x 1
        template <ForcefieldType T> cudaError_t get_bond_mae        (const IsomerBatch& B, CuArray<device_real_t>& bond_mae);       //M x 1
        template <ForcefieldType T> cudaError_t get_bond_rmse       (const IsomerBatch& B, CuArray<device_real_t>& bond_rmse);       //M x 1
        template <ForcefieldType T> cudaError_t get_bond_rrmse      (const IsomerBatch& B, CuArray<device_real_t>& bond_rrmse);      //M x 1
        template <ForcefieldType T> cudaError_t get_bond_max        (const IsomerBatch& B, CuArray<device_real_t>& bond_max);       //M x 1
        template <ForcefieldType T> cudaError_t get_angle_mean      (const IsomerBatch& B, CuArray<device_real_t>& angle_mean);     //M x 1
        template <ForcefieldType T> cudaError_t get_angle_mae       (const IsomerBatch& B, CuArray<device_real_t>& angle_mae);      //M x 1
        template <ForcefieldType T> cudaError_t get_angle_rmse      (const IsomerBatch& B, CuArray<device_real_t>& angle_rms);      //M x 1
        template <ForcefieldType T> cudaError_t get_angle_rrmse     (const IsomerBatch& B, CuArray<device_real_t>& angle_rrmse);     //M x 1
        template <ForcefieldType T> cudaError_t get_angle_max       (const IsomerBatch& B, CuArray<device_real_t>& angle_max);      //M x 1
        
        template <ForcefieldType T> cudaError_t get_flat_mean       (const IsomerBatch& B, CuArray<device_real_t>& flat_mean);      //M x 1
        template <ForcefieldType T> cudaError_t get_flat_max        (const IsomerBatch& B, CuArray<device_real_t>& flat_max);      //M x 1
        template <ForcefieldType T> cudaError_t get_flat_rmse       (const IsomerBatch& B, CuArray<device_real_t>& flat_rms);      //M x 1
        
        
        template <ForcefieldType T> cudaError_t get_gradient_mean   (const IsomerBatch& B, CuArray<device_real_t>& gradient_mean);  //M x 1
        template <ForcefieldType T> cudaError_t get_gradient_norm   (const IsomerBatch& B, CuArray<device_real_t>& gradient_norm);  //M x 1
        template <ForcefieldType T> cudaError_t get_gradient_rms    (const IsomerBatch& B, CuArray<device_real_t>& gradient_rms);   //M x 1
        template <ForcefieldType T> cudaError_t get_gradient_max    (const IsomerBatch& B, CuArray<device_real_t>& gradient_max);   //M x 1

    }

    namespace isomerspace_X0{
        //Returns total time spent in kernel in milliseconds. First call is not timed as that would require synchronous execution.
        std::chrono::microseconds time_spent(); 

        void reset_time(); //Resets time spent in kernel to 0.

        //Uses spherical projection and 'wrapping a sphere' technique to generate starting coordinates in 3D space for a batch of fullerenes.
        cudaError_t zero_order_geometry(IsomerBatch& B, const device_real_t scalerad, const LaunchCtx& ctx = LaunchCtx(), const LaunchPolicy policy = LaunchPolicy::SYNC);
    }

    namespace isomerspace_tutte{
        //Returns total time spent in kernel in milliseconds. First call is not timed as that would require synchronous execution.
        std::chrono::microseconds time_spent(); 

        void reset_time(); //Resets time spent in kernel to 0.
        //Computes tutte embedding for an entire batch of fullernes.
        cudaError_t tutte_layout(IsomerBatch& B, const size_t max_iterations = 10000000, const LaunchCtx& ctx = LaunchCtx(), const LaunchPolicy policy = LaunchPolicy::SYNC);
    }

    namespace isomerspace_dual{
        //Returns total time spent in kernel in milliseconds. First call is not timed as that would require synchronous execution.
        std::chrono::microseconds time_spent(); 

        void reset_time(); //Resets time spent in kernel to 0.
        //Computes the cubic neighbour list from the dual neighbour list
        cudaError_t dualise(IsomerBatch& B, const LaunchCtx& ctx = LaunchCtx(), const LaunchPolicy policy = LaunchPolicy::SYNC);
    }
}
#endif