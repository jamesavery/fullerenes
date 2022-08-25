#include "isomerspace_kernel.hh"
#include "cuda_execution.hh"
#include "cu_array.hh"
#include "chrono"

enum ForcefieldType {WIRZ, BUSTER, FLATNESS_ENABLED};

namespace gpu_kernels{
    namespace isomerspace_forcefield{
        int optimal_batch_size(const int N, const int device_id = 0);

        //Returns total time spent in kernel in milliseconds. First call is not timed as that would require synchronous execution.
        std::chrono::microseconds time_spent(); 
        //Uses forcefield optimization to relax the positions of atoms in all isomers of a batch.
        cudaError_t optimize_batch(IsomerBatch& B, const size_t iterations, const size_t max_iterations, const LaunchCtx& ctx = LaunchCtx(), const LaunchPolicy policy = LaunchPolicy::SYNC);
        cudaError_t test_fun(IsomerBatch& B, CuArray<device_real_t>& output);
        //Stat functions:
        cudaError_t get_bonds           (const IsomerBatch& B, CuArray<device_real_t>& bonds);          //N x M x 3
        cudaError_t get_angles          (const IsomerBatch& B, CuArray<device_real_t>& angles);         //N x M x 3
        cudaError_t get_dihedrals       (const IsomerBatch& B, CuArray<device_real_t>& dihedrals);      //N x M x 3

        cudaError_t get_energy          (const IsomerBatch& B, CuArray<device_real_t>& energy);         //M x 1
        cudaError_t get_bond_mean       (const IsomerBatch& B, CuArray<device_real_t>& bond_mean);      //M x 1
        cudaError_t get_angle_mean      (const IsomerBatch& B, CuArray<device_real_t>& angle_mean);     //M x 1
        cudaError_t get_dihedral_mean   (const IsomerBatch& B, CuArray<device_real_t>& dihedral_mean);  //M x 1

        cudaError_t get_bond_rms        (const IsomerBatch& B, CuArray<device_real_t>& bond_rms);       //M x 1
        cudaError_t get_angle_rms       (const IsomerBatch& B, CuArray<device_real_t>& angle_rms);      //M x 1
        cudaError_t get_dihedral_rms    (const IsomerBatch& B, CuArray<device_real_t>& dihedral_rms);   //M x 1
        cudaError_t get_bond_max        (const IsomerBatch& B, CuArray<device_real_t>& bond_max);       //M x 1
        cudaError_t get_angle_max       (const IsomerBatch& B, CuArray<device_real_t>& angle_max);      //M x 1
        cudaError_t get_dihedral_max    (const IsomerBatch& B, CuArray<device_real_t>& dihedral_max);   //M x 1
        cudaError_t get_energies        (const IsomerBatch& B, CuArray<device_real_t>& energies);       //M x 1
        
        cudaError_t get_gradient_mean   (const IsomerBatch& B, CuArray<device_real_t>& gradient_mean);  //M x 1
        cudaError_t get_gradient_norm   (const IsomerBatch& B, CuArray<device_real_t>& gradient_norm);  //M x 1
        cudaError_t get_gradient_rms    (const IsomerBatch& B, CuArray<device_real_t>& gradient_rms);   //M x 1
        cudaError_t get_gradient_max    (const IsomerBatch& B, CuArray<device_real_t>& gradient_max);   //M x 1

    }

    namespace isomerspace_X0{
        //Returns total time spent in kernel in milliseconds. First call is not timed as that would require synchronous execution.
        std::chrono::microseconds time_spent(); 
        //Uses spherical projection and 'wrapping a sphere' technique to generate starting coordinates in 3D space for a batch of fullerenes.
        cudaError_t zero_order_geometry(IsomerBatch& B, const device_real_t scalerad, const LaunchCtx& ctx = LaunchCtx(), const LaunchPolicy policy = LaunchPolicy::SYNC);
    }

    namespace isomerspace_tutte{
        //Returns total time spent in kernel in milliseconds. First call is not timed as that would require synchronous execution.
        std::chrono::microseconds time_spent(); 
        //Computes tutte embedding for an entire batch of fullernes.
        cudaError_t tutte_layout(IsomerBatch& B, const size_t max_iterations = 10000000, const LaunchCtx& ctx = LaunchCtx(), const LaunchPolicy policy = LaunchPolicy::SYNC);
    }

    namespace isomerspace_dual{
        //Returns total time spent in kernel in milliseconds. First call is not timed as that would require synchronous execution.
        std::chrono::microseconds time_spent(); 
        //Computes the cubic neighbour list from the dual neighbour list
        cudaError_t cubic_layout(IsomerBatch& B, const LaunchCtx& ctx = LaunchCtx(), const LaunchPolicy policy = LaunchPolicy::SYNC);
    }
}
