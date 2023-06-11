//Container for launch dimensions on all devices, can be stored statically in a function.
struct LaunchDims
{   
    //Figures out which device is being used and gets the appropriate grid dimensions.
    dim3 get_grid(){
        int device_id; cudaGetDevice(&device_id);
        return m_grid_dims[device_id];
    }

    //Figures out which device is being used and gets the appropriate block dimensions.
    dim3 get_block(){
        int device_id; cudaGetDevice(&device_id);
        return m_block_dims[device_id];
    }

    
    /**
    * @brief LaunchDims constructor
    * @param kernel pointer to the kernel function
    * @param block_size_ number of threads per block
    * @param smem_ amount of dynamic shared memory to allocate
    * @param num_tasks number of tasks of size block_size_ to be solved in the kernel launch
    */
    LaunchDims(const void* kernel, const int block_size_, const size_t smem_ = 0, const int num_tasks = -1) {
        cudaGetDeviceCount(&m_device_count);
        m_block_dims.resize(m_device_count);
        m_grid_dims.resize(m_device_count);
        m_properties.resize(m_device_count);
        //Fetch the multiprocessor count for each device.
        for (size_t i = 0; i < m_device_count; i++) {
            cudaGetDeviceProperties(&m_properties[i], i);
        }
        
        update_dims(kernel, block_size_, smem_, num_tasks);
    }
    //Only recomputes dimensions if a new block size or new amount of dynamic shared memory is specified.

    /**
     * @brief This function computes the optimal grid dimensions for a given kernel that maximizes occupancy.
     * @param kernel - The kernel to be launched
     * @param block_size_ - The size of the block to be used for the kernel launch
     * @param smem_ - The amount of shared memory per block to be allocated dynamically when the kernel is launched
     * @param num_tasks - The number of problems of size block_size_ to be solved by the kernel
     * @return cudaError_t - The error returned by the calculation if any
     */
    cudaError_t update_dims(const void* kernel, const int block_size_, const size_t smem_ = 0, const int num_tasks = -1){
        if(block_size_ == m_block_size && smem_ == m_smem) return cudaSuccess;
        m_block_size = block_size_; m_smem = smem_;
        int current_device; cudaGetDevice(&current_device);
        int temp_int; //CUDA api is pretty particular about integer types
        for (size_t i = 0; i < m_device_count; i++)
        {   
            cudaSetDevice(i);
            cudaOccupancyMaxActiveBlocksPerMultiprocessor(&temp_int,kernel,m_block_size,m_smem);
            int SMs = m_properties[i].multiProcessorCount;
            if( num_tasks == -1) {m_grid_dims[i].x = temp_int * SMs;}
            else{ 
                int num_blocks = 0;
                int best_blocks = temp_int*SMs;
                int best_remainder = temp_int*SMs;
                for (num_blocks = temp_int*SMs; num_blocks >= SMs ; num_blocks -= 1)
                {    
                    int remainder = (num_tasks / num_blocks + (num_tasks % num_blocks != 0)) * num_blocks - num_tasks;
                    if (remainder < best_remainder) {
                        best_blocks = num_blocks; 
                        best_remainder = remainder;}
                    if (best_remainder == 0) break;
                }
                m_grid_dims[i].x = best_blocks;
            }
            m_block_dims[i].x = m_block_size;
        }
        cudaSetDevice(current_device);
        cudaDeviceSynchronize();
        return cudaGetLastError();
    }

private:

    int m_block_size{};
    int m_smem{};
    int m_device_count{};
    int m_num_tasks{};
    std::vector<cudaDeviceProp> m_properties;
    std::vector<dim3> m_grid_dims;
    std::vector<dim3> m_block_dims;
};