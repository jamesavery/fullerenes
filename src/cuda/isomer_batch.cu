#include "fullerenes/isomer_batch.hh"

template<Device T>
IsomerBatch<T>::IsomerBatch(size_t n_atoms, size_t n_isomers, int device){
    this->n_atoms          = n_atoms;
    this->isomer_capacity  = n_isomers;
    this->n_faces          = n_atoms/2 + 2;
    this->m_device         = device;
    pointers =   {{"cubic_neighbours",(void**)&cubic_neighbours, sizeof(device_node_t)*n_atoms*3, true}, {"dual_neighbours", (void**)&dual_neighbours, sizeof(device_node_t) * (n_atoms/2 +2) * 6, true}, {"face_degrees", (void**)&face_degrees, sizeof(uint8_t)*(n_atoms/2 +2), true},{"X", (void**)&X, sizeof(device_real_t)*n_atoms*3, true}, {"xys", (void**)&xys, sizeof(device_real_t)*n_atoms*2, true}, {"statuses", (void**)&statuses, sizeof(IsomerStatus), false}, {"IDs", (void**)&IDs, sizeof(size_t), false}, {"iterations", (void**)&iterations, sizeof(size_t), false}};
    if (T == GPU){
        cudaSetDevice(m_device);
        for (size_t i = 0; i < pointers.size(); i++) {
            cudaMalloc(get<1>(pointers[i]), n_isomers * get<2>(pointers[i])); 
            cudaMemset(*get<1>(pointers[i]),0,n_isomers*get<2>(pointers[i]));
        }
    } else if(T == CPU){
        for (size_t i = 0; i < pointers.size(); i++) {
            //For asynchronous memory transfers host memory must be pinned. 
            cudaMallocHost(get<1>(pointers[i]), n_isomers * get<2>(pointers[i]));
            memset(*get<1>(pointers[i]),0, n_isomers*get<2>(pointers[i]));
        }
    }
    printLastCudaError("Failed to construct IsomerBatch");
    allocated = true;
}

template<Device T>
void IsomerBatch<T>::operator=(const IsomerBatch<T>& input){
    cudaSetDevice(m_device);
    pointers =   {{"cubic_neighbours",(void**)&cubic_neighbours, sizeof(device_node_t)*n_atoms*3, true}, {"dual_neighbours", (void**)&dual_neighbours, sizeof(device_node_t) * (n_atoms/2 +2) * 6, true}, {"face_degrees", (void**)&face_degrees, sizeof(uint8_t)*(n_atoms/2 +2), true},{"X", (void**)&X, sizeof(device_real_t)*n_atoms*3, true}, {"xys", (void**)&xys, sizeof(device_real_t)*n_atoms*2, true}, {"statuses", (void**)&statuses, sizeof(IsomerStatus), false}, {"IDs", (void**)&IDs, sizeof(size_t), false}, {"iterations", (void**)&iterations, sizeof(size_t), false}};
    if (allocated == true){
        for (size_t i = 0; i < pointers.size(); i++) {
            cudaFree(*get<1>(pointers[i]));
        }
        allocated = false;
    }
    cudaSetDevice(input.m_device);
    //Construct a tempory batch: allocates the needed amount of memory.
    this->m_device = input.m_device;
    this->isomer_capacity = input.isomer_capacity;
    this->n_atoms = input.n_atoms;
    this->n_faces = input.n_faces;
    
    
    //Copy contents of old batch into newly allocated memory.
    if (T == GPU){
        for (size_t i = 0; i < pointers.size(); i++) {
            cudaMalloc(get<1>(pointers[i]), isomer_capacity * get<2>(pointers[i]));
            cudaMemcpy(*get<1>(pointers[i]), *get<1>(input.pointers[i]), isomer_capacity * get<2>(pointers[i]), cudaMemcpyDeviceToDevice);
        }
    } else if(T == CPU){
        for (size_t i = 0; i < pointers.size(); i++) {
            //For asynchronous memory transfers host memory must be pinned. 
            cudaMallocHost(get<1>(pointers[i]), isomer_capacity * get<2>(pointers[i]));
	    //TODO: Isn't this a bug? Nothing is being copied!
        }
    }
    printLastCudaError("Failed to copy IsomerBatch");
}

template<Device T>
IsomerBatch<T>::~IsomerBatch(){
    if (allocated == true);
    {
    if (T == GPU){    
        cudaSetDevice(m_device);
        for (size_t i = 0; i < pointers.size(); i++) {
            cudaFree(*get<1>(pointers[i]));
        }
    } else if (T == CPU){
        for (size_t i = 0; i < pointers.size(); i++) {
            cudaFreeHost(*get<1>(pointers[i])); 
        }
    }
    }
    allocated = false;
}

template<Device T>
__global__
void reset_convergence_status_(IsomerBatch<T> B){
    for (int isomer_idx = blockIdx.x; isomer_idx < B.isomer_capacity; isomer_idx+=gridDim.x)
    {
        if (threadIdx.x == 0) B.statuses[isomer_idx] = B.statuses[isomer_idx] != IsomerStatus::EMPTY ? IsomerStatus::NOT_CONVERGED : IsomerStatus::EMPTY;
        if (threadIdx.x == 0) B.iterations[isomer_idx] = 0;
    }
    
}

template<Device T>
cudaError_t reset_convergence_statuses(IsomerBatch<T>& B, const LaunchCtx& ctx, const LaunchPolicy policy){
    cudaSetDevice(ctx.get_device_id());
    static LaunchDims dims((void*)reset_convergence_status_<GPU>, B.n_atoms);
    dims.update_dims((void*)reset_convergence_status_<GPU>,B.n_atoms,0,B.isomer_capacity);
    if(policy == LaunchPolicy::SYNC) ctx.wait();
    void* kargs[]{(void*)&B};
    cudaLaunchCooperativeKernel((void*)reset_convergence_status_<GPU>, dims.get_grid(), dims.get_block(), kargs, 0, ctx.stream);
    if(policy == LaunchPolicy::SYNC) ctx.wait();
    
    return cudaGetLastError();
}
