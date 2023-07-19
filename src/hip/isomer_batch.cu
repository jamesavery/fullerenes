#include "hip/hip_runtime.h"
#include "fullerenes/gpu/isomer_batch.hh"

IsomerBatch::IsomerBatch(size_t n_atoms, size_t n_isomers, Device buffer_type, int device){
    this->buffer_type      = buffer_type;
    this->n_atoms          = n_atoms;
    this->isomer_capacity  = n_isomers;
    this->n_faces          = n_atoms/2 + 2;
    this->m_device         = device;
    pointers =   {{"cubic_neighbours",(void**)&cubic_neighbours, sizeof(device_node_t)*n_atoms*3, true}, {"dual_neighbours", (void**)&dual_neighbours, sizeof(device_node_t) * (n_atoms/2 +2) * 6, true}, {"face_degrees", (void**)&face_degrees, sizeof(uint8_t)*(n_atoms/2 +2), true},{"X", (void**)&X, sizeof(device_real_t)*n_atoms*3, true}, {"xys", (void**)&xys, sizeof(device_hpreal_t)*n_atoms*2, true}, {"statuses", (void**)&statuses, sizeof(IsomerStatus), false}, {"IDs", (void**)&IDs, sizeof(size_t), false}, {"iterations", (void**)&iterations, sizeof(size_t), false}};
    if (buffer_type == GPU){
        hipSetDevice(m_device);
        for (size_t i = 0; i < pointers.size(); i++) {
            hipMalloc(get<1>(pointers[i]), n_isomers * get<2>(pointers[i])); 
            hipMemset(*get<1>(pointers[i]),0,n_isomers*get<2>(pointers[i]));
        }
    } else if(buffer_type == CPU){
        for (size_t i = 0; i < pointers.size(); i++) {
            //For asynchronous memory transfers host memory must be pinned. 
            hipHostMalloc(get<1>(pointers[i]), n_isomers * get<2>(pointers[i]));
            memset(*get<1>(pointers[i]),0, n_isomers*get<2>(pointers[i]));
        }
    }
    printLastCudaError("Failed to construct IsomerBatch");
    allocated = true;
}

void IsomerBatch::operator=(const IsomerBatch& input){
    hipSetDevice(m_device);
    pointers =   {{"cubic_neighbours",(void**)&cubic_neighbours, sizeof(device_node_t)*n_atoms*3, true}, {"dual_neighbours", (void**)&dual_neighbours, sizeof(device_node_t) * (n_atoms/2 +2) * 6, true}, {"face_degrees", (void**)&face_degrees, sizeof(uint8_t)*(n_atoms/2 +2), true},{"X", (void**)&X, sizeof(device_real_t)*n_atoms*3, true}, {"xys", (void**)&xys, sizeof(device_hpreal_t)*n_atoms*2, true}, {"statuses", (void**)&statuses, sizeof(IsomerStatus), false}, {"IDs", (void**)&IDs, sizeof(size_t), false}, {"iterations", (void**)&iterations, sizeof(size_t), false}};
    if (allocated == true){
        for (size_t i = 0; i < pointers.size(); i++) {
            hipFree(*get<1>(pointers[i]));
        }
        allocated = false;
    }
    hipSetDevice(input.m_device);
    //Construct a tempory batch: allocates the needed amount of memory.
    this->m_device = input.m_device;
    this->buffer_type = input.buffer_type;
    this->isomer_capacity = input.isomer_capacity;
    this->n_atoms = input.n_atoms;
    this->n_faces = input.n_faces;
    
    
    //Copy contents of old batch into newly allocated memory.
    if (buffer_type == GPU){
        for (size_t i = 0; i < pointers.size(); i++) {
            hipMalloc(get<1>(pointers[i]), isomer_capacity * get<2>(pointers[i]));
            hipMemcpy(*get<1>(pointers[i]), *get<1>(input.pointers[i]), isomer_capacity * get<2>(pointers[i]), hipMemcpyDeviceToDevice);
        }
    } else if(buffer_type == CPU){
        for (size_t i = 0; i < pointers.size(); i++) {
            //For asynchronous memory transfers host memory must be pinned. 
            hipHostMalloc(get<1>(pointers[i]), isomer_capacity * get<2>(pointers[i]));
        }
    }
    printLastCudaError("Failed to copy IsomerBatch");
}


IsomerBatch::~IsomerBatch(){
    if (allocated == true);
    {
    if (buffer_type == GPU){    
        hipSetDevice(m_device);
        for (size_t i = 0; i < pointers.size(); i++) {
            hipFree(*get<1>(pointers[i]));
        }
    } else{
        for (size_t i = 0; i < pointers.size(); i++) {
            hipHostFree(*get<1>(pointers[i])); 
        }
    }
    }
    allocated = false;
}

__global__
void reset_convergence_status_(IsomerBatch B){
    for (int isomer_idx = blockIdx.x; isomer_idx < B.isomer_capacity; isomer_idx+=gridDim.x)
    {
        if (threadIdx.x == 0) B.statuses[isomer_idx] = B.statuses[isomer_idx] != IsomerStatus::EMPTY ? IsomerStatus::NOT_CONVERGED : IsomerStatus::EMPTY;
        if (threadIdx.x == 0) B.iterations[isomer_idx] = 0;
    }
    
}

hipError_t reset_convergence_statuses(IsomerBatch& B, const LaunchCtx& ctx, const LaunchPolicy policy){
    hipSetDevice(ctx.get_device_id());
    static LaunchDims dims((void*)reset_convergence_status_, B.n_atoms);
    dims.update_dims((void*)reset_convergence_status_,B.n_atoms,0,B.isomer_capacity);
    if(policy == LaunchPolicy::SYNC) ctx.wait();
    void* kargs[]{(void*)&B};
    hipLaunchCooperativeKernel((void*)reset_convergence_status_, dims.get_grid(), dims.get_block(), kargs, 0, ctx.stream);
    if(policy == LaunchPolicy::SYNC) ctx.wait();
    
    return hipGetLastError();
}


GPUDataStruct::GPUDataStruct(size_t n_atoms, size_t n_isomers, Device buffer_type){
    this->buffer_type = buffer_type;
    this->n_atoms     = n_atoms;
    this->n_isomers   = n_isomers;
    this->n_faces     = n_atoms/2 + 2;
    if (buffer_type == GPU){
        for (size_t i = 0; i < pointers.size(); i++) {
            size_t num_elements = get<3>(pointers[i]) ? n_isomers * n_atoms: n_isomers;
            hipMalloc(get<1>(pointers[i]), num_elements* get<2>(pointers[i])); 
        }
    } else if(buffer_type == CPU){
        for (size_t i = 0; i < pointers.size(); i++) {
            size_t num_elements = get<3>(pointers[i]) ? n_isomers * n_atoms: n_isomers;
            //For asynchronous memory transfers host memory must be pinned. 
            hipHostMalloc(get<1>(pointers[i]), num_elements* get<2>(pointers[i]));
        }
    }
    printLastCudaError("Failed to construct GPUDataStruct");
    allocated = true;
    
}

GPUDataStruct::~GPUDataStruct(){
    if(allocated){
        if (buffer_type == GPU){    
            for (size_t i = 0; i < pointers.size(); i++) {
                hipFree(*get<1>(pointers[i]));
            }
        } else{
            for (size_t i = 0; i < pointers.size(); i++) {
                hipHostFree(*get<1>(pointers[i])); 
            }
        }
        allocated = false;
    }
}

void GPUDataStruct::allocate(GPUDataStruct& G, const size_t n_atoms, const size_t n_isomers, const Device buffer_type){
    if((!G.allocated)){
        G.buffer_type = buffer_type;
        G.n_atoms = n_atoms; 
        G.isomer_capacity = n_isomers;
        if (buffer_type == GPU){
            for (size_t i = 0; i < G.pointers.size(); i++) {
                size_t num_elements = get<3>(G.pointers[i]) ? n_isomers * n_atoms: n_isomers;
                hipMalloc(get<1>(G.pointers[i]), num_elements* get<2>(G.pointers[i])); 
            }
        }else{
            for (size_t i = 0; i < G.pointers.size(); i++) {
                size_t num_elements = get<3>(G.pointers[i]) ? n_isomers * n_atoms: n_isomers;
                //For asynchronous memory transfers host memory must be pinned. 
                hipHostMalloc(get<1>(G.pointers[i]), num_elements* get<2>(G.pointers[i]));
            }
        }        
        printLastCudaError("Failed to allocate struct");
        G.allocated = true;
    }
}

void GPUDataStruct::initialize(GPUDataStruct& G){
    if(G.buffer_type == GPU){
        for (size_t i = 0; i < G.pointers.size(); i++) {
            size_t num_elements = get<3>(G.pointers[i]) ? G.isomer_capacity * G.n_atoms: G.isomer_capacity;
            hipMemset(*get<1>(G.pointers[i]),0,num_elements*get<2>(G.pointers[i]));
        }
    } else {
        for (size_t i = 0; i < G.pointers.size(); i++) {
            size_t num_elements = get<3>(G.pointers[i]) ? G.isomer_capacity * G.n_atoms: G.isomer_capacity;
            memset(*get<1>(G.pointers[i]),0, num_elements*get<2>(G.pointers[i]));
        }
    }
}

void GPUDataStruct::copy(GPUDataStruct& destination, const GPUDataStruct& source, const hipStream_t stream){
    for (size_t i = 0; i < source.pointers.size(); i++)
    {
        size_t num_elements = get<3>(source.pointers[i]) ?  source.n_atoms * source.isomer_capacity : source.isomer_capacity;
        hipMemcpyAsync(*(get<1>(destination.pointers[i])) , *(get<1>(source.pointers[i])), get<2>(source.pointers[i])*num_elements, hipMemcpyKind(2*source.buffer_type +  destination.buffer_type), stream);
    }
    destination.n_isomers = source.n_isomers;
    printLastCudaError("Failed to copy struct");
}

void operator <<= (GPUDataStruct& destination, const GPUDataStruct& source){
    GPUDataStruct::copy(destination, source);
} 
