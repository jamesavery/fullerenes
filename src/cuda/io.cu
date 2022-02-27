#pragma once
#include "fullerenes/gpu/isomerspace_forcefield.hh"
#include "fullerenes/gpu/isomerspace_tutte.hh"
#include "auxiliary_cuda_functions.cu"
#include "fullerenes/gpu/gpudatastruct.hh"
#include "fullerenes/gpu/isomerspace_kernel.hh"

void GPUDataStruct::allocate(GPUDataStruct& G, size_t N, const size_t batch_size, const BufferType buffer_type){
    if((!G.allocated)){
        G.buffer_type = buffer_type;
        G.batch_size  = batch_size; 
        G.N           = N; 
        if (buffer_type == DEVICE_BUFFER){
            for (size_t i = 0; i < G.pointers.size(); i++) {
                size_t num_elements = get<3>(G.pointers[i]) ?  N*batch_size : batch_size;
                cudaMalloc(get<1>(G.pointers[i]), num_elements* get<2>(G.pointers[i])); 
            }
            printLastCudaError("Failed to allocate device struct");
        }else{
            for (size_t i = 0; i < G.pointers.size(); i++) {
                size_t num_elements = get<3>(G.pointers[i]) ?  N*batch_size : batch_size;
                *get<1>(G.pointers[i])= malloc(num_elements* get<2>(G.pointers[i])); 
            }
        }        
        G.allocated = true;
    }
}

void GPUDataStruct::free(GPUDataStruct& G){
    if(G.allocated){
        if (G.buffer_type == DEVICE_BUFFER){    
            for (size_t i = 0; i < G.pointers.size(); i++) {
                cudaFree(*get<1>(G.pointers[i]));
            }
            printLastCudaError("Failed to free device struct"); 
        } else{
            for (size_t i = 0; i < G.pointers.size(); i++) {
                std::free(*get<1>(G.pointers[i])); 
            }
        }
        G.allocated = false;
    }
}

void GPUDataStruct::copy(GPUDataStruct& destination, const GPUDataStruct& source){
    if(source.batch_size > 0){
    for (size_t i = 0; i < source.pointers.size(); i++)
    {
        size_t num_elements = get<3>(source.pointers[i]) ?  source.N*source.batch_size : source.batch_size;
        cudaMemcpy(*(get<1>(destination.pointers[i])) , *(get<1>(source.pointers[i])), get<2>(source.pointers[i])*num_elements, cudaMemcpyKind(2*source.buffer_type +  destination.buffer_type));
    }
    destination.N = source.N;
    destination.batch_size = source.batch_size;
    }
    else{
        std::cout << "WARNING: Call to copy made for 0 isomers \n";
    }
    printLastCudaError("Failed to copy struct");
}

void operator <<= (GPUDataStruct& destination, const GPUDataStruct& source){
    GPUDataStruct::copy(destination, source);
}

__global__
void clear_convergence_status(IsomerBatch G){
    if (threadIdx.x == 0) G.statuses[blockIdx.x] = G.statuses[blockIdx.x] != EMPTY ? NOT_CONVERGED : EMPTY;
    if (threadIdx.x == 0) G.iterations[blockIdx.x] = 0;
}

__global__
void input_type_conversion(IsomerBatch G){
    node_t input_neighbours[3]  = {reinterpret_cast<node_t*>(G.neighbours)[(threadIdx.x* + blockDim.x*blockIdx.x)*3], reinterpret_cast<node_t*>(G.neighbours)[(threadIdx.x* + blockDim.x*blockIdx.x)*3 + 1], reinterpret_cast<node_t*>(G.neighbours)[(threadIdx.x* + blockDim.x*blockIdx.x)*3 + 2]};
    double input_coordianates[3]= {reinterpret_cast<double*>(G.X)[(threadIdx.x* + blockDim.x*blockIdx.x)*3], reinterpret_cast<double*>(G.X)[(threadIdx.x* + blockDim.x*blockIdx.x)*3 + 1], reinterpret_cast<double*>(G.X)[(threadIdx.x* + blockDim.x*blockIdx.x)*3 + 2]};
    double input_xys[2]         = {reinterpret_cast<double*>(G.xys)[(threadIdx.x* + blockDim.x*blockIdx.x)*2], reinterpret_cast<double*>(G.xys)[(threadIdx.x* + blockDim.x*blockIdx.x)*2 + 1]};
    cg::sync(cg::this_grid());
    reinterpret_cast<device_node3*>(G.neighbours)[(threadIdx.x + blockDim.x*blockIdx.x)] = {input_neighbours[0], input_neighbours[1], input_neighbours[2]};
    reinterpret_cast<device_coord3d*>(G.X)[(threadIdx.x + blockDim.x*blockIdx.x)] = {input_coordianates[0], input_coordianates[1], input_coordianates[2]};
    reinterpret_cast<device_coord2d*>(G.xys)[(threadIdx.x + blockDim.x*blockIdx.x)] = {input_xys[0], input_xys[1]};
}

__global__
void output_type_conversion(IsomerBatch G){
    device_node3 output_neighbours = reinterpret_cast<device_node3*>(G.neighbours)[threadIdx.x + blockDim.x*blockIdx.x];
    device_coord3d output_coordinates = reinterpret_cast<device_coord3d*>(G.X)[threadIdx.x + blockDim.x*blockIdx.x];
    device_coord2d output_xys = reinterpret_cast<device_coord2d*>(G.xys)[threadIdx.x + blockDim.x*blockIdx.x];
    cg:.sync(cg::this_grid());
    reinterpret_cast<node_t*>(G.neighbours)[(threadIdx.x + blockDim.x*blockIdx.x)*3] = output_neighbours.x; reinterpret_cast<node_t*>(G.neighbours)[(threadIdx.x + blockDim.x*blockIdx.x)*3 + 1] = output_neighbours.y; reinterpret_cast<node_t*>(G.neighbours)[(threadIdx.x + blockDim.x*blockIdx.x)*3 + 2] = output_neighbours.z;
    reinterpret_cast<double*>(G.X)[(threadIdx.x + blockDim.x*blockIdx.x)*3] = output_coordinates.x; reinterpret_cast<double*>(G.X)[(threadIdx.x + blockDim.x*blockIdx.x)*3 + 1] = output_coordinates.y; reinterpret_cast<double*>(G.X)[(threadIdx.x + blockDim.x*blockIdx.x)*3 + 2] = output_coordinates.z;
    reinterpret_cast<double*>(G.xys)[(threadIdx.x + blockDim.x*blockIdx.x)*2] = output_xys.x; reinterpret_cast<double*>(G.xys)[(threadIdx.x + blockDim.x*blockIdx.x)*2 + 1] = output_xys.y;
}

void IsomerspaceKernel::kernel_to_kernel_copy(IsomerspaceKernel& S, IsomerspaceKernel& D){
    for (size_t i = 0; i < S.get_device_count(); i++)
    {
        cudaSetDevice(i);
        void* kernelArgs[] = {(void*)&S.d_batch[i]};
        safeCudaKernelCall((void*)clear_convergence_status, dim3(S.get_device_capacity(i), 1, 1), dim3(S.get_isomer_size(), 1, 1), kernelArgs, 0);
    }
    D.batch_size = S.batch_size;
    D.batch_sizes = S.batch_sizes;
    cudaDeviceSynchronize();
    printLastCudaError("Kernel to kernel copy failed: ");

    for (size_t i = 0; i < S.get_device_count(); i++){cudaSetDevice(i); D.d_batch[i] <<= S.d_batch[i];}
    //for (size_t i = 0; i < S.get_device_count(); i++){cudaSetDevice(i); D.h_batch[i] <<= D.d_batch[i];}
    //for (size_t i = 0; i < S.get_device_capacity(0); i++){cout << D.h_batch[0].statuses[i] << ", " << std::flush;}
}   

void IsomerspaceKernel::output_isomer(size_t i, size_t idx){
    IsomerBatch B    = h_batch[i];
    size_t offset    = idx*3*N;
    size_t c_offset  = idx*2*N;
    neighbours_t output(N); std::vector<coord3d> output_X(N); std::vector<coord2d> xys(N);
    for (size_t j = 0; j < N; j++) {
        output[j] = std::vector<node_t>(B.neighbours + offset + j*3, B.neighbours + offset + j*3 + 3);
        xys[j]      = {reinterpret_cast<device_coord2d*>(B.xys + c_offset)[j].x,reinterpret_cast<device_coord2d*>(B.xys + c_offset)[j].y};
        output_X[j] = {B.X[offset + j*3], B.X[offset + j*3 + 1], B.X[offset + j*3 + 2]};
    }
    Polyhedron P = Polyhedron(FullereneGraph(Graph(output,true)), output_X);
    P.layout2d = xys;
    output_queue.push({B.IDs[idx],P});
    B.statuses[idx]==CONVERGED ? converged_count++ : failed_count++;
}

void IsomerspaceKernel::insert_isomer(size_t i, size_t idx){

    IsomerBatch B = h_batch[i];
    assert(!insert_queue.empty());
    size_t ID        = insert_queue.front().first;
    Polyhedron P     = insert_queue.front().second;
    size_t offset  = idx*3*N;

    for (device_node_t u = 0; u < N; u++){
        for (int j = 0; j < 3; j++){
            device_node_t v = P.neighbours[u][j];
            size_t arc_index = u*3 + j + offset;
            
            
            B.neighbours  [arc_index]    = v;
            B.X           [arc_index]    = !P.points.empty() ? P.points[u][j] : 0.0;
        }   
        B.xys         [u*2 + idx*2*N] = !P.layout2d.empty() ? P.layout2d[u].first : 0.0;
        B.xys         [u*2 + idx*2*N + 1] = !P.layout2d.empty() ? P.layout2d[u].second : 0.0;
    }

    B.iterations[idx]   = 0;
    B.statuses[idx]    = NOT_CONVERGED;
    B.IDs[idx]         = ID;
    
    batch_size++;
    batch_sizes[i]++;
    insert_queue.pop();
}

void IsomerspaceKernel::update_batch(){
    while (batch_size < batch_capacity && !insert_queue.empty()){
        for (size_t i = 0; i < this->device_count; i++){
        if (batch_sizes[i] < device_capacities[i]){
            IsomerBatch B = h_batch[i];
            size_t idx       = index_queue[i].front();
            size_t offset  = idx*3*N;      //neighbour offset

            if ((B.statuses[idx] == CONVERGED) || (B.statuses[idx]==FAILED))
            {
                output_isomer(i,idx);
            }

            insert_isomer(i, idx);
            index_queue[i].pop();
            break;
        }
        }
    }
    if (batch_size > 0) for(int i = 0; i < device_count; i++) d_batch[i] <<= h_batch[i];
    if (insert_queue.empty()){
        for (size_t i = 0; i < device_count; i++)
        for (size_t j = 0; j < device_capacities[i]; j++){   
            if ((h_batch[i].statuses[j] == CONVERGED) || (h_batch[i].statuses[j]==FAILED)){
                output_isomer(i,j);
                h_batch[i].statuses[j] = EMPTY;
            }
        }
    }
}

void IsomerspaceKernel::output_batch_to_queue(){
    for(size_t i = 0; i < device_count; i++) {cudaSetDevice(i); h_batch[i] <<= d_batch[i];}
    for (size_t i = 0; i < device_count; i++)
        for (size_t j = 0; j < device_capacities[i]; j++){   
            //cout << h_batch[i].statuses[j] << ", " << std::flush;
            if ((h_batch[i].statuses[j] == CONVERGED) || (h_batch[i].statuses[j]==FAILED)){
                output_isomer(i,j);
                h_batch[i].statuses[j] = EMPTY;
            }
        }
}

void IsomerspaceKernel::insert_queued_isomers(){
    for (size_t i = 0; i < device_count; i++)
    for (size_t j = 0; j < device_capacities[i]; j++) if(!insert_queue.empty()) insert_isomer(i,j);
    for(size_t i = 0; i < device_count; i++) {cudaSetDevice(i); d_batch[i] <<= h_batch[i];}
}