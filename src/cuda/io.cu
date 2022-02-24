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
    for (size_t i = 0; i < destination.pointers.size(); i++)
    {
        size_t num_elements = get<3>(destination.pointers[i]) ?  destination.N*destination.batch_size : destination.batch_size;
        cudaMemcpy(*(get<1>(destination.pointers[i])) , *(get<1>(source.pointers[i])), get<2>(source.pointers[i])*num_elements, cudaMemcpyKind(2*source.buffer_type +  destination.buffer_type));
    }
    }
    else{
        std::cout << "WARNING: Call to copy made for 0 isomers \n";
    }
    printLastCudaError("Failed to copy struct");
}

void operator <<= (GPUDataStruct& destination, const GPUDataStruct& source){
    GPUDataStruct::copy(destination, source);
}


void operator <<= (IsomerspaceForcefield::IsomerBatch& a, const IsomerspaceForcefield::IsomerBatch& b){
    GPUDataStruct::copy(a,b);
    GPUDataStruct::copy(a.stats,b.stats);
}

void operator <<= (IsomerBatch& a, const IsomerBatch& b){
    GPUDataStruct::copy(a,b);
}

void IsomerspaceTutte::eject_isomer(size_t i, size_t idx){
    size_t n_offset  = idx*3*N;      //neighbour offset
    size_t c_offset  = idx*2*N;      //coords offset
    size_t f_offset  = idx*N;    //outer_face offset
    IsomerBatch B = h_batch[i];
    neighbours_t neighbours(N); std::vector<coord2d> xys(N);
    for (size_t j = 0; j < N; j++) {
        neighbours[j]       = std::vector<node_t>(B.neighbours + n_offset + j*3, B.neighbours + n_offset + j*3 + 3);
        xys[j]   = {reinterpret_cast<GPU_REAL2*>(B.xys + c_offset)[j].x,reinterpret_cast<GPU_REAL2*>(B.xys + c_offset)[j].y};
                            
    }
    FullereneGraph P = FullereneGraph(Graph(neighbours,true));
    P.layout2d       = xys;
    output_queue.push({B.IDs[idx],P});
    B.statuses[idx]==CONVERGED ? converged_count++ : failed_count++;
}

void IsomerspaceForcefield::eject_isomer(size_t i, size_t idx){
    IsomerBatchStats stats = h_batch[i].stats;
    size_t offset   = idx*3*N;
    neighbours_t output(N); std::vector<coord3d> output_X(N);
    for (size_t j = 0; j < N; j++) {
        output[j] = std::vector<node_t>(h_batch[i].neighbours + offset + j*3, h_batch[i].neighbours + offset + j*3 + 3);
        output_X[j] = {h_batch[i].X[offset + j*3], h_batch[i].X[offset + j*3 + 1], h_batch[i].X[offset + j*3 + 2]};
    }
    Polyhedron P = Polyhedron(FullereneGraph(Graph(output,true)), output_X);
    output_queue.push({stats.isomer_IDs[idx],P});
    stats.isomer_statuses[idx]==CONVERGED ? converged_count++ : failed_count++;
}

void IsomerspaceTutte::update_batch(){
    while (batch_size < batch_capacity && !insert_queue.empty()){
        for (size_t i = 0; i < this->device_count; i++)
        if (batch_sizes[i] < device_capacities[i]){
            IsomerBatch B = h_batch[i];
            size_t idx       = index_queue[i].front();
            size_t n_offset  = idx*3*N;      //neighbour offset
            size_t f_offset  = idx*N;    //outer_face offset
            if ((B.statuses[idx] == CONVERGED) || (B.statuses[idx]==FAILED))
            {
                eject_isomer(i,idx);
            }

            size_t ID        = insert_queue.front().first;
            FullereneGraph P = insert_queue.front().second;

            for (device_node_t u = 0; u < N; u++){
                for (int j = 0; j < 3; j++) B.neighbours[u*3 + j + n_offset]     = P.neighbours[u][j];
            }

            B.iterations[idx]   = 0;
            B.statuses[idx]    = NOT_CONVERGED;
            B.IDs[idx]         = ID;
            
            batch_size++;
            batch_sizes[i]++;
            insert_queue.pop();
            index_queue[i].pop();
            break;
    }

    }
    if (insert_queue.empty()){
        for (size_t i = 0; i < device_count; i++)
        for (size_t j = 0; j < device_capacities[i]; j++){   
            if ((h_batch[i].statuses[j] == CONVERGED) || (h_batch[i].statuses[j]==FAILED)){
                eject_isomer(i,j);
                h_batch[i].statuses[j] = EMPTY;
            }
        }
        
    }
}

void IsomerspaceForcefield::update_batch(){
    while (batch_size < batch_capacity && !insert_queue.empty()){
        for (size_t i = 0; i < this->device_count; i++)
        if (batch_sizes[i] < device_capacities[i]){
            IsomerBatchStats stats = h_batch[i].stats;
            size_t idx      = index_queue[i].front();
            size_t offset   = idx*3*N;
            if ((stats.isomer_statuses[idx] == CONVERGED) || (stats.isomer_statuses[idx]==FAILED))
            {
                eject_isomer(i,idx);
            }
            
            size_t ID       = insert_queue.front().first;
            Polyhedron P    = insert_queue.front().second;

            for (device_node_t u = 0; u < N; u++){
                for (int j = 0; j < 3; j++){
                    device_node_t v = P.neighbours[u][j];
                    size_t arc_index = u*3 + j + offset;
                    h_batch[i].neighbours  [arc_index] = v;
                    h_batch[i].next_on_face[arc_index] = P.next_on_face(u,v);
                    h_batch[i].prev_on_face[arc_index] = P.prev_on_face(u,v);
                    h_batch[i].face_right  [arc_index] = P.face_size(u,v);
                    h_batch[i].X           [arc_index] = P.points[u][j];
                }   
            }

            stats.iteration_counts[idx]   = 0;
            stats.isomer_statuses[idx]    = NOT_CONVERGED;
            stats.isomer_IDs[idx]         = ID;
            
            batch_size++;
            batch_sizes[i]++;
            insert_queue.pop();
            index_queue[i].pop();
            break;
        }
    }
    if (insert_queue.empty()){
        for (size_t i = 0; i < device_count; i++){
        IsomerBatchStats stats = h_batch[i].stats;
        for (size_t j = 0; j < device_capacities[i]; j++){   
            if ((stats.isomer_statuses[j] == CONVERGED) || (stats.isomer_statuses[j]==FAILED)){
                eject_isomer(i,j);
                stats.isomer_statuses[j] = EMPTY;
            }
        }
        }
    }
}

