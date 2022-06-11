#pragma once
#include "fullerenes/gpu/isomerspace_forcefield.hh"
#include "fullerenes/gpu/isomerspace_tutte.hh"
#include "fullerenes/gpu/gpudatastruct.hh"
#include "fullerenes/gpu/isomerspace_kernel.hh"




IsomerBatch::IsomerBatch(size_t n_atoms, size_t n_isomers, BufferType buffer_type){
    this->buffer_type      = buffer_type;
    this->n_atoms          = n_atoms;
    this->isomer_capacity  = n_isomers;
    this->n_faces          = n_atoms/2 + 2;
    pointers =   {{"cubic_neighbours",(void**)&cubic_neighbours, sizeof(device_node_t)*n_atoms*3, true}, {"dual_neighbours", (void**)&dual_neighbours, sizeof(device_node_t) * (n_atoms/2 +2) * 6, true}, {"face_degrees", (void**)&face_degrees, sizeof(uint8_t)*(n_atoms/2 +2), true},{"X", (void**)&X, sizeof(device_real_t)*n_atoms*3, true}, {"xys", (void**)&xys, sizeof(device_real_t)*n_atoms*2, true}, {"statuses", (void**)&statuses, sizeof(IsomerStatus), false}, {"IDs", (void**)&IDs, sizeof(size_t), false}, {"iterations", (void**)&iterations, sizeof(size_t), false}};
    if (buffer_type == DEVICE_BUFFER){
    for (size_t i = 0; i < pointers.size(); i++) {
        cudaMalloc(get<1>(pointers[i]), n_isomers * get<2>(pointers[i])); 
        cudaMemset(*get<1>(pointers[i]),0,n_isomers*get<2>(pointers[i]));
    }
    } else if(buffer_type == HOST_BUFFER){
        for (size_t i = 0; i < pointers.size(); i++) {
            //For asynchronous memory transfers host memory must be pinned. 
            cudaMallocHost(get<1>(pointers[i]), n_isomers * get<2>(pointers[i]));
            memset(*get<1>(pointers[i]),0, n_isomers*get<2>(pointers[i]));
        }
    }
    printLastCudaError("Failed to construct IsomerBatch");
    allocated = true;
}

GPUDataStruct::GPUDataStruct(size_t n_atoms, size_t n_isomers, BufferType buffer_type){
    this->buffer_type = buffer_type;
    this->n_atoms     = n_atoms;
    this->n_isomers   = n_isomers;
    this->n_faces     = n_atoms/2 + 2;
    if (buffer_type == DEVICE_BUFFER){
        for (size_t i = 0; i < pointers.size(); i++) {
            size_t num_elements = get<3>(pointers[i]) ? n_isomers * n_atoms: n_isomers;
            cudaMalloc(get<1>(pointers[i]), num_elements* get<2>(pointers[i])); 
        }
    } else if(buffer_type == HOST_BUFFER){
        for (size_t i = 0; i < pointers.size(); i++) {
            size_t num_elements = get<3>(pointers[i]) ? n_isomers * n_atoms: n_isomers;
            //For asynchronous memory transfers host memory must be pinned. 
            cudaMallocHost(get<1>(pointers[i]), num_elements* get<2>(pointers[i]));
        }
    }
    printLastCudaError("Failed to construct GPUDataStruct");
    allocated = true;
    
}

GPUDataStruct::~GPUDataStruct(){
    if(allocated){
        if (buffer_type == DEVICE_BUFFER){    
            for (size_t i = 0; i < pointers.size(); i++) {
                cudaFree(*get<1>(pointers[i]));
            }
        } else{
            for (size_t i = 0; i < pointers.size(); i++) {
                cudaFreeHost(*get<1>(pointers[i])); 
            }
        }
        allocated = false;
    }
}

void GPUDataStruct::allocate(GPUDataStruct& G, const size_t n_atoms, const size_t n_isomers, const BufferType buffer_type){
    if((!G.allocated)){
        G.buffer_type = buffer_type;
        G.n_atoms = n_atoms; 
        G.isomer_capacity = n_isomers;
        if (buffer_type == DEVICE_BUFFER){
            for (size_t i = 0; i < G.pointers.size(); i++) {
                size_t num_elements = get<3>(G.pointers[i]) ? n_isomers * n_atoms: n_isomers;
                cudaMalloc(get<1>(G.pointers[i]), num_elements* get<2>(G.pointers[i])); 
            }
        }else{
            for (size_t i = 0; i < G.pointers.size(); i++) {
                size_t num_elements = get<3>(G.pointers[i]) ? n_isomers * n_atoms: n_isomers;
                //For asynchronous memory transfers host memory must be pinned. 
                cudaMallocHost(get<1>(G.pointers[i]), num_elements* get<2>(G.pointers[i]));
            }
        }        
        printLastCudaError("Failed to allocate struct");
        G.allocated = true;
    }
}

void GPUDataStruct::initialize(GPUDataStruct& G){
    if(G.buffer_type == DEVICE_BUFFER){
        for (size_t i = 0; i < G.pointers.size(); i++) {
            size_t num_elements = get<3>(G.pointers[i]) ? G.isomer_capacity * G.n_atoms: G.isomer_capacity;
            cudaMemset(*get<1>(G.pointers[i]),0,num_elements*get<2>(G.pointers[i]));
        }
    } else {
        for (size_t i = 0; i < G.pointers.size(); i++) {
            size_t num_elements = get<3>(G.pointers[i]) ? G.isomer_capacity * G.n_atoms: G.isomer_capacity;
            memset(*get<1>(G.pointers[i]),0, num_elements*get<2>(G.pointers[i]));
        }
    }
}

void GPUDataStruct::copy(GPUDataStruct& destination, const GPUDataStruct& source, const cudaStream_t stream){
    for (size_t i = 0; i < source.pointers.size(); i++)
    {
        size_t num_elements = get<3>(source.pointers[i]) ?  source.n_atoms * source.isomer_capacity : source.isomer_capacity;
        cudaMemcpyAsync(*(get<1>(destination.pointers[i])) , *(get<1>(source.pointers[i])), get<2>(source.pointers[i])*num_elements, cudaMemcpyKind(2*source.buffer_type +  destination.buffer_type), stream);
    }
    destination.n_isomers = source.n_isomers;
    printLastCudaError("Failed to copy struct");
}

void operator <<= (GPUDataStruct& destination, const GPUDataStruct& source){
    GPUDataStruct::copy(destination, source);
} 

__device__ bool within_interval(size_t x, size_t a, size_t b){
    if (a < b) { return x >= a && x < b;} 
    else {return x >= a || x < b;}
    
}

__global__
void clear_convergence_status(IsomerBatch G){
    if (threadIdx.x == 0) G.statuses[blockIdx.x] = G.statuses[blockIdx.x] != EMPTY ? NOT_CONVERGED : EMPTY;
    if (threadIdx.x == 0) G.iterations[blockIdx.x] = 0;
}

__global__
void input_type_conversion(IsomerBatch G){ // TODO: Write nicer.
    node_t input_neighbours[3]  = {reinterpret_cast<node_t*>(G.cubic_neighbours)[(threadIdx.x* + blockDim.x*blockIdx.x)*3],   reinterpret_cast<node_t*>(G.cubic_neighbours) [(threadIdx.x* + blockDim.x*blockIdx.x)*3 + 1], reinterpret_cast<node_t*>(G.cubic_neighbours) [(threadIdx.x* + blockDim.x*blockIdx.x)*3 + 2]};
    double input_coordianates[3]= {reinterpret_cast<double*>(G.X)[(threadIdx.x* + blockDim.x*blockIdx.x)*3],            reinterpret_cast<double*>(G.X)          [(threadIdx.x* + blockDim.x*blockIdx.x)*3 + 1], reinterpret_cast<double*>(G.X)          [(threadIdx.x* + blockDim.x*blockIdx.x)*3 + 2]};
    double input_xys[2]         = {reinterpret_cast<double*>(G.xys)[(threadIdx.x* + blockDim.x*blockIdx.x)*2],          reinterpret_cast<double*>(G.xys)        [(threadIdx.x* + blockDim.x*blockIdx.x)*2 + 1]};
    cg::sync(cg::this_grid());
    size_t index = threadIdx.x + blockDim.x*blockIdx.x;
    reinterpret_cast<device_node3*>(G.cubic_neighbours)[index] = {device_node_t(input_neighbours[0]), device_node_t(input_neighbours[1]), device_node_t(input_neighbours[2])};
    reinterpret_cast<device_coord3d*>(G.X)       [index] = {device_real_t(input_coordianates[0]), device_real_t(input_coordianates[1]), device_real_t(input_coordianates[2])};
    reinterpret_cast<device_coord2d*>(G.xys)     [index] = {device_real_t(input_xys[0]), device_real_t(input_xys[1])}; 
}

__global__
void output_type_conversion(IsomerBatch G){
    device_node3 output_neighbours      = reinterpret_cast<device_node3*>  (G.cubic_neighbours)   [threadIdx.x + blockDim.x*blockIdx.x];
    device_coord3d output_coordinates   = reinterpret_cast<device_coord3d*>(G.X)            [threadIdx.x + blockDim.x*blockIdx.x];
    device_coord2d output_xys           = reinterpret_cast<device_coord2d*>(G.xys)          [threadIdx.x + blockDim.x*blockIdx.x];
    cg::sync(cg::this_grid());
    reinterpret_cast<node_t*>(G.cubic_neighbours) [(threadIdx.x + blockDim.x*blockIdx.x)*3] = output_neighbours.x;    reinterpret_cast<node_t*>(G.cubic_neighbours) [(threadIdx.x + blockDim.x*blockIdx.x)*3 + 1] = output_neighbours.y;    reinterpret_cast<node_t*>(G.cubic_neighbours) [(threadIdx.x + blockDim.x*blockIdx.x)*3 + 2] = output_neighbours.z;
    reinterpret_cast<double*>(G.X)          [(threadIdx.x + blockDim.x*blockIdx.x)*3] = output_coordinates.x;   reinterpret_cast<double*>(G.X)          [(threadIdx.x + blockDim.x*blockIdx.x)*3 + 1] = output_coordinates.y;   reinterpret_cast<double*>(G.X)          [(threadIdx.x + blockDim.x*blockIdx.x)*3 + 2] = output_coordinates.z;
    reinterpret_cast<double*>(G.xys)        [(threadIdx.x + blockDim.x*blockIdx.x)*2] = output_xys.x;           reinterpret_cast<double*>(G.xys)        [(threadIdx.x + blockDim.x*blockIdx.x)*2 + 1] = output_xys.y;
}
__device__ unsigned int queue_front = 0, queue_back = 0, queue_capacity = 0, queue_size = 0, num_queue_requests = 0;

__global__
void kernel_update_queue(IsomerBatch G, IsomerBatch queue_batch){
    __shared__ size_t queue_index;
    if ((threadIdx.x + blockIdx.x) == 0) {num_queue_requests = 0; queue_capacity = queue_batch.isomer_capacity;}
    GRID_SYNC
    if (G.statuses[blockIdx.x] != NOT_CONVERGED){
        if(threadIdx.x == 0){
            queue_index = atomicAdd(&queue_front, 1);
            queue_index = queue_index % queue_capacity;
            atomicAdd(&num_queue_requests, 1);
        } 
        cg::sync(cg::this_thread_block());
        assert(queue_index < queue_capacity);
        size_t queue_array_idx    = queue_index*blockDim.x+threadIdx.x;
        size_t global_idx         = blockDim.x*blockIdx.x + threadIdx.x;
        reinterpret_cast<device_coord2d*>(G.xys)[global_idx]         = reinterpret_cast<device_coord2d*>(queue_batch.xys)[queue_array_idx];  
        reinterpret_cast<device_coord3d*>(G.X)[global_idx]           = reinterpret_cast<device_coord3d*>(queue_batch.X)[queue_array_idx];  
        reinterpret_cast<device_node3*>(G.cubic_neighbours)[global_idx]    = reinterpret_cast<device_node3*>(queue_batch.cubic_neighbours)[queue_array_idx];  
        if (threadIdx.x == 0){
            G.IDs[blockIdx.x] = queue_batch.IDs[queue_index];
            G.iterations[blockIdx.x] = 0;
            G.statuses[blockIdx.x] = queue_batch.statuses[queue_index];
            queue_batch.statuses[queue_index] = EMPTY;
        }
    }
    GRID_SYNC
    if ((threadIdx.x + blockIdx.x) == 0) {
        //Main thread checks if the number of requests for new isomers is greater than the queue size, then either sets queuesize to 0 by subtracting queue size or subtracts the number of requests.
        bool enough_left_in_queue = queue_size >= num_queue_requests;
        queue_size -= enough_left_in_queue ? num_queue_requests : queue_size;
        queue_front = enough_left_in_queue ? (queue_back - queue_size + queue_capacity) % queue_capacity : queue_back;
    }
}

__global__
void kernel_initialize_queue(IsomerBatch queue_batch){
    if ((threadIdx.x + blockIdx.x) == 0)
    {
        queue_size = 0; queue_back = 0; queue_capacity = queue_batch.isomer_capacity; queue_front = 0;
    }
    
    for (size_t tid = blockDim.x * blockIdx.x + threadIdx.x; tid < queue_batch.isomer_capacity*queue_batch.n_atoms; tid+=blockDim.x*gridDim.x)
    {
        reinterpret_cast<device_coord3d*>(queue_batch.X)[tid]           = (device_coord3d){0.0,0.0,0.0};
        reinterpret_cast<device_node3*>(queue_batch.cubic_neighbours)[tid]    = (device_node3){0,0,0};
        reinterpret_cast<device_coord2d*>(queue_batch.xys)[tid]         = (device_coord2d){0.0,0.0};
    }
    for (size_t tid = blockDim.x * blockIdx.x + threadIdx.x; tid < queue_batch.isomer_capacity; tid+=blockDim.x*gridDim.x)
    {
        queue_batch.statuses[tid] = EMPTY;
        queue_batch.iterations[tid] = 0;
        queue_batch.IDs[tid] = 0;
    }
}

__global__
void kernel_push_batch(IsomerBatch input_batch, IsomerBatch queue_batch, device_real_t* global_reduction_array){
    extern __shared__ device_real_t sdata[];
    clear_cache(sdata, Block_Size_Pow_2 );
    size_t insert_size = size_t(global_reduction(sdata,global_reduction_array, device_real_t((input_batch.statuses[blockIdx.x] == NOT_CONVERGED) && (threadIdx.x == 0))));
    //Since queue is statically allocated, it needs to break if you try to exceed capacity.
    assert((queue_size + insert_size) <= queue_capacity);

    size_t queue_idx = (queue_back + blockIdx.x) % queue_capacity;
    reinterpret_cast<device_coord2d*>(queue_batch.xys)[queue_idx*blockDim.x + threadIdx.x]      = reinterpret_cast<device_coord2d*>(input_batch.xys)[blockIdx.x*blockDim.x + threadIdx.x];
    reinterpret_cast<device_coord3d*>(queue_batch.X)[queue_idx*blockDim.x + threadIdx.x]        = reinterpret_cast<device_coord3d*>(input_batch.X)[blockIdx.x*blockDim.x + threadIdx.x];
    reinterpret_cast<device_node3*>(queue_batch.cubic_neighbours)[queue_idx*blockDim.x + threadIdx.x] = reinterpret_cast<device_node3*>(input_batch.cubic_neighbours)[blockIdx.x*blockDim.x + threadIdx.x];
    if (threadIdx.x == 0)
    {
        queue_batch.IDs[queue_idx]          = input_batch.IDs[blockIdx.x];
        queue_batch.iterations[queue_idx]   = 0;
        queue_batch.statuses[queue_idx]     = input_batch.statuses[blockIdx.x];
    }
    GRID_SYNC
    if ((threadIdx.x + blockIdx.x) == 0) {
        queue_size += insert_size; queue_back = (queue_back + insert_size) % queue_capacity;
    }
}

void IsomerspaceForcefield::push_batch_from_kernel(const IsomerspaceKernel& input_kernel){
    synchronize();
    for (size_t i = 0; i < device_count; i++)
    {
        void* kernel_args_1[] = {(void*)&input_kernel.d_batch[i]};
        void* kernel_args_2[] = {(void*)&input_kernel.d_batch[i], &d_queues[i], &global_reduction_arrays[i]};
        safeCudaKernelCall((void*)clear_convergence_status, dim3(input_kernel.get_device_capacity(i), 1, 1), dim3(N, 1, 1), kernel_args_1, 0);
        safeCudaKernelCall((void*)kernel_push_batch,dim3(input_kernel.get_device_capacity(i), 1, 1), dim3(N,1,1), kernel_args_2, sizeof(device_real_t)*(N + Block_Size_Pow_2));
    }
    printLastCudaError("Pushing batch to queue Failed: ");
    for(size_t i = 0; i < device_count; i++) device_queue_size += input_kernel.d_batch[i].n_isomers;
    synchronize();
}

void IsomerspaceForcefield::update_device_batches(){
    synchronize();
    for (size_t i = 0; i < device_count; i++)
    {
        cudaStreamSynchronize(copy_stream[i]);
        void* kernelArgs[] = {(void*)&d_batch[i], &d_queues[i]};
        safeCudaKernelCall((void*)kernel_update_queue, dim3(device_capacities[i], 1, 1), dim3(N,1,1), kernelArgs, sizeof(size_t));
    }
    synchronize();
    printLastCudaError("Update device batch from queue Failed: ");
}

void IsomerspaceForcefield::move_to_output_buffer(){
    for (size_t i = 0; i < device_count; i++)
    {
        //Wait for GPU kernels to finish before moving data to second buffer.
        cudaStreamSynchronize(main_stream[i]);
        //Double buffering to allow for asynchronous GPU -> CPU copying, host data reodering and GPU execution.
        cudaSetDevice(i); GPUDataStruct::copy(d_output_batch[i], d_batch[i], copy_stream[i]);
    }
    printLastCudaError("Moving data to output buffer failed: ");
}



void IsomerspaceKernel::kernel_to_kernel_copy(IsomerspaceKernel& S, IsomerspaceKernel& D){
    //Wait for all previous kernel operations to finish before copying batches between kernels.
    S.synchronize(); D.synchronize();
    for (size_t i = 0; i < S.get_device_count(); i++)
    {
        cudaSetDevice(i);
        void* kernelArgs[] = {(void*)&S.d_batch[i]};
        //Convergence is kernel dependent, i.e. X0 may have converged but we need to clear this when we transfer the batch to the forcefield kernel.
        safeCudaKernelCall((void*)clear_convergence_status, dim3(S.get_device_capacity(i), 1, 1), dim3(S.get_isomer_size(), 1, 1), kernelArgs, 0);
    }
    D.batch_size = S.batch_size;
    D.batch_sizes = S.batch_sizes;
    printLastCudaError("Kernel to kernel copy failed: ");

    for (size_t i = 0; i < S.get_device_count(); i++){cudaSetDevice(i); D.d_batch[i] <<= S.d_batch[i];}
    //Wait for transfer to finish
    S.synchronize(); D.synchronize();
}

void IsomerspaceKernel::copy_metadata(){
    for (size_t i = 0; i < device_count; i++)
    {
        cudaSetDevice(i);
        cudaMemcpy(d_batch[i].statuses,     h_batch[i].statuses,    sizeof(IsomerStatus)* device_capacities[i], cudaMemcpyHostToDevice);
        cudaMemcpy(d_batch[i].iterations,   h_batch[i].iterations,  sizeof(size_t)* device_capacities[i],       cudaMemcpyHostToDevice);
        cudaMemcpy(d_batch[i].IDs,          h_batch[i].IDs,         sizeof(size_t)* device_capacities[i],       cudaMemcpyHostToDevice);
    }
}

void IsomerspaceKernel::output_isomer(size_t i, size_t idx){
    IsomerBatch& B   = h_batch[i];
    size_t offset    = idx*3*N;
    size_t c_offset  = idx*2*N;
    neighbours_t output(N); std::vector<coord3d> output_X(N); std::vector<coord2d> xys(N);

    for (size_t j = 0; j < N; j++) {
        output[j]   = std::vector<node_t>(B.cubic_neighbours + offset + j*3, B.cubic_neighbours + offset + j*3 + 3);
        xys[j]      = {reinterpret_cast<device_coord2d*>(B.xys + c_offset)[j].x,reinterpret_cast<device_coord2d*>(B.xys + c_offset)[j].y};
        output_X[j] = {B.X[offset + j*3], B.X[offset + j*3 + 1], B.X[offset + j*3 + 2]};
    }
    Polyhedron P(FullereneGraph(Graph(output,true)), output_X);
    P.layout2d = xys;
    output_queue.push({B.IDs[idx],P});
    B.statuses[idx]==CONVERGED ? converged_count++ : failed_count++;
    B.n_isomers--; if(queue_mode == DEVICE_QUEUE) device_queue_size--;
}

void IsomerspaceKernel::insert_isomer(size_t i, size_t idx){

    IsomerBatch& B = h_batch[i];
    assert(!insert_queue.empty());
    cudaSetDevice(i);
    size_t ID        = insert_queue.front().first;
    Polyhedron &P     = insert_queue.front().second;

    bool is_points_empty    = P.points.empty();
    bool is_layout2d_empty  = P.layout2d.empty();

    size_t offset  = idx*3*N;
    
    for (device_node_t u = 0; u < N; u++){
        for (int j = 0; j < 3; j++){
            device_node_t v = P.neighbours[u][j];
            size_t arc_index = u*3 + j + offset;
            
            B.cubic_neighbours  [arc_index]    = v;
            if(!is_points_empty){ B.X           [arc_index] = P.points[u][j];}
        }   
        if(!is_layout2d_empty){
        B.xys         [u*2 + idx*2*N] = P.layout2d[u].first;
        B.xys         [u*2 + idx*2*N + 1] = P.layout2d[u].second;
        }
    }
    
    B.iterations[idx]   = 0;
    B.statuses[idx]    = NOT_CONVERGED;
    B.IDs[idx]         = ID;
    
    batch_size++;
    batch_sizes[i]++;
    B.n_isomers++;
    insert_queue.pop();
}

void IsomerspaceKernel::update_batch(){
    while (batch_size < batch_capacity && !insert_queue.empty()){
        for (size_t i = 0; i < this->device_count; i++){
        if (batch_sizes[i] < device_capacities[i]){
            IsomerBatch B = h_batch[i];
            size_t idx       = index_queue[i].front();
	    //            size_t offset  = idx*3*N;      //neighbour offset; TODO: DENNE BRUGES ALDRIG. Er det med vilje?

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
    for(size_t i = 0; i < device_count; i++) d_batch[i] <<= h_batch[i];
}

void IsomerspaceKernel::set_all_converged(){
    for (size_t i = 0; i < device_count; i++)
    {
        cudaSetDevice(i);

        cudaMemcpy(h_batch[i].statuses, d_batch[i].statuses, sizeof(IsomerStatus)*device_capacities[i], cudaMemcpyDeviceToHost);
        for (size_t j = 0; j < device_capacities[i]; j++)
        {
            h_batch[i].statuses[j] = (h_batch[i].statuses[j] != EMPTY)  ? CONVERGED : EMPTY;
        }
        cudaMemcpy(d_batch[i].statuses, h_batch[i].statuses, sizeof(IsomerStatus)*device_capacities[i], cudaMemcpyHostToDevice);
    }
    
}

void IsomerspaceKernel::output_batch_to_queue(){
    
    for(size_t i = 0; i < device_count; i++) {
        cudaSetDevice(i); 
        //This copy operation depends on output residing in the output batch
        GPUDataStruct::copy(h_batch[i],d_output_batch[i],copy_stream[i]); 
        //We need to wait for all queued copy stream operations to finish before reading from the host memory. 
        cudaStreamSynchronize(copy_stream[i]);
    }
    
    for (size_t i = 0; i < device_count; i++)
    for (size_t j = 0; j < device_capacities[i]; j++){   
        if ((h_batch[i].statuses[j] == CONVERGED) || (h_batch[i].statuses[j]==FAILED)){
            output_isomer(i,j);
            h_batch[i].statuses[j] = EMPTY;
        }
    }
}

void IsomerspaceKernel::insert_queued_isomers(){
    batch_size = 0;
    batch_sizes = {};
    for (size_t i = 0; i < device_count; i++) {
        h_batch[i].n_isomers = 0;
        for (size_t j = 0; j < device_capacities[i]; j++) if(!insert_queue.empty()) insert_isomer(i,j); else h_batch[i].statuses[j] = EMPTY;
            printLastCudaError("Failed to insert isomers: ");
    }
    for (size_t i = 0; i < device_count; i++) d_batch[i] <<= h_batch[i];
}
