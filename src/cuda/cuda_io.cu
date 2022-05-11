#include "fullerenes/gpu/cuda_io.hh"
#include "type_traits"

namespace cuda_io{

    cudaError_t output_to_queue(std::queue<std::tuple<Polyhedron, size_t, IsomerStatus>>& queue, IsomerBatch& batch, const bool copy_2d_layout){
        //Batch needs to exist on the host. For performance reasons we don't want to create a new batch here and copy to that, cudaMalloc is expensive.
        printLastCudaError();
        if (batch.buffer_type != HOST_BUFFER) assert(false); 

        size_t N = batch.n_atoms;
        for (size_t isomer_idx = 0; isomer_idx < batch.isomer_capacity; isomer_idx++){   
            //Only insert the isomer if it has finished (either CONVERGED or FAILED)
            if(!(batch.statuses[isomer_idx] == CONVERGED || batch.statuses[isomer_idx] == FAILED)) continue;
            
            //Graphs always have a neighbour array.
            neighbours_t out_neighbours(N);
 
            std::vector<coord2d> output_2D;

            for (size_t i = 0; i < N; i++){
                //Fill in cubic neighbours
                out_neighbours[i] = {batch.neighbours[isomer_idx*N*3 + i*3], batch.neighbours[isomer_idx*N*3 + i*3 + 1], batch.neighbours[isomer_idx*N*3 + i*3 + 2]};
                
            }

            std::vector<coord3d> output_X(N);
            for (size_t i = 0; i < N; i++){
                for (int j = 0; j < 3; ++j){
                    output_X[i][j] = batch.X[isomer_idx*N*3 + i*3 + j];
                }
            }

            //If 2D layout is true, allocate memory and copy 2D layout. If this is not needed disable it for performance gains.
            if (copy_2d_layout){
                output_2D = std::vector<coord2d>(N);
                for (size_t i = 0; i < N; i++){
                    output_2D[i] = {batch.xys[isomer_idx*N*2 + i*2], batch.xys[isomer_idx*N*2 + i*2 + 1]};
                }
                Polyhedron P(FullereneGraph(Graph(out_neighbours, true)),output_X);
                P.layout2d = output_2D;
                queue.push({P,batch.IDs[isomer_idx], batch.statuses[isomer_idx]});    
            }else{
                queue.push({Polyhedron(Graph(out_neighbours, true),output_X),batch.IDs[isomer_idx], batch.statuses[isomer_idx]});
            }
        }
        return cudaGetLastError();
    }


    cudaError_t copy(   IsomerBatch& destination, //Copy data to this batch
                        const IsomerBatch& source, //Copy data from this batch
                        const LaunchCtx& ctx, //Optional: specify which launch context to perform the copy operation in.
                        const LaunchPolicy policy, //Optional: specifies whether to synchronize the stream before and after copying)
                        const std::pair<int,int>& lhs_range, //Optional: provide a range of indices to assign similar to slices in numpy eg. {0,5} = [0:5]
                        const std::pair<int,int>& rhs_range //Optional: provide a range of indices to copy from similar to slices in numpy eg. {0,5} = [0:5]
                        ){
        //Iterate over the data fields of the IsomerBatch (pseudo reflection) and copy the contents of each using the provided stream.
        if(policy == LaunchPolicy::SYNC) {ctx.wait();}
        for (size_t i = 0; i < source.pointers.size(); i++)
        {

            int num_isomers = (lhs_range.second > -1 && lhs_range.first > -1) ? lhs_range.second - lhs_range.first : destination.isomer_capacity;
            int num_elements = get<3>(source.pointers[i]) ?  source.n_atoms * num_isomers : num_isomers;
            int lhs_offset = lhs_range.first > 0 ? lhs_range.first * max((size_t)1, source.n_atoms*get<3>(source.pointers[i])) * get<2>(source.pointers[i]) : 0;
            int rhs_offset = rhs_range.first > 0 ? rhs_range.first * max((size_t)1, source.n_atoms*get<3>(source.pointers[i])) * get<2>(source.pointers[i]) : 0;
            char* lhs_ptr = (char*)(*(get<1>(destination.pointers[i])));
            char* rhs_ptr = (char*)(*(get<1>(source.pointers[i])));
            cudaMemcpyAsync(lhs_ptr + lhs_offset, rhs_ptr + rhs_offset, get<2>(source.pointers[i])*num_elements, cudaMemcpyKind(2*source.buffer_type +  destination.buffer_type), ctx.stream);
        }
        destination.n_isomers = source.n_isomers;
        printLastCudaError("Failed to copy struct");
        if(policy == LaunchPolicy::SYNC) {ctx.wait();}
        return cudaGetLastError();
    }

    cudaError_t free(IsomerBatch& batch){
        for (int i = 0; i < batch.pointers.size(); i++){
            if(batch.buffer_type == DEVICE_BUFFER) {cudaFree(*get<1>(batch.pointers[i]));} else {cudaFreeHost(*get<1>(batch.pointers[i]));}
        }
        return cudaGetLastError();
    }
    
    cudaError_t resize(IsomerBatch& batch, const size_t new_capacity, const LaunchCtx& ctx, const LaunchPolicy policy, int front, int back){
        //Construct a tempory batch: allocates the needed amount of memory.
        IsomerBatch temp_batch = IsomerBatch(batch.n_atoms, new_capacity, batch.buffer_type);
        //Copy contents of old batch into newly allocated memory.
        IsomerBatch::copy(temp_batch, batch, ctx.stream);
        for (int i = 0; i < batch.pointers.size(); i++)
        {
            void* temp_ptr = *get<1>(batch.pointers[i]);
            printLastCudaError("Free failed");
            //Reassign pointers of the input batch, to the new memory
            *get<1>(batch.pointers[i]) = *get<1>(temp_batch.pointers[i]);
            //Assign old pointers to temporary object, let destructor take care of cleanup.
            *get<1>(temp_batch.pointers[i]) = temp_ptr;
        }
        batch.isomer_capacity = temp_batch.isomer_capacity;
        return cudaGetLastError();
    }

}


bool IsomerBatch::operator==(const IsomerBatch& b){
    bool passed = true;
    
    if (! (buffer_type == b.buffer_type && isomer_capacity == b.isomer_capacity && n_atoms == b.n_atoms) ) {
        return false;
    }else if(buffer_type == HOST_BUFFER){
        for(int i = 0; i < isomer_capacity; ++i){
            passed &= statuses[i] == b.statuses[i];
            passed &= IDs[i] == b.IDs[i];
            passed &= iterations[i] = b.iterations[i];
        }
        for(int i = 0; i < isomer_capacity * n_atoms * 3; ++i){
            passed &= X[i] == b.X[i];
            passed &= xys[i] == b.xys[i];
            passed &= neighbours[i] == b.neighbours[i];
        }
        return passed;
    } else{
        std::cout << "== operator only supported for HOST_BUFFER" << std::endl;
        return false;
    }
    return false;
}



std::ostream& operator << (std::ostream& os, const IsomerBatch& a){
    os << "Batch Dimensions: (" << a.isomer_capacity << ", " << a.n_atoms << ")\n";
    os << "ID List: " << "\n[";
    for (size_t i = 0; i < a.isomer_capacity - 1; i++){os << a.IDs[i] << ", ";}
    os << a.IDs[a.isomer_capacity-1] << "]\n"; 
    
    os << "Status List: " << "\n[";
    for (size_t i = 0; i < a.isomer_capacity - 1; i++){os << a.statuses[i] << ", ";}
    os << a.statuses[a.isomer_capacity-1] << "]\n"; 

    os << "X Lists: " << "\n[";
    for (size_t i = 0; i < a.isomer_capacity; i++){
        os << "[";
        for (size_t j = 0; j < a.n_atoms*3 - 1; j++)
        {
            os << a.X[i*a.n_atoms*3 + j] << ", ";
        }
        if(i != (a.isomer_capacity - 1)) {
            os << a.X[(i+1)*a.n_atoms*3 - 1] << "], ";
        } else{
            os << a.X[(i+1)*a.n_atoms*3 - 1] << "]]\n";
        }
    }

    os << "2D Layouts: " << "\n[";
    for (size_t i = 0; i < a.isomer_capacity; i++){
        os << "[";
        for (size_t j = 0; j < a.n_atoms*2 - 1; j++)
        {
            os << a.xys[i*a.n_atoms*2 + j] << ", ";
        }
        if(i != (a.isomer_capacity - 1)) {
            os << a.xys[(i+1)*a.n_atoms*2 - 1] << "], ";
        } else{
            os << a.xys[(i+1)*a.n_atoms*2 - 1] << "]]\n";
        }
    }

    os << "Neighbour Lists: " << "\n[";
    for (size_t i = 0; i < a.isomer_capacity; i++){
        os << "[";
        for (size_t j = 0; j < a.n_atoms*3 - 1; j++)
        {
            os << a.neighbours[i*a.n_atoms*3 + j] << ", ";
        }
        if(i != (a.isomer_capacity - 1)) {
            os << a.neighbours[(i+1)*a.n_atoms*3 - 1] << "], ";
        } else{
            os << a.neighbours[(i+1)*a.n_atoms*3 - 1] << "]]\n";
        }
    }


    return os;
}

 