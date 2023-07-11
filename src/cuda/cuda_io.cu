#include "fullerenes/gpu/cuda_io.hh"
#include "fullerenes/gpu/misc_cuda.cuh"
#include "type_traits"
#include <numeric>
namespace cuda_io{
    
    template cudaError_t output_to_queue<HOST_BUFFER>(std::queue<std::tuple<Polyhedron, size_t, IsomerStatus>>& queue, IsomerBatch<HOST_BUFFER>& batch, const bool copy_2d_layout);
    template cudaError_t reset_convergence_statuses<DEVICE_BUFFER>(IsomerBatch<DEVICE_BUFFER>& B, const LaunchCtx& ctx, const LaunchPolicy policy);
    template double average_iterations<HOST_BUFFER>(const IsomerBatch<HOST_BUFFER>& input);
    template double average_iterations<DEVICE_BUFFER>(const IsomerBatch<DEVICE_BUFFER>& input);
    template void is_close<HOST_BUFFER>(const IsomerBatch<HOST_BUFFER>& a, const IsomerBatch<HOST_BUFFER>& b, float tol, bool verbose);
    template void sort<HOST_BUFFER>(IsomerBatch<HOST_BUFFER>& input_batch, const BatchMember key, const SortOrder order);
    template void sort<DEVICE_BUFFER>(IsomerBatch<DEVICE_BUFFER>& input_batch, const BatchMember key, const SortOrder order);
    template int count_batch_status<HOST_BUFFER>(const IsomerBatch<HOST_BUFFER>& input, const IsomerStatus status);
    template int count_batch_status<DEVICE_BUFFER>(const IsomerBatch<DEVICE_BUFFER>& input, const IsomerStatus status);


    template <BufferType T>
    cudaError_t output_to_queue(std::queue<std::tuple<Polyhedron, size_t, IsomerStatus>>& queue, IsomerBatch<T>& batch, const bool copy_2d_layout){
        //Batch needs to exist on the host. For performance reasons we don't want to create a new batch here and copy to that, cudaMalloc is expensive.

        size_t N = batch.n_atoms;
        for (size_t isomer_idx = 0; isomer_idx < batch.isomer_capacity; isomer_idx++){   
            //Only insert the isomer if it has finished (either IsomerStatus::CONVERGED or IsomerStatus::FAILED)
            if(batch.statuses[isomer_idx] == IsomerStatus::EMPTY) continue;
            
            //Graphs always have a neighbour array.
            neighbours_t out_neighbours(N);

            std::vector<coord2d> output_2D;

            for (size_t i = 0; i < N; i++){
                //Fill in cubic cubic_neighbours
                out_neighbours[i] = {batch.cubic_neighbours[isomer_idx*N*3 + i*3], batch.cubic_neighbours[isomer_idx*N*3 + i*3 + 1], batch.cubic_neighbours[isomer_idx*N*3 + i*3 + 2]};
                
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
    
    template <BufferType T, BufferType U>
    cudaError_t copy(   IsomerBatch<T>& destination, //Copy data to this batch
                        const IsomerBatch<U>& source, //Copy data from this batch
                        const LaunchCtx& ctx, //Optional: specify which launch context to perform the copy operation in.
                        const LaunchPolicy policy, //Optional: specifies whether to synchronize the stream before and after copying)
                        const std::pair<int,int>& lhs_range, //Optional: provide a range of indices to assign similar to slices in numpy eg. {0,5} = [0:5]
                        const std::pair<int,int>& rhs_range //Optional: provide a range of indices to copy from similar to slices in numpy eg. {0,5} = [0:5]
                        ){
        assert(lhs_range.second - lhs_range.first == rhs_range.second - rhs_range.first);
        if ( lhs_range.second == -1 && lhs_range.first == -1 ) assert(source.isomer_capacity <= destination.isomer_capacity);
        
        if(U == DEVICE_BUFFER) {cudaSetDevice(source.get_device_id()); }
        else{ cudaSetDevice(destination.get_device_id());}
        
        //Iterate over the data fields of the IsomerBatch (pseudo reflection) and copy the contents of each using the provided stream.
        if(policy == LaunchPolicy::SYNC) {ctx.wait();}
        source.get_print_mode() ? destination.set_print_verbose() : destination.set_print_simple();
        for (size_t i = 0; i < source.pointers.size(); i++)
        {
            int num_isomers = (lhs_range.second > -1 && lhs_range.first > -1) ? lhs_range.second - lhs_range.first : source.isomer_capacity;
            int lhs_offset = lhs_range.first > 0 ? lhs_range.first * get<2>(source.pointers[i]) : 0;
            int rhs_offset = rhs_range.first > 0 ? rhs_range.first * get<2>(source.pointers[i]) : 0;
            char* lhs_ptr = (char*)(*(get<1>(destination.pointers[i])));
            char* rhs_ptr = (char*)(*(get<1>(source.pointers[i])));

            if( (source.get_device_id() != destination.get_device_id()) && ( (T==DEVICE_BUFFER) && (U == DEVICE_BUFFER))){
                cudaMemcpyPeerAsync(lhs_ptr + lhs_offset, destination.get_device_id(), rhs_ptr + rhs_offset, source.get_device_id(), num_isomers * get<2>(source.pointers[i]), ctx.stream);
            } else{
                cudaMemcpyAsync(lhs_ptr + lhs_offset, rhs_ptr + rhs_offset, num_isomers * get<2>(source.pointers[i]), cudaMemcpyKind(2*U +  T), ctx.stream);
            }
        }
        destination.n_isomers = source.n_isomers;
        printLastCudaError("cuda_io::Failed to copy struct");
        if(policy == LaunchPolicy::SYNC) {ctx.wait();}
        return cudaGetLastError();
    }

    template <BufferType T>
    cudaError_t free(IsomerBatch<T>& batch){
        cudaSetDevice(batch.get_device_id());
        for (int i = 0; i < batch.pointers.size(); i++){
            if(T == DEVICE_BUFFER) {cudaFree(*get<1>(batch.pointers[i]));} else {cudaFreeHost(*get<1>(batch.pointers[i]));}
        }
        return cudaGetLastError();
    }

    template <BufferType T>    
    cudaError_t resize(IsomerBatch<T>& batch, const size_t new_capacity, const LaunchCtx& ctx, const LaunchPolicy policy, int front, int back){
        cudaSetDevice(batch.get_device_id());
        //Construct a tempory batch: allocates the needed amount of memory.
        IsomerBatch<T> temp_batch = IsomerBatch<T>(batch.n_atoms, new_capacity, batch.get_device_id());
        //Copy contents of old batch into newly allocated memory.
        if (new_capacity < batch.isomer_capacity){
            cuda_io::copy(temp_batch, batch, ctx, LaunchPolicy::SYNC, {0,new_capacity}, {0,new_capacity});
        } else{
            cuda_io::copy(temp_batch, batch, ctx);
        }
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

    template <typename T>
    std::tuple<int,float> compare_isomer_meta(T* a, T* b, int n_isomers, int n_elements_per_isomers){
        int n_correct = 0;

        for (size_t i = 0; i < n_isomers; ++i){
            for (size_t j = 0; j < n_elements_per_isomers; j++)
            {
                auto idx = i*n_elements_per_isomers + j;
                if(a[idx] == b[idx]){ 
                    ++n_correct;
                }
            }
        }
        return {n_correct, (float)(100*n_correct)/(float)(n_isomers*n_elements_per_isomers)};
    }

    std::tuple<int, float, float> compare_isomer_arrays(float* a, float* b, int n_isomers, int n_elements_per_isomer, float rtol, bool verbose, float zero_threshold){
        float average_rdifference = (float)0.0;
        float average_adifference = (float)0.0;
        int n_correct = 0;
        int n_valid = 0;
        std::string fail_string{"\nFailed in: ["};
        for (int i = 0; i < n_isomers; ++i){
            bool isomer_passed = true;
            for (int j = 0; j < n_elements_per_isomer; ++j){   
                auto idx = i*n_elements_per_isomer + j;
                float rdiff = std::abs((a[idx] - b[idx])/(b[idx] + (std::abs(b[idx]) < zero_threshold) ));
                float adiff = std::abs(a[idx] - b[idx]);
                bool is_invalid = (std::isnan(a[idx]) || std::isnan(b[idx])) || (std::isinf(a[idx]) || std::isinf(b[idx]));
                if(rdiff <= rtol || a[idx] == b[idx] || (std::isnan(a[idx]) && std::isnan(b[idx])) || (std::isinf(a[idx]) && std::isinf(b[idx])) ){
                    ++n_correct;
                } else{
                    isomer_passed = false;
                }
                if (!is_invalid)
                {
                    average_rdifference += rdiff;
                    average_adifference += adiff;
                    ++n_valid;
                }
                
            }
            if(!isomer_passed && verbose) fail_string += to_string(i) + ", ";
        }
        if (n_correct < (n_isomers*n_elements_per_isomer) && verbose) {std::cout << fail_string << std::endl;}
        return {n_correct, average_adifference / (float)(n_valid), average_rdifference / (float)(n_valid)};
    }

    template <BufferType T>
    void is_close(const IsomerBatch<T>& a, const IsomerBatch<T>& b, float tol, bool verbose){

        if (! (a.isomer_capacity == b.isomer_capacity && a.n_atoms == b.n_atoms) ) {
            std::cout << "IsomerBatches are of different types and or dimensions \n";
        }{
            auto [n_correct_X,          accuracy0, avg_diff0] = compare_isomer_arrays(a.X, b.X, a.isomer_capacity, a.n_atoms*3, tol, verbose);
	    //            auto [n_correct_xys,        accuracy1, avg_diff1] = compare_isomer_arrays(a.xys, b.xys, a.isomer_capacity, a.n_atoms*2, tol, verbose);
            auto [n_correct_neighbours, accuracy2] = compare_isomer_meta(a.cubic_neighbours, b.cubic_neighbours, a.isomer_capacity, a.n_atoms*3);
            auto [n_correct_IDs,        accuracy3] = compare_isomer_meta(a.IDs, b.IDs, a.isomer_capacity, 1);
            auto [n_correct_statuses,   accuracy4] = compare_isomer_meta(a.statuses, b.statuses, a.isomer_capacity, 1);
            auto [n_correct_iterations, accuracy5] = compare_isomer_meta(a.iterations, b.iterations, a.isomer_capacity, 1);
            
            std::cout << "X                 | " << n_correct_X << "/" << a.isomer_capacity*a.n_atoms*3 << " | Acc:" << accuracy0 << "%\t | Avg Diff: " << avg_diff0 << "\n";
	    //            std::cout << "xys               | " << n_correct_xys << "/" << a.isomer_capacity*a.n_atoms*2 << " | Acc:" << accuracy1 << "%\t | Avg Diff: " << avg_diff1 << "\n";
            std::cout << "cubic_neighbours  | " << n_correct_neighbours << "/" << a.isomer_capacity*a.n_atoms*3 << " | Acc:" << accuracy2 << "%\t\n";
            std::cout << "IDs               | " << n_correct_IDs << "/" << a.isomer_capacity << " | Acc:" << accuracy3 << "%\t\n";
            std::cout << "statuses          | " << n_correct_statuses << "/" << a.isomer_capacity << " | Acc:" << accuracy4 << "%\t\n";
            std::cout << "iterations        | " << n_correct_iterations << "/" << a.isomer_capacity << " | Acc:" << accuracy5 << "%\t\n";
            
        }
    }

    template <BufferType T>
    int count_batch_status(const IsomerBatch<T>& input, const IsomerStatus status){
        auto fun = [&](const IsomerBatch<HOST_BUFFER>& batch){
            int res = 0;
            for (size_t i = 0; i < batch.isomer_capacity; i++)
            {
                if(batch.statuses[i] == status) ++res;
            }
            return res;
        }; 
        if (T == DEVICE_BUFFER){
            IsomerBatch<HOST_BUFFER> temp(input.n_atoms, input.isomer_capacity);
            cuda_io::copy(temp, input);
            return fun(temp);
        } else if(T == HOST_BUFFER){
            IsomerBatch<HOST_BUFFER> temp(input.n_atoms, input.isomer_capacity);
            cuda_io::copy(temp, input);
            return fun(temp);
        }
        return 0;
    }

    template <BufferType T>
    double average_iterations(const IsomerBatch<T>& input){
        auto fun = [](const IsomerBatch<HOST_BUFFER>& batch){
            u_int64_t result = 0;
            auto non_empty_isomers = 0;
            for (size_t i = 0; i < batch.isomer_capacity; i++) {
                result += batch.statuses[i] != IsomerStatus::EMPTY ? batch.iterations[i] : 0;
                if(batch.statuses[i] != IsomerStatus::EMPTY) ++non_empty_isomers;
            }
            return (double)result / (double)non_empty_isomers;
        };
        if (T == DEVICE_BUFFER){
            IsomerBatch<HOST_BUFFER> temp(input.n_atoms, input.isomer_capacity);
            cuda_io::copy(temp, input);
            return fun(temp);
        } else if(T == HOST_BUFFER){
            IsomerBatch<HOST_BUFFER> temp(input.n_atoms, input.isomer_capacity);
            cuda_io::copy(temp, input);
            return fun(temp);
        }
        return 0;
    }

    
    //Sorts the batch by whatever metadata key you pass, and if the key is identical for two isomers they are sorted by ID.
    template <BufferType T>
    void sort(IsomerBatch<T>& input_batch, const BatchMember key, const SortOrder order){
        IsomerBatch<HOST_BUFFER> temp(input_batch.n_atoms, input_batch.isomer_capacity);

        auto sort_fun = [&](IsomerBatch<HOST_BUFFER>& batch){
            std::map<int,int> index_map{};
            std::vector<int> lookup_indices(batch.isomer_capacity);
            std::iota(lookup_indices.begin(), lookup_indices.end(),0);
            switch (key)
            {
            case IDS:
                if (order == ASCENDING) std::sort(lookup_indices.begin(), lookup_indices.end(), [&](int a, int b) -> bool {return batch.IDs[a] < batch.IDs[b];});
                else if (order == DESCENDING) std::sort(lookup_indices.begin(), lookup_indices.end(), [&](int a, int b) -> bool {return batch.IDs[a] > batch.IDs[b];});
                break;

            case ITERATIONS:
                if (order == ASCENDING) std::sort(lookup_indices.begin(), lookup_indices.end(), [&](int a, int b) -> bool {return batch.iterations[a] == batch.iterations[b] ? batch.IDs[a] < batch.IDs[b] : batch.iterations[a] < batch.iterations[b];});
                if (order == DESCENDING) std::sort(lookup_indices.begin(), lookup_indices.end(), [&](int a, int b) -> bool {return batch.iterations[a] == batch.iterations[b] ? batch.IDs[a] > batch.IDs[b] : batch.iterations[a] > batch.iterations[b];});
                break;

            case STATUSES:
                if (order == ASCENDING) std::sort(lookup_indices.begin(), lookup_indices.end(), [&](int a, int b) -> bool {return batch.statuses[a] == batch.statuses[b] ? batch.IDs[a] < batch.IDs[b] : batch.statuses[a] < batch.statuses[b];});
                if (order == DESCENDING) std::sort(lookup_indices.begin(), lookup_indices.end(), [&](int a, int b) -> bool {return batch.statuses[a] == batch.statuses[b] ? batch.IDs[a] > batch.IDs[b] : batch.statuses[a] > batch.statuses[b];});
                break;

            default:
                if (order == ASCENDING) std::sort(lookup_indices.begin(), lookup_indices.end(), [&](int a, int b) -> bool {return batch.IDs[a] < batch.IDs[b];});
                if (order == DESCENDING) std::sort(lookup_indices.begin(), lookup_indices.end(), [&](int a, int b) -> bool {return batch.IDs[a] > batch.IDs[b];});
                break;
            }

            for (int i = 0; i < batch.isomer_capacity; i++){
                int offset      = lookup_indices[i] * batch.n_atoms*3;
                int face_offset = lookup_indices[i] * batch.n_faces;
                for (int j = 0; j < batch.n_atoms * 3; j++){
                    temp.X[i*batch.n_atoms*3 + j] = batch.X[offset + j];
                    temp.cubic_neighbours[i*batch.n_atoms*3 + j] = batch.cubic_neighbours[offset + j];
                }
                for (int j = 0; j < batch.n_faces; ++j){
                    temp.face_degrees[i*batch.n_faces + j] = batch.face_degrees[face_offset + j];
                    for (int k = 0; k < 6; k++) temp.dual_neighbours[(i*batch.n_faces + j)*6 + k] = batch.dual_neighbours[(face_offset + j)*6 + k];
                }
                offset = lookup_indices[i] * batch.n_atoms*2;
                for (size_t j = 0; j < batch.n_atoms*2; j++){
                    temp.xys[i*batch.n_atoms*2 + j] = batch.xys[offset +j];
                }
                temp.statuses[i] = batch.statuses[lookup_indices[i]];
                temp.iterations[i] = batch.iterations[lookup_indices[i]];
                temp.IDs[i] = batch.IDs[lookup_indices[i]];
            }
        };
        if (T == DEVICE_BUFFER){
            IsomerBatch<HOST_BUFFER> temp2(input_batch.n_atoms, input_batch.isomer_capacity);
            copy(temp2, input_batch);
            sort_fun(temp2);

        } else if(T == HOST_BUFFER){
            IsomerBatch<HOST_BUFFER> temp2(input_batch.n_atoms, input_batch.isomer_capacity);
            copy(temp2, input_batch);
            sort_fun(temp2);
        }
        copy(input_batch,temp);   
    }

    template <typename T>
    T mean(std::vector<T>& input){
      T result = input[0];
      for(size_t i=1;i<input.size();i++) result += input[i];
      
      return result/input.size();
    }
  
    //Standard deviation of a time vector. Used in benchmarking.
    std::chrono::nanoseconds sdev(std::vector<std::chrono::nanoseconds>& input){
    auto mean_ = mean(input);
    chrono::nanoseconds result = std::chrono::nanoseconds(0);
    for (int i = 0; i < input.size(); i++){
            result += (input[i] - mean_) *  (input[i] - mean_).count();
    }
    return std::chrono::nanoseconds((int)std::sqrt( (result / input.size()).count()));
    }

    std::chrono::microseconds sdev(std::vector<std::chrono::microseconds>& input){
    auto mean_ = mean(input);
    chrono::microseconds result = std::chrono::microseconds(0);
    for (int i = 0; i < input.size(); i++){
            result += (input[i] - mean_) *  (input[i] - mean_).count();
    }
    return std::chrono::microseconds((int)std::sqrt( (result / input.size()).count()));
    }

    double sdev(std::vector<double>& input){
        auto mean_ = mean(input);
        double result = 0;
        for (int i = 0; i < input.size(); i++){
                result += (input[i] - mean_) *  (input[i] - mean_);
        }
        return std::sqrt( (result / input.size()));
    }

    __global__
    void reset_convergence_status_(IsomerBatch<DEVICE_BUFFER> B){
        for (int isomer_idx = blockIdx.x; isomer_idx < B.isomer_capacity; isomer_idx+=gridDim.x)
        {
	  if (threadIdx.x == 0) B.statuses[isomer_idx] = (B.statuses[isomer_idx] == IsomerStatus::EMPTY)? IsomerStatus::EMPTY : IsomerStatus::NOT_CONVERGED;
	  if (threadIdx.x == 0) B.iterations[isomer_idx] = 0;
        }
        
    }
    template <BufferType T>
    cudaError_t reset_convergence_statuses(IsomerBatch<T>& B, const LaunchCtx& ctx, const LaunchPolicy policy){
        cudaSetDevice(ctx.get_device_id());
        static LaunchDims dims((void*)reset_convergence_status_, B.n_atoms);
        dims.update_dims((void*)reset_convergence_status_,B.n_atoms,0,B.isomer_capacity);
        if(policy == LaunchPolicy::SYNC) ctx.wait();
        void* kargs[]{(void*)&B};
        cudaLaunchCooperativeKernel((void*)reset_convergence_status_, dims.get_grid(), dims.get_block(), kargs, 0, ctx.stream);
        if(policy == LaunchPolicy::SYNC) ctx.wait();
        
        return cudaGetLastError();
    }

}

template <BufferType T>
void IsomerBatch<T>::append(const Graph& graph, const size_t id){
    if (m_size == isomer_capacity) {std::cout << "IsomerBatch is full" << std::endl; assert(false);}
    auto Nf = this->n_atoms / 2 + 2;
    if (graph.neighbours.size() != Nf) {std::cout << "Graph has wrong number of faces" + std::to_string(graph.neighbours.size()) + " " + std::to_string(Nf) << std::endl; assert(false);}
    size_t face_offset = m_size * Nf;
    for(node_t u=0;u< graph.neighbours.size();u++){
        face_degrees[face_offset + u] = graph.neighbours[u].size();
        for(int j=0;j<graph.neighbours[u].size();j++){
            dual_neighbours[6*(face_offset+u)+j] = graph.neighbours[u][j];
        }
    }
    statuses[m_size] = IsomerStatus::NOT_CONVERGED;
    iterations[m_size] = 0;
    IDs[m_size] = id;
    m_size++;
}

template <BufferType T>
void IsomerBatch<T>::append(const PlanarGraph& G, const size_t id, const bool insert_2d){ 
    //Extract the graph information (neighbours) from the PlanarGraph object and insert it at the appropriate location in the queue.
    size_t offset = m_size* n_atoms;
    for(node_t u=0;u<n_atoms;u++){
        for(int j=0;j<3;j++){
            cubic_neighbours[3*(offset+u)+j] = G.neighbours[u][j];
        }
        if(insert_2d){
            xys[2*(offset+u) + 0] = G.layout2d[u].first;
            xys[2*(offset+u) + 1] = G.layout2d[u].second; 
        }
    }
    statuses[m_size] = IsomerStatus::NOT_CONVERGED;
    iterations[m_size] = 0;
    IDs[m_size] = id;
    m_size++;
}

template <BufferType T>
void IsomerBatch<T>::append(const Polyhedron& graph, const size_t id){
    if (m_size == isomer_capacity) throw std::runtime_error("IsomerBatch is full");
    if (graph.neighbours.size() != n_atoms) throw std::runtime_error("Graph has wrong number of faces");
    size_t offset = m_size * n_atoms;
    for(node_t u=0;u< graph.neighbours.size();u++){
        for(int j=0;j<graph.neighbours[u].size();j++){
            cubic_neighbours[3*(offset+u)+j] = graph.neighbours[u][j];
            X[3*(offset+u)+j] = graph.points[u][j];
        }
    }
    statuses[m_size] = IsomerStatus::NOT_CONVERGED;
    iterations[m_size] = 0;
    IDs[m_size] = id;
    m_size++;
}

template <BufferType T>
void IsomerBatch<T>::clear(const LaunchCtx& ctx, const LaunchPolicy policy){
    //Set statuses to empty and iterations to 0
    if (T == DEVICE_BUFFER){
        cudaSetDevice(ctx.get_device_id());
        if(policy == LaunchPolicy::SYNC) ctx.wait();
        cudaMemsetAsync((void*)statuses, int(IsomerStatus::EMPTY), isomer_capacity * sizeof(IsomerStatus), ctx.stream);
        cudaMemsetAsync((void*)iterations, 0, isomer_capacity * sizeof(size_t), ctx.stream);
        if(policy == LaunchPolicy::SYNC) ctx.wait();
    } else if (T == HOST_BUFFER){
        std::fill(statuses, statuses + isomer_capacity, IsomerStatus::EMPTY);
        std::fill(iterations, iterations + isomer_capacity, 0);
    }
    m_size = 0;
}

template <BufferType T>
void IsomerBatch<T>::shrink_to_fit(){
    auto first_empty = -1;
    auto sort_fun = [&](IsomerBatch<HOST_BUFFER>& batch){
        cuda_io::sort(batch, STATUSES, DESCENDING);
        for (size_t i = 0; i < batch.isomer_capacity; i++){
            if (batch.statuses[i] == IsomerStatus::EMPTY){
                first_empty = i; break;
            }
        }
        if(first_empty == -1) first_empty = batch.isomer_capacity;
    };


    if (T == DEVICE_BUFFER){
        IsomerBatch<HOST_BUFFER> temp(n_atoms, isomer_capacity);
        cuda_io::copy(temp, *this);
        sort_fun(temp);
        cuda_io::copy(*this, temp);
        cuda_io::resize(*this, first_empty);
    } else if (T == HOST_BUFFER){
        IsomerBatch<HOST_BUFFER> temp(n_atoms, isomer_capacity);
        cuda_io::copy(temp, *this);
        sort_fun(temp);
        cuda_io::copy(*this, temp);
        cuda_io::resize(*this, first_empty);
    }
    
}
template <BufferType T>
std::optional<Polyhedron> IsomerBatch<T>::get_isomer(const size_t index) const {
    auto N = n_atoms;
    neighbours_t out_neighbours(N);
    if (index > isomer_capacity){
        std::cout << "Index out of range";
        return std::nullopt;
    }
    std::vector<coord2d> output_2D;

    for (size_t i = 0; i < N; i++){
        //Fill in cubic cubic_neighbours
        out_neighbours[i] = {cubic_neighbours[index*N*3 + i*3], cubic_neighbours[index*N*3 + i*3 + 1], cubic_neighbours[index*N*3 + i*3 + 2]};
        
    }

    std::vector<coord3d> output_X(N);
    for (size_t i = 0; i < N; i++){
        for (int j = 0; j < 3; ++j){
            output_X[i][j] = X[index*N*3 + i*3 + j];
        }
    }
    return Polyhedron(Graph(out_neighbours, true),output_X);
}

template <BufferType T>
std::optional<Polyhedron> IsomerBatch<T>::get_isomer_by_id(const size_t ID) const {
    auto N = n_atoms;
    neighbours_t out_neighbours(N, std::vector<int>(3, -1));
    int index = -1;
    for (size_t i = 0; i < isomer_capacity; i++) {if(IDs[i] == ID) index = i;}
    if (index == -1) {
        std::cout << "ID not in batch." << endl;
        return std::nullopt;
    }


    std::vector<coord2d> output_2D;

    for (size_t i = 0; i < N; i++){
        //Fill in cubic cubic_neighbours
        out_neighbours[i] = {cubic_neighbours[index*N*3 + i*3], cubic_neighbours[index*N*3 + i*3 + 1], cubic_neighbours[index*N*3 + i*3 + 2]};
        
    }

    std::vector<coord3d> output_X(N);
    for (size_t i = 0; i < N; i++){
        for (int j = 0; j < 3; ++j){
            output_X[i][j] = X[index*N*3 + i*3 + j];
        }
    }
    return Polyhedron(Graph(out_neighbours, true),output_X);
}

template <BufferType T>
std::vector<size_t> IsomerBatch<T>::find_ids(const IsomerStatus status){
    std::vector<size_t> result;
    for (size_t i = 0; i < isomer_capacity; i++){
        if(statuses[i] == status) result.push_back(IDs[i]);
    }
    std::sort(result.begin(), result.end());
    return result;
}

template <BufferType T>
bool IsomerBatch<T>::operator==(const IsomerBatch<T>& b){
    bool passed = true;
    auto device_real_t_equals = [](const float a, const float b){
        return (a == b) || (std::isnan(a) && std::isnan(b)) || (std::isinf(a) && std::isinf(b));
    };
    
    if (! (isomer_capacity == b.isomer_capacity && n_atoms == b.n_atoms) ) {

        return false;
    }else if(T == HOST_BUFFER){

        for(int i = 0; i < isomer_capacity; ++i){
            passed = passed && statuses[i] == b.statuses[i];
            passed = passed && IDs[i] == b.IDs[i];
            passed = passed && iterations[i] == b.iterations[i];
        }
        for(int i = 0; i < isomer_capacity * n_atoms * 3; ++i){
            passed = passed && device_real_t_equals(X[i],b.X[i]);
            passed = passed && cubic_neighbours[i] == b.cubic_neighbours[i];
        }
        for(int i = 0; i < isomer_capacity * n_atoms * 2; ++i){
            passed = passed && device_real_t_equals(xys[i],b.xys[i]);}

        //for(int i = 0; i < isomer_capacity * (n_atoms/2 + 1)*6; ++i){
        //    passed = passed && dual_neighbours[i] == b.dual_neighbours[i];
        //}
        return passed;
    } else{
        std::cout << "== operator only supported for HOST_BUFFER" << std::endl;
        return false;
    }
    return false;
}

template <BufferType T>
void IsomerBatch<T>::print(const BatchMember param, const std::pair<int,int>& range){
    auto print_fun = [&](const IsomerBatch<HOST_BUFFER>& a){
        std::ostream& os = std::cout;
        auto start = range.first < 0 ? 0 : range.first;
        auto end = range.second < 0 ? a.isomer_capacity : range.second;
        switch (param)
        {
        case IDS:
            os << "ID List: " << "\n[";
            for (int i = start; i < end - 1; i++){os << a.IDs[i] << ", ";}
            os << a.IDs[end-1] << "]\n"; 
            break;
        case STATUSES:
            os << "Status List: " << "\n[";
            for (int i = start; i < end - 1; i++){os << static_cast<int>(a.statuses[i]) << ", ";}
            os << static_cast<int>(a.statuses[end-1]) << "]\n"; 
            break;
        case ITERATIONS:
            os << "Iteration List: " << "\n[";
            for (int i = start; i < end - 1; i++){os << static_cast<int>(a.iterations[i]) << ", ";}
            os << a.iterations[end-1] << "]\n"; 
            break;
        case COORDS3D:
            os << "X Lists: " << "\n[";
            for (int i = start; i < end; i++){
                os << "[";
                for (int j = 0; j < a.n_atoms*3 - 1; j++)
                {
                    os << a.X[i*a.n_atoms*3 + j] << ", ";
                }
                if(i != (end - 1)) {
                    os << a.X[(i+1)*a.n_atoms*3 - 1] << "], ";
                } else{
                    os << a.X[(i+1)*a.n_atoms*3 - 1] << "]]\n";
                }
            }
            break;
        case COORDS2D:
            os << "2D Layouts: " << "\n[";
            for (int i = start; i < end; i++){
                os << "[";
                for (int j = 0; j < a.n_atoms*2 - 1; j++)
                {
                    os << a.xys[i*a.n_atoms*2 + j] << ", ";
                }
                if(i != (end - 1)) {
                    os << a.xys[(i+1)*a.n_atoms*2 - 1] << "], ";
                } else{
                    os << a.xys[(i+1)*a.n_atoms*2 - 1] << "]]\n";
                }
            }
            break;
        case CUBIC_NEIGHBOURS:
            os << "Cubic Neighbour Lists: " << "\n[";
            for (int i = start; i < end; i++){
                os << "[";
                for (int j = 0; j < a.n_atoms*3 - 1; j++)
                {
                    os << a.cubic_neighbours[i*a.n_atoms*3 + j] << ", ";
                }
                if(i != (end - 1)) {
                    os << a.cubic_neighbours[(i+1)*a.n_atoms*3 - 1] << "], ";
                } else{
                    os << a.cubic_neighbours[(i+1)*a.n_atoms*3 - 1] << "]]\n";
                }
            }
            break;
        case DUAL_NEIGHBOURS:
            os << "Dual Neighbour Lists: " << "\n[";
            for (int i = start; i < end; i++){
                os << "[";
                for (int j = 0; j  < a.n_faces; j++){
                    if (a.face_degrees[i*a.n_faces + j] <= 0) continue;
                    os << "[";
                    for (int k = 0; k < a.face_degrees[i*a.n_faces + j] - 1; k++)
                    {
                        os << a.dual_neighbours[(i*a.n_faces + j)*6 + k] << ",";
                    }
                    os << a.dual_neighbours[(i*a.n_faces + j)*6 + a.face_degrees[i*a.n_faces + j] - 1];
                    if (j != a.n_faces - 1)  os << "],";
                    else os << "]";
                }
                
                if(i != (end - 1)) {

                    os << "], ";
                } else{
                    os << "]]\n";
                }
            }
            break;

        default:
            break;
        }
    };
    if(T == DEVICE_BUFFER){
        IsomerBatch<HOST_BUFFER> temp(n_atoms, isomer_capacity);
        cuda_io::copy(temp, *this);
        print_fun(temp);
    } else if(T == HOST_BUFFER){
        IsomerBatch<HOST_BUFFER> temp(n_atoms, isomer_capacity);
        cuda_io::copy(temp, *this);
        print_fun(temp);
    }
}

/* template <BufferType T>
std::ostream& operator << (std::ostream& os, const IsomerBatch<T>& input){
    auto print_fun = [&](const IsomerBatch<HOST_BUFFER>& a){
        os << "Batch Dimensions: (" << a.isomer_capacity << ", " << a.n_atoms << ")\n";
        os << "ID List: " << "\n[";
        for (size_t i = 0; i < a.isomer_capacity - 1; i++){os << a.IDs[i] << ", ";}
        os << a.IDs[a.isomer_capacity-1] << "]\n"; 
        
        os << "Status List: " << "\n[";
        for (size_t i = 0; i < a.isomer_capacity - 1; i++){os << static_cast<int>(a.statuses[i]) << ", ";}
        os << static_cast<int>(a.statuses[a.isomer_capacity-1]) << "]\n"; 
        if (a.verbose){
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

            os << "Cubic Neighbour Lists: " << "\n[";
            for (size_t i = 0; i < a.isomer_capacity; i++){
                os << "[";
                for (size_t j = 0; j < a.n_atoms*3 - 1; j++)
                {
                    os << a.cubic_neighbours[i*a.n_atoms*3 + j] << ", ";
                }
                if(i != (a.isomer_capacity - 1)) {
                    os << a.cubic_neighbours[(i+1)*a.n_atoms*3 - 1] << "], ";
                } else{
                    os << a.cubic_neighbours[(i+1)*a.n_atoms*3 - 1] << "]]\n";
                }
            }

            os << "Dual Neighbour Lists: " << "\n[";
            for (size_t i = 0; i < a.isomer_capacity; i++){
                os << "[";
                for (int j = 0; j  < a.n_faces; j++){
                    if (a.face_degrees[i*a.n_faces + j] <= 0) continue;
                    os << "[";
                    for (size_t k = 0; k < a.face_degrees[i*a.n_faces + j] - 1; k++)
                    {
                        os << a.dual_neighbours[(i*a.n_faces + j)*6 + k] << ",";
                    }
                    os << a.dual_neighbours[(i*a.n_faces + j)*6 + a.face_degrees[i*a.n_faces + j] - 1];
                    if (j != a.n_faces - 1)  os << "],";
                    else os << "]";
                }
                
                if(i != (a.isomer_capacity - 1)) {

                    os << "], ";
                } else{
                    os << "]]\n";
                }
            }
        }
    };
    if (T == DEVICE_BUFFER){
        IsomerBatch<T> output(input.n_atoms, input.isomer_capacity);
        cuda_io::copy(output, input);
        print_fun(output);
    }else{
        print_fun(input);
    }

    return os;
} */



 
void template_instantiaons(){
    IsomerBatch<HOST_BUFFER> a;
    IsomerBatch<DEVICE_BUFFER> b;

    IsomerBatch<HOST_BUFFER> c(10,10);
    IsomerBatch<DEVICE_BUFFER> d(10,10, 0);

    a.set_print_simple(); b.set_print_simple();
    a.set_print_verbose(); b.set_print_verbose();
    a.get_print_mode(); b.get_print_mode();
    a.print(COORDS3D); b.print(COORDS3D);
    a.get_device_id(); b.get_device_id();
    a.size(); b.size();
    a.capacity(); b.capacity();
    a.get_isomer(0); b.get_isomer(0);
    a.get_isomer_by_id(0); b.get_isomer_by_id(0);
    a.find_ids(IsomerStatus::CONVERGED); b.find_ids(IsomerStatus::CONVERGED);
    a.shrink_to_fit(); b.shrink_to_fit();
    a.append(Graph(), 0); b.append(Graph(), 0);
    a.append(PlanarGraph(), 0); b.append(PlanarGraph(), 0);
    a.append(Polyhedron(), 0); b.append(Polyhedron(), 0);
    a.clear(); b.clear();
    a == a; b == b;
    a != a; b != b;
    //std::cout << c << d << std::endl;

}